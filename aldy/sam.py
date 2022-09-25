# 786
# Aldy source: sam.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple, Dict, List, Optional
from collections import defaultdict, Counter
from statistics import mean
import pysam
import os
import os.path
import gzip
import tarfile
import pickle
import tempfile

from .indelpost import Variant, VariantAlignment
from .common import log, GRange, AldyException, script_path, Timing, chr_prefix
from .gene import Gene, CNConfigType
from .coverage import Coverage
from .profile import Profile


class Sample:
    """Parse read alignments in a SAM/BAM/CRAM/VCF/dump format"""

    def __init__(
        self,
        gene: Gene,
        profile: Optional[Profile],
        path: str,
        reference: Optional[str] = None,
        debug: Optional[str] = None,
    ):
        """
        :param gene: Gene instance.
        :param profile: Profile instance. `None` if loading a dump file.
        :param path: Path to a SAM/BAM/CRAM/VCF/dump file.
        :param reference: Reference genome path for reading CRAM files.
            Default: None.
        :param debug: When set, create a `{debug}.dump` file for debug purposes.
            Default: None.

        :raise: :py:class:`aldy.common.AldyException` if the sample is invalid
            (e.g., the coverage of the copy-number neutral region is too low).
        """

        self.name = os.path.basename(path).split(".")[0]
        """Sample name."""

        self.gene = gene
        """Gene instance."""

        self.profile = profile
        """Profile instance."""

        self.path = path
        """File path."""

        self._dump_cn = defaultdict(int)
        """dict[int, int]: Coverage of the copy-number neutral region"""

        self._dump_reads: List[Tuple[Tuple, List]] = []
        """list[tuple[tuple, list]]: Read information."""

        self._indel_sites = {
            (pos, op): [0, 0] for pos, op in gene.mutations if op[:3] in ["ins", "del"]
        }
        self._indel_sites_eqs = {}

        self._multi_sites = {
            m.pos: m.op
            for _, a in gene.alleles.items()
            for m in a.func_muts
            if ">" in m.op and len(m.op) > 3
        }
        """Multi-substitutions (e.g., `A.C>T.G`)."""

        self.phaseable = {
            pos: i for i, pos in enumerate(sorted({pos for pos, _ in gene.mutations}))
        }
        """Locations that should be phased."""

        self.phases: Dict[str, Dict[int, str]] = {}
        """Phasing information."""

        self._fusion_counter: Dict = {}
        """Fusion read coverage (for long reads)."""

        self.is_long_read = False
        """Set if long-read data is used."""

        with Timing("[sam] Read SAM"):
            self.kind, _ = detect_genome(path)
            self.genome = gene.genome

            if self.kind == "vcf":
                try:
                    norm, muts = self._load_vcf(path)
                except ValueError:
                    raise AldyException(f"VCF {path} is not indexed")
            elif self.kind == "dump":
                norm, muts = self._load_dump(path)
            else:
                if self.profile and self.profile.sam_long_reads:
                    self.is_long_read = True
                    norm, muts = self._load_long_sam(path, reference, debug)
                else:
                    norm, muts = self._load_sam(path, reference, debug)
                if self.profile and self.profile.cn_region:
                    self._dump_cn = self._load_cn_region(
                        path, reference, self.profile.cn_region
                    )
            self._make_coverage(norm, muts)
            if self.kind == "sam" and debug:
                self._dump_alignments(f"{debug}.{gene.name}", norm, muts)

        assert self.profile, "profile not set"
        if self.profile.cn_region:
            self.coverage._normalize_coverage()
        log.debug("[sam] avg_coverage= {:.1f}x", self.coverage.average_coverage())
        if self.profile.cn_region and self.coverage.diploid_avg_coverage() < 2:
            raise AldyException(
                "The average coverage of the sample is too low ({:.1f}).".format(
                    self.coverage.diploid_avg_coverage()
                )
            )

    def _load_sam(self, sam_path: str, reference=None, debug=None):
        """
        Load the read, mutation and coverage data from a SAM/BAM/CRAM file.
        :raise: :py:class:`aldy.common.AldyException` if the file lacks an index.
        """

        log.debug("[sam] path= {}", os.path.abspath(sam_path))
        if reference:
            log.debug("[sam] reference= {}", os.path.abspath(reference))
        assert self.profile, "profile not set"

        norm: dict = defaultdict(list)
        muts: dict = defaultdict(list)

        with pysam.AlignmentFile(  # type: ignore
            sam_path, reference_filename=reference
        ) as sam:
            # Check do we have proper index to speed up the queries
            try:
                has_index = sam.check_index()
            except AttributeError:
                # SAM files do not have an index. BAMs might also lack it
                has_index = False
            except ValueError:
                raise AldyException(f"Cannot check index of {sam_path}")

            self._prefix = chr_prefix(
                self.gene.chr, [x["SN"] for x in sam.header["SQ"]]
            )

            with tempfile.TemporaryDirectory() as tmp:
                self._realign_indels(tmp, sam, reference)
            if has_index:
                iter = sam.fetch(
                    region=self.gene.get_wide_region().samtools(prefix=self._prefix)
                )
            else:
                log.warn("SAM/BAM index not found. Reading will be slow.")
                iter = sam.fetch()
            # Fetch the reads
            for read in iter:
                if not read.cigartuples:  # only valid alignments
                    continue
                if read.is_supplementary:  # avoid supplementary alignments
                    continue
                if "H" in read.cigarstring:  # avoid hard-clipped reads
                    continue
                if not read.query_sequence:
                    continue
                # ensure that it is a proper gene read
                if not _in_region(self.gene.get_wide_region(), read, self._prefix):
                    continue

                # Handle 10X Genomics tags
                if read.has_tag("BX"):
                    fragment = read.get_tag("BX")
                    if read.has_tag("XC"):
                        fragment = f"{fragment}:{read.get_tag('XC')}"
                    elif read.has_tag("MI"):
                        fragment = f"{fragment}:{read.get_tag('MI')}"
                    self.is_long_read = True
                else:
                    fragment = read.query_name
                r = self._parse_read(
                    fragment,
                    read.reference_start,
                    read.cigartuples,
                    read.query_sequence,
                    norm,
                    muts,
                    read.mapping_quality,
                    read.query_qualities,
                )
                if r and debug:
                    self._dump_reads.append(r)
        return norm, muts

    def _load_vcf(self, vcf_path: str):
        """Load the read, mutation and coverage data from a VCF file."""

        log.debug("[vcf] path= {}", os.path.abspath(vcf_path))

        norm = {
            p: [(40, 40)] * 20
            for p in range(
                self.gene.get_wide_region().start - 500,
                self.gene.get_wide_region().end + 1,
            )
        }
        muts: dict = defaultdict(list)

        def get_mut(pos, ref, alt):
            off = 0
            while off < len(ref) and off < len(alt) and ref[off] == alt[off]:
                off += 1
            if len(ref) - off == 1 and len(alt) - off == 1:
                if alt[off] == self.gene[off + pos]:
                    return off + pos, "_"
                return off + pos, f"{self.gene[off + pos]}>{alt[off]}"
            elif len(ref) > len(alt) and len(alt) - off == 0:
                return off + pos, f"del{self.gene[off + pos : pos + len(ref)]}"
            elif len(ref) < len(alt) and len(ref) - off == 0:
                return off + pos, f"ins{alt[off:]}"
            else:
                log.trace(f"[sam] ignoring {pos}: {ref}->{alt}")
                return pos, None

        with pysam.VariantFile(vcf_path) as vcf:  # type: ignore
            self._prefix = chr_prefix(self.gene.chr, list(vcf.header.contigs))

            samples = list(vcf.header.samples)
            self.name = sample = samples[0]
            if len(samples) > 1:
                log.warn("WARNING: Multiple VCF samples found; using the first one.")
            log.info("Using VCF sample {}", sample)
            for read in vcf.fetch(
                region=self.gene.get_wide_region().samtools(prefix=self._prefix)
            ):
                g = sorted(y for y in read.samples[sample]["GT"] if y is not None)
                if len(g) != 2 or self.gene[read.pos - 1] == "N":
                    continue  # ignore polyploid and incomplete cases
                dump_arr = {}
                if len(read.ref) == 1 and read.ref != self.gene[read.pos - 1]:
                    hgvs = [(read.pos - 1, f"{self.gene[read.pos - 1]}>{read.ref}")]
                else:
                    hgvs = [(read.pos - 1, "_")]
                hgvs += [get_mut(read.pos - 1, read.ref, a) for a in read.alleles[1:]]
                for gt in g:
                    pos, op = hgvs[gt]
                    if op == "_":
                        continue
                    muts[pos, op] += [(40, 40)] * 10
                    norm[pos] = norm[pos][:-10]
                    dump_arr[pos] = op

                # Handle multi-SNPs
                for pos, op in self._multi_sites.items():
                    if pos not in dump_arr:
                        continue
                    l, r = op.split(">")
                    if all(
                        dump_arr.get(pos + p, "-") == f"{l[p]}>{r[p]}"
                        for p in range(len(l))
                        if l[p] != "."
                    ):
                        for p in range(len(l)):
                            if l[p] != ".":
                                np = pos + p, f"{l[p]}>{r[p]}"
                                muts[np] = muts[np][:-10]
                                if p:
                                    norm[pos + p] += [(40, 40)] * 10
                        muts[pos, op] += [(40, 40)] * 10
        return norm, muts

    def _load_dump(self, dump_path: str):
        """Load Aldy dump data."""

        log.debug("[dump] path= {}", os.path.abspath(dump_path))
        if dump_path.endswith(".tar.gz"):
            tar = tarfile.open(dump_path, "r:gz")

            f = [i for i in tar.getnames() if i.endswith(f".{self.gene.name}.dump")]
            if not f:
                raise AldyException("Invalid dump file")
            log.debug("Found {} in the archive", f[0])
            data = tar.extractfile(f[0])
            assert data, "Malformed dump file"
            fd = gzip.open(data)
        else:
            fd = gzip.open(dump_path, "rb")

        log.warn("Loading debug dump from {}", dump_path)
        (
            self.name,
            self.profile,
            self._dump_cn,
            norm,
            muts,
            phases,
            self._fusion_counter,
            self._indel_sites,
        ) = pickle.load(
            fd
        )  # type: ignore
        self.phases = {f"r{i}": v for i, v in enumerate(phases)}
        norm = {p: [q for q, n in c.items() for _ in range(n)] for p, c in norm.items()}
        muts = {p: [q for q, n in c.items() for _ in range(n)] for p, c in muts.items()}
        return norm, muts

    def _dump_alignments(self, debug: str, norm, muts):
        """Pickle alignment data for debug purposes."""
        with open(f"{debug}.genome", "w") as fd:
            print(self.gene.genome, file=fd)
        with gzip.open(f"{debug}.dump", "wb") as fd:  # type: ignore
            pickle.dump(
                (
                    self.name,
                    self.profile,
                    self._dump_cn,
                    {p: Counter(q) for p, q in norm.items()},
                    {p: Counter(q) for p, q in muts.items()},
                    [v for v in self.phases.values() if len(v) > 1],
                    self._fusion_counter,
                    self._indel_sites,  # TODO: remove
                ),
                fd,  # type: ignore
            )

    def _realign_indels(self, tmp, sam, reference, long_reads=False):
        """Realign reads around database indels via indelpost module."""

        assert self.profile, "profile not loaded"

        rname = f"{self._prefix}{self.gene.chr}"
        if not reference:
            reference = f"{tmp}/ref.fa"
            with open(reference, "w") as f:
                print(f">{rname}", file=f)
                print("N" * self.gene._lookup_range[0], end="", file=f)
                print(self.gene._lookup_seq, end="", file=f)
                sz = sam.get_reference_length(rname)
                print("N" * (sz - self.gene._lookup_range[1]), file=f)
            with open(f"{reference}.fai", "w") as f:
                print(rname, sz, len(rname) + 2, sz, sz + 1, sep="\t", file=f)
        ref = pysam.FastaFile(reference)  # type: ignore

        prev_indel = None
        for pos, op in sorted(self._indel_sites, key=lambda x: (x[0], -len(x[1]))):
            p = pos
            if op.startswith("del"):
                if "ins" in op:
                    pd, pi = op[3:].split("ins")
                    o1, o2 = pd, pi
                else:
                    p -= 1
                    o = self.gene[p]
                    o1, o2 = o + op[3:], o
            else:
                o = self.gene[p]
                o1, o2 = o, o + op[3:]

            v = Variant(rname, p + 1, o1, o2, ref)  # type: ignore

            if long_reads:
                # Speed-up: just generate equivalent indels, no need for the realignment
                for ev in v.generate_equivalents():
                    np, no = ev.pos - 1, ""
                    if len(ev.ref) < len(ev.alt) and ev.alt.startswith(ev.ref):
                        np += len(ev.ref)
                        no = "ins" + ev.alt[len(ev.ref) :]
                    elif len(ev.ref) > len(ev.alt) and ev.ref.startswith(ev.alt):
                        np += len(ev.alt)
                        no = "del" + ev.ref[len(ev.alt) :]
                    if no:
                        self._indel_sites_eqs[np, no] = (pos, op)
                continue

            valn = VariantAlignment(  # type: ignore
                v,
                sam,
                mapping_quality_threshold=self.profile.min_mapq,
                base_quality_threshold=self.profile.min_quality,
                # needed to account for indel and database errors
                exact_match_for_shiftable=True,
            )

            phased = valn.phase()
            if len(phased.ref) - len(phased.alt) != len(v.ref) - len(v.alt):
                continue  # HACK: this indicates a subsumed indel

            self._indel_sites[pos, op] = list(valn.count_alleles())
            on_target = {r.query_name for r in valn.fetch_reads("target")}
            if (
                prev_indel
                and pos == prev_indel[0]
                and op.startswith("ins")
                and prev_indel[1].startswith(op)  # handle indel repeats
            ):
                # do not consider previously used reads
                off, on = self._indel_sites[pos, op]
                x = len(on_target & prev_indel[2])
                self._indel_sites[pos, op] = [off + x, on - x]
            if self._indel_sites[pos, op][1]:
                log.debug(
                    "[indel] {}:{} -> {}", pos + 1, op, self._indel_sites[pos, op]
                )
                prev_indel = (pos, op, on_target)
        # assert False

    def _load_cn_region(self, path, reference, cn_region):
        """
        Load the copy-number neutral coverage data from a SAM/BAM/CRAM file.

        :param cn_region: Copy-number neutral region to be used for coverage rescaling.
        """
        self._dump_cn = defaultdict(int)
        with pysam.AlignmentFile(  # type: ignore
            path, reference_filename=reference
        ) as sam:
            # Check do we have proper index to speed up the queries
            try:
                has_index = sam.check_index()
            except AttributeError:
                # SAM files do not have an index. BAMs might also lack it
                has_index = False
            except ValueError:
                raise AldyException(f"Cannot check index of {self.path}")
            self._prefix = chr_prefix(
                self.gene.chr, [x["SN"] for x in sam.header["SQ"]]
            )
            # Set it to _fetched_ if a CN-neutral region is user-provided.
            # Then read the CN-neutral region.
            if has_index:
                iter = sam.fetch(region=cn_region.samtools(prefix=self._prefix))
            else:
                log.warn("SAM/BAM index not found. Reading will be slow.")
                iter = sam.fetch()
            for read in iter:
                if _in_region(cn_region, read, self._prefix):
                    start = read.reference_start
                    if read.cigartuples is None:
                        continue
                    if read.is_supplementary:
                        continue
                    for op, size in read.cigartuples:
                        if op in [0, 7, 8, 2]:
                            for i in range(size):
                                self._dump_cn[start + i] += 1
                            start += size
        return self._dump_cn

    def _make_coverage(self, norm, muts):
        """Populate coverage data."""

        coverage: Dict[int, Dict[str, List]] = dict()
        for pos, cov in norm.items():
            if len(cov) == 0:
                continue
            coverage.setdefault(pos, {})["_"] = cov
        bounds = min(self.gene.chr_to_ref), max(self.gene.chr_to_ref)
        for (pos, mut), cov in muts.items():
            if pos not in coverage:
                coverage[pos] = {}
            if not bounds[0] <= pos <= bounds[1] and mut[:3] != "ins":
                mut = "_"  # ignore mutations outside of the region of interest
            coverage.setdefault(pos, {}).setdefault(mut, []).extend(cov)
        for pos, op in self._multi_sites.items():
            if pos in coverage and op in coverage[pos]:
                log.debug(f"[sam] multi-SNP {pos}{op}: {len(coverage[pos][op])} reads")
        assert self.profile, "profile not set"
        self.coverage = Coverage(
            self.gene,
            self.profile,
            self,
            {p: {m: v for m, v in coverage[p].items() if len(v) > 0} for p in coverage},
            self._indel_sites,
            self._dump_cn,
        )
        """Sample coverage data."""

        return self.coverage

    def _parse_read(
        self,
        fragment,
        ref_start,
        cigar,
        seq,
        norm,
        muts,
        mq=None,
        qual=None,
    ):
        """
        Parse a read.

        :param gene: Gene instance.
        :param fragment: Read name.
        :param ref_start: Start position in the reference genome.
        :param cigar: List of CIGAR operations in tuple format (operation, size).
        :param seq: Read sequence.
        :param norm: Positions within the read that have not been mutated.
        :param muts: Positions within the read that have been mutated.
        :param mq: Mapping quality.
        :param seq: Read qualities (normalized).

        .. note:: `norm` and `muts` are modified.

        .. note:: Qualities will be binned.
        """

        def bin_quality(q):
            """
            Quality score binning.
            See https://www.illumina.com/content/dam/illumina-marketing/docum ents/products/technotes/technote_understanding_quality_scores.pdf  # noqa
            """
            if q < 2:
                return int(q)
            if q < 10:
                return 6
            if q < 20:
                return 15
            if q < 29:
                return 25
            if q < 39:
                return 35
            return 40

        phase = self.phases.setdefault(fragment, {})
        dump_arr = []
        start, s_start = ref_start, 0
        prev_q = 10
        for op, size in cigar:
            if op == 2:  # Deletion
                mut = (start, "del" + self.gene[start : start + size])
                for i in range(size):
                    muts[start + i, "-"].append((bin_quality(mq), bin_quality(prev_q)))
                dump_arr.append(mut)
                if start in self.phaseable:
                    phase[start] = mut[1]
                if self._indel_sites_eqs and mut in self._indel_sites_eqs:
                    self._indel_sites[self._indel_sites_eqs[mut]][1] += 1
                start += size
            elif op == 1:  # Insertion
                mut = (start, "ins" + seq[s_start : s_start + size])
                q = mean(qual[s_start : s_start + size]) if qual else prev_q
                muts[mut].append((bin_quality(mq), bin_quality(q)))
                prev_q = q
                dump_arr.append(mut)
                if start in self.phaseable:
                    phase[start] = mut[1]
                if self._indel_sites_eqs and mut in self._indel_sites_eqs:
                    self._indel_sites[self._indel_sites_eqs[mut]][1] += 1
                s_start += size
            elif op == 4:  # Soft-clip
                s_start += size
            elif op in [0, 7, 8]:  # M, X and =
                for i in range(size):
                    q = qual[s_start + i] if qual else prev_q
                    if (
                        start + i in self.gene
                        and self.gene[start + i] != seq[s_start + i]
                    ):
                        mut = (start + i, f"{self.gene[start + i]}>{seq[s_start + i]}")
                        dump_arr.append(mut)
                        muts[mut].append((bin_quality(mq), bin_quality(q)))
                        if start + i in self.phaseable:
                            phase[start + i] = mut[1]
                    else:  # We ignore all mutations outside the RefSeq region
                        norm[start + i].append((bin_quality(mq), bin_quality(q)))
                        if start + i in self.phaseable:
                            phase[start + i] = "_"
                    prev_q = q
                start += size
                s_start += size

        dump_arr_pos = {p for p, _ in dump_arr}
        for pos, op in self._multi_sites.items():
            if pos not in dump_arr_pos:
                continue
            l, r = op.split(">")
            if all(
                (pos + p, f"{l[p]}>{r[p]}") in dump_arr
                for p in range(len(l))
                if l[p] != "."
            ):
                items = []
                for p in range(len(l)):
                    if l[p] != ".":
                        items.append(muts[pos + p, f"{l[p]}>{r[p]}"].pop())
                        if p:  # no idea why...
                            norm[pos + p].append(items[-1])
                # TODO: use sth else instead of mean?
                muts[pos, op].append(
                    (mean(mq for mq, _ in items), mean(q for _, q in items))
                )

        for pos, op in self._indel_sites:
            if ref_start <= pos < start:
                self._indel_sites[pos, op][0] += 1
        read_pos = (ref_start, start, len(seq))
        return read_pos, dump_arr

    def _load_long_sam(self, sam_path: str, reference=None, debug=None):
        """
        Load the read, mutation and coverage data from a long-read SAM/BAM/CRAM file.
        This function remaps long reads for selected genes.
        Tested on PacBio HiFi data only.

        .. note:: Currently only supports CYP2D6.
        """

        log.debug("[sam] pacbio_path= {}", os.path.abspath(sam_path))
        assert self.genome, "Genome not provided"
        assert self.genome == "hg38", "Only hg38 supported at this moment"
        assert self.profile, "profile not set"

        norm: dict = defaultdict(list)
        muts: dict = defaultdict(list)

        self._index = None
        ref = script_path(f"aldy.resources.genes/{self.gene.name.lower()}.fa.gz")
        if os.path.exists(ref):
            import mappy

            self._index = mappy.Aligner(ref, preset=self.profile.sam_mappy_preset)

            seq_name = self._index.seq_names[0]
            self._seq_range = list(map(int, seq_name.split(":")[1].split("-")))
            self._index_gene = []
            for _, g in enumerate(self.gene.regions):
                st = min(r.start for r in g.values())
                ed = max(r.end for r in g.values())
                idx = mappy.Aligner(
                    seq=self._index.seq(seq_name)[
                        st - self._seq_range[0] : ed - self._seq_range[0]
                    ],
                    preset=self.profile.sam_mappy_preset,
                )
                self._index_gene.append((idx, st))

            self._fusion_signatures = {}
            self._fusion_counter = {}
            for n, cn in self.gene.cn_configs.items():
                if cn.kind == CNConfigType.LEFT_FUSION:
                    s, pos = 0, next(i for i, j in cn.cn[0].items() if j == 1)
                elif cn.kind == CNConfigType.RIGHT_FUSION:
                    s, pos = 1, next(i for i, j in cn.cn[0].items() if j == 0)
                else:
                    continue
                sig = (
                    (1 - s, list(cn.cn[0])[list(cn.cn[0]).index(pos) - 1]),
                    (s, pos),
                )
                self._fusion_counter[n] = [0, 0]
                if (
                    sig[0][1] in self.gene.unique_regions
                    or sig[1][1] in self.gene.unique_regions
                ):
                    self._fusion_signatures[sig] = n

            log.debug("[sam] PacBio remapping enabled")

        with pysam.AlignmentFile(  # type: ignore
            sam_path, reference_filename=reference
        ) as sam, Timing("[sam] Remap"):
            # Assumes SAM index exists
            self._prefix = chr_prefix(
                self.gene.chr, [x["SN"] for x in sam.header["SQ"]]
            )

            with tempfile.TemporaryDirectory() as tmp:
                self._realign_indels(tmp, sam, reference, True)
                for po, (off, on) in self._indel_sites.items():
                    self._indel_sites[po] = [off - on, on]

            wide = self.gene.get_wide_region()
            end = wide.end
            if self.gene.name == "CYP2D6":
                end = {"hg19": 42_551_000, "hg38": 42_155_000}[self.genome]
            rng = GRange(wide.chr, wide.start, end)
            counter = 0
            iter = sam.fetch(region=rng.samtools(prefix=self._prefix))
            for read in iter:
                if not read.cigartuples or "H" in read.cigarstring:
                    continue
                if read.reference_start >= rng.end or read.reference_end < rng.start:
                    continue

                seq = read.query_sequence
                qual = read.query_qualities if read.query_qualities else [40] * len(seq)
                if self._index:
                    # Split reads
                    pieces = []
                    s_start = 0
                    gene_bnd = self.gene.regions[0]["up"].end - self._seq_range[0]
                    pseudo_bnd = self.gene.regions[1]["up"].end - self._seq_range[0]

                    while s_start < len(seq):
                        s = seq[s_start:]
                        q = qual[s_start:]
                        h = self._map(self._index, s)
                        if (
                            h
                            and h.r_st <= gene_bnd < h.r_en
                            and gene_bnd - h.r_st > 100
                        ):
                            s = seq[s_start : s_start + (gene_bnd - h.r_st)]
                            q = qual[s_start : s_start + (gene_bnd - h.r_st)]
                            h = self._map(self._index, s)
                        elif (
                            h
                            and h.r_st <= pseudo_bnd < h.r_en
                            and pseudo_bnd - h.r_st > 100
                        ):
                            s = seq[s_start : s_start + (pseudo_bnd - h.r_st)]
                            q = qual[s_start : s_start + (pseudo_bnd - h.r_st)]
                            h = self._map(self._index, s)
                        if not h:
                            s_start += 100
                        else:
                            pcs, ed = self._split_read(s, q, h)
                            pieces += pcs
                            s_start += ed
                    for ref_start, seq, qual, cigar in pieces:
                        r = self._parse_read(
                            counter,
                            ref_start,
                            cigar,
                            seq,
                            norm,
                            muts,
                            read.mapping_quality,
                            qual,
                        )
                        counter += 1
                        if r and debug:
                            self._dump_reads.append(r)
                else:
                    r = self._parse_read(
                        counter,
                        read.reference_start,
                        read.cigartuples,
                        seq,
                        norm,
                        muts,
                        read.mapping_quality,
                        read.query_qualities,
                    )
                    counter += 1
                    if r and debug:
                        self._dump_reads.append(r)
        return norm, muts

    def _map(self, idx, seq):
        """
        Map a sequence to a mappy reference index.
        :returns: The first primary hit or `None` if no primary hits are found.
        """

        for hit in idx.map(seq):
            if hit.is_primary:
                return hit
        return None

    def _split_read(self, seq, qual, hit):
        """Split a long read."""

        def _rev_cigar(c):
            """Reverse CIGAR tuples for pysam/mappy compatibility."""
            return [(op, size) for (size, op) in c]

        regs = self._get_gene_regions(
            hit.r_st + self._seq_range[0] - 1, _rev_cigar(hit.cigar)
        )
        if regs:
            splits = []
            for fs, fn in self._fusion_signatures.items():
                brk = [r for r in regs if r[0][1] == fs[1][1]]
                if not brk:
                    continue
                self._fusion_counter[fn][1] += 1
                brk = brk[0]
                lg, brk, rg = fs[1][0], brk[1] + brk[2], fs[0][0]
                lh = self._map(self._index_gene[lg][0], seq[:brk])
                rh = self._map(self._index_gene[rg][0], seq[brk:])
                if lh and rh:
                    unmap = len(seq) - (lh.q_en - lh.q_st) - (rh.q_en - rh.q_st)
                    unmap += lh.NM + rh.NM
                    splits.append([fn, unmap, brk, lh, rh, lg, rg])
            splits.sort(key=lambda x: x[1])
            if splits and splits[0][1] <= hit.NM:
                fn, _, brk, lh, rh, lg, rg = splits[0]
                self._fusion_counter[fn][0] += 1
                return (
                    [
                        (
                            lh.r_st + self._index_gene[lg][1] - 1,
                            seq[lh.q_st : lh.q_en],
                            qual[lh.q_st : lh.q_en],
                            _rev_cigar(lh.cigar),
                        ),
                        (
                            rh.r_st + self._index_gene[rg][1] - 1,
                            seq[brk + rh.q_st : brk + rh.q_en],
                            qual[brk + rh.q_st : brk + rh.q_en],
                            _rev_cigar(rh.cigar),
                        ),
                    ],
                    max(lh.q_en, brk + rh.q_en),
                )
        return [
            (
                hit.r_st + self._seq_range[0] - 1,
                seq[hit.q_st : hit.q_en],
                qual[hit.q_st : hit.q_en],
                _rev_cigar(hit.cigar),
            )
        ], hit.q_en

    def _get_gene_regions(self, r_start, cigar):
        """Get gene regions that cover the given CIGAR range."""

        regs = []
        start, s_start = r_start, 0
        for op, size in cigar:
            if op == 1 or op == 4:
                s_start += size
            elif op == 2 or op == 3:
                start += size
            elif op == 0 or op == 7 or op == 8:
                for i in range(size):
                    rg = self.gene.region_at(start + i)
                    if rg:
                        if not regs or regs[-1][0] != rg:
                            regs.append([rg, s_start + i, 0])
                        regs[-1][2] += 1
                start += size
                s_start += size
        return regs


def detect_genome(sam_path: str) -> Tuple[str, Optional[str]]:
    """Detect file type and the reference genome."""

    try:
        pysam.set_verbosity(0)  # type: ignore
        with pysam.AlignmentFile(sam_path) as sam:  # type: ignore
            try:
                sam.check_index()
            except AttributeError:
                pass  # SAM files do not have an index. BAMs might also lack it
            except ValueError:
                raise AldyException(
                    f"File {sam_path} has no index (it must be indexed)"
                )

            #: str: Check if we need to append 'chr' or not
            _prefix = chr_prefix("1", [x["SN"] for x in sam.header["SQ"]])
            # Detect the genome
            hg19_lens = {"1": 249250621, "10": 135534747, "22": 51304566}
            hg38_lens = {"1": 248956422, "10": 133797422, "22": 50818468}
            chrs = {x["SN"]: int(x["LN"]) for x in sam.header["SQ"]}
            is_hg19, is_hg38 = True, True
            for c in "1 10 22".split():
                if _prefix + c in chrs:
                    is_hg19 &= chrs[_prefix + c] == hg19_lens[c]
                    is_hg38 &= chrs[_prefix + c] == hg38_lens[c]
            if is_hg19:
                return "sam", "hg19"
            elif is_hg38:
                return "sam", "hg38"
            else:
                return "sam", None
    except (ValueError, OSError):
        try:
            with pysam.VariantFile(sam_path):  # type: ignore
                return "vcf", None
        except (ValueError, OSError):
            if sam_path.endswith(".tar.gz"):
                tar = tarfile.open(sam_path, "r:gz")
                f = [i for i in tar.getnames() if i.endswith(".genome")]
                if not f:
                    raise AldyException("Invalid dump file")
                data = tar.extractfile(f[0])
                if data:
                    genome = data.read().decode("utf-8").strip()
                    return "dump", genome
    return "", None


def _in_region(region: GRange, read, prefix: str) -> bool:
    """Check if a read intersects a given gene region."""

    if read.reference_id == -1 or read.reference_name != prefix + region.chr:
        return False
    if read.reference_end is None:
        return False

    a = (read.reference_start, read.reference_end)
    b = (region.start, region.end)
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]
