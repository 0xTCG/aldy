# 786
# Aldy source: sam.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import collections
from typing import Tuple, Dict, List, Optional, Any, Set

import pysam
import os
import os.path
import yaml
import gzip
import tarfile
import pickle


from dataclasses import dataclass
from collections import defaultdict, Counter
from natsort import natsorted
from statistics import mean
from .common import log, GRange, AldyException, script_path, Timing, chr_prefix
from .gene import Gene, Mutation, CNConfigType
from .coverage import Coverage
from .profile import Profile


@dataclass
class Sample:
    """
    Interface for reading SAM/BAM/CRAM files that parses and stores read alignments.
    """

    coverage: Coverage
    name: str
    _prefix: str
    _insertion_sites: Set[Mutation]
    _insertion_counts: Dict[Tuple[int, int], int]
    _insertion_reads: Dict[Mutation, Dict[Any, int]]
    _multi_sites: Dict[int, str]

    # Dump data
    _dump_cn: Dict[int, int]  # pos -> coverage
    _dump_reads: List[Tuple[Tuple[int, int, int], List[Mutation]]]

    def __init__(
        self,
        gene: Gene,
        profile: Profile,
        sam_path: Optional[str] = None,
        reference: Optional[str] = None,
        debug: Optional[str] = None,
        vcf_path: Optional[str] = None,
    ):
        """
        Initialize a :obj:`Sample` object.

        :param sam_path: Path to a SAM/BAM/CRAM file.
        :param gene: Gene instance.
        :param profile:
            Profile data (:obj:`Profile`) loaded from a YML profile (e.g. 'prgnseq-v1')
            or a profile SAM/BAM file.
        :param reference:
            Reference genome for reading CRAM files.
            Default is None.
        :param debug:
            If set, create a "`debug`.dump" file for debug purposes.
            Default is None.

        :raise: :obj:`aldy.common.AldyException` if the average coverage of the
                copy-number neutral region is too low (less than 2).
        """

        if sam_path:
            self.name = os.path.basename(sam_path).split(".")[0]
        else:
            self.name = ""

        # Get the list of indel sites that should be corrected
        # TODO: currently uses only functional indels;
        #       other indels should be corrected as well
        self._dump_cn = defaultdict(int)
        self._dump_reads = []
        self._insertion_sites = {
            m for a in gene.alleles.values() for m in a.func_muts if m.op[:3] == "ins"
        }
        self._insertion_counts: Dict[Tuple[int, int], int] = defaultdict(int)
        self._insertion_reads: Dict[Mutation, Dict[Any, int]] = {
            m: defaultdict(int) for m in self._insertion_sites
        }
        self._multi_sites = {
            m.pos: m.op
            for _, a in gene.alleles.items()
            for m in a.func_muts
            if ">" in m.op and len(m.op) > 3
        }
        self.gene = gene

        # Phasing information!
        self.grid_columns = sorted({pos for pos, _ in gene.mutations})
        self.grid_columns = {pos: i for i, pos in enumerate(self.grid_columns)}
        self.phases = {}
        self.moved = {}

        self.profile = profile
        self._fusion_counter = {}

        with Timing("[sam] Read SAM"):
            is_sam = False
            group_indels = False
            assert vcf_path or sam_path
            if vcf_path:
                self.path = vcf_path
                try:
                    norm, muts = self._load_vcf(vcf_path, gene)
                except ValueError:
                    raise AldyException(f"VCF {vcf_path} is not indexed")
            elif sam_path and sam_path.endswith(".dump"):
                self.path = sam_path
                norm, muts = self._load_dump(sam_path, gene.name)
                group_indels = True
            elif sam_path and sam_path.endswith(".tar.gz"):
                self.path = sam_path
                norm, muts = self._load_dump(sam_path, gene.name)
            else:
                assert sam_path
                self.path = sam_path
                is_sam = True
                if self.profile.name.startswith("pacbio"):
                    _, norm, muts = self._load_pacbio_sam(
                        sam_path, gene, reference, debug
                    )
                else:
                    _, norm, muts = self._load_sam(sam_path, gene, reference, debug)
                    group_indels = True
                if self.profile.cn_region:
                    self._dump_cn = self._load_cn_region(
                        self.profile.cn_region, sam_path, reference
                    )
            self._make_coverage(gene, norm, muts, group_indels=group_indels)
            if is_sam and debug:
                self._dump_alignments(f"{debug}.{gene.name}", norm, muts)

        if self.profile.cn_region:
            self.coverage._normalize_coverage()
        log.debug("[sam] avg_coverage= {:.1f}x", self.coverage.average_coverage())
        if self.profile.cn_region and self.coverage.diploid_avg_coverage() < 2:
            raise AldyException(
                "The average coverage of the sample is too low ({:.1f}).".format(
                    self.coverage.diploid_avg_coverage()
                )
            )

    def _load_sam(self, sam_path: str, gene: Gene, reference=None, debug=None):
        """
        Load the read, mutation and coverage data from a SAM/BAM/CRAM file.

        :param sam_path: Path to a SAM/BAM/CRAM file.
        :param gene: Gene instance.
        :param reference:
            Reference genome for reading CRAM files.
            Default is None.
        :param debug:
            If set, create a "`debug`.dump" file for debug purposes.
            Default is None.

        :raise: :obj:`aldy.common.AldyException` if a BAM/CRAM file lacks an index.
        """

        log.debug("[sam] path= {}", os.path.abspath(sam_path))
        if reference:
            log.debug("[sam] reference= {}", os.path.abspath(reference))

        norm: dict = defaultdict(list)
        muts: dict = defaultdict(list)

        with pysam.AlignmentFile(sam_path, reference_filename=reference) as sam:
            # Check do we have proper index to speed up the queries
            try:
                has_index = sam.check_index()
            except AttributeError:
                # SAM files do not have an index. BAMs might also lack it
                has_index = False
            except ValueError:
                raise AldyException(f"Cannot check index of {sam_path}")

            self._prefix = chr_prefix(gene.chr, [x["SN"] for x in sam.header["SQ"]])

            if has_index:
                iter = sam.fetch(
                    region=gene.get_wide_region().samtools(prefix=self._prefix)
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
                if not _in_region(gene.get_wide_region(), read, self._prefix):
                    continue

                if read.has_tag("BX"):
                    fragment = read.get_tag("BX")
                    if read.has_tag("XC"):
                        fragment = f"{fragment}:{read.get_tag('XC')}"
                    elif read.has_tag("MI"):
                        fragment = f"{fragment}:{read.get_tag('MI')}"
                else:
                    fragment = read.query_name
                r = self._parse_read(
                    gene,
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
        return has_index, norm, muts

    def _load_cn_region(self, cn_region, path, reference=None):
        """
        Load copy-number-neutral coverage from a SAM/BAM file.

        :param cn_region:
            Copy-number neutral region to be used for coverage rescaling.
            If None, profile loading and coverage rescaling will be skipped
            (and Aldy will require a ``--cn`` parameter to be user-provided).
        """
        self._dump_cn = collections.defaultdict(int)
        with pysam.AlignmentFile(path, reference_filename=reference) as sam:
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

    def _load_vcf(self, vcf_path: str, gene: Gene):
        """
        Load the read, mutation and coverage data from a VCF file.
        """

        log.debug("[vcf] path= {}", os.path.abspath(vcf_path))

        norm = {
            p: [(40, 40)] * 20
            for p in range(
                gene.get_wide_region().start - 500, gene.get_wide_region().end + 1
            )
        }
        muts: dict = defaultdict(list)

        def get_mut(pos, ref, alt):
            off = 0
            while off < len(ref) and off < len(alt) and ref[off] == alt[off]:
                off += 1
            if len(ref) - off == 1 and len(alt) - off == 1:
                if alt[off] == gene[off + pos]:
                    return off + pos, "_"
                return off + pos, f"{gene[off + pos]}>{alt[off]}"
            elif len(ref) > len(alt) and len(alt) - off == 0:
                return off + pos, f"del{gene[off + pos : pos + len(ref)]}"
            elif len(ref) < len(alt) and len(ref) - off == 0:
                return off + pos, f"ins{alt[off:]}"
            else:
                log.trace(f"[sam] ignoring {pos}: {ref}->{alt}")
                return pos, None

        with pysam.VariantFile(vcf_path) as vcf:
            self._prefix = chr_prefix(gene.chr, list(vcf.header.contigs))

            samples = list(vcf.header.samples)
            self.name = sample = samples[0]
            if len(samples) > 1:
                log.warn("WARNING: Multiple VCF samples found; using the first one.")
            log.info("Using VCF sample {}", sample)
            for read in vcf.fetch(
                region=gene.get_wide_region().samtools(prefix=self._prefix)
            ):
                g = sorted(y for y in read.samples[sample]["GT"] if y is not None)
                if len(g) != 2 or gene[read.pos - 1] == "N":
                    continue  # ignore polyploid and incomplete cases
                dump_arr = {}
                if len(read.ref) == 1 and read.ref != gene[read.pos - 1]:
                    hgvs = [(read.pos - 1, f"{gene[read.pos - 1]}>{read.ref}")]
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
                                muts[pos + p, f"{l[p]}>{r[p]}"] -= 10
                                if p:
                                    norm[pos + p] += [(40, 40)] * 10
                        muts[pos, op] += 10
        return norm, muts

    def _load_dump(self, dump_path: str, gene: str):
        log.debug("[dump] path= {}", os.path.abspath(dump_path))

        if dump_path.endswith(".tar.gz"):
            tar = tarfile.open(dump_path, "r:gz")

            f = [i for i in tar.getnames() if i.endswith(f".{gene}.dump")]
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
            self._dump_cn,
            norm,
            muts,
            phases,
            self._fusion_counter,
            self._insertion_reads,
        ) = pickle.load(fd)
        self.phases = {f"r{i}": v for i, v in enumerate(phases)}
        norm = {p: [q for q, n in c.items() for _ in range(n)] for p, c in norm.items()}
        muts = {p: [q for q, n in c.items() for _ in range(n)] for p, c in muts.items()}
        for (start, op), quals in muts.items():
            if op.startswith("ins"):
                self._insertion_counts[start, len(op) - 3] += len(quals)

        return norm, muts

    def _make_coverage(self, gene, norm, muts, group_indels=True):
        # Establish the coverage dictionary
        coverage: Dict[int, Dict[str, List]] = dict()
        for pos, cov in norm.items():
            if len(cov) == 0:
                continue
            coverage.setdefault(pos, {})["_"] = cov
        bounds = min(gene.chr_to_ref), max(gene.chr_to_ref)
        for m, cov in muts.items():
            pos, mut = m
            if pos not in coverage:
                coverage[pos] = {}
            if not bounds[0] <= pos <= bounds[1] and mut[:3] != "ins":
                mut = "_"  # ignore mutations outside of the region of interest
            coverage.setdefault(pos, {}).setdefault(mut, []).extend(cov)
        if group_indels:
            self._group_indels(gene, coverage)
        for pos, op in self._multi_sites.items():
            if pos in coverage and op in coverage[pos]:
                log.debug(
                    f"[sam] multi-SNP {pos}{op} with {len(coverage[pos][op])} reads"
                )
        for mut in self._insertion_sites:
            if mut.pos in coverage and mut.op in coverage[mut.pos]:
                total = sum(
                    len(v) for o, v in coverage[mut.pos].items() if o[:3] != "ins"
                )
                if total > 0:
                    corrected = self._correct_ins_coverage(mut, total)
                    diff = corrected - len(coverage[mut.pos][mut.op])
                    coverage[mut.pos][mut.op].extend(
                        [(mean(mq for mq, _ in coverage[mut.pos][mut.op]), 10)] * diff
                    )

        self.coverage = Coverage(
            self.gene,
            self.profile,
            self,
            {p: {m: v for m, v in coverage[p].items() if len(v) > 0} for p in coverage},
            self._dump_cn,
        )

    def _group_indels(self, gene, coverage):
        # Group ambiguous deletions
        indels = defaultdict(set)
        for m in gene.mutations:
            if "del" in m[1]:
                indels[m[0]].add(m[1])
        for pos in coverage:
            for mut in coverage[pos]:
                if mut[:3] != "del" or "N" in mut[3:]:
                    continue
                potential = [pos]
                sz = len(mut[3:])
                deleted = gene[pos : pos + sz]
                for p in range(pos - sz, -1, -sz):
                    if gene[p : p + sz] != deleted:
                        break
                    potential.append(p)
                for p in range(pos + sz, max(coverage), sz):
                    if gene[p : p + sz] != deleted:
                        break
                    potential.append(p)
                potential = [
                    p
                    for p in set(potential) & set(indels)
                    for m in indels[p]
                    if m == mut
                ]
                if len(potential) == 1 and potential[0] != pos:
                    new_pos = potential[0]
                    log.trace(
                        f"[sam] relocate {len(coverage[pos][mut])} "
                        + f"from {pos+1}{mut} to {new_pos+1}{mut}"
                    )
                    coverage[new_pos].setdefault(mut, []).extend(coverage[pos][mut])
                    coverage[pos][mut] = []
                    self.moved[pos, mut] = (new_pos, mut)
        # Group ambiguous insertions
        indels = defaultdict(set)
        for m in gene.mutations:
            if "ins" in m[1]:
                indels[m[0]].add(m[1])
        for pos in coverage:
            for mut in coverage[pos]:
                if mut[:3] != "ins":
                    continue
                inserted, sz = mut[3:], len(mut[3:])
                potential = []
                for p in range(pos, -1, -sz):
                    potential.append(p)
                    if gene[p : p + sz] != inserted:
                        break
                for p in range(pos, max(coverage), sz):
                    potential.append(p)
                    if gene[p : p + sz] != inserted:
                        break
                potential = [
                    p
                    for p in set(potential) & set(indels)
                    for m in indels[p]
                    if m == mut
                ]
                if len(potential) == 1 and potential[0] != pos:
                    new_pos = potential[0]
                    log.debug(
                        f"[sam] relocate {len(coverage[pos][mut])} "
                        + f"from {pos+1}{mut} to {new_pos+1}{mut}"
                    )
                    coverage[new_pos].setdefault(mut, []).extend(coverage[pos][mut])
                    coverage[pos][mut] = []

                    L = len(mut) - 3
                    self._insertion_counts[new_pos, L] += self._insertion_counts[pos, L]
                    self._insertion_counts[pos, L] = 0
                    self.moved[pos, mut] = (new_pos, mut)

    def _parse_read(
        self,
        gene: Gene,
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
        Parse a :obj:`pysam.AlignedSegment` read.

        :param gene: Gene instance.
        :param ref_start: Start position in the reference genome.
        :param cigar: List of CIGAR operations in tuple format (operation, size).
        :param seq: Read sequence.
        :param norm: Positions within the read that have not been mutated.
        :param muts: Positions within the read that have been mutated.

        .. note:: `norm` and `muts` are modified.
        """

        def bin_quality(q):
            """Inspired by https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf"""
            if q < 2:
                return q
            if 2 <= q < 10:
                return 6
            if 10 <= q < 20:
                return 15
            if 20 <= q < 29:
                return 25
            if 30 <= q < 39:
                return 35
            return 40

        phase = self.phases.setdefault(fragment, {})
        dump_arr = []
        start, s_start = ref_start, 0
        prev_q = 10
        for op, size in cigar:
            if op == 2:  # Deletion
                mut = (start, "del" + gene[start : start + size])
                muts[mut].append((bin_quality(mq), bin_quality(prev_q)))
                dump_arr.append(mut)
                if start in self.grid_columns:
                    phase[start] = mut[1]
                start += size
            elif op == 1:  # Insertion
                mut = (start, "ins" + seq[s_start : s_start + size])
                q = mean(qual[s_start : s_start + size]) if qual else prev_q
                muts[mut].append((bin_quality(mq), bin_quality(q)))
                prev_q = q
                # HACK: just store the length due to seq. errors
                self._insertion_counts[start, size] += 1
                dump_arr.append(mut)
                if start in self.grid_columns:
                    phase[start] = mut[1]
                s_start += size
            elif op == 4:  # Soft-clip
                s_start += size
            elif op in [0, 7, 8]:  # M, X and =
                for i in range(size):
                    q = qual[s_start + i] if qual else prev_q
                    if start + i in gene and gene[start + i] != seq[s_start + i]:
                        mut = (start + i, f"{gene[start + i]}>{seq[s_start + i]}")
                        dump_arr.append(mut)
                        muts[mut].append((bin_quality(mq), bin_quality(q)))
                        if start + i in self.grid_columns:
                            phase[start + i] = mut[1]
                    else:  # We ignore all mutations outside the RefSeq region
                        norm[start + i].append((bin_quality(mq), bin_quality(q)))
                        if start + i in self.grid_columns:
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
        read_pos = (ref_start, start, len(seq))
        for (pos, _), data in self._insertion_reads.items():
            if pos >= ref_start - 100 and pos < start + 100:
                data[read_pos] += 1
        return read_pos, dump_arr

    def _correct_ins_coverage(self, mut: Mutation, total):
        """
        Fix low coverage of large tandem insertions.

        Targets the cases where reference looks like
            ...X...
        and the donor genome looks like
            ...XX... (i.e. X is a tandem insertion).

        Any read that covers only one tandem (e.g. read X...) will get perfectly aligned
        (without trigerring an insertion tag) to the tail of the tandem insertion.
        This results in under-estimation of the insertion abundance will (as aligner
        could not assign `I` CIGAR to the correct insertion).

        This function attempts to correct this bias.
        """

        MARGIN_RATIO = 0.20
        INSERT_LEN_RESCALE = 10.0
        INSERT_LEN_AMPLIFY = 3

        orig_cov = self._insertion_counts[mut.pos, len(mut.op) - 3]
        new_total = 0.0
        for (orig_start, orig_end, read_len), cov in self._insertion_reads[mut].items():
            ins_len = len(mut.op) - 3
            # HACK: This is horrible way of detecting *40 (min 20% on the sides? max??)
            min_margin = MARGIN_RATIO * max(1, ins_len / INSERT_LEN_RESCALE) * read_len
            if orig_end >= mut.pos + max(
                INSERT_LEN_AMPLIFY * ins_len, min_margin
            ) and orig_start <= mut.pos - max(
                (INSERT_LEN_AMPLIFY - 1) * ins_len, min_margin
            ):
                new_total += cov
        new_total += orig_cov
        new_cov = int(total * (orig_cov / new_total))
        if new_cov > orig_cov:
            log.debug(f"[sam] rescale {mut} from {orig_cov} to {new_cov}")
            return new_cov
        return orig_cov

    def _load_pacbio_sam(self, sam_path: str, gene: Gene, reference=None, debug=None):
        """
        Load the read, mutation and coverage data from a PacBio SAM/BAM/CRAM file.
        This function remaps long PacBio reads for selected genes.
        Tested on PacBio HIFI data only.

        :param sam_path: Path to a PacBio SAM/BAM/CRAM file.
        :param gene: Gene instance.
        :param reference:
            Reference genome for reading CRAM files.
            Default is None.
        :param debug:
            If set, create a "`debug`.dump" file for debug purposes.
            Default is None.

        :raise: :obj:`aldy.common.AldyException` if a BAM/CRAM file lacks an index.
        """
        log.debug("[sam] pacbio_path= {}", os.path.abspath(sam_path))

        norm: dict = defaultdict(list)
        muts: dict = defaultdict(list)

        self._index = None
        ref = script_path(f"aldy.resources.genes/{gene.name.lower()}.fa.gz")
        if os.path.exists(ref):
            import mappy

            self._index = mappy.Aligner(ref, preset="map-hifi")
            seq_name = self._index.seq_names[0]
            self._seq_range = list(map(int, seq_name.split(":")[1].split("-")))
            self._index_gene = []
            for _, g in enumerate(gene.regions):
                st = min(r.start for r in g.values())
                ed = max(r.end for r in g.values())
                idx = mappy.Aligner(
                    seq=self._index.seq(seq_name)[
                        st - self._seq_range[0] : ed - self._seq_range[0]
                    ],
                    preset="map-hifi",
                )
                self._index_gene.append((idx, st))

            self._fusion_signatures = {}
            self._fusion_counter = {}
            for n, cn in gene.cn_configs.items():
                if cn.kind == CNConfigType.LEFT_FUSION:
                    sig, pos = 0, next(i for i, j in cn.cn[0].items() if j == 1)
                elif cn.kind == CNConfigType.RIGHT_FUSION:
                    sig, pos = 1, next(i for i, j in cn.cn[0].items() if j == 0)
                else:
                    continue
                sig = (
                    (1 - sig, list(cn.cn[0])[list(cn.cn[0]).index(pos) - 1]),
                    (sig, pos),
                )
                self._fusion_counter[n] = [0, 0]
                if sig[0][1] in gene.unique_regions or sig[1][1] in gene.unique_regions:
                    self._fusion_signatures[sig] = n

            log.debug("[sam] PacBio remapping enabled")

        with pysam.AlignmentFile(sam_path, reference_filename=reference) as sam, Timing(
            "[sam] Remap"
        ):
            # Assumes SAM index exists
            self._prefix = chr_prefix(gene.chr, [x["SN"] for x in sam.header["SQ"]])

            rng = [*gene.get_wide_region()]
            if gene.name == "CYP2D6":
                rng[2] = 42_155_000  # include CYP2D8 [hg38]
            rng = GRange(*rng)
            counter = 0
            iter = sam.fetch(region=rng.samtools(prefix=self._prefix))
            for read in iter:
                if not read.cigartuples or "H" in read.cigarstring:
                    continue
                if read.reference_start >= rng.end or read.reference_end < rng.start:
                    continue

                seq = read.query_sequence
                qual = read.query_qualities
                if self._index:
                    pieces = []
                    s_start = 0
                    while s_start < len(seq):
                        s = seq[s_start:]
                        q = qual[s_start:]
                        h = self._map(self._index, s)
                        if h and h.r_st <= 134_991 < h.r_en and 134_991 - h.r_st > 100:
                            s = seq[s_start : s_start + (134_991 - h.r_st)]
                            q = qual[s_start : s_start + (134_991 - h.r_st)]
                            h = self._map(self._index, s)
                        elif (
                            h and h.r_st <= 148_730 < h.r_en and 148_730 - h.r_st > 100
                        ):
                            s = seq[s_start : s_start + (148_730 - h.r_st)]
                            q = qual[s_start : s_start + (148_730 - h.r_st)]
                            h = self._map(self._index, s)
                        if not h:
                            s_start += 100
                        else:
                            pcs, ed = self._split_read(gene, s, q, h)
                            pieces += pcs
                            s_start += ed
                    for (ref_start, seq, qual, cigar) in pieces:
                        r = self._parse_read(
                            gene,
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
                        gene,
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
        return True, norm, muts

    def _map(self, idx, seq):
        for hit in idx.map(seq):
            if hit.is_primary:
                return hit
        return None

    def _split_read(self, gene, seq, qual, hit):
        def _rev_cigar(c):
            """Reverse CIGAR tuples for pysam/mappy compatibility."""
            return [(op, size) for (size, op) in c]

        regs = self._get_gene_regions(
            gene, hit.r_st + self._seq_range[0] - 1, _rev_cigar(hit.cigar)
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

    def _get_gene_regions(self, gene, r_start, cigar):
        regs = []
        start, s_start = r_start, 0
        for op, size in cigar:
            if op == 1 or op == 4:
                s_start += size
            elif op == 2 or op == 3:
                start += size
            elif op == 0 or op == 7 or op == 8:
                for i in range(size):
                    rg = gene.region_at(start + i)
                    if rg:
                        if not regs or regs[-1][0] != rg:
                            regs.append([rg, s_start + i, 0])
                        regs[-1][2] += 1
                start += size
                s_start += size
        return regs

    # ----------------------------------------------------------------------------------
    # Coverage-rescaling functions
    # ----------------------------------------------------------------------------------

    def _dump_alignments(self, debug: str, norm, muts):
        with open(f"{debug}.genome", "w") as fd:
            print(self.gene.genome, file=fd)
        with gzip.open(f"{debug}.dump", "wb") as fd:
            pickle.dump(
                (
                    self.name,
                    self._dump_cn,
                    {p: Counter(q) for p, q in norm.items()},
                    {p: Counter(q) for p, q in muts.items()},
                    [v for v in self.phases.values() if len(v) > 1],
                    self._fusion_counter,
                    self._insertion_reads,
                ),
                fd,
            )


def detect_genome(sam_path: str) -> Tuple[str, Optional[str]]:
    try:
        pysam.set_verbosity(0)
        with pysam.AlignmentFile(sam_path) as sam:
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
            with pysam.VariantFile(sam_path):
                return "vcf", None
        except (ValueError, OSError):
            if sam_path.endswith(".tar.gz"):
                tar = tarfile.open(sam_path, "r:gz")
                f = [i for i in tar.getnames() if i.endswith(f".genome")]
                if not f:
                    raise AldyException("Invalid dump file")
                data = tar.extractfile(f[0])
                if data:
                    genome = data.read().decode("utf-8").strip()
                    return "dump", genome
    return "", None


def _chr_prefix(ch: str, chrs: List[str]) -> str:  # assumes ch is not prefixed
    """
    Check if `ch` (*without any chr prefix*) should be prefixed with "chr" or not.

    Returns:
        str: Prefix to be prepended to chromosome (empty if needed).

    Params:
        ch (str): chromosome name
        chrs (list[str]): list of chromosome names in the alignment file
    """
    if ch not in chrs and "chr" + ch in chrs:
        return "chr"
    return ""


def _in_region(region: GRange, read: pysam.AlignedSegment, prefix: str) -> bool:
    """
    Returns ``True`` if a read is located within a given gene region.

    Notes:
        The region is padded with 500bp on the left side.
    """

    if read.reference_id == -1 or read.reference_name != prefix + region.chr:
        return False
    if read.reference_end is None:
        return False

    a = (read.reference_start, read.reference_end)
    b = (region.start, region.end)
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]
