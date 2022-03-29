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
import struct
import tarfile

from dataclasses import dataclass
from collections import defaultdict
from natsort import natsorted
from .common import log, GRange, AldyException, script_path, Timing
from .gene import Gene, Mutation, CNConfigType
from .coverage import Coverage


class Profile:
    def __init__(self, name, cn_region, data, neutral_value=0, cn_solution=None):
        self.name = name
        self.cn_region = cn_region
        self.data = data
        self.cn_solution = cn_solution
        self.neutral_value = neutral_value


@dataclass
class Sample:
    """
    Interface for reading SAM/BAM/CRAM files that parses and stores read alignments.
    """

    coverage: Coverage
    sample_name: str
    min_cov: float
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
        sam_path: Optional[str] = None,
        threshold: float = 0.5,
        profile: Optional[Profile] = None,
        reference: Optional[str] = None,
        debug: Optional[str] = None,
        vcf_path: Optional[str] = None,
        min_cov: float = 1.0,
    ):
        """
        Initialize a :obj:`Sample` object.

        :param sam_path: Path to a SAM/BAM/CRAM file.
        :param gene: Gene instance.
        :param threshold:
            Threshold for filtering out low quality mutations. Ranges from 0 to 1.
            Check :obj:`aldy.coverage.Coverage` for more information.
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

        # Get the list of indel sites that should be corrected
        # TODO: currently uses only functional indels;
        #       other indels should be corrected as well
        self.min_cov = min_cov
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
        self.grid = []

        self.profile = profile
        self._fusion_counter = None

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
            elif sam_path and sam_path.endswith(".tar.gz"):
                self.path = sam_path
                norm, muts = self._load_dump(sam_path, gene.name)
            else:
                assert sam_path
                self.path = sam_path
                is_sam = True
                if self.profile and self.profile.name.startswith("pacbio"):
                    _, norm, muts = self._load_pacbio_sam(
                        sam_path, gene, reference, debug
                    )
                else:
                    _, norm, muts = self._load_sam(sam_path, gene, reference, debug)
                    group_indels = True
                if self.profile:
                    self._dump_cn = self._load_cn_region(self.profile.cn_region)
            self._make_coverage(gene, norm, muts, threshold, group_indels=group_indels)
            if is_sam and debug:
                self._dump_alignments(f"{debug}.{gene.name}")

        if self.profile and self.profile.cn_region:
            self.coverage._normalize_coverage()
        log.debug("[sam] avg_coverage= {:.1f}x", self.coverage.average_coverage())
        if self.profile and self.coverage.diploid_avg_coverage() < 2:
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
        :param threshold:
            Threshold for filtering out low quality mutations. Ranges from 0 to 1.
            Check :obj:`aldy.coverage.Coverage` for more information.
        :param reference:
            Reference genome for reading CRAM files.
            Default is None.
        :param debug:
            If set, create a "`debug`.dump" file for debug purposes.
            Default is None.

        :raise: :obj:`aldy.common.AldyException` if a BAM/CRAM file lacks an index.
        """

        log.debug("[sam] path= {}", os.path.abspath(sam_path))
        self.sample_name = os.path.basename(sam_path).split(".")[0]

        norm: dict = defaultdict(int)
        muts: dict = defaultdict(int)

        with pysam.AlignmentFile(sam_path, reference_filename=reference) as sam:
            # Check do we have proper index to speed up the queries
            try:
                has_index = sam.check_index()
            except AttributeError:
                # SAM files do not have an index. BAMs might also lack it
                has_index = False
            except ValueError:
                raise AldyException(f"Cannot check index of {sam_path}")

            self._prefix = _chr_prefix(gene.chr, [x["SN"] for x in sam.header["SQ"]])

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
                # ensure that it is a proper gene read
                if not _in_region(gene.get_wide_region(), read, self._prefix):
                    continue

                r = self._parse_read(
                    gene,
                    read.reference_start,
                    read.cigartuples,
                    read.query_sequence,
                    norm,
                    muts,
                )
                if r and debug:
                    self._dump_reads.append(r)
        return has_index, norm, muts

    def _load_cn_region(self, cn_region):
        """
        Load copy-number-neutral coverage from a SAM/BAM file.

        :param cn_region:
            Copy-number neutral region to be used for coverage rescaling.
            If None, profile loading and coverage rescaling will be skipped
            (and Aldy will require a ``--cn`` parameter to be user-provided).
        """
        self._dump_cn = collections.defaultdict(int)
        with pysam.AlignmentFile(self.path) as sam:
            # Check do we have proper index to speed up the queries
            try:
                has_index = sam.check_index()
            except AttributeError:
                # SAM files do not have an index. BAMs might also lack it
                has_index = False
            except ValueError:
                raise AldyException(f"Cannot check index of {self.path}")
            self._prefix = _chr_prefix(
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
            p: 20
            for p in range(
                gene.get_wide_region().start - 500, gene.get_wide_region().end + 1
            )
        }
        muts: dict = defaultdict(int)

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
            self._prefix = _chr_prefix(gene.chr, list(vcf.header.contigs))

            samples = list(vcf.header.samples)
            self.sample_name = sample = samples[0]
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
                    muts[pos, op] += 10
                    norm[pos] -= 10
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
                                    norm[pos + p] += 10
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
            fd = gzip.open(tar.extractfile(f[0]))
        else:
            fd = gzip.open(dump_path, "rb")

        self.sample_name = "DUMP"
        self.sample_name = os.path.basename(dump_path)
        norm: dict = defaultdict(int)
        muts: dict = defaultdict(int)

        log.warn("Loading debug dump from {}", dump_path)
        l, h, i = (
            struct.calcsize("<l"),
            struct.calcsize("<h"),
            struct.calcsize("<i"),
        )
        cn_len = struct.unpack("<l", fd.read(l))[0]
        for _ in range(cn_len):
            i, v = struct.unpack("<ll", fd.read(l + l))
            self._dump_cn[i] = v

        ld = struct.unpack("<l", fd.read(l))[0]
        for _ in range(ld):
            ref_start, ref_end, read_len, num_mutations = struct.unpack(
                "<llhh", fd.read(l + l + h + h)
            )
            ref_end += ref_start
            for j in range(ref_start, ref_end):
                norm[j] += 1
            mut_set = set()
            for _ in range(num_mutations):
                mut_start, op_len = struct.unpack("<hh", fd.read(h + h))
                mut_start += ref_start
                op = fd.read(op_len).decode("ascii")
                muts[mut_start, op] += 1
                mut_set.add((mut_start, op))
                if op[:3] == "del":
                    for j in range(0, len(op) - 3):
                        norm[mut_start + j] -= 1
                elif op[:3] == "ins":
                    self._insertion_counts[mut_start, len(op) - 3] += 1
                else:
                    norm[mut_start] -= 1
            mut_set_pos = {p for p, _ in mut_set}
            for pos, op in self._multi_sites.items():
                if pos not in mut_set_pos:
                    continue
                ll, r = op.split(">")
                if all(
                    (pos + p, f"{ll[p]}>{r[p]}") in mut_set
                    for p in range(len(ll))
                    if ll[p] != "."
                ):
                    for p in range(len(ll)):
                        if ll[p] != ".":
                            muts[pos + p, f"{ll[p]}>{r[p]}"] -= 1
                            if p:
                                norm[pos + p] += 1
                    muts[pos, op] += 1
            for data in self._insertion_reads.values():
                data[ref_start, ref_end, read_len] += 1

        frags = []
        ld = struct.unpack("<l", fd.read(l))[0]
        for _ in range(ld):
            ln = struct.unpack("<l", fd.read(l))[0]
            frags.append([])
            for _ in range(ln):
                mut_start, op_len = struct.unpack("<ll", fd.read(l + l))
                op = fd.read(op_len).decode("ascii")
                frags[-1].append((mut_start, op))
        fd.close()

        return norm, muts, frags

    def _make_coverage(self, gene, norm, muts, threshold, group_indels=True):
        # Establish the coverage dictionary
        coverage: Dict[int, Dict[str, int]] = dict()
        for pos, cov in norm.items():
            if cov == 0:
                continue
            if pos not in coverage:
                coverage[pos] = {}
            coverage[pos]["_"] = cov
        bounds = min(gene.chr_to_ref), max(gene.chr_to_ref)
        for m, cov in muts.items():
            pos, mut = m
            if pos not in coverage:
                coverage[pos] = {}
            if not bounds[0] <= pos <= bounds[1] and mut[:3] != "ins":
                mut = "_"  # ignore mutations outside of the region of interest
            if mut not in coverage[pos]:
                coverage[pos][mut] = 0
            coverage[pos][mut] += cov
        if group_indels:
            self._group_indels(gene, coverage)
        for pos, op in self._multi_sites.items():
            if pos in coverage and op in coverage[pos]:
                log.debug(f"[sam] multi-SNP {pos}{op} with {coverage[pos][op]} reads")
        for mut in self._insertion_sites:
            if mut.pos in coverage and mut.op in coverage[mut.pos]:
                total = sum(v for o, v in coverage[mut.pos].items() if o[:3] != "ins")
                if total > 0:
                    coverage[mut.pos][mut.op] = self._correct_ins_coverage(mut, total)

        self.coverage = Coverage(
            self.gene,
            self.profile,
            self,
            {p: {m: v for m, v in coverage[p].items() if v > 0} for p in coverage},
            self._dump_cn,
            threshold,
            self.min_cov,
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
                        f"[sam] relocate {coverage[pos][mut]} "
                        + f"from {pos+1}{mut} to {new_pos+1}{mut}"
                    )
                    if mut not in coverage[new_pos]:
                        coverage[new_pos][mut] = 0
                    coverage[new_pos][mut] += coverage[pos][mut]
                    coverage[pos][mut] = 0
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
                    log.trace(
                        f"[sam] relocate {coverage[pos][mut]} "
                        + f"from {pos+1}{mut} to {new_pos+1}{mut}"
                    )
                    if mut not in coverage[new_pos]:
                        coverage[new_pos][mut] = 0
                    coverage[new_pos][mut] += coverage[pos][mut]

                    L = len(mut) - 3
                    self._insertion_counts[new_pos, L] += self._insertion_counts[pos, L]
                    self._insertion_counts[pos, L] = 0

                    coverage[pos][mut] = 0

    def _parse_read(self, gene: Gene, ref_start, cigar, seq, norm, muts):
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

        grid = {}
        regions = []
        dump_arr = []
        start, s_start = ref_start, 0
        for op, size in cigar:
            if op == 2:  # Deletion
                mut = (start, "del" + gene[start : start + size])
                muts[mut] += 1
                dump_arr.append(mut)
                if gene.region_at(start) and (
                    not regions or gene.region_at(start) != regions[-1]
                ):
                    regions.append(gene.region_at(start))
                if start in self.grid_columns:
                    grid[start] = mut
                start += size
            elif op == 1:  # Insertion
                mut = (start, "ins" + seq[s_start : s_start + size])
                muts[mut] += 1
                # HACK: just store the length due to seq. errors
                self._insertion_counts[start, size] += 1
                dump_arr.append(mut)
                s_start += size
            elif op == 4:  # Soft-clip
                s_start += size
            elif op in [0, 7, 8]:  # M, X and =
                for i in range(size):
                    if gene.region_at(start + i) and (
                        not regions or gene.region_at(start + i) != regions[-1]
                    ):
                        regions.append(gene.region_at(start + i))
                    if start + i in gene and gene[start + i] != seq[s_start + i]:
                        mut = (start + i, f"{gene[start + i]}>{seq[s_start + i]}")
                        dump_arr.append(mut)
                        muts[mut] += 1
                        if start + i in self.grid_columns:
                            grid[start + i] = mut
                    else:  # We ignore all mutations outside the RefSeq region
                        norm[start + i] += 1
                        if start + i in self.grid_columns:
                            grid[start + i] = (start + i, "_")
                start += size
                s_start += size

        self.grid.append(grid)
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
                for p in range(len(l)):
                    if l[p] != ".":
                        muts[pos + p, f"{l[p]}>{r[p]}"] -= 1
                        if p:
                            norm[pos + p] += 1
                muts[pos, op] += 1

        read_pos = (ref_start, start, len(seq))
        for data in self._insertion_reads.values():
            data[read_pos] += 1  # type: ignore
        return read_pos, dump_arr

    def _correct_ins_coverage(self, mut: Mutation, total) -> int:
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
        :param threshold:
            Threshold for filtering out low quality mutations. Ranges from 0 to 1.
            Check :obj:`aldy.coverage.Coverage` for more information.
        :param reference:
            Reference genome for reading CRAM files.
            Default is None.
        :param debug:
            If set, create a "`debug`.dump" file for debug purposes.
            Default is None.

        :raise: :obj:`aldy.common.AldyException` if a BAM/CRAM file lacks an index.
        """
        log.debug("[sam] pacbio_path= {}", os.path.abspath(sam_path))
        self.sample_name = os.path.basename(sam_path).split(".")[0]

        norm: dict = defaultdict(int)
        muts: dict = defaultdict(int)

        self._index = None
        if gene.name == "CYP2D6":
            import mappy

            ref = script_path(f"aldy.resources.genes/{gene.name}.fa.gz")
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
            self._prefix = _chr_prefix(gene.chr, [x["SN"] for x in sam.header["SQ"]])

            rng = [*gene.get_wide_region()]
            if gene.name == "CYP2D6":
                rng[2] = 42_155_000  # include CYP2D8
            rng = GRange(*rng)
            iter = sam.fetch(region=rng.samtools(prefix=self._prefix))
            for read in iter:
                if not read.cigartuples or "H" in read.cigarstring:
                    continue
                if read.reference_start >= rng.end or read.reference_end < rng.start:
                    continue

                seq = read.query_sequence
                if self._index:
                    pieces = []
                    s_start = 0
                    while s_start < len(seq):
                        s = seq[s_start:]
                        h = self._map(self._index, s)
                        if h and h.r_st <= 134_991 < h.r_en and 134_991 - h.r_st > 100:
                            s = seq[s_start : s_start + (134_991 - h.r_st)]
                            h = self._map(self._index, s)
                        elif (
                            h and h.r_st <= 148_730 < h.r_en and 148_730 - h.r_st > 100
                        ):
                            s = seq[s_start : s_start + (148_730 - h.r_st)]
                            h = self._map(self._index, s)
                        if not h:
                            s_start += 100
                        else:
                            pcs, ed = self._split_read(gene, s, h)
                            pieces += pcs
                            s_start += ed
                    for (ref_start, seq, cigar) in pieces:
                        r = self._parse_read(gene, ref_start, cigar, seq, norm, muts)
                        if r and debug:
                            self._dump_reads.append(r)
                else:
                    r = self._parse_read(
                        gene, read.reference_start, read.cigartuple, seq, norm, muts
                    )
                    if r and debug:
                        self._dump_reads.append(r)
        return True, norm, muts

    def _map(self, idx, seq):
        for hit in idx.map(seq):
            if hit.is_primary:
                return hit
        return None

    def _split_read(self, gene, seq, hit):
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
                            _rev_cigar(lh.cigar),
                        ),
                        (
                            rh.r_st + self._index_gene[rg][1] - 1,
                            seq[brk + rh.q_st : brk + rh.q_en],
                            _rev_cigar(rh.cigar),
                        ),
                    ],
                    max(lh.q_en, brk + rh.q_en),
                )
        return [
            (
                hit.r_st + self._seq_range[0] - 1,
                seq[hit.q_st : hit.q_en],
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
                    # if not rg and 42_149_886 <= start + i <= 42_155_001:
                    #     rg = (2, '2D8')
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

    def load_profile(gene, profile, cn_region=None):
        prof = None
        is_yml = True
        if os.path.exists(profile) and os.path.isfile(profile):
            ext = os.path.splitext(profile)
            if ext[-1] in [".bam", ".sam", ".cram"]:
                assert gene and cn_region
                regions = {
                    (gene.name, r, gi): rng
                    for gi, gr in enumerate(gene.regions)
                    for r, rng in gr.items()
                }
                prof = load_sam_profile(
                    profile_path,
                    regions=regions,
                    genome=gene.genome,
                    cn_region=cn_region,
                )
                is_yml = False
            else:
                with open(profile_path) as f:
                    prof = yaml.safe_load(f)
        else:
            profile_path = script_path(
                "aldy.resources.profiles/{}.yml".format(profile.lower())
            )
            with open(profile_path) as f:
                prof = yaml.safe_load(f)
        if not prof:
            raise AldyException(f"Could not load profile from {profile}")
        if is_yml and profile != "illumina" and cn_region:
            raise AldyException("-n is only valid with illumina or BAM profile")
        if is_yml and profile == "illumina" and cn_region:
            prof["neutral"]["value"] = cn_region.end - cn_region.start
        if "neutral" not in prof or "value" not in prof["neutral"]:
            raise AldyException("Profile missing neutral region")
        if gene.name not in prof:
            raise AldyException(f"Profile missing {gene.name}")
        return Profile(
            profile,
            GRange(*prof["neutral"][gene.genome]),
            prof,
            prof["neutral"]["value"],
        )

    def _dump_alignments(self, debug: str):
        with gzip.open(f"{debug}.dump", "wb") as fd:
            fd.write(struct.pack("<l", len(self._dump_cn)))
            for i, v in self._dump_cn.items():
                fd.write(struct.pack("<ll", i, v))

            fd.write(struct.pack("<l", len(self._dump_reads)))
            for (s, e, l), m in self._dump_reads:
                fd.write(struct.pack("<llhh", s, e - s, l, len(m)))
                for p, md in m:
                    fd.write(struct.pack("<hh", p - s, len(md)))
                    fd.write(md.encode("ascii"))

            # fd.write(struct.pack("<l", len(self._dump_phase)))
            # for r in self._dump_phase:
            #     fd.write(struct.pack("<l", len(r)))
            #     for p, md in r:
            #         fd.write(struct.pack("<ll", p, len(md)))
            #         fd.write(md.encode("ascii"))


def load_sam_profile(
    sam_path: str,
    factor: float = 2.0,
    regions: Dict[Tuple[str, str, int], GRange] = dict(),
    cn_region: Optional[GRange] = None,
    genome: Optional[str] = "hg19",
) -> Dict[str, Dict[str, List[float]]]:
    """
    Load the profile information from a SAM/BAM file.

    Returns:
        list[str, str, int, float]: list of tuples
        ``(gene_name, chromosome, loci, coverage)``.

    Params:
        factor (float):
            Scaling factor. Default is 2.0 (for two copies).
        regions (list[:obj:`GRange`], optional):
            List of regions to be extracted.

    Notes:
        Profiles that were used in Aldy paper:

            1. PGRNseq-v1/v3: NA17642 was used for all genes
                (n.b. PGXT147 with rescale 2.52444127771 was used for CYP2B6 beta).
            2. PGRNseq-v2: NA19789.bam was used for all genes.
            3. Illumina: by definition contains all ones (uniform coverage profile).
    """

    if not genome:
        genome = "hg19"
    if len(regions) == 0:
        import pkg_resources

        gene_regions = {}
        for g in sorted(pkg_resources.resource_listdir("aldy.resources", "genes")):
            if g[-4:] != ".yml":
                continue
            gg = Gene(script_path(f"aldy.resources.genes/{g}"), genome=genome)
            for gi, gr in enumerate(gg.regions):
                for r, rng in gr.items():
                    gene_regions[gg.name, r, gi] = rng
    else:
        gene_regions = regions

    default_cn_neutral_region = {
        "hg19": GRange("22", 42547463, 42548249),
        "hg38": GRange("22", 42151472, 42152258),
    }
    gene_regions["neutral", "value", 0] = (
        cn_region if cn_region else default_cn_neutral_region[genome]
    )

    chr_regions: Dict[str, Tuple[int, int]] = {}
    for c, s, e in gene_regions.values():
        if c not in chr_regions:
            chr_regions[c] = (s, e)
        else:
            chr_regions[c] = (min(s, chr_regions[c][0]), max(e, chr_regions[c][1]))

    cov: dict = defaultdict(lambda: defaultdict(int))
    for c, (s, e) in natsorted(chr_regions.items()):
        if sam_path == "<illumina>":
            continue
        with pysam.AlignmentFile(sam_path) as sam:
            region = GRange(c, s, e).samtools(
                pad_left=1000,
                pad_right=1000,
                prefix=_chr_prefix(c, [x["SN"] for x in sam.header["SQ"]]),
            )
            log.info("Scanning {}...", region)
            try:
                for read in sam.fetch(region=region):
                    start, s_start = read.reference_start, 0
                    if not read.cigartuples:
                        continue

                    for op, size in read.cigartuples:
                        if op == 2:
                            for i in range(size):
                                cov[c][start + i] += 1
                            start += size
                        elif op == 1:
                            s_start += size
                        elif op == 4:
                            s_start += size
                        elif op in [0, 7, 8]:
                            for i in range(size):
                                cov[c][start + i] += 1
                            start += size
                            s_start += size
            except ValueError:
                log.warn("Cannot fetch {}", region)

    d: Any = {}
    for (g, r, ri), (c, s, e) in gene_regions.items():
        if g not in d:
            d[g] = {}
        if r not in d[g]:
            d[g][r] = [0]
        if ri >= len(d[g][r]):
            d[g][r].append(0)
        if sam_path == "<illumina>":
            d[g][r][ri] = e - s
        elif g == "neutral":
            d[g][r][ri] = sum(cov[c][i] for i in range(s, e))
        else:
            d[g][r][ri] = sum(cov[c][i] / factor for i in range(s, e))
    d["neutral"]["value"] = d["neutral"]["value"][0]
    d["neutral"][genome] = [*gene_regions["neutral", "value", 0]]
    return d


def detect_genome(sam_path: str) -> Tuple[str, Optional[str]]:
    try:
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
            _prefix = _chr_prefix("1", [x["SN"] for x in sam.header["SQ"]])
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
            pass
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
