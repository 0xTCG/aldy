# 786
# Aldy source: sam.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple, Dict, List, Optional, Any, Set

import pysam
import os
import os.path
import re
import yaml
import gzip
import struct
import tarfile

from dataclasses import dataclass
from collections import defaultdict
from natsort import natsorted
from .common import log, GRange, AldyException, script_path, Timing
from .coverage import Coverage
from .gene import Gene, Mutation


CIGAR2CODE: Dict[int, int] = {ord(y): x for x, y in enumerate("MIDNSHP=XB")}
"""CIGAR characters for an ASCII ordinal of CIGAR character"""

CIGAR_REGEX = re.compile(r"(\d+)([MIDNSHP=XB])")
"""Regex that matches valid CIGAR strings"""


DEFAULT_CN_NEUTRAL_REGION = {
    "hg19": GRange("22", 42547463, 42548249),
    "hg38": GRange("22", 42151472, 42152258),
}
"""Default copy-number neutral region (exon 4-6 of the CYP2D8 gene)"""


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

    # Phasing data

    def __init__(
        self,
        gene: Gene,
        sam_path: Optional[str] = None,
        threshold: float = 0.5,
        profile: Optional[str] = None,
        reference: Optional[str] = None,
        cn_region: Optional[GRange] = None,
        debug: Optional[str] = None,
        vcf_path: Optional[str] = None,
        min_cov: float = 1.0,
        phase: bool = False,
    ):
        """
        Initialize a :obj:`Sample` object.

        :param sam_path: Path to a SAM/BAM/CRAM file.
        :param gene: Gene instance.
        :param threshold:
            Threshold for filtering out low quality mutations. Ranges from 0 to 1.
            Check :obj:`aldy.coverage.Coverage` for more information.
        :param profile:
            Profile specification (e.g. 'prgnseq-v1') or a profile SAM/BAM file.
        :param reference:
            Reference genome for reading CRAM files.
            Default is None.
        :param cn_region:
            Copy-number neutral region to be used for coverage rescaling.
            If None, profile loading and coverage rescaling will be skipped
            (and Aldy will require a ``--cn`` parameter to be user-provided).
            Default is ``DEFAULT_CN_NEUTRAL_REGION``.
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
        self._dump_phase = []
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

        with Timing("Read SAM"):
            has_index = False
            frags = []
            is_sam = False
            assert vcf_path or (sam_path and profile)
            if vcf_path:
                try:
                    norm, muts = self._load_vcf(vcf_path, gene)
                except ValueError:
                    raise AldyException(f"VCF {vcf_path} is not indexed")
            elif sam_path and sam_path.endswith(".dump"):
                norm, muts, frags = self._load_dump(sam_path, gene.name)
            elif sam_path and sam_path.endswith(".tar.gz"):
                norm, muts, frags = self._load_dump(sam_path, gene.name)
            else:
                assert sam_path
                is_sam = True
                has_index, norm, muts = self._load_sam(
                    sam_path, gene, reference, cn_region, debug
                )
            self._make_coverage(gene, norm, muts, threshold)
            if is_sam and (debug or phase):
                if has_index:
                    with Timing("Read phase"):
                        frags = self._get_phases(sam_path, gene, reference)
                if debug:
                    self._dump_alignments(f"{debug}.{gene.name}")
            if phase and frags:
                self.coverage.fragments = frags

        if cn_region:
            assert profile
            self.detect_cn(gene, profile, cn_region)
        log.debug("[sam] avg_coverage= {:.1f}x", self.coverage.average_coverage())
        if cn_region and self.coverage.diploid_avg_coverage() < 2:
            raise AldyException(
                "The average coverage of the sample is too low ({:.1f}).".format(
                    self.coverage.diploid_avg_coverage()
                )
            )

    def _load_sam(
        self,
        sam_path: str,
        gene: Gene,
        reference: Optional[str] = None,
        cn_region: Optional[GRange] = None,
        debug: Optional[str] = None,
    ):
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
        :param cn_region:
            Copy-number neutral region to be used for coverage rescaling.
            If None, profile loading and coverage rescaling will be skipped
            (and Aldy will require a ``--cn`` parameter to be user-provided).
            Default is ``DEFAULT_CN_NEUTRAL_REGION`` (CYP2D8).
        :param debug:
            If set, create a "`debug`.dump" file for debug purposes.
            Default is None.

        :raise: :obj:`aldy.common.AldyException` if a BAM/CRAM file lacks an index.
        """

        log.debug("[sam] path= {}", os.path.abspath(sam_path))
        self.sample_name = os.path.basename(sam_path).split(".")[0]

        norm: dict = defaultdict(int)
        muts: dict = defaultdict(int)

        # Try to read CN-neutral region if there is an index (fast)
        # If not, the main loop will catch it later on
        def read_cn_read(read):
            """
            Check is pysam.AlignmentSegment valid CN-neutral read,
            and if so, parse it
            """
            start = read.reference_start
            if read.cigartuples is None or (
                (read.flag & 1) and (read.flag & 0x40) == 0
            ):
                return
            for op, size in read.cigartuples:
                if op in [0, 7, 8, 2]:
                    for i in range(size):
                        self._dump_cn[start + i] += 1
                    start += size

        has_index = True
        with pysam.AlignmentFile(sam_path, reference_filename=reference) as sam:
            # Check do we have proper index to speed up the queries
            try:
                sam.check_index()
            except AttributeError:
                has_index = (
                    False  # SAM files do not have an index. BAMs might also lack it
                )
            except ValueError:
                raise AldyException(f"File {sam_path} has no index")

            self._prefix = _chr_prefix(gene.chr, [x["SN"] for x in sam.header["SQ"]])

            # Set it to _fetched_ if a CN-neutral region is user-provided.
            # Then read the CN-neutral region.
            is_cn_region_fetched = cn_region is None
            if cn_region and sam.has_index():
                log.debug("[sam] has_index= True")
                for read in sam.fetch(region=cn_region.samtools(prefix=self._prefix)):
                    read_cn_read(read)
                is_cn_region_fetched = True

            # Fetch the reads
            for read in sam.fetch(
                region=gene.get_wide_region().samtools(prefix=self._prefix)
            ):
                # If we haven't obtained CN-neutral region so far, do it now
                if (
                    not is_cn_region_fetched
                    and cn_region
                    and _in_region(cn_region, read, self._prefix)
                ):  # type: ignore
                    read_cn_read(read)

                r = self._parse_read(read, gene, norm, muts)
                if r and debug:
                    self._dump_reads.append(r)
        return has_index, norm, muts

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

    def _make_coverage(self, gene, norm, muts, threshold):
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
            {p: {m: v for m, v in coverage[p].items() if v > 0} for p in coverage},
            threshold,
            self._dump_cn,
            self.sample_name,
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
                    log.debug(
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
                    log.debug(
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

    def _parse_read(self, read, gene: Gene, norm, muts):
        """
        Parse a :obj:`pysam.AlignedSegment` read.

        :param read: pysam read.
        :param gene: Gene instance.
        :param norm: Positions within the read that have not been mutated.
        :param muts: Positions within the read that have been mutated.

        .. note:: `norm` and `muts` are modified.
        """

        if not read.cigartuples:  # only valid alignments
            return None
        if read.flag & 0x800:  # avoid supplementary alignments
            return None
        if "H" in read.cigarstring:  # avoid hard-clipped reads
            return None
        if not _in_region(
            gene.get_wide_region(), read, self._prefix
        ):  # ensure that it is a proper gene read
            return None

        dump_arr = []
        start, s_start = read.reference_start, 0
        for op, size in read.cigartuples:
            if op == 2:  # Deletion
                mut = (start, "del" + gene[start : start + size])
                muts[mut] += 1
                dump_arr.append(mut)
                start += size
            elif op == 1:  # Insertion
                mut = (start, "ins" + read.query_sequence[s_start : s_start + size])
                muts[mut] += 1
                # HACK: just store the length due to seq. errors
                self._insertion_counts[start, size] += 1
                dump_arr.append(mut)
                s_start += size
            elif op == 4:  # Soft-clip
                s_start += size
            elif op in [0, 7, 8]:  # M, X and =
                for i in range(size):
                    if (
                        start + i in gene
                        and gene[start + i] != read.query_sequence[s_start + i]
                    ):
                        mut = (
                            start + i,
                            f"{gene[start + i]}>{read.query_sequence[s_start + i]}",
                        )
                        dump_arr.append(mut)
                        muts[mut] += 1
                    else:  # We ignore all mutations outside the RefSeq region
                        norm[start + i] += 1
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
                for p in range(len(l)):
                    if l[p] != ".":
                        muts[pos + p, f"{l[p]}>{r[p]}"] -= 1
                        if p:
                            norm[pos + p] += 1
                muts[pos, op] += 1

        read_pos = (read.reference_start, start, len(read.query_sequence))
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

    def _get_phases(self, sam_path, gene, reference, max_insert=1000):
        index = sorted(p for p, m in self.coverage._coverage.items() if len(m) > 1)
        fragments = []

        def add(lo, extra, *args):
            frag = {}
            for r in args:
                start, s_start = r.reference_start, 0
                for op, size in r.cigartuples:
                    if op == 2:  # Deletion
                        start += size
                    elif op == 1:  # Insertion
                        s_start += size
                    elif op == 4:  # Soft-clip
                        s_start += size
                    elif op in [0, 7, 8]:  # M, X and =
                        for i in range(size):
                            while lo < len(index) and index[lo] < start + i:
                                lo += 1
                            if lo < len(index) and index[lo] == start + i:
                                if r.query_qualities[s_start + i] >= 10:
                                    frag.setdefault(start + i, set()).add(
                                        r.query_sequence[s_start + i]
                                    )
                        start += size
                        s_start += size
            read = []
            for pos, snps in frag.items():
                if len(snps) == 1:  # avoid contradicting pairs
                    read.append((pos, snps.pop()))
            if len(read) > 1:  # links at least 2 mutations
                fragments.append((extra, read))

        with pysam.AlignmentFile(sam_path, reference_filename=reference) as sam:
            sam.check_index()

            vlo = 0
            seen = {}
            # Fetch the reads
            for r in sam.fetch(
                region=gene.get_wide_region().samtools(prefix=self._prefix)
            ):
                if not r.cigartuples:  # only valid alignments
                    continue
                if r.is_supplementary or r.is_secondary or r.is_duplicate:
                    continue
                if "H" in r.cigarstring:  # avoid hard-clipped reads
                    continue
                if not _in_region(gene.get_wide_region(), r, self._prefix):
                    continue  # ensure that it is a proper gene read

                # any SNPs here?
                while vlo < len(index) and index[vlo] < r.reference_start:
                    vlo += 1
                if vlo == len(index) or index[vlo] >= r.reference_end:
                    continue

                name = r.query_name
                if len(name) > 2 and (name[-2] == "#" or name[-2] == "/"):
                    name = name[:-2]

                # 10X data support
                extra = name
                if r.has_tag("BX"):
                    extra = r.get_tag("BX")
                if extra and r.has_tag("XC"):
                    extra = f"{extra}:{r.get_tag('XC')}"
                elif extra and r.has_tag("MI"):
                    extra = f"{extra}:{r.get_tag('MI')}"

                if (
                    r.mate_is_unmapped
                    or r.reference_id != r.next_reference_id
                    or r.template_length > max_insert
                ):
                    add(vlo, extra, r)
                elif name in seen:  # Paired reads
                    add(vlo, extra, seen[name][0], r)
                    del seen[name]
                elif (
                    r.next_reference_start < r.reference_start
                ):  # already seen but not added
                    add(vlo, extra, r)
                else:
                    seen[name] = (r, vlo, extra)
            for name, (r, vlo, extra) in seen.items():
                add(vlo, extra, r)

        fragments.sort(key=lambda x: x[0])
        new_fragments = []
        prev_fid, prev_read = "", []
        for fid, read in fragments:
            if fid != prev_fid:
                if prev_read:
                    new_fragments.append(prev_read)
                prev_fid, prev_read = fid, read
            else:
                prev_read += read
        if prev_read:
            new_fragments.append(prev_read)
        for r in new_fragments:
            r.sort()

        self._dump_phase = new_fragments
        return new_fragments

    # ----------------------------------------------------------------------------------
    # Coverage-rescaling functions
    # ----------------------------------------------------------------------------------

    def detect_cn(self, gene: Gene, profile: str, cn_region: GRange) -> None:
        """
        Rescale the ``self.coverage`` to fit the sequencing profile.

        :param gene: Gene instance.
        :param profile: Profile identifier. Can be a SAM/BAM file as well.
        :param cn_region: Coordinates of the copy-number neutral region.

        .. notes:: This function assumes that ``self.coverage`` is set and modifies it.
        """

        if os.path.exists(profile) and os.path.isfile(profile):
            ext = os.path.splitext(profile)
            if ext[-1] in [".bam", ".sam"]:
                prof = self._load_profile(
                    profile,
                    is_bam=True,
                    gene=gene,
                    cn_region=cn_region,
                )
            else:
                prof = self._load_profile(profile)
        else:
            profile_path = script_path(
                "aldy.resources.profiles/{}.yml".format(profile.lower())
            )
            prof = self._load_profile(profile_path)
            if "neutral" not in prof or "cn" not in prof["neutral"]:
                raise AldyException("Profile missing neutral region")
            if gene.name not in prof:
                raise AldyException(f"Profile missing {gene.name}")
        self.coverage._normalize_coverage(
            prof[gene.name], gene.regions, cn_region, prof["neutral"]["cn"][0]
        )

    def _load_profile(
        self,
        profile_path: str,
        is_bam: bool = False,
        gene: Optional[Gene] = None,
        cn_region: Optional[GRange] = None,
    ) -> Dict[str, Dict[str, List[float]]]:
        """
        Load a coverage profile.

        :param profile_path:
            Path to a profile file.
            SAM/BAM files are also accepted if `is_bam` is set
            (profile will be dynamically calculated in that case).
        :param is_bam: A flag indicating if the `profile_path` is SAM/BAM or not.
        :param gene_region: Region to be extracted from the profile.
        :param cn_region: Copy-number neutral region to be extracted from the profile.

        :return: A profile dictionary where keys are chromosome IDs (e.g. '7' for chr7)
                 and values are dictionaries that map the genomic locus to the
                 corresponding profile coverage.
                 This is a :obj:`collections.defaultdict` that uses 0 as a placeholder
                 for the missing loci.

        :raise: :obj:`aldy.common.AldyException` if ``is_bam`` is set
                but ``gene_region`` and ``cn_region`` are not.
        """

        if is_bam:
            assert gene and cn_region
            regions = {}
            for gi, gr in enumerate(gene.regions):
                for r, rng in gr.items():
                    regions[gene.name, r, gi] = rng
            regions["neutral", "cn", 0] = cn_region
            return load_sam_profile(profile_path, regions=regions)
        else:
            with open(profile_path) as f:
                return yaml.safe_load(f)

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

            fd.write(struct.pack("<l", len(self._dump_phase)))
            for r in self._dump_phase:
                fd.write(struct.pack("<l", len(r)))
                for p, md in r:
                    fd.write(struct.pack("<ll", p, len(md)))
                    fd.write(md.encode("ascii"))


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
    gene_regions["neutral", "cn", 0] = (
        cn_region if cn_region else DEFAULT_CN_NEUTRAL_REGION[genome]
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
            d[g][r][ri] = sum(1.0 for i in range(s, e))
        else:
            d[g][r][ri] = sum(cov[c][i] * (factor / 2.0) for i in range(s, e))
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

    a = (read.reference_start, read.reference_end)
    b = (region.start, region.end)
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]
