# 786
# Aldy source: sam.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple, Dict, List, Optional, Any

import pysam
import os
import os.path
import re
import gzip
import struct
import collections

from natsort import natsorted
from .common import log, GRange, AldyException, script_path
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


class Sample:
    """
    Interface for reading SAM/BAM/CRAM files that parses and stores read alignments.
    """

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

        if vcf_path:
            self.load_vcf(vcf_path, gene, debug)
        else:
            assert sam_path and profile
            self.load_aligned(sam_path, gene, threshold, reference, cn_region, debug)
            if cn_region:
                self.detect_cn(gene, profile, cn_region)
            log.debug("[sam] avg_coverage= {:.1f}x", self.coverage.average_coverage())
            if cn_region and self.coverage.diploid_avg_coverage() < 2:
                raise AldyException(
                    "The average coverage of the sample is too low ({:.1f}).".format(
                        self.coverage.diploid_avg_coverage()
                    )
                )

    def load_aligned(
        self,
        sam_path: str,
        gene: Gene,
        threshold: float,
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

        # dict of int: int: the coverage of the CN-neutral region
        cnv_coverage: Dict[int, int] = collections.defaultdict(int)
        # Get the list of indel sites that should be corrected
        # TODO: currently uses only functional indels; other indels should be
        # corrected as well
        _insertion_sites = {
            m for _, a in gene.alleles.items() for m in a.func_muts if m.op[:3] == "ins"
        }
        _insertion_reads: Dict[Tuple, int] = collections.defaultdict(int)
        _multi_sites = {
            m.pos: m.op
            for _, a in gene.alleles.items()
            for m in a.func_muts
            if ">" in m.op and len(m.op) > 3
        }
        muts: dict = collections.defaultdict(int)
        norm: dict = collections.defaultdict(int)

        if sam_path[-5:] == ".dump":
            self._load_dump(sam_path)
        else:
            with pysam.AlignmentFile(sam_path, reference_filename=reference) as sam:
                # Check do we have proper index to speed up the queries
                try:
                    sam.check_index()
                except AttributeError:
                    pass  # SAM files do not have an index. BAMs might also lack it
                except ValueError:
                    raise AldyException(
                        f"File {sam_path} has no index (it must be indexed)"
                    )

                #: str: Check if we need to append 'chr' or not
                self._prefix = _chr_prefix(
                    gene.chr, [x["SN"] for x in sam.header["SQ"]]
                )

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
                                cnv_coverage[start + i] += 1
                            start += size

                # Set it to _fetched_ if a CN-neutral region is user-provided.
                # Then read the CN-neutral region.
                is_cn_region_fetched = cn_region is None
                if cn_region and sam.has_index():
                    log.debug("[sam] has_index= True")
                    for read in sam.fetch(
                        region=cn_region.samtools(prefix=self._prefix)
                    ):
                        read_cn_read(read)
                    is_cn_region_fetched = True

                # Fetch the reads
                total = 0
                dump_data = []
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

                    r = self._parse_read(
                        read, gene, norm, muts, _multi_sites, bool(debug)
                    )
                    if r:
                        total += 1
                        if _insertion_sites:
                            _insertion_reads[r[0]] += 1
                        if debug:
                            dump_data.append(r)
                if debug:
                    self._dump()
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
        self._group_deletions(gene, coverage)
        for pos, op in _multi_sites.items():
            if pos in coverage and op in coverage[pos]:
                log.debug(f"[sam] multi-SNP {pos}{op} with {coverage[pos][op]} reads")
        for mut in _insertion_sites:
            if mut.pos in coverage and mut.op in coverage[mut.pos]:
                total = sum(coverage[mut.pos].values())
                if total > 0:
                    coverage[mut.pos][mut.op] = self._correct_ins_coverage(
                        mut, coverage[mut.pos][mut.op], total, _insertion_reads
                    )

        #: dict of int: (dict of str: int): coverage dictionary
        #: keys are loci in the reference genome, while values are dictionaries
        #: that describe the coverage of each mutation
        # (_ stands for non-mutated nucleotide)
        self.coverage = Coverage(
            {p: {m: v for m, v in coverage[p].items() if v > 0} for p in coverage},
            threshold,
            cnv_coverage,
            os.path.basename(sam_path).split(".")[0],
        )

    def _group_deletions(self, gene, coverage):
        # Group ambiguous deletions
        indels = collections.defaultdict(set)
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
        indels = collections.defaultdict(set)
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
                    coverage[pos][mut] = 0

    def _parse_read(
        self, read, gene: Gene, norm, muts, multi_sites=None, dump=False,
    ) -> Optional[Tuple[Tuple[int, int, int], Any]]:
        """
        Parse a :obj:`pysam.AlignedSegment` read.

        :param read: pysam read.
        :param gene: Gene instance.
        :param norm: Positions within the read that have not been mutated.
        :param muts: Positions within the read that have been mutated.

        .. note:: `norm` and `muts` are modified.
        """

        if not _in_region(
            gene.get_wide_region(), read, self._prefix
        ):  # ensure that it is a proper gene read
            return None
        if not read.cigartuples:  # only valid alignments
            return None
        if read.flag & 0x800:  # avoid supplementary alignments
            return None
        if "H" in read.cigarstring:  # avoid hard-clipped reads
            return None

        dump_arr = {}
        start, s_start = read.reference_start, 0
        for op, size in read.cigartuples:
            if op == 2:  # Deletion
                mut = (
                    start,
                    "del" + gene[start : start + size],
                )
                muts[mut] += 1
                dump_arr[start] = mut[1]
                start += size
            elif op == 1:  # Insertion
                mut = (start, "ins" + read.query_sequence[s_start : s_start + size])
                muts[mut] += 1
                dump_arr[start] = mut[1]
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
                        dump_arr[start + i] = mut[1]
                        muts[mut] += 1
                    else:  # We ignore all mutations outside the RefSeq region
                        norm[start + i] += 1
                start += size
                s_start += size

        if multi_sites and set(multi_sites.keys()) & set(dump_arr.keys()):
            for pos, op in multi_sites.items():
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
                            muts[pos + p, dump_arr[pos + p]] -= 1
                            if p:
                                norm[pos + p] += 1
                    muts[pos, op] += 1
        return (
            (read.reference_start, start, len(read.query_sequence)),
            list(dump_arr.items()),
        )

    def _correct_ins_coverage(self, mut: Mutation, orig_cov, total, reads) -> int:
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

        new_total = 0.0
        for (orig_start, orig_end, read_len), cov in reads.items():
            ins_len = len(mut.op) - 3
            # HACK: This is rather hackish way of detecting *40
            # min 20% on the sides, max ... how much?
            min_margin = MARGIN_RATIO * max(1, ins_len / INSERT_LEN_RESCALE) * read_len
            if orig_end >= mut.pos + max(
                INSERT_LEN_AMPLIFY * ins_len, min_margin
            ) and orig_start <= mut.pos - max(
                (INSERT_LEN_AMPLIFY - 1) * ins_len, min_margin
            ):
                new_total += cov
        if new_total >= orig_cov:
            new_cov = int(total * (orig_cov / new_total))
            if new_cov > orig_cov:
                log.debug(f"[sam] rescale {mut} from {orig_cov} to {new_cov} ")
                return new_cov
        return orig_cov

    def load_vcf(self, vcf_path: str, gene: Gene, debug: Optional[str] = None):
        """
        Load the read, mutation and coverage data from a VCF file.
        """

        log.debug("[vcf] path= {}", os.path.abspath(vcf_path))

        _multi_sites = {
            m.pos: m.op
            for an, a in gene.alleles.items()
            for m in a.func_muts
            if ">" in m.op and len(m.op) > 3
        }
        muts: dict = collections.defaultdict(int)
        norm = {
            p: 20
            for p in range(
                gene.get_wide_region().start - 500, gene.get_wide_region().end + 1
            )
        }

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
            #: str: Check if we need to append 'chr' or not
            self._prefix = _chr_prefix(gene.chr, list(vcf.header.contigs))

            samples = list(vcf.header.samples)
            sample = samples[0]
            if len(samples) > 1:
                log.warn(
                    "WARNING: Multiple VCF samples found; only using the first one."
                )
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
                    hgvs = [(read.pos - 1, f"_")]
                hgvs += [get_mut(read.pos - 1, read.ref, a) for a in read.alleles[1:]]
                for gt in g:
                    pos, op = hgvs[gt]
                    if op == "_":
                        continue
                    muts[pos, op] += 10
                    norm[pos] -= 10
                    dump_arr[pos] = op
                if _multi_sites and set(_multi_sites.keys()) & set(dump_arr.keys()):
                    for pos, op in _multi_sites.items():
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
                                    muts[pos + p, dump_arr[pos + p]] -= 1
                            muts[pos, op] += 1

        # Establish the coverage dictionary
        coverage: Dict[int, Dict[str, int]] = dict()
        for pos, cov in norm.items():
            if cov == 0:
                continue
            if pos not in coverage:
                coverage[pos] = {}
            coverage[pos]["_"] = cov
        for m, cov in muts.items():
            pos, mut = m
            if pos not in coverage:
                coverage[pos] = {}
            if not min(gene.chr_to_ref) <= pos <= max(gene.chr_to_ref):
                mut = "_"  # ignore mutations outside of the region of interest
            coverage[pos][mut] = cov
        self._group_deletions(gene, coverage)
        for pos, op in _multi_sites.items():
            if pos in coverage and op in coverage[pos]:
                log.debug(f"[sam] multi-SNP {pos}{op} with {coverage[pos][op]} reads")
        self.coverage = Coverage(
            {p: {m: v for m, v in coverage[p].items() if v > 0} for p in coverage},
            0.5,
            {},
            sample,
        )

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
                    gene_region=gene.get_wide_region(),
                    cn_region=cn_region,
                )
            else:
                prof = self._load_profile(profile)
        else:
            profile_path = script_path(
                "aldy.resources.profiles/{}.profile".format(profile.lower())
            )
            prof = self._load_profile(profile_path)
        self.coverage._normalize_coverage(prof, gene.regions, cn_region)

    def _load_profile(
        self,
        profile_path: str,
        is_bam: bool = False,
        gene_region: Optional[GRange] = None,
        cn_region: Optional[GRange] = None,
    ) -> Dict[str, Dict[int, float]]:
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

        profile: Dict[str, Dict[int, float]] = collections.defaultdict(
            lambda: collections.defaultdict(int)
        )
        if is_bam:
            if gene_region and cn_region:
                ptr = Sample.load_sam_profile(
                    profile_path, regions=[gene_region, cn_region]
                )
                for _, c, p, v in ptr:
                    profile[c][p] = v
            else:
                raise AldyException("Region parameters must be provided")
        else:

            def read_file(f, profile):
                for line in f:
                    if line[0] == "#":
                        continue  # skip comments
                    _, ch, pos, val = line.strip().split()
                    profile[ch.decode("utf-8")][int(pos)] = float(val)

            try:
                with gzip.open(profile_path) as f:
                    read_file(f, profile)
            except OSError:
                with open(profile_path) as f:
                    read_file(f, profile)
        return profile

    def _dump(self):
        with gzip.open(f"{debug}.dump", "wb") as fd:
            if len(cnv_coverage) == 0:
                m, M = 0, -1
            else:
                m, M = min(cnv_coverage.keys()), max(cnv_coverage.keys())
            fd.write(struct.pack("<ll", m, M))
            for i in range(m, M + 1):
                fd.write(struct.pack("<i", cnv_coverage[i]))
            fd.write(struct.pack("<l", len(dump_data)))
            for (s, e, l), m in dump_data:
                fd.write(struct.pack("<llhh", s, e - s, l, len(m)))
                for p, md in m:
                    fd.write(struct.pack("<hh", p - s, len(md)))
                    fd.write(md.encode("ascii"))

    def _load_dump(self, path):
        with gzip.open(path, "rb") as fd:
            log.warn("Loading debug dump from {}", path)
            l, h, i = (
                struct.calcsize("<l"),
                struct.calcsize("<h"),
                struct.calcsize("<i"),
            )
            m, M = struct.unpack("<ll", fd.read(l + l))
            for j in range(m, M + 1):
                (cnv_coverage[j],) = struct.unpack("<i", fd.read(i))
            (ld,) = struct.unpack("<l", fd.read(l))
            for _ in range(ld):
                ref_start, ref_end, read_len, num_mutations = struct.unpack(
                    "<llhh", fd.read(l + l + h + h)
                )
                ref_end += ref_start
                for j in range(ref_start, ref_end):
                    norm[j] += 1
                insertions = set()
                for _ in range(num_mutations):
                    mut_start, op_len = struct.unpack("<hh", fd.read(h + h))
                    mut_start += ref_start
                    op = fd.read(op_len).decode("ascii")
                    muts[mut_start, op] += 1
                    if op[:3] == "del":
                        for j in range(0, len(op) - 3):
                            norm[mut_start + j] -= 1
                    elif op[:3] == "ins":
                        insertions.add((mut_start, len(op) - 3))
                    else:
                        norm[mut_start] -= 1
                for ins in _tandem_sites:
                    ins_len = len(ins.op) - 3
                    if (ins.pos, ins_len) in insertions:
                        _tandem_sites[ins] = (
                            _tandem_sites[ins][0],
                            _tandem_sites[ins][1] + 1,
                        )
                    else:
                        _tandem_sites[ins][0][(ref_start, ref_end, read_len)] += 1

    # ----------------------------------------------------------------------------------

    @staticmethod
    def detect_genome(sam_path: str):
        def detect_vcf(sam_path):
            with pysam.VariantFile(sam_path):
                return True

        def detect_sam(sam_path):
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
                    return "hg19"
                elif is_hg38:
                    return "hg38"
                else:
                    return None

        try:
            return "sam", detect_sam(sam_path)
        except ValueError:
            if detect_vcf(sam_path):
                return "vcf", None
        return "", None

    @staticmethod
    def load_sam_profile(
        sam_path: str,
        factor: float = 2.0,
        regions: Optional[List[GRange]] = None,
        cn_region: Optional[GRange] = None,
    ) -> List[Tuple[str, str, int, float]]:
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
        if regions is None:
            import pkg_resources
            import yaml

            gene_regions = []
            for g in pkg_resources.resource_listdir("aldy.resources", "genes"):
                if g[-4:] != ".yml":
                    continue
                with open(script_path(f"aldy.resources.genes/{g}")) as f:
                    yml = yaml.safe_load(f)
                    ref = yml["reference"]["mappings"]["hg19"][0]
                    reg = yml["structure"]["regions"]["hg19"]
                    rmi = min(min(r) for _, r in reg.items())
                    rma = max(max(r) for _, r in reg.items())
                    gene_regions.append((yml["name"], GRange(ref, rmi, rma)))
            gene_regions = natsorted(gene_regions, key=lambda x: x[1])
        else:
            gene_regions = [(str(i), r) for i, r in enumerate(sorted(regions))]
        gene_regions.append(
            ("CN", cn_region if cn_region else DEFAULT_CN_NEUTRAL_REGION["hg19"])
        )
        result: List[Tuple[str, str, int, float]] = []
        for gene, location in gene_regions:
            with pysam.AlignmentFile(sam_path) as sam:
                prefix = _chr_prefix(location.chr, sam.header["SQ"])
                region = location.samtools(pad_left=1000, pad_right=1000, prefix=prefix)
                cov: dict = collections.defaultdict(int)
                log.info("Generating profile for {} ({})", gene, region)
                try:
                    for read in sam.fetch(region=region):
                        start, s_start = read.reference_start, 0
                        if not read.cigartuples:
                            continue

                        for op, size in read.cigartuples:
                            if op == 2:
                                for i in range(size):
                                    cov[start + i] += 1
                                start += size
                            elif op == 1:
                                s_start += size
                            elif op == 4:
                                s_start += size
                            elif op in [0, 7, 8]:
                                for i in range(size):
                                    cov[start + i] += 1
                                start += size
                                s_start += size
                    result += [
                        (gene, location.chr, p, cov[p] * (factor / 2.0))
                        for p in range(location.start, location.end)
                    ]
                except ValueError:
                    log.warn("Cannot fetch gene {} ({})", gene, region)
        return result


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

    return (
        read.reference_id != -1
        and read.reference_name == prefix + region.chr
        and region.start - 500 <= read.reference_start <= region.end
    )
