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

from .common import log, GRange, AldyException, script_path
from .coverage import Coverage
from .gene import Gene, Mutation


CIGAR2CODE: Dict[int, int] = {ord(y): x for x, y in enumerate("MIDNSHP=XB")}
"""str: CIGAR characters for an ASCII ordinal of CIGAR character"""

CIGAR_REGEX = re.compile(r"(\d+)([MIDNSHP=XB])")
"""regex: Regex that matches valid CIGAR strings"""


DEFAULT_CN_NEUTRAL_REGION = GRange("22", 42547463, 42548249)
"""obj:`aldy.common.GRange` Default copy-number neutral region
   (exon 4-6 of the CYP2D8 gene)"""


class Sample:
    """
    Interface for reading SAM/BAM/CRAM files that parses
    and stores read alignments.

    Attributes:
        coverage (:obj:`aldy.coverage.Coverage`):
            The coverage data for the sample.
            Consult the :obj:`aldy.coverage.Coverage`.
        reads (dict[str, tuple[:obj:`pysam.AlignedSegment`,
               :obj:`pysam.AlignedSegment`]]):
            Dictionary that keeps the paired-end reads together in a tuple.
            Key is the read name. Used only if ``phase`` is set.

    Static methods:
        ``load_sam_profile`` (documentation below).
    """

    def __init__(
        self,
        sam_path: str,
        gene: Gene,
        threshold: float,
        profile: str,
        phase: bool = False,
        reference: Optional[str] = None,
        cn_region: Optional[GRange] = DEFAULT_CN_NEUTRAL_REGION,
        debug: Optional[str] = None,
    ) -> None:
        """
        Initialize a :obj:`Sample` object.

        Args:
            sam_path (str):
                Path to a SAM/BAM/CRAM file.
            gene (:obj:`aldy.gene.Gene`):
                Gene instance.
            threshold (float):
                Threshold for filtering out low quality mutations. Ranges from 0 to 1.
                Check :obj:`aldy.coverage.Coverage` for more information.
            profile (str, optional):
                Profile specification (e.g. 'prgnseq-v1') or a profile SAM/BAM file.
            phase (bool):
                Construct basic rudimentary phasing of the reads to aid the genotyping.
                DEPRECATED: currently does nothing.
                Default is ``False``.
            reference (str, optional):
                Reference genome for reading CRAM files.
                Default is None.
            cn_region (:obj:`aldy.common.GRange`, optional):
                Copy-number neutral region to be used for coverage rescaling.
                If None, profile loading and coverage rescaling will be skipped
                (and Aldy will require a ``--cn`` parameter to be user-provided).
                Default is ``DEFAULT_CN_NEUTRAL_REGION``.
            debug (str, optional):
                If set, create a "`debug`.dump" file for debug purposes.
                Default is None.

        Raises:
            :obj:`aldy.common.AldyException` if the average coverage of the
            copy-number neutral region is too low (less than 2).
        """

        self.load_aligned(sam_path, gene, threshold, phase, reference, cn_region, debug)
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
        phase: bool = False,
        reference: Optional[str] = None,
        cn_region: Optional[GRange] = None,
        debug: Optional[str] = None,
    ) -> None:
        """
        Load the read, mutation and coverage data from a SAM/BAM/CRAM file.

        Args:
            sam_path (str):
                Path to a SAM/BAM/CRAM file.
            gene (:obj:`aldy.gene.Gene`):
                Gene instance.
            threshold (float):
                Threshold for filtering out low quality mutations. Ranges from 0 to 1.
                Check :obj:`aldy.coverage.Coverage` for more information.
            phase (bool):
                Construct basic rudimentary phasing of the reads to aid the genotyping.
                DEPRECATED: currently does nothing.
                Default is False.
            reference (str, optional):
                Reference genome for reading CRAM files.
                Default is None.
            cn_region (:obj:`aldy.common.GRange`, optional):
                Copy-number neutral region to be used for coverage rescaling.
                If None, profile loading and coverage rescaling will be skipped
                (and Aldy will require a ``--cn`` parameter to be user-provided).
                Default is ``DEFAULT_CN_NEUTRAL_REGION`` (CYP2D8).
            debug (str, optional):
                If set, create a "`debug`.dump" file for debug purposes.
                Default is None.

        Raises:
            :obj:`aldy.common.AldyException` if a BAM/CRAM file lacks an index.
        """

        # store paired-end reads (only if phasing is on)
        self.reads: Dict[str, Tuple[pysam.AlignedSegment, pysam.AlignedSegment]] = {}

        log.debug("[sam] path= {}", os.path.abspath(sam_path))

        # dict of int: int: the coverage of the CN-neutral region
        cnv_coverage: Dict[int, int] = collections.defaultdict(int)
        # Get the list of indel sites that should be corrected
        # TODO: currently uses only functional indels; other indels should be
        # corrected as well
        _indel_sites: Dict[Mutation, tuple] = {
            m: (collections.defaultdict(int), 0)
            for an, a in gene.alleles.items()
            for m in a.func_muts
            if m.op[:3] == "ins"
        }
        _multi_sites: Dict[Mutation, tuple] = {
            m.pos: m.op
            for an, a in gene.alleles.items()
            for m in a.func_muts
            if ">" in m.op and len(m.op) > 3
        }
        muts: dict = collections.defaultdict(int)
        norm: dict = collections.defaultdict(int)

        if sam_path[-5:] == ".dump":
            if phase:
                raise AldyException("Debug dumps do not work with --phase parameter")
            with gzip.open(sam_path, "rb") as fd:
                log.warn("Loading debug dump from {}", sam_path)
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
                    for __ in range(num_mutations):
                        mut_start, op_len = struct.unpack("<hh", fd.read(h + h))
                        mut_start += ref_start
                        op = fd.read(op_len).decode("ascii")
                        muts[mut_start, op] += 1
                        if op[:3] == "DEL":
                            for j in range(0, len(op) - 4):
                                norm[mut_start + j] -= 1
                        elif op[:3] == "INS":
                            insertions.add((mut_start, len(op) - 4))
                        else:
                            norm[mut_start] -= 1
                    for ins in _indel_sites:
                        ins_len = len(ins.op) - 4
                        if (ins.pos, ins_len) in insertions:
                            _indel_sites[ins] = (
                                _indel_sites[ins][0],
                                _indel_sites[ins][1] + 1,
                            )
                        else:
                            _indel_sites[ins][0][(ref_start, ref_end, read_len)] += 1
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
                        read,
                        gene,
                        norm,
                        muts,
                        _indel_sites,
                        _multi_sites,
                        debug is not None,
                    )
                    if r:
                        total += 1
                        if debug:
                            dump_data.append(r)
                        if phase:
                            if read.query_name in self.reads:
                                assert self.reads[read.query_name][1] is None
                                read1 = self.reads[read.query_name][0]
                                if read.reference_start < read1.reference_start:
                                    read, read1 = read1, read
                                self.reads[read.query_name] = (read1, read)
                            else:
                                self.reads[read.query_name] = (read, None)

                if debug:
                    if phase:
                        raise AldyException(
                            "Debug dumps do not work with --phase parameter"
                        )
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
        for mut, ins_cov in _indel_sites.items():
            if mut.pos in coverage and mut.op in coverage[mut.pos]:
                self._correct_ins_coverage(mut, ins_cov)

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
                deleted = "".join(gene[i] for i in range(pos, pos + sz))
                for p in range(pos - sz, -1, -sz):
                    if "".join(gene[i] for i in range(p, p + sz)) != deleted:
                        break
                    potential.append(p)
                for p in range(pos + sz, max(coverage), sz):
                    if "".join(gene[i] for i in range(p, p + sz)) != deleted:
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
                    log.debug(f"[sam] relocate {pos}{mut} to {new_pos}{mut}")
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
                    if "".join(gene[i] for i in range(p, p + sz)) != inserted:
                        break
                for p in range(pos, max(coverage), sz):
                    potential.append(p)
                    if "".join(gene[i] for i in range(p, p + sz)) != inserted:
                        break
                potential = [
                    p
                    for p in set(potential) & set(indels)
                    for m in indels[p]
                    if m == mut
                ]
                if len(potential) == 1 and potential[0] != pos:
                    new_pos = potential[0]
                    log.debug(f"[sam] relocate {pos}{mut} to {new_pos}{mut}")
                    if mut not in coverage[new_pos]:
                        coverage[new_pos][mut] = 0
                    coverage[new_pos][mut] += coverage[pos][mut]
                    coverage[pos][mut] = 0

    def _parse_read(
        self,
        read: pysam.AlignedSegment,
        gene: Gene,
        norm: dict,
        muts: dict,
        indel_sites=None,
        multi_sites=None,
        dump=False,
    ) -> Optional[Tuple[Tuple[int, int, int], Any]]:
        """
        Parse a :obj:`pysam.AlignedSegment` read.

        Params:
            read (:obj:`pysam.AlignedSegment`)
            gene (:obj:`aldy.gene.Gene`)

        Returns:
            optional: None if parsing was not successful.
            Otherwise, returns the tuple consisting of:

                - start and end positions of the read, and
                - the list of mutations found in the read.

        Params that are mutated:
            norm (:obj:`collections.defaultdict(int)`):
                dict of positions within the read that have not been mutated.
            muts (:obj:`collections.defaultdict(int)`):
                dict of positions within the read that have been mutated.
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

        insertions = set()
        dump_arr = {}
        start, s_start = read.reference_start, 0
        for op, size in read.cigartuples:
            if op == 2:  # Deletion
                mut = (
                    start,
                    "del" + "".join(gene[i] for i in range(start, start + size)),
                )
                muts[mut] += 1
                dump_arr[start] = mut[1]
                start += size
            elif op == 1:  # Insertion
                mut = (start, "ins" + read.query_sequence[s_start : s_start + size])
                muts[mut] += 1
                dump_arr[start] = mut[1]
                insertions.add((start, size))
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
                    muts[pos, op] += 1
        if indel_sites:
            for ins in indel_sites:
                ins_len = len(ins.op) - 4
                if (ins.pos, ins_len) in insertions:
                    indel_sites[ins] = (indel_sites[ins][0], indel_sites[ins][1] + 1)
                else:
                    indel_sites[ins][0][
                        (read.reference_start, start, len(read.query_sequence))
                    ] += 1
        return (
            (read.reference_start, start, len(read.query_sequence)),
            list(dump_arr.items()),
        )

    def _correct_ins_coverage(self, mut: Mutation, j) -> None:
        """
        Fix low coverage of large tandem insertions.

        Targets the cases where reference looks like
            ...X...
        and the donor genome looks like
            ...XX... (i.e. X is a tandem insertion).
        Any read that covers only one tandem (e.g. read X...)
        will get perfectly aligned (without trigerring an insertion tag)
        to the tail of the tandem insertion.
        This results in under-estimation of the insertion abundance will
        (as aligner could not assign `I` CIGAR to the correct insertion).
        This function attempts to correct this bias.

        Notes:
            This function modifies ``self.coverage``.
        """

        MARGIN_RATIO = 0.20
        INSERT_LEN_RESCALE = 10.0
        INSERT_LEN_AMPLIFY = 3

        total = self.coverage.total(mut.pos)
        if total == 0:
            return
        current_ratio = float(self.coverage[mut]) / total

        j0 = 0
        for orig_start, start, read_len in j[0]:
            ins_len = len(mut.op) - 4
            # HACK: This is rather hackish way of detecting *40
            # min 20% on the sides, max ... how much?
            min_margin = MARGIN_RATIO * max(1, ins_len / INSERT_LEN_RESCALE) * read_len
            if start < mut.pos + max(
                INSERT_LEN_AMPLIFY * ins_len, min_margin
            ) or orig_start > mut.pos - max(
                (INSERT_LEN_AMPLIFY - 1) * ins_len, min_margin
            ):
                continue
            j0 += j[0][(orig_start, start, read_len)]

        new_ratio = j[1] / float(j0 + j[1])
        new_coverage = int(self.coverage[mut] * (new_ratio / current_ratio))
        log.debug(
            f"[sam] rescale {mut} from {self.coverage[mut]} to {new_coverage} "
            + f"(indel:total = {j[1]}:{j0 + j[1]})",
        )
        self.coverage._coverage[mut.pos][mut.op] = new_coverage

    # ----------------------------------------------------------------------------------
    # Coverage-rescaling functions
    # ----------------------------------------------------------------------------------

    def detect_cn(self, gene: Gene, profile: str, cn_region: GRange) -> None:
        """
        Rescale the ``self.coverage`` to fit the sequencing profile.

        Params:
            gene (:obj:`aldy.gene.Gene`):
                Gene instance.
            profile (str):
                Profile identifier (e.g. 'pgrnseq-v1'). Can be a SAM/BAM file as well.
            cn_region (:obj:`aldy.common.GRange`):
                Coordinates of the copy-number neutral region.

        Notes:
            This function assumes that ``self.coverage`` is set.
            It modifies ``self.coverage``.
        """

        if os.path.exists(profile):
            ext = os.path.splitext(profile)
            if ext[-1] in [".bam", ".sam"]:
                prof = self._load_profile(
                    profile, is_bam=True, gene_region=gene.region, cn_region=cn_region
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

        Returns:
            defaultdict[str, dict[int, float]]: A profile dictionary where
            keys are chromosome IDs (e.g. '7' for chr7)
            and values are dictionaries that map the genomic locus to the
            corresponding profile coverage.
            This is a :obj:`collections.defaultdict` that uses 0 as a placeholder
            for the missing loci.

        Args:
            profile_path:
                Path to a profile file.
                SAM/BAM files are also accepted if `is_bam` is set
                (profile will be dynamically calculated in that case).
            is_bam (bool):
                A flag indicating if the `profile_path` is SAM/BAM or not.
            gene_region (:obj:`aldy.common.GRange`, optional):
                Region to be extracted from the profile.
                Default is None.
            cn_region (:obj:`aldy.common.GRange`, optional):
                Copy-number neutral region to be extracted from the profile.
                Default is None.

        Raises:
            :obj:`aldy.common.AldyException` if ``is_bam`` is set but ``gene_region``
            and ``cn_region`` are not.
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
            with open(profile_path) as f:
                for line in f:
                    if line[0] == "#":
                        continue  # skip comments
                    ch, pos, val = line.strip().split()[1:]
                    profile[ch][int(pos)] = float(val)
        return profile

    # ----------------------------------------------------------------------------------

    @staticmethod
    def detect_genome(sam_path: str):
        with pysam.AlignmentFile(sam_path) as sam:
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

                1. PGRNseq-v1: PGXT104 was used for all genes
                   (n.b. PGXT147 with rescale 2.52444127771 was used for CYP2B6 beta).
                2. PGRNseq-v2: NA19789.bam was used for all genes.
                3. Illumina: by definition contains all ones (uniform coverage profile).
        """
        if regions is None:
            gene_regions = sorted(
                [  # paper gene coordinates in hg19
                    # TODO: Auto populate
                    ("CYP3A5", GRange("7", 99245000, 99278000)),
                    ("CYP3A4", GRange("7", 99354000, 99465000)),
                    ("CYP2C19", GRange("10", 96445000, 96615000)),
                    ("CYP2C9", GRange("10", 96691000, 96754000)),
                    ("CYP2C8", GRange("10", 96796000, 96830000)),
                    ("CYP4F2", GRange("19", 15619000, 16009500)),
                    ("CYP2A6", GRange("19", 41347500, 41400000)),
                    ("CYP2D6", GRange("22", 42518900, 42553000)),
                    ("TPMT", GRange("6", 18126541, 18157374)),
                    ("DPYD", GRange("1", 97541298, 98388615)),
                ],
                key=lambda x: x[1],
            )
        else:
            gene_regions = [(str(i), r) for i, r in enumerate(sorted(regions))]
        if cn_region:
            gene_regions.append(("CN", cn_region))
        result: List[Tuple[str, str, int, float]] = []
        for gene, location in gene_regions:
            with pysam.AlignmentFile(sam_path) as sam:
                prefix = _chr_prefix(location.chr, sam.header["SQ"])
                region = location.samtools(pad_left=1, prefix=prefix)
                cov: dict = collections.defaultdict(
                    lambda: collections.defaultdict(int)
                )
                log.info("Generating profile for {} ({})", gene, region)
                try:
                    for read in sam.fetch(region=region):
                        start, s_start = read.reference_start, 0
                        if not read.cigartuples:
                            continue

                        for op, size in read.cigartuples:
                            if op == 2:
                                for i in range(size):
                                    cov[location.chr][start + i] += 1
                                start += size
                            elif op == 1:
                                s_start += size
                            elif op == 4:
                                s_start += size
                            elif op in [0, 7, 8]:
                                for i in range(size):
                                    cov[location.chr][start + i] += 1
                                start += size
                                s_start += size
                    result += [
                        (gene, c, p, cov[c][p] * (factor / 2.0))
                        for c in sorted(cov.keys())
                        for p in sorted(cov[c].keys())
                        if location.start - 500 <= p <= location.end + 500
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
