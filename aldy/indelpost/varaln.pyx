#cython: embedsignature=True
#cython: profile=False

cimport cython
import random
import numpy as np
from functools import partial

from .pileup import (
    update_pileup,
    retarget,
    update_read_info,
    check_overhangs,
    filter_spurious_overhangs,
)
from .gappedaln import find_by_normalization
from .softclip import find_by_softclip_split
from .localn import find_by_smith_waterman_realn, make_aligner

from .alleles import phase_nearby_variants
from .contig import compare_contigs
from .utilities import (
    get_local_reference,
    relative_aln_pos,
    split_cigar,
    most_common_gap_ptrn,
    get_gap_ptrn2,
    most_common
)
from .pileup cimport make_pileup
from .utilities cimport split
from .contig cimport Contig, FailedContig

from .variant cimport Variant, NullVariant
from .local_reference cimport UnsplicedLocalReference

from pysam.libcalignmentfile cimport AlignmentFile

random.seed(123)

cdef class VariantAlignment:
    """This class accepts the target indel as :class:`~indelpost.Variant`
    and the BAM file as `pysam.AlignmentFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile>`__ .
    Upon creation, the target indel is searched through realignment.
    Evaluation of the equality between :class:`~indelpost.VariantAlignment` objects triggers local-phasing around the indels
    to test the alignment identity.


    Parameters
    ----------
    variant : Variant
        :class:`~indelpost.Variant` object representing the target indel.

    bam : pysam.AlignmentFile
        BAM file supplied as
        `pysam.AlignmentFile <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile>`__ object.

    window : integer
        analyzes region input_indel_pos +/- window (defalt 50).

    exclude_duplicates : bool
        True (default) to exclude reads marked duplicate by
        `MarkDuplicate <https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates>`__.

    downsample_threshold : integer
        downsamples to the threshold if the covearge at the input locus exceeds the threshold (default to 1000).

    mapping_quality_threshold : integer
        reads with a mapping quality (MAPQ) < the threshold (default 1) will be filtered.

    base_quality_threshold : integer
        non-reference base-calls with a quality score < the threshold (default 20) will be masked as "N" in the analysis.

    retarget_search_window : integer
        attempts to re-target in input_indel_pos +/- window (defalt 30) if the exact match is not found after realignment.

    retarget_similarity_cutoff : float
        re-targets to the indel sequence with a highest
        `Ratcliff/Obershelp similarity score <https://docs.python.org/3/library/difflib.html#difflib.get_close_matches>`__ > cutoff (default 0.7).

    exact_match_for_shiftable : bool
        True to require exact (equivalent)  match for indels that are equivalently alignable by shifting the position (default True)

    match_score : integer
        score (default 3) for matched bases in Smith-Waterman local alignment.

    mismatch_penalty : interger
        penalty (default 2) for mismatched bases in Smith-Waterman local alignment.

    gap_open_penalty : integer
        penalty (default 3) to create a gap in Smith-Waterman local alignment.

    gap_extension_penalty : integer
        penalty (default 1) to expend a created gap by one base in Smith-Waterman local alignment.

    auto_adjust_extension_penalty : bool
        True (default) to auto-adjust the gap open and extend penalties to find the input indel by realignment.

    no_realignment : bool
        True to only analyzed gap-aligned indels (default False)
    """
    def __cinit__(
        self,
        Variant target,
        AlignmentFile bam,
        int window=50,
        bint exclude_duplicates=True,
        int retarget_search_window=30,
        float retarget_similarity_cutoff=0.7,
        bint exact_match_for_shiftable=True,
        int mapping_quality_threshold=1,
        int downsample_threshold=1000,
        int base_quality_threshold=20,
        int match_score=3,
        int mismatch_penalty=2,
        int gap_open_penalty=3,
        int gap_extension_penalty=1,
        bint auto_adjust_extension_penalty=True,
        bint no_realignment=False,
    ):

        self.target, second_target = target, target

        is_complex_input = False
        if not target.is_non_complex_indel() and target.is_indel:
            is_complex_input = True

            if auto_adjust_extension_penalty:
                decomposed_variants = target.decompose_complex_variant(match_score, mismatch_penalty)
            else:
                decomposed_variants = target.decompose_complex_variant(match_score, mismatch_penalty, gap_open_penalty, gap_extension_penalty)

            decomposed_indels = [i for i in decomposed_variants if i.is_indel]
            decomposed_indels.sort(key=lambda x : len(x.indel_seq))
            self.__target = decomposed_indels[-1]
            self.target = self.__target
            if len(decomposed_indels) > 1:
                second_target = decomposed_indels[-2]

            #self.__target = max(decomposed_indels, key=lambda l : len(l.indel_seq))
            #self.target = self.__target
        else:
            self.__target = target.normalize()

        self.bam = bam
        self.window = window
        self.exclude_duplicates = exclude_duplicates
        self.retarget_window = retarget_search_window
        self.retarget_cutoff = retarget_similarity_cutoff
        self.exact_match_for_shiftable = exact_match_for_shiftable
        self.mapqthresh = mapping_quality_threshold
        self.downsamplethresh = downsample_threshold
        self.basequalthresh = base_quality_threshold
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_open_penalty = gap_open_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.auto_adjust_extension_penalty = auto_adjust_extension_penalty
        self.no_realignment = no_realignment
        self.is_complex_input = is_complex_input
        self.second_target = second_target
        self.unspliced_local_reference = UnsplicedLocalReference(
                                            self.__target.chrom,
                                            self.__target.pos,
                                            self.__target.reference.get_reference_length(self.__target.chrom),
                                            self.window,
                                            self.__target.reference
                                         )
        self.__pileup, self.contig = self.__parse_pileup()

    cdef __parse_pileup(self, Contig contig=None, bint retargeted=False, bint skip_read_end_check=False):
        """Dictize reads for target indel and make target indel template by consensus

        bam (pysam.AlignmentFile)
        target_indel (Variant): indel of interest
        smith_waterman (bool): True to perform SW realn
        template (bool): None if no indel template supplied
        """

        cdef list pileup

        read_end_evidence_only = False

        # equivalence search
        if retargeted:
            pileup = self.__pileup
            sample_factor = self.__sample_factor
        else:
            pileup, self.__sample_factor = make_pileup(
                self.__target,
                self.bam,
                self.unspliced_local_reference,
                exclude_duplicates=self.exclude_duplicates,
                window=self.window,
                downsamplethresh=self.downsamplethresh,
                basequalthresh=self.basequalthresh,
            )

            (
            self.__target,
            pileup,
            exptension_penalty_used,
            self._observed_pos,
            read_end_evidence_only
            ) = find_by_normalization(
                self.__target,
                pileup,
                self.window,
                self.match_score,
                self.mismatch_penalty,
                self.gap_open_penalty,
                self.gap_extension_penalty,
                self.basequalthresh,
            )

            if skip_read_end_check:
                read_end_evidence_only = False

            if self.target != self.__target:
                    self.__target, pileup = update_pileup(
                        pileup,
                        self.__target,
                        self.window,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                        self.basequalthresh,
                        bypass_search=True
                    )

            contig = Contig(
                self.__target,
                preprocess_for_contig_construction(
                    self.__target,
                    self.target,
                    pileup,
                    self.unspliced_local_reference,
                    self.window,
                    self.match_score,
                    self.mismatch_penalty,
                    self.gap_open_penalty,
                    exptension_penalty_used,
                ),
                self.unspliced_local_reference,
                self.basequalthresh,
                self.mapqthresh,
            )

            self.is_spurious_overhang = False
            if contig.failed and not self.no_realignment:
                within = self.retarget_window

                grid = generate_grid(self.auto_adjust_extension_penalty,
                                     self.gap_open_penalty,
                                     self.gap_extension_penalty,
                                     self.__target,
                       )

                ans = check_overhangs(pileup)

                if ans:
                    intron, overhangs = ans[0], ans[1]
                    non_spurious_overhangs = filter_spurious_overhangs(
                        self.__target,
                        intron,
                        overhangs,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                    )
                    if not non_spurious_overhangs:
                        contig = Contig(self.__target, [], self.unspliced_local_reference, self.basequalthresh, self.mapqthresh)
                        self.is_spurious_overhang = True
                        return pileup, contig
                    else:
                        res = grid_search(
                            self.__target,
                            non_spurious_overhangs,
                            self.window,
                            self.mapqthresh,
                            within,
                            self.retarget_cutoff,
                            self.match_score,
                            self.mismatch_penalty,
                            grid,
                            self.unspliced_local_reference,
                            self.exact_match_for_shiftable,
                        )

                        if res:
                            self.gap_open_penalty, self.gap_extension_penalty = res[2], res[3]
                        else:
                            contig = Contig(self.__target, [], self.unspliced_local_reference, self.basequalthresh, self.mapqthresh)
                            self.is_spurious_overhang = True
                            return pileup, contig
                else:
                    res = grid_search(
                        self.__target,
                        pileup,
                        self.window,
                        self.mapqthresh,
                        within,
                        self.retarget_cutoff,
                        self.match_score,
                        self.mismatch_penalty,
                        grid,
                        self.unspliced_local_reference,
                        self.exact_match_for_shiftable,
                    )

                    if res:
                        self.gap_open_penalty, self.gap_extension_penalty = res[2], res[3]

                # if retargeted successfully -> make contig based on the retarget
                if res:
                    self.__target, retarget_reads = res[0], res[1]

                    self.__target, self.__pileup = update_pileup(
                        pileup,
                        self.__target,
                        self.window,
                        self.match_score,
                        self.mismatch_penalty,
                        self.gap_open_penalty,
                        self.gap_extension_penalty,
                        self.basequalthresh,
                        bypass_search=True,
                    )

                    contig = Contig(
                        self.__target,
                        preprocess_for_contig_construction(
                            self.__target,
                            self.__target,
                            self.__pileup,
                            self.unspliced_local_reference,
                            self.window,
                            self.match_score,
                            self.mismatch_penalty,
                            self.gap_open_penalty,
                            self.gap_extension_penalty,
                        ),
                        self.unspliced_local_reference,
                        self.basequalthresh,
                        self.mapqthresh
                    )

                    # 2nd-pass using the retarget
                    return self.__parse_pileup(contig=contig, retargeted=True)

                # no target in this pileup
                else:
                    if self.is_complex_input:
                        try:
                            self.__target = self.second_target
                            self.target = self.second_target
                            self.is_complex_input = False
                            return self.__parse_pileup(contig=None, retargeted=False, skip_read_end_check=True)
                        except:
                             pileup, contig
                    else:
                        return pileup, contig

        #_target = [read for read in pileup if read["is_target"]]
        #_nontarget = [read for read in pileup if not read["is_target"]]
        #_target = sorted(_target, key=partial(centrality, target_pos=self.__target.pos))
        #_nontarget = [read for read in pileup if not read["is_target"]]
        #_nontarget = sorted(_nontarget, key=partial(centrality, target_pos=self.__target.pos))

        #print(_target[0])
        #print(_nontarget[0])

        # soft-clip realn & SW realn
        if contig.qc_passed and not self.no_realignment:

            orig_contig = contig

            # realign reads that are not through retarget path
            if not retargeted:
               cutoff = 1.0
               within = 30

               target = [read for read in pileup if read["is_target"]]
               nontarget = [read for read in pileup if not read["is_target"]]

               grid = generate_grid(self.auto_adjust_extension_penalty,
                                    self.gap_open_penalty,
                                    self.gap_extension_penalty,
                                    self.__target,
                      )

               res = grid_search(
                       self.__target,
                       nontarget,
                       self.window,
                       self.mapqthresh,
                       within,
                       cutoff,
                       self.match_score,
                       self.mismatch_penalty,
                       grid,
                       self.unspliced_local_reference,
                       self.exact_match_for_shiftable,
                   )

               if res:
                   nontarget = [read for read in nontarget if read not in res[1]]

                   pileup = target + res[1] + nontarget
                   self.gap_open_penalty, self.gap_extension_penalty = res[2], res[3]

                   self.__target, pileup = update_pileup(
                       pileup,
                       self.__target,
                       self.window,
                       self.match_score,
                       self.mismatch_penalty,
                       self.gap_open_penalty,
                       self.gap_extension_penalty,
                       self.basequalthresh,
                       bypass_search=True,
                   )

                   #equivalent but different position
                   if self.__target == res[0]:
                       self.__target = res[0]

               else:
                  pileup = target + nontarget

            if self.__target.count_repeats() == 0:
                pileup = find_by_softclip_split(self.__target, contig, pileup)


            if read_end_evidence_only:
                target_pileup = [read for read in pileup if read["is_target"]]

            pileup = find_by_smith_waterman_realn(
                self.__target,
                contig,
                pileup,
                self.match_score,
                self.mismatch_penalty,
                self.gap_open_penalty,
                self.gap_extension_penalty,
                self.basequalthresh
            )

            if read_end_evidence_only:
                newly_identified = [read for read in pileup if read["is_target"] and not read in target_pileup]
                if newly_identified:
                    expected_ptrn = most_common_gap_ptrn(newly_identified)

                    indels = []
                    contig_seq = contig.get_contig_seq()
                    alinger = make_aligner(contig_seq, self.match_score, self.mismatch_penalty)
                    for new_one in newly_identified:
                        if "N" not in new_one["cigar_string"] and is_perfect_match(alinger, contig_seq, new_one["read_seq"]):
                            indels += [i[-1] for i in new_one["I"]] + [d[-1] for d in new_one["D"]]

                    if indels:
                        try:
                            self.__target = most_common(indels)
                        except:
                            target_pos = self.__target.pos
                            indels.sort(key=lambda x : abs(x.pos - target_pos))

                        return self.__parse_pileup(contig=None, retargeted=False, skip_read_end_check=True)

            contig = Contig(
                self.__target,
                preprocess_for_contig_construction(
                    self.__target,
                    self.target,
                    pileup,
                    self.unspliced_local_reference,
                    self.window,
                    self.match_score,
                    self.mismatch_penalty,
                    self.gap_open_penalty,
                    self.gap_extension_penalty,
                ),
                self.unspliced_local_reference,
                self.basequalthresh,
                self.mapqthresh,
            )

            contig = compare_contigs(orig_contig, contig, self.__target.pos)

        return pileup, contig


    def __eq__(self, other):
        # alignment equivalence
        my_contig, other_contig = self.contig, other.contig

        if my_contig.failed or other_contig.failed:
            return False

        # check eq in phased form
        phasing_mode = "local"
        my_phased, other_phased = self.phase(how=phasing_mode), other.phase(how=phasing_mode)

        return (my_phased == other_phased)


    def __hash__(self):
        hashable = self.phase(how="local")
        return hash(hashable)


    def get_contig(self):
        """returns :class:`~indelpost.Contig` object built from reads supporting the target indel.
        :class:`~indelpost.FailedContig` is returned if contig assembly is not successful.
        """
        contig = self.contig
        if contig and not contig.failed:
            return contig
        else:
            failed = FailedContig()
            alt_cnt = self.count_alleles()[1]
            if alt_cnt:

                dirty_target_pileup = [read["is_dirty"] for read in self.__pileup if read["is_target"]]
                if sum(dirty_target_pileup) == len(dirty_target_pileup):
                    failed.is_low_quality = True
                else:
                    failed.failed_anyway = True
            else:
                failed.target_not_found = True

            return failed


    def get_target_indel(self):
        """returns :class:`~indelpost.Variant` object representing the actual target indel analyzed.
        The target indel may be different from the input indel due to the internal realignment.
        :class:`~indelpost.NullVariant` is returned when no target is found.
        """
        alt_cnt = self.count_alleles()[1]
        if alt_cnt:
            return self.__target
        else:
            return NullVariant(self.__target.chrom, self.__target.pos, self.__target.reference)


    def fetch_reads(self, how='target'):
        """returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__
        of `pysam.AlignedSegment <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment>`__ objects,
        depending on the strategy specified.

        Parameters
        ----------
            how : string
                specifies the read fetch strategy. "target" (default) to retrieve reads supporting the target.
                "non_target" for reads w/o target. "covering" for reads covering the locus w/ and w/o the target.
        """

        if how == "target":
            return [read["read"] for read in self.__pileup if read["is_target"]]
        elif how == "non_target":
            pos, indel_len = self._observed_pos, len(self.target.indel_seq)
            r_pos = max(v.pos for v in self.target.generate_equivalents())
            margin = r_pos - pos
            del_len = indel_len if self.target.is_del else 0
            targets = [read["read_name"] for read in self.__pileup if read["is_target"]]
            return [read["read"] for read in self.__pileup if count_as_non_target(read, pos, del_len, margin) and not read["read_name"] in targets]
        elif how == "covering":
            return [read["read"] for read in self.__pileup if read["is_covering"]]
        else:
            raise Exception("fetch stragety must be either of target, non_target, covering")


    def count_alleles(
        self, fwrv=False, by_fragment=False, three_class=False, estimated_count=False, quality_window=None, quality_threshold=None
    ):
        """returns a `tuple <https://docs.python.org/3/library/stdtypes.html#tuple>`__ of
        read counts:  (#non_target reads, #target reads).

        Parameters
        ----------
            fwrv : bool
                breaks down to the forward and reverse read counts.
                ( (non_target-fw, non_target-rv), (target-fw, target-rv) ).
            by_fragment : bool
                counts by fragment. Overlapping fw and rv reads are counted as one.
            three_class : bool
                breaks down to (reference, non_reference_non_target, target) counts.
            estimated_count : bool
                True to return estimated count when the coverage is higher than :attr:`~indelpost.VariantAlignment.downsample_threshold`.
            quality_window : integer
                specifies the range of base call quality filter. indel pos +/- quality_window.
            quality_threshold : integer
                filters reads with the median base call quality in indel pos +/- quality_window < quality_threshold.
        """

        cdef dict read

        pos, indel_len = self._observed_pos, len(self.target.indel_seq)

        r_pos = max(v.pos for v in self.__target.generate_equivalents())

        margin = r_pos - pos

        del_len = indel_len if self.target.is_del else 0

        reads = self.__pileup
        if quality_window and quality_threshold:
            reads = [
                read
                for read in reads
                if is_quality_read(read, pos, quality_window, quality_threshold)
            ]

        if three_class:
            for read in reads:
                read["is_locally_ref"] = is_locally_ref(read, pos)

        fw_target = [
            read["read_name"]
            for read in reads
            if read["is_target"] and not read["is_reverse"]
        ]
        rv_target = [
            read["read_name"]
            for read in reads
            if read["is_target"] and read["is_reverse"]
        ]

        fw_non_target = [
            read["read_name"]
            for read in reads
            if count_as_non_target(read, pos, del_len, margin) and not read["is_reverse"]
        ]
        rv_non_target = [
            read["read_name"]
            for read in reads
            if count_as_non_target(read, pos, del_len, margin) and read["is_reverse"]
        ]

        fw_target = set(fw_target)
        rw_target = set(rv_target)
        fwrv_target = fw_target | rw_target

        fw_non_target = set(fw_non_target) - fwrv_target
        rv_non_target = set(rv_non_target) - fwrv_target

        est = self.__sample_factor if estimated_count else 1

        if three_class:
            fw_ref = [
                read["read_name"]
                for read in reads if read["is_locally_ref"] and not read["is_reverse"] and read["read_name"] in fw_non_target
            ]

            fw_ref = set(fw_ref)
            fw_non_ref_non_target = fw_non_target - fw_ref

            rv_ref = [
                read["read_name"]
                for read in reads if read["is_locally_ref"] and read["is_reverse"] and read["read_name"] in rv_non_target
            ]
            rv_ref = set(rv_ref)
            rv_non_ref_non_target = rv_non_target - rv_ref

            if fwrv:
                return (
                    (
                        int(len(fw_ref) * est),
                        int(len(rv_ref) * est),
                    ),
                    (
                        int(len(fw_non_ref_non_target) * est),
                        int(len(rv_non_ref_non_target) * est),
                    ),
                    (
                        int(len(fw_target) * est),
                        int(len(rv_target) * est),
                    ),
                )

            else:
                if by_fragment:
                    fwrv_ref = len(fw_ref | rv_ref)
                    fwrv_non_ref_non_target = len(fw_non_ref_non_target | rv_non_ref_non_target)
                    fwrv_target = len(fw_target | rw_target)
                else:
                    fwrv_ref = len(fw_ref) + len(rv_ref)
                    fwrv_non_ref_non_target =  len(fw_non_ref_non_target) + len(rv_non_ref_non_target)
                    fwrv_target = len(fw_target) + len(rw_target)

                return (int(fwrv_ref), int(fwrv_non_ref_non_target), int(fwrv_target))

        if fwrv:
            return (
                (
                    int(len(fw_non_target) * est),
                    int(len(rv_non_target) * est),
                ),
                (
                    int(len(fw_target) * est),
                    int(len(rv_target) * est),
                ),
            )
        else:
            if by_fragment:
                fwrv_non_target = len(fw_non_target | rv_non_target)
                fwrv_target = len(fw_target | rw_target)
            else:
                fwrv_non_target = len(fw_non_target) + len(rv_non_target)
                fwrv_target = len(fw_target) + len(rw_target)

            return (
                int(fwrv_non_target * est),
                int(fwrv_target * est),
            )


    def phase(
        self,
        how="local",
        local_threshold=20,
        longest_common_substring_threshold=15,
        indel_repeat_threshold=None,
        mutation_density_threshold=0.05,
    ):
        """returns a :class:`~indelpost.Variant` object represeting a phased target indel.
        :class:`~indelpost.NullVariant` is returned if no target is found.
        Phasing is limited to the exon where the target indel is located for RNA-Seq.
        Refer to :meth:`~indelpost.Contig.get_phasables()` to retrieve phasable variants across exons.

        Parameters
        ----------
        how : string
            - "local" (default) phasing by recursive elimination of longest common substrings (LCS) and locality  scoring (see local_threshold).
            - "greedy" phase all phasable events.
            - "complex" phasing by the "local" option + exclusivity check. The exclusivity check removes phasable events that are also observed in non-target reads such as nearby homozygous germline events.
        local_threshold : integer
            local (and complex) phasing method checks if phasable SNVs are locally clustered around the indel.
            For i-th base apart from indel, the score gains one if mismatch and, if match, loses
            1*min(j/local_thresh, 1) where j = i for short indels (< 10-nt) and 0.6*i for longer indels.
        longest_common_substring_threshold : integer
            removes common substrings between the reference and contig assembled from target reads that are longer than longest_common_substring_threshold (default 15).
        indel_repeat_threshold : integer
            do not phase indels that repeat more than indel_repeat_threshold (default None)
        mutation_density_threshold : float
            do not phase if the pileup contains too many mutations (possibly error-prone dirty region).
            In non-target reads, #non-ref bases/#all bases > mutation_density_threshold (default 0.05)
        """
        if how == "complex":
            hard, to_complex = False, True
        elif how == "greedy":
            hard, to_complex = True, False
        elif how == "local":
            hard, to_complex = False, False
        else:
            raise Exception("phasing stragety must be either of local, greedy, complex")

        if indel_repeat_threshold is None:
            indel_repeat_threshold = np.inf

        return phase_nearby_variants(
            self.__target,
            self.contig,
            self.__pileup,
            self.basequalthresh,
            local_threshold,
            longest_common_substring_threshold,
            indel_repeat_threshold,
            mutation_density_threshold,
            hard,
            to_complex
        )


def is_quality_read(read, pos, qualitywindow, qualitythresh):

    try:
        lt_qual, rt_qual = read["lt_qual"], read["rt_qual"]
    except:
        lt_qual, rt_qual = split(
            read["read_qual"],
            read["cigar_string"],
            pos,
            read["read_start"],
            is_for_ref=False,
            reverse=False,
        )

    if lt_qual and rt_qual:
        lt_median = np.median(lt_qual[-min(len(lt_qual), qualitywindow) :])
        rt_median = np.median(rt_qual[: min(len(rt_qual), qualitywindow)])

        return lt_median > qualitythresh and rt_median > qualitythresh


def is_locally_ref(read, pos):
    if read["is_reference_seq"]:
        return True


    try:
        lt_seq, rt_seq = read["lt_seq"], read["rt_seq"]
    except:
        lt_seq, rt_seq = split(read["read_seq"],
            read["cigar_string"],
            pos,
            read["read_start"],
            is_for_ref=False,
            reverse=False,
        )

    try:
        lt_ref, rt_ref = read["lt_ref"], read["rt_ref"]
    except:
        lt_ref, rt_ref = split(read["ref_seq"],
            read["cigar_string"],
            pos,
            read["aln_start"],
            is_for_ref=True,
            reverse=False,
        )



    lt_seq_len = len(lt_seq)
    lt_ref_len = len(lt_ref)
    #clipped
    if not lt_ref_len:
        return False

    lt_len = min(5, lt_seq_len, lt_ref_len)

    rt_seq_len = len(rt_seq)
    rt_ref_len = len(rt_ref)
    if not rt_ref_len:
        return False

    rt_len = min(5, rt_seq_len, rt_ref_len)

    if lt_seq[-lt_len:] == lt_ref[-lt_len:] and rt_seq[:rt_len] == rt_ref[:rt_len]:
        return True
    else:
        return False


cdef bint count_as_non_target(dict read, int pos, int del_len, int margin):
    if read["is_target"]:
        return False

    cdef int aln_start = read["aln_start"]
    cdef int aln_end = read["aln_end"]

    # undetermined reads
    if read.get("undetermined", False):
        return False

    if read["is_covering"]:
        covering_subread = read["covering_subread"]
        if covering_subread[1] <= pos + margin:
            return False

        if pos < aln_start or aln_end < pos:
            return False
    else:
        if aln_end < pos:
            return False

        if del_len:
            if pos + del_len < aln_start:
                return False
        else:
            return False

    return True

def centrality(read, target_pos):
    relative_pos = relative_aln_pos(read["ref_seq"], read["cigar_list"], read["aln_start"], target_pos)
    return abs(0.5 - relative_pos)

cdef list preprocess_for_contig_construction(
    Variant target,
    Variant orig_target,
    list pileup,
    UnsplicedLocalReference unspl_loc_ref,
    int window,
    int match_score,
    int mismatch_penalty,
    int gap_open_penalty,
    int gap_extension_penalty,
):

    cdef dict read
    cdef int clips, nonclips

    if not pileup:
        return pileup

    targetpileup = [read for read in pileup if read["is_target"] and not read["is_dirty"]]

    if not targetpileup:
        return targetpileup

    clipped_targetpileup = [
        read for read in targetpileup if "S" in read["cigar_string"]
    ]
    nonclipped_targetpileup = [
        read for read in targetpileup if not "S" in read["cigar_string"]
        and (read.get("lt_cigar", None) and read.get("rt_cigar", None))
    ]

    clips = len(clipped_targetpileup)
    nonclips = len(nonclipped_targetpileup)

    if target == orig_target and nonclips > 9:
        random.seed(123)
        targetpileup = random.sample(nonclipped_targetpileup, 10)
        targetpileup = [right_aligner(read, target) for read in targetpileup]
    else:
        targetpileup = sorted(targetpileup, key=partial(centrality, target_pos=target.pos))

        unspl_ref_seq, unspl_lt_len = get_local_reference(orig_target, pileup, window, unspl_loc_ref, unspliced=True)
        unspl_aligner = make_aligner(unspl_ref_seq, match_score, mismatch_penalty)
        unspl_start = orig_target.pos + 1 - unspl_lt_len

        is_gapped_aln = False
        targetpileup = [
            update_spliced_read_info(
                read,
                target,
                orig_target,
                is_gapped_aln,
                window,
                match_score,
                mismatch_penalty,
                gap_open_penalty,
                gap_extension_penalty,
                unspl_loc_ref,
            )
            if "N" in read["cigar_string"] else
            update_read_info(
                read,
                target,
                is_gapped_aln,
                gap_open_penalty,
                gap_extension_penalty,
                unspl_aligner,
                unspl_ref_seq,
                unspl_start,
            )
            for read in targetpileup
        ]

        targetpileup = [read for read in targetpileup if read is not None and (read.get("lt_cigar", None) and read.get("rt_cigar", None))]

        _targetpileup = [
            read for read in targetpileup if read.get("cigar_updated", False)
        ]

        if _targetpileup:
            targetpileup = _targetpileup
        else:
            return targetpileup

    return targetpileup


def update_spliced_read_info(
    read,
    target,
    orig_target,
    is_gapped_aln,
    window,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    unspl_loc_ref,
):
    ref_seq, lt_len = get_local_reference(orig_target, [read], window, unspl_loc_ref)
    aligner = make_aligner(ref_seq, match_score, mismatch_penalty)
    ref_start = orig_target.pos + 1 - lt_len

    read = update_read_info(
        read,
        target,
        is_gapped_aln,
        gap_open_penalty,
        gap_extension_penalty,
        aligner,
        ref_seq,
        ref_start
    )

    return right_aligner(read, target)

def right_aligner(read, target):
    """Right align indels around splice site"""

    if (
        "N" not in read["cigar_string"]
        or (
            "I" not in read["cigar_string"]
            and "D" not in read["cigar_string"]
        )
    ):
        return read

    cigar_lst = read["cigar_list"]

    query_pos = 0
    ref_pos = read["aln_start"]
    new_cigar = []
    prev_event = "A"
    skip_next = False
    for i, c in enumerate(cigar_lst):
        event, event_len = c[-1], int(c[:-1])

        if event_len < 0:
            return None

        query_move = 0 if event in ("D", "N", "H", "P") else event_len
        ref_move = 0 if event in ("I", "H", "P") else event_len

        if event in ("I", "D") and prev_event == "N":
            try:
                nxt_c = cigar_lst[i + 1]
                nxt_event, nxt_event_len = nxt_c[-1], int(nxt_c[:-1])
                if nxt_event != "M":
                    raise Exception
            except:
                return None

            chrom, reference = target.chrom, target.reference
            padding_base = reference.fetch(chrom, ref_pos - 2, ref_pos - 1)
            if event == "I":
                ins_seq = read["read_seq"][query_pos: query_pos + event_len]
                ref = padding_base
                alt = padding_base + ins_seq
            else:
                del_seq = reference.fetch(chrom, ref_pos - 1, ref_pos - 1 + event_len)
                ref = padding_base + del_seq
                alt = padding_base

            right_aligned_vars = Variant(
                                    chrom,
                                    ref_pos - 1,
                                    ref,
                                    alt,
                                    reference,
                                    skip_validation=True
                                ).generate_equivalents()

            diff = max(v.pos for v in right_aligned_vars) - ref_pos + 1
            if diff > 0:
                new_cigar += [str(diff) + "M", str(event_len) + event, str(nxt_event_len - diff) + "M"]
            else:
                return None

            ref_pos += query_move + nxt_event_len
            query_pos += ref_move + nxt_event_len
            skip_next = True

        else:
            if skip_next:
                skip_next = False
            else:
                query_pos += query_move
                ref_pos += ref_move
                new_cigar.append(c)

        prev_event = event

    read["cigar_list"] = new_cigar
    read["cigar_string"] = "".join(new_cigar)

    try:
        if target in right_aligned_vars:
            rt_aln_pos = target.pos + diff
            read["lt_cigar"], read["rt_cigar"] = split_cigar(read["cigar_string"], rt_aln_pos, read["read_start"])
            read["lt_flank"], read["rt_flank"] = split(
                                                      read["read_seq"],
                                                      read["cigar_string"],
                                                      rt_aln_pos,
                                                      read["read_start"],
                                                      is_for_ref=False,
                                                      reverse=False
                                                 )
            read["lt_qual"], read["rt_qual"] = split(
                                                    read["read_qual"],
                                                    read["cigar_string"],
                                                    rt_aln_pos,
                                                    read["read_start"],
                                                    is_for_ref=False,
                                                    reverse=False
                                               )
            read["lt_ref"], read["rt_ref"] = split(
                                                    read["ref_seq"],
                                                    read["cigar_string"],
                                                    rt_aln_pos, read["aln_start"],
                                                    is_for_ref=True,
                                                    reverse=False
                                             )
            read["target_right_shifted"] = rt_aln_pos

            indel_len = len(target.indel_seq)
            if target.is_ins:
                read["rt_flank"] = read["rt_flank"][indel_len :]
                read["rt_qual"] = read["rt_qual"][indel_len :]
            else:
                read["rt_ref"] = read["rt_ref"][indel_len :]
        else:
            read["lt_cigar"], read["rt_cigar"] = split_cigar(read["cigar_string"], target.pos, read["read_start"])
    except:
        pass

    return read


def generate_grid (auto_adjust_extension_penalty,
                   gap_open_penalty,
                   gap_extension_penalty,
                   target,
):
    if auto_adjust_extension_penalty:
        if (gap_open_penalty, gap_extension_penalty) != (3, 1):
            if len(target.indel_seq) < 20:
                grid = [(gap_open_penalty, gap_extension_penalty),
                            (3, 1), (3, 0), (5, 1), (5, 0), (4, 1), (4, 0)
                       ]
            else:
                grid = [(gap_open_penalty, gap_extension_penalty),
                            (3, 0), (3, 1), (5, 1), (5, 0), (4, 1), (4, 0)
                       ]
        else:
            if len(target.indel_seq) < 20:
                grid = [(3, 1), (3, 0), (5, 1), (5, 0), (4, 1), (4, 0)]
            else:
                grid = [(3, 0), (3, 1), (5, 1), (5, 0), (4, 1), (4, 0)]
    else:
        grid = [(gap_open_penalty, gap_extension_penalty)]

    return grid


def grid_search(
    target,
    pileup,
    window,
    mapq_thresh,
    within,
    retarget_cutoff,
    match_score,
    mismatch_penalty,
    grid,
    unspl_loc_ref,
    exact_match_for_shiftable,
):
    # grid = [(gap.open, gap.ext)]
    h = 0
    responses, scores, hs = [], [], []
    while h < len(grid):
        res = retarget(
            target,
            pileup,
            window,
            mapq_thresh,
            within,
            retarget_cutoff,
            match_score,
            mismatch_penalty,
            grid[h][0],
            grid[h][1],
            unspl_loc_ref,
            exact_match_for_shiftable,
        )

        if res:
            score = res[2]
            responses.append(res)
            #scores.append(score)
            hs.append(h)

            # exact match
            if score == 1.0:
                scores.append(score * len(res[1]))
                # cnt exact hit?
                #break
            else:
                scores.append(score)
        h += 1

    if responses:
        idx = scores.index(max(scores))
        best_res = responses[idx]
        best_params = grid[hs[idx]]

        is_gapped_aln=False # to be removed

        candidate = best_res[0]

        gap_open_penalty, gap_extension_penalty = best_params[0], best_params[1]

        updated_reads = []
        for read, aligner, ref_seq, ref_start in zip(best_res[1], best_res[5], best_res[3], best_res[4]):

            updated_read = update_read_info(
                                read,
                                candidate,
                                is_gapped_aln,
                                gap_open_penalty,
                                gap_extension_penalty,
                                aligner,
                                ref_seq,
                                ref_start,
                            )

            updated_reads.append(updated_read)


        return candidate, updated_reads, gap_open_penalty, gap_extension_penalty
    else:
        return None


def is_perfect_match(aligner, contig_seq, read_seq):
    aligner.setRead(read_seq)
    _aln = aligner.align(gap_open=len(read_seq), gap_extension=len(read_seq))
    _contig = contig_seq[_aln.reference_start : _aln.reference_end]
    _read = read_seq[_aln.read_start : _aln.read_end]

    return _contig == _read


