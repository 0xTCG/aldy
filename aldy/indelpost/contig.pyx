#cython: embedsignature=True
#cython: profile=False

import random
import numpy as np
from collections import OrderedDict, namedtuple

from .utilities import *

from .variant cimport Variant
from .local_reference cimport UnsplicedLocalReference

from .consensus import make_consensus


random.seed(123)


cdef class Contig:
    """This class represents a consensus contig assembled from a subset (not all) of reads supporting the target indel.
    """
    def __cinit__(self, Variant target, list pileup, UnsplicedLocalReference unspl_loc_ref, int basequalthresh, int mapqthresh, double low_consensus_thresh=0.7, int donwsample_lim=100):
        self.target = target
        self.pileup = pileup

        self.targetpileup = self.__preprocess(mapqthresh, donwsample_lim)

        if self.targetpileup:
            consensus = make_consensus(self.target, self.targetpileup, basequalthresh)
            if consensus:

                self.splice_pattern = get_local_reference(self.target, consensus[2], 50, unspl_loc_ref, unspliced=False, splice_pattern_only=True)

                rt_aln_consensus = False
                rt_aligned_indel_seq = consensus[3]
                if rt_aligned_indel_seq and len(rt_aligned_indel_seq) == len(consensus[2]):
                    if len(set(rt_aligned_indel_seq)) == 1:
                        rt_aln_consensus = rt_aligned_indel_seq[0]

                self.__make_contig(consensus[0], consensus[1], rt_aln_consensus, basequalthresh)

                self.failed = False
            else:
                self.qc_passed = False
                self.failed = True
        else:
            self.qc_passed = False
            self.failed = True


    def __preprocess(self, mapqthresh, donwsample_lim):
        targetpileup = [read for read in self.pileup if read is not None and read["is_target"]]

        self.mapq = 0

        #self.splice_pattern = get_local_reference(self.target, targetpileup, window=50, unspliced=False, splice_pattern_only=True)

        #self.is_target_right_aligned = sum(read.get("target_right_aligned", 0) for read in targetpileup)

        if not targetpileup:
            return targetpileup

        if len(targetpileup) > donwsample_lim:
            targetpileup = random.sample(targetpileup, donwsample_lim)

        self.mapq = np.percentile([read["mapq"] for read in targetpileup], 50)
        self.low_qual_mapping_rate = sum(read["mapq"] < mapqthresh for read in targetpileup) / len(targetpileup)

        return targetpileup


    def __make_contig(self, lt_consensus, rt_consensus, rt_aln_consensus, basequalthresh):
        self.__index_by_genome_coord(lt_consensus[0], rt_consensus[0])

        self.lt_reference_seq = ""
        self.lt_target_block_reference_seq = ""
        self.lt_consensus_seq = ""
        self.lt_target_block_consensus_seq = ""
        self.lt_consensus_scores = []
        self.lt_target_block_consensus_scores = []

        self.indel_seq = ""

        self.rt_reference_seq = ""
        self.rt_target_block_reference_seq = ""
        self.rt_consensus_seq = ""
        self.rt_target_block_consensus_seq = ""
        self.rt_consensus_scores = []
        self.rt_target_block_consensus_scores = []

        exon_start, exon_end = -np.inf, np.inf
        if self.splice_pattern:
            for exon in self.splice_pattern:
                if exon[0] <= self.target.pos <= exon[1]:
                    exon_start, exon_end = exon[0], exon[1]

        for k, v in self.contig_dict.items():
            #if k < self.target.pos:
            if k < self.lt_end_pos:
                self.lt_reference_seq += v[0]
                self.lt_consensus_seq += v[1]
                self.lt_consensus_scores.extend([v[2]] * len(v[1]))
                if exon_start <= k:
                    self.lt_target_block_reference_seq += v[0]
                    self.lt_target_block_consensus_seq += v[1]
                    self.lt_target_block_consensus_scores.extend([v[2]] * len(v[1]))

            #elif k == self.target.pos:
            elif k == self.lt_end_pos:
                self.lt_reference_seq += v[0][0]
                self.lt_target_block_reference_seq += v[0][0]

                self.lt_consensus_seq += v[1][0]
                self.lt_target_block_consensus_seq += v[1][0]

                self.lt_consensus_scores.append(v[2])
                self.lt_target_block_consensus_scores.extend([v[2]])

                if rt_aln_consensus:
                    self.indel_seq = rt_aln_consensus
                else:
                    self.indel_seq = self.target.indel_seq

            #elif k > self.target.pos:
            elif k > self.lt_end_pos:
                self.rt_reference_seq += v[0]
                self.rt_consensus_seq += v[1]
                self.rt_consensus_scores.extend([v[2]] * len(v[1]))
                if k <= exon_end:
                    self.rt_target_block_reference_seq += v[0]
                    self.rt_target_block_consensus_seq += v[1]
                    self.rt_target_block_consensus_scores.extend([v[2]] * len(v[1]))

        self.start = lt_consensus[1]
        self.end = rt_consensus[1]

        self.__profile_non_target_variants()

        self.qc_passed = self.__qc()


    def __index_by_genome_coord(self, lt_index, rt_index):
        self.lt_genomic_index = lt_index
        self.rt_genomic_index = rt_index

        lt_end_pos = next(iter(lt_index))

        #hotfix for right-aln-case
        self.lt_end_pos = lt_end_pos

        # the target may be of low quality ("N")
        if "N" in rt_index[lt_end_pos][1]:
            rt_index[lt_end_pos] = (
                                     rt_index[lt_end_pos][0],
                                     self.target.alt,
                                     rt_index[lt_end_pos][2],
                                     rt_index[lt_end_pos][3]
            )

        genome_indexed_contig = lt_index
        genome_indexed_contig.update(rt_index)
        self.contig_dict = OrderedDict(sorted(genome_indexed_contig.items()))

        target_pos_alleles = self.contig_dict[lt_end_pos]
        ref, alt = target_pos_alleles[0], target_pos_alleles[1]
        if len(ref) < len(alt):
            the_shorter, the_longer = ref, alt
        else:
            the_shorter, the_longer = alt, ref

        self.is_non_complex_at_target_pos = the_longer[: len(the_shorter)] == the_shorter
        self.target_ref = ref[1 :]
        self.target_alt = alt[1 :]


    def __profile_non_target_variants(self):
        non_target_variants = [
            Variant(self.target.chrom, k, v[0], v[1], self.target.reference, skip_validation=True)
            for k, v in self.contig_dict.items()
            if v[0] and v[0] != v[1] and k != self.target.pos
        ]
        self.non_target_indels = [var for var in non_target_variants if var.is_indel]
        self.mismatches = [var for var in non_target_variants if not var.is_indel]

        self.gaps = [
            str(len(var.indel_seq)) + var.variant_type for var in self.non_target_indels
        ]
        self.gaps.append(str(len(self.target.indel_seq)) + self.target.variant_type)


    def __qc(self):

        lt_n, lt_len = self.lt_consensus_seq.count("N"), len(self.lt_consensus_seq)
        rt_n, rt_len = self.rt_consensus_seq.count("N"), len(self.rt_consensus_seq)

        qc_stats ={}

        qc_stats["low_qual_base_frac"] = low_qual_fraction(self.targetpileup)

        qc_stats["clip_rate"] = sum(True for k, v in self.contig_dict.items() if not v[0]) / len(self.contig_dict)

        lt_n_proportion = lt_n / lt_len
        rt_n_proportion = rt_n / rt_len
        qc_stats["n_rate"] = (lt_n + rt_n) / (lt_len + rt_len)

        low_consensus_rate_lt = (
            sum(score < self.low_consensus_thresh for score in self.lt_consensus_scores) / lt_len
        )
        low_consensus_rate_rt = (
            sum(score < self.low_consensus_thresh for score in self.rt_consensus_scores) / rt_len
        )

        qc_stats["low_consensus_rate"] = (low_consensus_rate_lt * lt_len + low_consensus_rate_rt * rt_len) / (lt_len + rt_len)

        self.qc_stats = qc_stats
        if qc_stats["low_qual_base_frac"] > 0.2:
            return False
        #elif qc_stats["clip_rate"] > 0.1: #abolish clip rate filter
        #    return False
        elif qc_stats["n_rate"] > 0.1:
            return False
        elif low_consensus_rate_lt > 0.2 or low_consensus_rate_rt > 0.2:
            return False
        else:
            return True


    def _get_splice_patterns(self):
        spls = self.splice_pattern
        if spls:
            intervals = []
            i, last_idx = 0, len(spls) - 1
            while i < last_idx:
                start = spls[i][1] + 1
                end = spls[i + 1][0] - 1
                intervals.append((start, end))
                i += 1

            return intervals


    def get_alignment(self):
        """shows the contig alignment as `namedtuple <https://docs.python.org/3/library/collections.html#collections.namedtuple>`__.
        ContigAlignment.chrom returns chromosome. ContigAlignment.aln returns an alignment dictionary as
        `OrderedDict <https://docs.python.org/3/library/collections.html#collections.OrderedDict>`__.
        The dictionary key is the position and the value is a `tuple <https://docs.python.org/3/library/stdtypes.html#tuple>`__ as (REF, ALT).
        ContigAlignment.spliced_intervals returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__
        of spliced intervals involved in the contig. `None <https://docs.python.org/3/c-api/none.html>`__ if the region of interest is not spliced.
        """
        ContigAlignment = namedtuple("ContigAlignment", "chrom aln spliced_intervals")

        data = []
        for k, v in self.contig_dict.items():
            if v[1] and v[0]:
                data.append((k, (v[0], v[1])))

        caln = ContigAlignment(chrom=self.target.chrom, aln=OrderedDict(data), spliced_intervals=self._get_splice_patterns())

        return caln


    def get_phasables(self):
        """returns a `list <https://docs.python.org/3/library/stdtypes.html#list>`__  of :class:`~indelpost.Variant` objects phasable with the target indel.
        """
        phasables = []
        chrom = self.target.chrom
        reference = self.target.reference
        for k, v in self.contig_dict.items():
            if v[1] and v[0] and v[1] != v[0]:
                phasables.append(Variant(chrom, k, v[0], v[1], reference, skip_validation=True))

        return phasables


    def get_reference_seq(self, split=False):
        """returns the reference contig sequence.

        Parameters
        ----------
            split : bool
                splits the reference at the target indel position.
        """
        if self.failed:
            return None

        if split:
            if self.is_non_complex_at_target_pos:
                if self.target.is_del:
                    return self.lt_reference_seq, self.indel_seq, self.rt_reference_seq
                else:
                    return  self.lt_reference_seq, "", self.rt_reference_seq
            else:
                return self.lt_reference_seq, self.target_ref, self.rt_reference_seq
        else:
            if self.target.is_non_complex_indel:
                refseq = (
                    self.lt_reference_seq + self.indel_seq + self.rt_reference_seq
                    if self.target.is_del
                    else self.lt_reference_seq + self.rt_reference_seq
                )
            else:
                return self.lt_reference_seq + self.target_ref + self.rt_reference_seq

            return refseq


    def get_contig_seq(self, split=False):
        """returns the contig sequence.

        Parameters
        ----------
            split : bool
                splits the contig sequence at the target indel position.
        """
        if self.failed:
            return None

        if split:
            if self.is_non_complex_at_target_pos:
                if self.target.is_ins:
                    return self.lt_consensus_seq, self.indel_seq, self.rt_consensus_seq
                else:
                    return self.lt_consensus_seq, "", self.rt_consensus_seq
            else:
                return self.lt_consensus_seq, self.target_alt, self.rt_consensus_seq
        else:
            if self.target.is_non_complex_indel:
                conseq = (
                    self.lt_consensus_seq + self.indel_seq + self.rt_consensus_seq
                    if self.target.is_ins
                    else self.lt_consensus_seq + self.rt_consensus_seq
                )
            else:
                conseq = self.lt_consensus_seq + self.target_alt + self.rt_consensus_seq

            return conseq


cdef class FailedContig:
    """This class is returned when contig assembly has failed. Boolean expression evaluates to `False <https://docs.python.org/3/library/stdtypes.html#boolean-values>`__.

    Parameters
    ----------
    target_not_found : bool
        True when contig was not constructed because the target indel was not found.

    is_low_quality : bool
        True when contig was not constructed because the target indel was only suppported by low quality reads.
        Reads are defined low quality if 15% or more bases have quality score < :attr:`~indelpost.VariantAlignment.base_quality_threshold`.

    failed_anyway : bool
        True when contig was not constructed due to the other reasons.
    """
    def __cinit__(self):
        self.target_not_found = False
        self.is_low_quality = False
        self.failed_anyway = False

    def __bool__(self):
        return False


def compare_contigs(orig_contig, new_contig, target_pos):
    if new_contig.failed:
        return orig_contig

    orig_len = len(orig_contig.get_reference_seq())
    orig_clip_rate = orig_contig.qc_stats["clip_rate"]

    new_len = len(new_contig.get_reference_seq())
    new_clip_rate = new_contig.qc_stats["clip_rate"]

    orig_score = contig_centerness_score(orig_contig, target_pos)
    new_score = contig_centerness_score(new_contig, target_pos)

    if new_clip_rate > 0.1:
        return orig_contig

    condition1 = (new_len <= orig_len)
    condition2 = (new_clip_rate > orig_clip_rate and new_clip_rate >=0.03)
    condition3 = (orig_score < new_score)

    if sum([condition1, condition2, condition3]) >= 2:
        return orig_contig
    else:
        return new_contig


def contig_centerness_score(contig, target_pos):
    lt_cnt, rt_cnt = 0, 0
    for k, v in contig.contig_dict.items():
        if v[0] and v[1]:
            if k <= target_pos:
                lt_cnt += 1
            else:
                rt_cnt += 1

    return 0.5 - min(lt_cnt, rt_cnt) /(lt_cnt + rt_cnt)

