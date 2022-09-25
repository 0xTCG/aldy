# cython: embedsignature=True
from .variant cimport Variant

cdef class FailedContig:
    cdef public bint target_not_found
    cdef public bint is_low_quality
    cdef public bint failed_anyway

cdef class Contig:
    cdef Variant target
    cdef list pileup, targetpileup
    cdef public bint qc_passed, failed
    cdef int donwsample_lim
    cdef object lt_genomic_index, rt_genomic_index
    cdef int start, end
    cdef str lt_reference_seq, rt_reference_seq
    cdef str lt_target_block_reference_seq, rt_target_block_reference_seq
    cdef bint is_non_complex_at_target_pos
    cdef str target_ref, target_alt
    cdef double low_consensus_thresh
    cdef public int is_target_right_aligned
    cdef public object contig_dict
    cdef public str lt_consensus_seq, rt_consensus_seq, indel_seq
    cdef public str lt_target_block_consensus_seq, rt_target_block_consensus_seq
    cdef public list lt_consensus_scores, rt_consensus_scores, mismatches, non_target_indels, gaps
    cdef public list lt_target_block_consensus_scores, rt_target_block_consensus_scores
    cdef public int mapq, lt_end_pos
    cdef public double low_qual_mapping_rate
    cdef public dict qc_stats
    cdef public tuple splice_pattern
