# cython: embedsignature=True

#from pysam.libcbcf cimport VariantRecord, VariantFile
from pysam.libcalignmentfile cimport AlignmentFile
from .variant cimport Variant
from .contig cimport Contig
from .local_reference cimport UnsplicedLocalReference

cdef class VariantAlignment:
    cdef Variant target, __target, second_target
    cdef readonly AlignmentFile bam
    cdef int window, retarget_window, mapqthresh
    cdef int downsamplethresh, basequalthresh, match_score, _observed_pos
    cdef int mismatch_penalty, gap_open_penalty, gap_extension_penalty
    cdef float retarget_cutoff, __sample_factor,
    cdef bint exclude_duplicates, exact_match_for_shiftable, auto_adjust_extension_penalty, no_realignment
    cdef list __pileup
    cdef readonly is_spurious_overhang, is_complex_input
    cdef readonly Contig contig
    cdef readonly UnsplicedLocalReference unspliced_local_reference

    cdef __parse_pileup(self, Contig contig=*, bint retargeted=*, bint skip_read_end_check=*)
