from .variant cimport Variant
from cpython cimport array
import array

cdef tuple split(object data, str cigarstring, int target_pos, int string_pos, bint is_for_ref, bint reverse)

cdef tuple locate_indels(str cigarstring, int aln_start_pos)

cdef list get_spliced_subreads(str cigarstring, int read_start_pos, int read_end_pos)

cdef int count_lowqual_non_ref_bases(str read_seq, str ref_seq, array.array quals, list cigar_list, int basequalthresh)
