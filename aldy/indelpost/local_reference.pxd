from pysam.libcfaidx cimport FastaFile


cdef class UnsplicedLocalReference:
    cdef:
        str chrom, ref, alt, unspliced_local_reference
        int pos, window, local_ref_start, ref_len
        FastaFile reference
        readonly int left_len
    
