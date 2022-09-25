from .variant cimport Variant
from .local_reference cimport UnsplicedLocalReference
from pysam.libcalignmentfile cimport AlignmentFile

cdef tuple make_pileup(
    Variant target,
    AlignmentFile bam,
    UnsplicedLocalReference unspl_loc_ref,
    bint exclude_duplicates,
    int window,
    int downsamplethresh,
    int basequalthresh
)
