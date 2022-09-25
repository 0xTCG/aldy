
#cython: boundscheck=False, wraparound=False
from typing import (
    NamedTuple,
    Union
)

from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cython.operator cimport postincrement as inc
from cpython.bytes import PyBytes_Check
from libc.stdint cimport int32_t, uint32_t, uint16_t, int8_t, uint8_t

"""
What is a CIGAR?
http://genome.sph.umich.edu/wiki/SAM#What_is_a_CIGAR.3F
"""

cdef int8_t* DNA_BASE_LUT = [
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
]

cdef inline void dnaToInt8(const char *c_str, int8_t *arr, int len):
    for i in range(0, len):
        arr[i] = DNA_BASE_LUT[c_str[i]]

cdef inline str _str(s):
    if isinstance(s, bytes):
        # encode to the specific encoding used inside of the module
        return (<bytes>s).decode('utf8')
    else:
        return s

cdef extern from "Python.h":
    cdef int PyBytes_AsStringAndSize(object, char **, Py_ssize_t *)
    cdef char* PyUnicode_AsUTF8AndSize(object, Py_ssize_t *)

cdef inline char* obj_to_cstr_len(object o1, Py_ssize_t *length):
    cdef char* c_str1
    if PyBytes_Check(o1):
        if PyBytes_AsStringAndSize(o1, &c_str1, length) == -1:
            raise TypeError("obj_to_cstr_len: PyBytes_AsStringAndSize error")
        return c_str1
    else:
        c_str1 = PyUnicode_AsUTF8AndSize(o1, length)
        if c_str1 == NULL:
            raise OSError("obj_to_cstr_len: PyUnicode_AsUTF8AndSize error")
    return c_str1

cdef extern from "ssw.h":
    # leave out a few members
    ctypedef struct s_profile:
        const int8_t* read
        const int8_t* mat
        int32_t readLen
        int32_t n
        uint8_t bias

    ctypedef struct s_align:
        uint16_t score1
        uint16_t score2
        int32_t ref_begin1
        int32_t ref_end1
        int32_t read_begin1
        int32_t read_end1
        int32_t ref_end2
        uint32_t* cigar
        int32_t cigarLen

    s_profile* ssw_init(const int8_t*, const int32_t, const int8_t*, const int32_t, const int8_t)
    void init_destroy (s_profile*)
    s_align* ssw_align (const s_profile*, const int8_t*, int32_t, const uint8_t, const uint8_t, const uint8_t, const uint16_t, const int32_t, const int32_t)
    void align_destroy (s_align*)
    char cigar_int_to_op (uint32_t)
    uint32_t cigar_int_to_len(uint32_t)
    uint32_t to_cigar_int(uint32_t, char)

Alignment = NamedTuple("Alignment", [
        ('CIGAR', str),
        ('optimal_score', int),
        ('sub_optimal_score', int),
        ('reference_start', int),
        ('reference_end', int),
        ('read_start', int),
        ('read_end', int)
    ]
)

STR_T = Union[str, bytes]


cdef class SSW:

    cdef int8_t* score_matrix
    cdef s_profile* profile

    cdef object read
    cdef int8_t* read_arr
    cdef Py_ssize_t read_length

    cdef object reference
    cdef int8_t* ref_arr
    cdef Py_ssize_t ref_length

    def __cinit__(self, int match_score=2, int mismatch_penalty=2):
        self.score_matrix = NULL
        self.profile = NULL
        self.read_arr = NULL
        self.ref_arr = NULL
    # end

    def __init__(self,  int match_score=2,
                        int mismatch_penalty=2):
        """ Requires a

        Args:
            match_score (int): for scoring matches
            mismatch_penalty (int): for scoring mismatches
        """
        self.score_matrix = <int8_t*> PyMem_Malloc(25*sizeof(int8_t))
        self.buildDNAScoreMatrix(   <uint8_t>match_score,
                                    <uint8_t> mismatch_penalty,
                                    self.score_matrix)
        self.read = None
        self.reference = None
    # end def

    def __dealloc__(self):
        PyMem_Free(self.score_matrix)

        if self.profile != NULL:
            init_destroy(self.profile)
            self.profile = NULL

        if self.read_arr != NULL:
            PyMem_Free(self.read_arr)

        if self.ref_arr != NULL:
            PyMem_Free(self.ref_arr)
    # end def

    def setRead(self, read: STR_T):
        """ Set the query read string

        Args:
            read:  String-like (str or bytestring) that represents the read.
                    Must be set
        """
        cdef Py_ssize_t read_length
        cdef const char* read_cstr = obj_to_cstr_len(read, &read_length)
        cdef int8_t* read_arr = <int8_t*> PyMem_Malloc(read_length*sizeof(char))

        if self.profile != NULL:
            init_destroy(self.profile)
            self.profile = NULL
        if self.read_arr != NULL:
            PyMem_Free(self.read_arr)
            self.read_arr = NULL

        self.read = read
        self.read_arr = read_arr
        self.read_length = read_length
        dnaToInt8(read_cstr, read_arr, read_length)

        self.profile = ssw_init(read_arr,
                                <int32_t> read_length,
                                self.score_matrix,
                                5,
                                2 # don't know best score size
                                )
    # end def

    def setReference(self, reference: STR_T):
        """Set the query reference string

        Args:
            reference:  String-like (str or bytestring) that represents the
                reference sequence must be set
        """
        cdef Py_ssize_t ref_length
        cdef const char* ref_cstr = obj_to_cstr_len(reference, &ref_length)
        cdef int8_t* ref_arr = <int8_t*> PyMem_Malloc(ref_length*sizeof(char))
        dnaToInt8(ref_cstr, ref_arr, ref_length)
        self.reference = reference
        if self.ref_arr != NULL:
            PyMem_Free(self.ref_arr)
            self.ref_arr = NULL
        self.ref_arr = ref_arr
        self.ref_length = ref_length
    # end def

    cdef s_align* align_c(self,
        int gap_open,
        int gap_extension,
        Py_ssize_t start_idx,
        int32_t mod_ref_length) except NULL:
        """C version of the alignment code
        """
        cdef Py_ssize_t read_length
        cdef const char* read_cstr = obj_to_cstr_len(self.read, &read_length)
        cdef s_align* result = NULL
        cdef int32_t mask_len = read_length // 2

        mask_len = 15 if mask_len < 15 else mask_len

        if self.profile != NULL:
            result = ssw_align ( self.profile,
                                &self.ref_arr[start_idx],
                                mod_ref_length,
                                gap_open,
                                gap_extension,
                                1, 0, 0, mask_len)
        else:
            raise ValueError("Must set profile first")
        if result == NULL:
            raise ValueError("Problem Running alignment, see stdout")
        return result
    # end def

    def align(self,
        int gap_open = 3,
        int gap_extension = 1,
        Py_ssize_t start_idx = 0,
        Py_ssize_t end_idx = 0) -> Alignment:
        '''Align a read to the reference with optional index offseting

        returns a dictionary no matter what as align_c can't return
        NULL

        Args:
            gap_open (int):         penalty for gap_open. default 3
            gap_extension (int):    penalty for gap_extension. default 1
            start_idx (Py_ssize_t): index to start search. default 0
            end_idx (Py_ssize_t):   index to end search (trying to avoid a target region).
                                    default 0 means use whole reference length

        Returns:
            Alignment with keys `CIGAR`,        <for depicting alignment>
                                `optimal_score`,
                                `sub-optimal_score`,
                                `reference_start`, <index into reference>
                                `reference_end`,   <index into reference>
                                `read_start`,  <index into read>
                                `read_end`     <index into read>

        Raises
            ValueError
        '''
        cdef Py_ssize_t c
        cdef char letter
        cdef int letter_int
        cdef uint32_t length
        cdef int32_t search_length
        cdef Py_ssize_t end_idx_final

        if start_idx < 0 or end_idx < 0:
            raise ValueError("negative indexing not supported")
        if end_idx > self.ref_length or start_idx > self.ref_length:
            err = "start_idx: {} or end_idx: {} can't be greater than ref_length: {}".format(
                                                start_idx,
                                                end_idx,
                                                self.ref_length)
            raise ValueError(err)
        if end_idx == 0:
            end_idx_final = self.ref_length
        else:
            end_idx_final = end_idx
        search_length = end_idx_final - start_idx

        cigar = None

        if self.reference is None:
            raise ValueError("call setReference first")

        cdef s_align* result =  self.align_c(gap_open, gap_extension, start_idx, search_length)
        if result.cigar != NULL:
            cigar = ""
            for c in range(result.cigarLen):
                letter = cigar_int_to_op(result.cigar[c])
                letter_int = letter
                length = cigar_int_to_len(result.cigar[c])
                cigar += "%d%s" % (<int>length, chr(letter_int))
        out = Alignment(
                cigar,
                result.score1,
                result.score2,
                result.ref_begin1,
                result.ref_end1,
                result.read_begin1,
                result.read_end1
        )
        #print("RAW BEGIN")
        #self.printResult_c(result)
        #print("RAW END")
        align_destroy(result)
        return out
    # end def

    cdef int buildDNAScoreMatrix(self,
        const uint8_t match_score,
        const uint8_t mismatch_penalty,
        int8_t* matrix) except -1:
        """
        mismatch_penalty should be positive

        The score matrix looks like
                            A,  C,  G,  T,  N
        score_matrix  = {   2, -2, -2, -2,  0, // A
                           -2,  2, -2, -2,  0, // C
                           -2, -2,  2, -2,  0, // G
                           -2, -2, -2,  2,  0, // T
                            0,  0,  0,  0,  0  // N
                        }
        """
        cdef Py_ssize_t i, j
        cdef Py_ssize_t idx = 0;
        for i in range(4):
            for j in range(4):
                if i == j:
                    matrix[idx] =  <int8_t> match_score
                else:
                    matrix[idx] = <int8_t> (-mismatch_penalty)
                inc(idx)
            matrix[idx] = 0;
            inc(idx)
        for i in range(5):
            matrix[inc(idx)] = 0
        return 0
    # end def
# end class

def force_align( read: STR_T,
                reference: STR_T,
                force_overhang: bool = False,
                aligner: SSW = None) -> Alignment:
    '''Enforces no gaps by raising the ``gap_open`` penalty

    Args:
        read:
        reference:
        force_overhang: Make sure only one end overhangs
        aligner: pass an existing :class:`SSW` object
    Raises:
        ValueError for no solution found
    '''
    a: SSW = SSW() if aligner is None else aligner
    a.setRead(read)
    a.setReference(reference)
    len_x: int = len(read)
    # set the gap_open penalty high to drop all splits in alignment
    # pick the first max hit
    res: Alignment = a.align(gap_open=len_x)
    if res.optimal_score < 4:
        raise ValueError("No solution found")
    if force_overhang:
        # read must align to either the beginning or end of the reference string
        if (res.reference_start != 0 or
            res.reference_end != len(reference) - 1):
            raise ValueError("Read does not align to one overhang")
    return res
# end def

def format_force_align(  read: STR_T,
                        reference: STR_T,
                        alignment: Alignment,
                        do_print: bool = False):
    '''Does not truncate strings

    Args:
        read:
        reference:
        alignment:
        do_print: default is False
    '''
    start_ref: int = alignment.reference_start
    start_read: int = alignment.read_start
    buffer_ref: str = ''
    buffer_read: str = ''
    if start_ref < start_read:
        buffer_ref = ' '*(start_read - start_ref)
    else:
        buffer_read = ' '*(start_ref - start_read)
    ref_out = buffer_ref + _str(reference)
    read_out = buffer_read + _str(read)
    if do_print:
        print(ref_out)
        print(read_out)
    return ref_out, read_out
# end def
