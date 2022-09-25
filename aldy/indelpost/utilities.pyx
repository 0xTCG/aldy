#cython: profile=False

import re
from cpython cimport array
import array

import numpy as np
from collections import namedtuple
from operator import mul
from functools import reduce
from pysam.libcbcf cimport VariantRecord, VariantRecordFilter, VariantFile

from .variant cimport Variant
from .local_reference cimport UnsplicedLocalReference

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


def most_common(lst):
    alst = list(set(lst))
    alst.sort()
    return max(alst, key=lst.count)


def get_gap_ptrn(read):
    return "".join([c for c in read["cigar_list"] if "D" in c or "I" in c])

def get_gap_ptrn2(read):
    ptrn = ""
    pos = read["aln_start"]
    for cigar in read["cigar_list"]:
        event, event_len = cigar[-1], int(cigar[:-1])
        if event in ("M", "X", "="):
            pos += event_len
        elif event in ["I", "D", "N"]:
            ptrn += "{}@{}".format(cigar, pos-1)
            if event == "D":
                pos += event_len
    return ptrn

def most_common_gap_pattern(targetpileup):
    ptrns = [get_gap_ptrn(read) for read in targetpileup]
    return most_common(ptrns)

def most_common_gap_ptrn(targetpileup):
    ptrns = [get_gap_ptrn2(read) for read in targetpileup]
    return most_common(ptrns)

cpdef list to_flat_list(list lst_of_lst):
    cdef list lst
    return [i for lst in lst_of_lst for i in lst]


cpdef list to_flat_vcf_records(VariantRecord record):

    cdef str alt

    VcfRec = namedtuple(
        "VcfRec", "chrom pos id ref alt qual filter info format samples orig"
    )

    if not record.alts:
        return []

    flat_record = [
        VcfRec(
            chrom=record.chrom,
            pos=record.pos,
            id=record.id,
            ref=record.ref,
            alt=alt,
            qual=record.qual,
            filter=record.filter,
            info=record.info,
            format=record.format,
            samples=record.samples,
            orig=record,
        )
        for alt in record.alts
    ]

    return flat_record


cpdef dict to_dict(object record):
    cdef str k

    d = {}
    for k, v in record.items():
        if isinstance(v, tuple):
            d[k] = ",".join([str(i) for i in v])
        else:
            d[k] = v

    if d:
        return d


cpdef bint match_indels(Variant query, Variant subject, str matchby, bint indel_only):
    if matchby != "normalization" and indel_only and not query.is_indel:
        return False

    if matchby == "normalization":
        return query == subject

    elif matchby == "locus":
        if query.chrom != subject.chrom:
            return False

        query.normalize(inplace=True)
        subject.normalize(inplace=True)

        return query.pos == subject.pos

    elif matchby == "exact":
        return (
            (query.chrom == subject.chrom)
            and (query.pos == subject.pos)
            and (query.ref == subject.ref)
            and (query.alt == subject.alt)
        )


cpdef double linguistic_complexity(str seq):
    cdef int i, j, n
    n = len(seq)
    if n <= 1:
        return float(n)
    else:
        usage = []
        for i in range(1, n):
            i_mer = [seq[j : j + i] for j in range(n - i + 1)]
            usage.append(len(set(i_mer)) / min(4 ** i, n - i + 1))

        return reduce(mul, usage)


cpdef double low_qual_fraction(list pileup):
    cdef dict read
    cdef int pileup_vol = 1
    cdef int low_qual_vol = 0

    for read in pileup:
        pileup_vol += len(read["read_seq"])
        low_qual_vol += read["low_qual_base_num"]

    return low_qual_vol / pileup_vol


def to_minimal_repeat_unit(seq):
    """Find repeat unit in indel sequence
    """
    mid = int(len(seq) / 2)
    min_unit = seq

    j = 1
    found = False

    while j <= mid and not found:
        tandems = [seq[i : i + j] for i in range(0, len(seq), j)]
        if len(set(tandems)) == 1:
            found = True
            min_unit = list(set(tandems))[0]
        j += 1

    return min_unit


def repeat_counter(query_seq, flank_seq):
    """
    """
    qlen, flen = len(query_seq), len(flank_seq)
    count = 0

    if flen < qlen:
        return count

    for i in range(0, flen, qlen):
        if flank_seq[i  : i + qlen] == query_seq:
            count += 1
        else:
            break

    return count


cdef int count_lowqual_non_ref_bases(
    str read_seq,
    str ref_seq,
    array.array quals,
    list cigar_list,
    int basequalthresh,
):
    cdef int i = 0, j = 0, k = 0, cnt = 0, event_len = 0
    cdef str event = "", cigar = "" ,

    for cigar in cigar_list:
        event, event_len = cigar[-1], int(cigar[:-1])

        if event in ("M", "=", "X"):
            while k < event_len:
                if read_seq[i] != ref_seq[j] and quals[i] < basequalthresh:
                    cnt += 1
                i += 1
                j += 1
                k += 1
            k = 0
        elif event in ("I", "S"):
            while k < event_len:
                if quals[i] < basequalthresh:
                    cnt += 1
                i += 1
                k += 1
            k = 0
        elif event == "D":
            j += event_len

    return cnt


cpdef list get_mapped_subreads(str cigarstring, int aln_start_pos, int aln_end_pos):

    cdef int event_len, current_pos
    cdef str cigar, event
    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef list res = []

    current_pos = aln_start_pos
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        if event in ("M", "X", "="):
            res.append((current_pos, (current_pos + event_len - 1)))
            current_pos += event_len
        elif event in ("I", "S", "H", "P"):
            pass
        else:
            current_pos += event_len

    return res


cdef list get_spliced_subreads(str cigarstring, int read_start_pos, int read_end_pos):

    cdef int i = 0
    cdef int event_len
    cdef str cigar, event, prev_event
    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef list pos_lst = [read_start_pos]
    cdef list res = []

    if not "N" in cigarstring:
        return [(read_start_pos, read_end_pos)]

    prev_event = "A"
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        if event == "N":
            pos_lst.append(read_start_pos - 1)
        elif prev_event == "N":
            pos_lst.append(read_start_pos)

        if event in ("I", "H", "P"):
            pass
        else:
            read_start_pos += event_len

        prev_event = event

    if prev_event != "N":
        pos_lst.append(read_end_pos)

    while i < len(pos_lst):
        res.append(pos_lst[i : i+2])
        i += 2

    return res


cpdef int get_end_pos(int read_start_pos, str lt_flank, str cigarstring):
    read_start_pos -= 1

    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef str cigar, event
    cdef int event_len, i = 0, flank_len = len(lt_flank)

    while flank_len > 0:
        cigar = cigar_lst[i]
        event, event_len = cigar[-1], int(cigar[:-1])

        if event == "D" or event == "N":
            read_start_pos += event_len
        elif event == "I":
            flank_len -= event_len
        elif event == "H" or event == "P":
            pass
        else:
            flank_len -= event_len
            read_start_pos += event_len

        i += 1

    return read_start_pos + flank_len


cdef tuple locate_indels(str cigarstring, int aln_start_pos):
    aln_start_pos -= 1

    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef str cigar, event
    cdef int event_len

    cdef list ins = []
    cdef list dels = []
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])
        if event == "I":
            ins.append((aln_start_pos, event_len))
        elif event == "D":
            dels.append((aln_start_pos, event_len))
            aln_start_pos += event_len
        elif event == "H" or event == "P":
            pass
        else:
            aln_start_pos += event_len

    return ins, dels


cpdef tuple split_cigar(str cigarstring, int target_pos, int start):

    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef str cigar, event
    cdef int event_len, move, diff
    cdef list lt_lst = []
    cdef list rt_lst = cigar_lst

    start -= 1
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        move = 0 if event in ("I", "H", "P") else event_len
        start += move
        rt_lst = rt_lst[1 :]

        if target_pos <= start:
            diff = start - target_pos
            lt_cigar = str(event_len - diff) + event
            lt_lst.append(lt_cigar)

            if diff:
                    rt_lst = [str(diff) + event] + rt_lst

            return lt_lst, rt_lst
        else:
            lt_lst.append(cigar)


def merge_consecutive_gaps(cigar_lst):

    merged_lst = []
    while cigar_lst:
        c = cigar_lst[0]
        cigar_lst = cigar_lst[1:]

        if "I" in c or "D" in c:
            i = 0
            is_gap = True
            while i < len(cigar_lst) and is_gap:
                tmp = cigar_lst[i]
                is_gap = True if "I" in tmp or "D" in tmp else False
                i += 1

            if i - 1:
                c += "".join(cigar_lst[: i - 1])
                cigar_lst = cigar_lst[i - 1 :]

        merged_lst.append(c)

    return merged_lst


def make_insertion_first(cigarstring):

    cigar_lst = cigar_ptrn.findall(cigarstring)

    merged_cigar_lst = merge_consecutive_gaps(cigar_lst)
    new_cigar = []
    for c in merged_cigar_lst:
        if "I" in c and "D" in c:
            c_lst = cigar_ptrn.findall(c)
            if "D" in c_lst[0]:
                swapped = c_lst[::-1]
                new_cigar.append("".join(swapped))
            else:
                new_cigar.append("".join(c_lst))
        else:
            new_cigar.append(c)

    return "".join(new_cigar)


def relative_aln_pos(ref_seq, cigar_lst, aln_start, target_pos, include_clip=False):

    current_pos = aln_start - 1
    ref_seq_pos = 0
    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        event = "M" if include_clip and event == "S" else event

        if event == "M" or event == "D":
            current_pos += event_len
            ref_seq_pos += event_len
        elif event in ("I", "H", "P"):
            pass
        else:
            current_pos += event_len

        if current_pos >= target_pos:
            break

    ref_seq_pos += (target_pos - current_pos)

    return ref_seq_pos / len(ref_seq)


cdef tuple split(
    object data,
    str cigarstring,
    int target_pos,
    int string_pos,
    bint is_for_ref,
    bint reverse,
):

    cdef list cigar_lst = cigar_ptrn.findall(cigarstring)
    cdef int _size = len(cigar_lst)

    cdef str cigar, event
    cdef int event_len, d_move, g_move

    cdef double [:] data_moves = np.zeros((_size,))
    cdef double [:] genome_moves = np.zeros((_size,))

    cdef int i = 0, j = 0

    for cigar in cigar_lst:
        event, event_len = cigar[-1], int(cigar[:-1])

        if event == "N":
            d_move = 0
            g_move = event_len
        elif event == "I":
            g_move = 0
            d_move = 0 if is_for_ref else event_len
        elif event == "D":
            g_move = event_len
            d_move = event_len if is_for_ref else 0
        elif event == "H" or event == "P":
            d_move = 0
            g_move = 0
        else:
            g_move, d_move = event_len, event_len

        data_moves[i] = d_move
        genome_moves[i] = g_move
        i += 1

    if reverse:
        string_pos += 1
        data = data[::-1]
        data_moves = data_moves[::-1]
        genome_moves = genome_moves[::-1]
    else:
        string_pos -= 1

    for d_move, g_move in zip(data_moves, genome_moves):
        if reverse:
            if target_pos < string_pos:
                string_pos -= g_move
            else:
                break
        else:
            if string_pos < target_pos:
                string_pos += g_move
            else:
                break
        j += d_move

    diff = string_pos - (target_pos + 1)if reverse else target_pos - string_pos
    if reverse:
        lt = data[j + diff :]
        lt = lt[::-1]
        rt = data[: j + diff]
        rt = rt[::-1]
    else:
        lt = data[: j + diff]
        rt = data[j + diff :]

    return lt, rt


cpdef tuple get_local_reference(
    Variant target,
    list pileup,
    int window,
    UnsplicedLocalReference unspl_loc_ref,
    bint unspliced=False,
    bint splice_pattern_only=False,
):

    cdef str span
    cdef tuple ptrn

    chrom, pos, reference = target.chrom, target.pos, target.reference

    if unspliced:
        splice_patterns = None
    else:
        splice_patterns = [read["splice_pattern"] for read in pileup if read["splice_pattern"] != ("", "")]

    ref_len = reference.get_reference_length(chrom)

    cdef list spl_ptrn = []

    if splice_patterns:
        lt_patterns = [ptrn[0] for ptrn in splice_patterns if ptrn[0]]
        if lt_patterns:
            lt_pattern = most_common(lt_patterns)
            lt_spl_pos = []
            for span in lt_pattern.split(":"):
                lt_spl_pos += [int(i) for i in span.split("-")]
        else:
            lt_spl_pos = []

        rt_patterns = [ptrn[1] for ptrn in splice_patterns if ptrn[1]]
        if rt_patterns:
            rt_pattern = most_common(rt_patterns)
            rt_spl_pos = []
            for span in rt_pattern.split(":"):
                rt_spl_pos += [int(i) for i in span.split("-")]
        else:
            rt_spl_pos = []

        spl_pos = lt_spl_pos + rt_spl_pos
        last_idx = len(spl_pos) - 1

        left_len = 0
        first_pass = False
        local_reference = ""
        for i, x in enumerate(spl_pos):
            if i == 0:
                lt_end = max(0, x - window * 2)
                local_reference += reference.fetch(chrom, lt_end, x - 1)
                rt_end = x - 1
                if x + 1 < rt_end:
                    spl_ptrn.append((x + 1, rt_end))
                else:
                    spl_ptrn.append((lt_end, rt_end))
            elif i % 2 == 1 and i != last_idx:
                local_reference += reference.fetch(chrom, x, spl_pos[i+1] - 1)
                rt_end = spl_pos[i+1] - 1
                spl_ptrn.append((x + 1, rt_end))
            elif i % 2 == 0:
                pass
            elif i == last_idx:
                rt_end = min(x + window * 2, ref_len)
                local_reference += reference.fetch(chrom, x, rt_end)
                spl_ptrn.append((x + 1, rt_end))

            if pos <= rt_end and not first_pass:
                left_len = len(local_reference) - (rt_end - pos)
                first_pass = True

    else:
        local_reference = unspl_loc_ref.fetch_ref_seq(pos, window)
        #left_len = unspl_loc_ref.left_len
        #local_reference = reference.fetch(chrom, max(0, pos - window * 3), min(pos + window * 3, ref_len))
        left_len = pos - max(0, pos - window * 3)

    if splice_pattern_only:
        return tuple(spl_ptrn)

    return local_reference, left_len
