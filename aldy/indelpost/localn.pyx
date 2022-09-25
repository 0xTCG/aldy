#cython: profile=False

import re
import numpy as np

from .sswpy import SSW
from .consensus import is_compatible
from .utilities import get_mapped_subreads, get_end_pos, make_insertion_first, to_minimal_repeat_unit
from .utilities cimport split, count_lowqual_non_ref_bases

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


def find_by_smith_waterman_realn(
    target_indel,
    contig,
    pileup,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    basequalthresh,
    mapq_lim=1,
):
    """Annotate if reads contain target indel

    Args:
        target_indel (Variant)
        template (dict): indel template
        pileup (list): list of dictized reads (dict)
        mapq (int): minimum mapping quality
    Returns:
        pileup (list): annotated pileup
    """
    indel_type, indel_seq = target_indel.variant_type, target_indel.indel_seq

    mut_ref_lt, mut_ref_mid, mut_ref_rt = contig.get_contig_seq(split=True)
    ref_ref = contig.get_reference_seq()
    mut_ref = mut_ref_lt + mut_ref_mid + mut_ref_rt

    mut_aligner = make_aligner(mut_ref, match_score, mismatch_penalty)
    ref_aligner = make_aligner(ref_ref, match_score, mismatch_penalty)

    pileup = [findall_mismatches(read) for read in pileup]

    pileup = [
        is_target_by_ssw(
            read,
            target_indel,
            contig,
            mut_ref_lt,
            mut_ref_mid,
            mut_ref_rt,
            mut_aligner,
            ref_aligner,
            match_score,
            mismatch_penalty,
            gap_open_penalty,
            gap_extension_penalty,
            indel_type,
            basequalthresh,
            mapq_lim,
        )
        for read in pileup
    ]

    return pileup


def findall_mismatches(read, end_trim=0):

    if read["is_reference_seq"]:
        read["mismatches"] = []
        return read

    aln_start, aln_end = read["aln_start"], read["aln_end"]
    mapped_subpreads = get_mapped_subreads(read["cigar_string"], aln_start, aln_end)

    mismatches = []
    for subread in mapped_subpreads:
        start, end = subread[0], subread[1]
        span = end - start + 1

        # trim clipped segments
        cigarstring = read["cigar_string"]
        if "S" in cigarstring:
            cigarlst = read["cigar_list"]
            read_seq = read["read_seq"]
            quals = read["read_qual"]

            if "S" in cigarlst[0]:
                cigarlst = cigarlst[1:]
                read_seq = read_seq[read["start_offset"] :]
                quals = quals[read["start_offset"] :]

            if "S" in cigarlst[-1]:
                cigarlst = cigarlst[:-1]
                read_seq = read_seq[: -read["end_offset"]]
                quals = quals[: -read["end_offset"]]

            cigarstring = "".join(cigarlst)
        else:
            read_seq = read["read_seq"]
            quals = read["read_qual"]

        lt_seq, rt_seq = split(
            read_seq, cigarstring, start, aln_start, is_for_ref=False, reverse=False,
        )
        lt_qual, rt_qual = split(
            quals, cigarstring, start, aln_start, is_for_ref=False, reverse=False,
        )
        lt_ref, rt_ref = split(
            read["ref_seq"],
            cigarstring,
            start,
            aln_start,
            is_for_ref=True,
            reverse=False,
        )

        mapped_seq = lt_seq[-1] + rt_seq[: span - 1]
        mapped_qual = [lt_qual[-1]] + list(rt_qual[: span - 1])
        mapped_ref = lt_ref[-1] + rt_ref[: span - 1]

        pos = start  # pos for first elem of the rt side
        for r, a, q in zip(mapped_ref, mapped_seq, mapped_qual):
            if r != a:
                if aln_start + end_trim <= pos <= aln_end - end_trim:
                    mismatches.append((pos, r.upper(), a, q))

            pos += 1

    read["mismatches"] = mismatches

    return read


def is_worth_realn(read, target_indel, qual_lim=23):

    if read["covering_subread"]:
        is_covered = True
        covering_start, covering_end = (
            read["covering_subread"][0],
            read["covering_subread"][1],
        )
    else:
        is_covered = False
        if target_indel.is_ins:
            return False
        else:
            covering_start = target_indel.pos
            covering_end = covering_start + len(target_indel.ref)

    dist_to_left_end = target_indel.pos - read["aln_start"]
    dist_to_right_end = read["aln_end"] - target_indel.pos
    if dist_to_left_end < 0:
        is_lefty = True
    elif dist_to_right_end < 0:
        is_lefty = False
    else:
        is_lefty = (dist_to_left_end <= dist_to_right_end)

    start_cigar, end_cigar = read["cigar_list"][0], read["cigar_list"][-1]

    # start clipped
    if is_lefty and covering_start < read["aln_start"] <= covering_end and int(start_cigar[:-1]) > 2:
        return True

    # end clipped
    if not is_lefty and  covering_start <= read["aln_end"] < covering_end and int(end_cigar[:-1]) > 2:
        return True

    mismatches = [
        var
        for var in read["mismatches"]
        if covering_start <= var[0] <= covering_end and var[3] > qual_lim
    ]

    shiftable_pos = [v.pos for v in target_indel.generate_equivalents()]
    lt_pos, rt_pos = min(shiftable_pos), max(shiftable_pos)

    if lt_pos < rt_pos:
        if is_lefty:
            if lt_pos < read["aln_start"]:
                lt_end_read = read["read_seq"][: (rt_pos - read["aln_start"])]
                lt_end_ref = read["ref_seq"][: (rt_pos - read["aln_start"])]
                if lt_end_read == lt_end_ref:
                    return False
        else:
            if read["aln_end"] <= rt_pos:
                rt_end_read = read["read_seq"][-(read["aln_end"] - lt_pos) :]
                rt_end_ref = read["ref_seq"][-(read["aln_end"] - lt_pos) :]
                if rt_end_read == rt_end_ref:
                    return False


    if mismatches:
        if is_lefty:
            lt_most_pos = min(var[0] for var in mismatches)
            with_end_mismatches = abs(lt_most_pos - read["aln_start"]) < 4
        else:
            rt_most_pos = max(var[0] for var in mismatches)
            with_end_mismatches = abs(rt_most_pos - read["aln_end"]) < 4

        if with_end_mismatches:
            return True
        else:
            if is_covered:
                return True
            else:
                return False

    indels = [
        var for var in read["I"] + read["D"] if covering_start <= var[0] <= covering_end
    ]

    if indels:
        return True

    return False

def is_target_by_ssw(
    read,
    target_indel,
    contig,
    mut_ref_lt,
    mut_ref_mid,
    mut_ref_rt,
    mut_aligner,
    ref_aligner,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    indel_type,
    basequalthresh,
    mapq_lim,
    mapped_base_cnt_thresh=40,
    allow_mismatches=10,
):

    # already found
    if read["is_target"]:
        return read

    if read["is_reference_seq"] or read["mapq"] <= mapq_lim or not is_worth_realn(read, target_indel):
        read["is_target"] = False
        return read

    read_seq = read["read_seq"]

    ref_aln = align(ref_aligner, read_seq, gap_open_penalty, gap_extension_penalty)
    # forced to align by setting gap_open_penalty=len(read_seq)
    mut_aln = align(mut_aligner, read_seq, len(read_seq), gap_extension_penalty)

    if mut_aln.optimal_score <= ref_aln.optimal_score:
        read["is_target"] = False
        return read

    mut_aln_cigar = mut_aln.CIGAR
    mut_aln_ref_start = mut_aln.reference_start
    mut_aln_ref_end = mut_aln.reference_end
    mut_aln_read_start = mut_aln.read_start
    mut_aln_read_end = mut_aln.read_end
    mut_cigar_list = cigar_ptrn.findall(mut_aln.CIGAR)

    target_covered = is_covering_target(
                    read["read_name"],
                    read_seq,
                    target_indel.indel_seq,
                    mut_ref_lt,
                    mut_ref_mid,
                    mut_ref_rt,
                    mut_aln_cigar,
                    len(read_seq),
                    mut_aln_ref_start,
                    mut_aln_ref_end,
                    mut_aln_read_start,
                    mut_aln_read_end,
                    target_indel.count_repeats(),
                   )

    if target_covered == 1:
        read["is_target"] = True
        return read
    elif target_covered == -1:
        read["undetermined"] = True
        return read
    else:
        return read

cdef int is_covering_target(
    str readname,
    str read_seq,
    str indel_seq,
    str mut_ref_lt,
    str mut_ref_mid,
    str mut_ref_rt,
    str mut_aln_cigar,
    int read_seq_len,
    int ref_aln_start,
    int ref_aln_end,
    int read_aln_start,
    int read_aln_end,
    int n_repeats,
):
    cdef list mut_cigar_list = cigar_ptrn.findall(mut_aln_cigar)

    if len(mut_cigar_list) > 1:
        return False

    cdef str repeat_unit = to_minimal_repeat_unit(indel_seq)
    cdef int mapped_len = int(mut_cigar_list[0][:-1])
    cdef int mut_ref_lt_len = len(mut_ref_lt)
    cdef int mut_ref_mid_len = len(mut_ref_mid)
    cdef int mut_ref_rt_len = len(mut_ref_rt)
    cdef int lt_consumption
    cdef int rt_consumption
    cdef int total_consumption = read_aln_end - read_aln_start + 1
    cdef bint is_lt_indel_end_reached
    cdef bint is_rt_indel_end_reached
    cdef bint is_lt_read_consumed = (read_aln_start == 0)
    cdef bint is_rt_read_consumed = (read_aln_end == read_seq_len - 1)

    # not covering indel
    if ref_aln_end < mut_ref_lt_len:
        return False
    elif mut_ref_lt_len + mut_ref_mid_len <= ref_aln_start:
        return False

    # ins
    if mut_ref_mid_len:
        if ref_aln_start < mut_ref_lt_len:
            lt_consumption = mut_ref_lt_len - ref_aln_start

            if total_consumption > lt_consumption + mut_ref_mid_len:
                lt_read = read_seq[read_aln_start : read_aln_start + lt_consumption]
                mid_read = read_seq[read_aln_start + lt_consumption : read_aln_start + lt_consumption + mut_ref_mid_len]
                rt_read = read_seq[read_aln_start + lt_consumption + mut_ref_mid_len : read_aln_end + 1]

                lt_check = is_compatible_repeats(lt_read, repeat_unit, n_repeats, is_left=True)
                rt_check = is_compatible_repeats(rt_read, repeat_unit, n_repeats, is_left=False)

                if lt_check and rt_check:
                    #True
                    return 1
                else:
                    #False (undetermined)
                    return -1
            else:
                ### abolish heuristic for reapeat
                #if n_repeats:
                #    return False
                if is_rt_read_consumed:
                        #True
                        return 1
                else:
                    # aligned from left and ended in inserted seq
                    rt_consumption = total_consumption - lt_consumption
                    consumed_mid_seq = mut_ref_mid[:rt_consumption:]

                    if consumed_mid_seq == read_seq[-rt_consumption:]:
                        #True
                        return 1
                    else:
                        #False
                        return 0
        else:
            # no left side alignment

            ### abolish heuristic for repeat
            # instead check the exact match
            #if n_repeats:
            #    return False

            if is_lt_read_consumed:
                #True
                return 1
            else:
                # aligned from right and ended in inserted seq
                lt_consumption = mut_ref_lt_len + mut_ref_mid_len - ref_aln_start
                consumed_mid_seq = mut_ref_mid[-lt_consumption:]

                # require exact match
                if consumed_mid_seq == read_seq[:lt_consumption]:
                    #True
                    return 1
                else:
                    #False
                    return 0

    else:
        lt_consumption = mut_ref_lt_len - ref_aln_start
        rt_consumption = total_consumption - lt_consumption

        lt_read = read_seq[read_aln_start : read_aln_start + lt_consumption]
        rt_read = read_seq[read_aln_start + lt_consumption  : read_aln_end]

        lt_check = is_compatible_repeats(lt_read, repeat_unit, n_repeats, is_left=True)
        rt_check = is_compatible_repeats(rt_read, repeat_unit, n_repeats, is_left=False)

        if lt_check and rt_check:
            pass
        else:
            #False (undetermined)
            return -1

        if lt_consumption <= rt_consumption:
            if is_lt_read_consumed:
                #True
                return 1
            else:
                if lt_consumption > 2:
                    #True
                    return 1
                else:
                    #False
                    return 0
        else:
            if is_rt_read_consumed:
                #True
                return 1
            else:
                if rt_consumption > 2:
                    #True
                    return 1
                else:
                    #False
                    return 0


def is_compatible_repeats(seq, repeat_unit, expected_n_repeats, is_left):
    unit_len = len(repeat_unit)


    if is_left:
        seq = seq[::-1]
        repeat_unit = repeat_unit[::-1]

    cnt = 0
    is_done = False
    while not is_done:
        if not seq:
            is_done = True
        elif repeat_unit == seq[: unit_len]:
            seq = seq[unit_len :]
            cnt += 1
        else:
            is_done = True

    is_compatible = True
    if not seq:
        is_compatible = False
    else:
        if cnt and cnt != expected_n_repeats:
            is_compatible = False

    return is_compatible




def make_aligner(ref_seq, match_score, mismatch_penalty):
    aligner = SSW(match_score=match_score, mismatch_penalty=mismatch_penalty)
    aligner.setReference(ref_seq)
    return aligner


def align(aligner, read_seq, gap_open_penalty, gap_extension_penalty):
    aligner.setRead(read_seq)
    return aligner.align(gap_open=gap_open_penalty, gap_extension=gap_extension_penalty)


def parse_read_by_mut_aln(mut_aln, contig, read, indel_type):
    """Decompose read based on the SSW alignment result for evaluation

    Args:
        mut_aln (named tuple): Alignment(CIGAR, optimal_score, sub_optimal_score, reference_start, reference_end, read_start, read_end)
        template (dict): indel template
        read (dict): dictized read
        indel_type (str): "I" for insertion or "D" for deletion
    Returns:
        read (dict): dictized read with read seq decomposed
    """
    lt_len, indel_len, rt_len = (
        len(contig.lt_consensus_seq),
        len(contig.indel_seq),
        len(contig.rt_consensus_seq),
    )

    read_seq = read["read_seq"]
    read_qual = read["read_qual"]
    ref_start, ref_end = mut_aln.reference_start, mut_aln.reference_end
    aln_start, aln_end = mut_aln.read_start, mut_aln.read_end

    lt_flank, mid_seq, rt_flank = "", "", ""
    lt_qual, rt_qual = [], []

    if ref_start <= lt_len:
        lt_diff = lt_len - ref_start
        lt_flank = read_seq[aln_start : aln_start + lt_diff]
        lt_qual = read_qual[aln_start : aln_start + lt_diff]
        if indel_type == "I":
            end_point = min(aln_start + lt_diff + indel_len, aln_end)
            mid_seq = read_seq[aln_start + lt_diff : end_point]
        else:
            rt_flank = read_seq[aln_start + lt_diff :]
            rt_qual = read_qual[aln_start + lt_diff :]
            del_pos = get_end_pos(
                read["read_start"] + aln_start, lt_flank, read["cigar_string"]
            )

            lt_ref, rt_ref = split(
                read["ref_seq"],
                read["cigar_string"],
                del_pos,
                read["aln_start"],
                is_for_ref=True,
                reverse=False,
            )

            read["del_pos"] = del_pos
            read["del_seq"] = rt_ref[:indel_len]

    if lt_len + indel_len <= ref_end and indel_type == "I":
        rt_diff = ref_end - (lt_len + indel_len)
        rt_flank = read_seq[aln_end - rt_diff : aln_end]
        rt_qual = read_qual[aln_end - rt_diff : aln_end]
        end_point = max(aln_start, aln_end - rt_diff - indel_len)
        mid_seq = read_seq[end_point : aln_end - rt_diff]

    read["lt_flank"] = lt_flank
    read["lt_qual"] = lt_qual
    read["indel_seq"] = mid_seq
    read["rt_flank"] = rt_flank
    read["rt_qual"] = rt_qual

    return read


def findall_indels(ref_aln, genome_aln_pos, ref_seq, read_seq, report_snvs=False, basequals=None):

    genome_aln_pos -= 1
    ref_idx = ref_aln.reference_start
    read_idx = ref_aln.read_start

    lt_clipped = read_seq[:read_idx]

    indels, snvs = [], []

    for token in cigar_ptrn.findall(make_insertion_first(ref_aln.CIGAR)):
        event, event_len = token[-1], int(token[:-1])

        if event == "I" or event == "D":
            indel = {}

            indel["pos"] = genome_aln_pos
            indel["lt_ref"] = ref_seq[:ref_idx]
            indel["lt_flank"] = read_seq[:read_idx]

            if basequals:
                indel["lt_qual"] = basequals[:read_idx]

            if event == "I":
                indel["indel_type"] = "I"
                indel["indel_seq"] = read_seq[read_idx : read_idx + event_len]
                indel["rt_ref"] = ref_seq[ref_idx:]
                indel["rt_flank"] = read_seq[read_idx + event_len :]
                indel["ref_idx"] = ref_idx
                indel["read_idx"] = read_idx

                if basequals:
                    indel["rt_qual"] = basequals[read_idx + event_len :]

                read_idx += event_len
            else:
                indel["indel_type"] = "D"
                indel["indel_seq"] = ""
                indel["del_seq"] = ref_seq[ref_idx : ref_idx + event_len]
                indel["rt_ref"] = ref_seq[ref_idx + event_len :]
                indel["rt_flank"] = read_seq[read_idx:]
                indel["ref_idx"] = ref_idx
                indel["read_idx"] = read_idx

                if basequals:
                    indel["rt_qual"] = basequals[read_idx:]

                ref_idx += event_len
                genome_aln_pos += event_len

            indels.append(indel)

        else:
            if report_snvs:
                i = 0
                while i < event_len:
                     if ref_seq[ref_idx + i : ref_idx + i + 1] != read_seq[read_idx + i : read_idx + i + 1]:
                         snv = {}
                         snv["pos"] = genome_aln_pos + i + 1
                         snv["ref"] = ref_seq[ref_idx + i : ref_idx + i + 1]
                         snv["alt"] = read_seq[read_idx + i : read_idx + i + 1]

                         snvs.append(snv)

                     i += 1

            ref_idx += event_len
            read_idx += event_len
            genome_aln_pos += event_len


    rt_clipped = read_seq[read_idx:]
    for indel in indels:
        indel["lt_clipped"] = lt_clipped
        indel["rt_clipped"] = rt_clipped

    if report_snvs:
        return indels, snvs
    else:
        return indels
