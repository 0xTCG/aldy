#cython: profile=False

from .variant import Variant
from .utilities cimport split
from .utilities import split_cigar, get_local_reference, relative_aln_pos, most_common
from .localn import make_aligner, align, findall_indels, findall_mismatches, is_worth_realn

from .utilities cimport split


def find_by_normalization(
    target,
    pileup,
    window,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    basequalthresh=24,
    is_first_pass=True,
):
    """Search indels equivalent to the target indel

    Args:
        target (Variant): target indel
        pileup (list): a list of dictized read (dict)
    Returns:
        annoated pileup (list): a list of dictized read (dict)
    """

    pos = target.pos

    if target.is_indel:
        pileup = [is_target_by_normalization(read, target) for read in pileup]
        _pos = [read.get("observed_pos", pos) for read in pileup if read["is_target"]]
        pos = most_common(_pos) if _pos else pos
        are_read_ends = [read.get("is_read_end", False) for read in pileup if read["is_target"] and not read["is_dirty"]]
        if are_read_ends:
            read_end_evidence_only = all(are_read_ends)
        else:
            read_end_evidence_only = False
    else:
        is_single = (target.variant_type == "S")
        alt_bases = target.alt
        pileup = [is_substitute_target(read, pos, alt_bases, is_single) for read in pileup]
        read_end_evidence_only = False

    return target, pileup, gap_extension_penalty, pos, read_end_evidence_only

    ##target, extension_penalty_used = seek_larger_gapped_aln(
    #    target,
    #    pileup,
    #    window,
    #    match_score,
    #    mismatch_penalty,
    #    gap_open_penalty,
    #    gap_extension_penalty,
    #    basequalthresh,
    #    is_first_pass,
    #)

    #if is_first_pass and extension_penalty_used == 255:
    #    return find_by_normalization(
    #        target,
    #        pileup,
    #        window,
    #        match_score,
    #        mismatch_penalty,
    #        gap_open_penalty,
    #        gap_extension_penalty,
    #        basequalthresh=24,
    #        is_first_pass=False,
    #    )
    #else:
    #    return target, pileup, extension_penalty_used




def is_target_by_normalization(read, target):
    """Check if read contains an indel equivalent to target

    read (dict): dictized read
    target (Variant)
    """
    if read.get("is_target", False):
        return read
    else:
        read["is_target"] = False

    # trivial case
    if read["is_reference_seq"]:
        read["is_worth_realn"] = False
        read["mismatches"] = []
        return read

    for indel in read[target.variant_type]:
        if target == indel[-1]:

            pos = target.pos

            if avoid_left_aln(read, target):
                pos = indel[0] # pos observed in alignment
                read["observed_pos"] = pos

            read["is_target"] = True

            # trim clipped bases
            lt_offset = read["start_offset"]
            read["lt_flank"] = indel[1]
            read["lt_ref"] = indel[4]
            read["lt_qual"] = indel[6]

            read["indel_seq"] = indel[2]

            rt_offset = read["end_offset"]
            read["rt_flank"] = indel[3]
            read["rt_ref"] = indel[5]
            read["rt_qual"] = indel[7]

            read["lt_cigar"], read["rt_cigar"] = split_cigar(
                read["cigar_string"], pos, read["read_start"]
            )

            if min(len(read["lt_ref"]), len(read["rt_ref"])) / len(read["read_seq"]) < 0.15 or "H" in read["cigar_string"]:
                read["is_read_end"] = True
            else:
                read["is_read_end"] = False

    findall_mismatches(read)
    read["is_worth_realn"] = is_worth_realn(read, target, qual_lim=23)

    return read


def avoid_left_aln(read, target):
    #check if normalized pos should be used

    pos = target.pos #normalized
    if "N" in read["cigar_string"]:
        if  read["covering_subread"]:
            return pos < read["covering_subread"][0]
        else:
            return True

    return False


def is_substitute_target(read, target_pos, alt_bases, is_single):
    """this is experimental
    """
    if read.get("is_target", False):
        return read
    else:
        read["is_target"] = False

    if read["is_reference_seq"] or not read["is_covering"]:
        return read

    lt_end, rt_end = split(
                        read["read_seq"],
                        read["cigar_string"],
                        target_pos,
                        read["read_start"],
                        False,
                        False,
                     )

    lt_base = lt_end[-1] if len(lt_end) else ""
    if is_single:
        read["is_target"] = (lt_base == alt_bases)
    else:
        rt_fetch_len = min(len(alt_bases) - 1, len(rt_end))
        rt_bases = rt_end[: rt_fetch_len]
        read["is_target"] = match_mnv(alt_bases, lt_base, rt_bases)

    return read


def match_mnv(truth, lt_base, rt_bases):

    if lt_base and rt_bases:
        return truth == (lt_base + rt_bases)
    elif lt_base:
        return truth[0] == lt_base
    else:
        return truth[-len(rt_bases)] == rt_bases


def get_most_centered_read(target, pileup, target_annotated=True):

    most_centered_read = None
    center_score = 0

    if target_annotated:
        targetpileup = [
            read for read in pileup if read["is_target"] and not read["is_dirty"]
        ]
    else:
        targetpileup = [read for read in pileup if not read["is_dirty"]]

    if targetpileup:
        dist2center = [
            0.5
            - relative_aln_pos(
                read["ref_seq"], read["cigar_list"], read["aln_start"], target.pos,
            )
            for read in targetpileup
        ]

        abs_dist2center = [abs(i) for i in dist2center]
        most_central = min(abs_dist2center)
        most_centered_read = targetpileup[abs_dist2center.index(most_central)]
        center_score = dist2center[abs_dist2center.index(most_central)]

    return most_centered_read, center_score


def get_closest_gap(center_score, read_end, target, pileup):
    pos_look_up = {}
    read_look_up = {}
    for read in pileup:
        if (
            not read["is_reference_seq"]
            and read["is_covering"]
            and (read["D"] or read["I"])
        ):

            gaps = []
            if center_score >= 0:
                if (
                    read["aln_start"] < target.pos - len(read_end)
                    and read["is_covering"]
                ):
                    gaps = [
                        i[-1] for i in read["D"] + read["I"] if i[-1] != target
                    ]
            else:
                if read["aln_end"] > target.pos + len(read_end) and read["is_covering"]:
                    gaps = [
                        i[-1] for i in read["D"] + read["I"] if i[-1] != target
                    ]

            for i in gaps:
                if i in pos_look_up:
                    read_look_up[i].append(read)
                else:
                    pos_look_up[i] = abs(i.pos - target.pos)
                    read_look_up[i] = [read]
    if pos_look_up:
        closest_gap = min(pos_look_up, key=pos_look_up.get)
        closest_gap_reads = read_look_up[closest_gap]
        central_closest_gap_read, score = get_most_centered_read(
            closest_gap, closest_gap_reads, target_annotated=False,
        )

        if closest_gap and central_closest_gap_read:
            return closest_gap, central_closest_gap_read
        else:
            return None
    else:
        return None


def seek_larger_gapped_aln(
    target,
    pileup,
    window,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    basequalthresh,
    is_first_pass,
):

    read, center_score = get_most_centered_read(target, pileup)

    if not read:
        return target, gap_extension_penalty
    else:
        read_seq, ref_seq, cigarstring = (
            read["read"].query_alignment_sequence,
            read["ref_seq"],
            read["cigar_string"],
        )
        lt_read, rt_read = split(
            read_seq,
            cigarstring,
            target.pos,
            read["aln_start"],
            is_for_ref=True,
            reverse=False,
        )

        lt_ref, rt_ref = split(
            ref_seq,
            cigarstring,
            target.pos,
            read["aln_start"],
            is_for_ref=True,
            reverse=False,
        )

        lt_qual, rt_qual = split(
            read["read_qual"],
            cigarstring,
            target.pos,
            read["read_start"],
            is_for_ref=False,
            reverse=False,
        )

        if target.is_ins:
            rt_read = rt_read[len(target.indel_seq):]
        elif target.is_del:
            rt_ref = rt_ref[len(target.indel_seq):]

        with_end_mut = False
        if center_score >= 0:
            if lt_read != lt_ref and min(lt_qual) > basequalthresh:
                with_end_mut = True
        else:
            if rt_read != rt_ref and min(rt_qual) > basequalthresh:
                with_end_mut = True

    if is_first_pass and with_end_mut:

        read_end = lt_read if center_score >= 0 else rt_read

        if len(read_end) / len(read["read_seq"]) < 0.25:
            res = get_closest_gap(center_score, read_end, target, pileup)
            if res:
                closet_gap, closest_gap_read = res[0], res[1]
                subject_aligned_seq = closest_gap_read["read"].query_alignment_sequence
                query_aligned_seq = read["read"].query_alignment_sequence

                diff = len(query_aligned_seq) - len(subject_aligned_seq)
                if diff > 0:
                    if center_score >= 0:
                        query_aligned_seq = query_aligned_seq[:-diff]
                    else:
                        query_aligned_seq = query_aligned_seq[diff:]

                if read_end in query_aligned_seq and len(query_aligned_seq) > 30:
                    if query_aligned_seq in subject_aligned_seq:
                        return closet_gap, 255

    if "N" in read["cigar_string"]:
        ref_seq, lt_len = get_local_reference(target, [read], window)
    else:
        ref_seq, lt_len = get_local_reference(target, [read], window, unspliced=True)

    orig_gap_extension_penalty = gap_extension_penalty
    gap_extension_penalty = (
        0 if abs(center_score) > 0.35 and with_end_mut else gap_extension_penalty
    )
    aln = align(
        make_aligner(ref_seq, match_score, mismatch_penalty),
        read_seq,
        gap_open_penalty,
        gap_extension_penalty,
    )

    genome_aln_pos = target.pos + 1 - lt_len + aln.reference_start

    indels, mismatches = findall_indels(aln, genome_aln_pos, ref_seq, read_seq, report_snvs=True)
    if len(indels) != 1 or mismatches:
        return target, orig_gap_extension_penalty

    closest = min([abs(target.pos - indel["pos"]) for indel in indels])
    if "N" in read["cigar_string"] and closest > 3:
        return target, gap_extension_penalty

    candidates = [
        indel for indel in indels if abs(target.pos - indel["pos"]) == closest
    ]

    if candidates:
        candidate = candidates[0]
        if candidate["indel_type"] == "I":
            ref = candidate["lt_ref"][-1]
            alt = ref + candidate["indel_seq"]
        else:
            alt = candidate["lt_ref"][-1]
            ref = alt + candidate["del_seq"]

        target = Variant(target.chrom, candidate["pos"], ref, alt, target.reference)

    return target, gap_extension_penalty
