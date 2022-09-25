#!/usr/bin/env python3

from .utilities cimport split
from .contig cimport Contig
from .variant cimport Variant
from .utilities import get_end_pos
from .consensus import is_compatible


def find_by_softclip_split(target, contig, pileup):
    """Annotate if reads contain target indel

    Args:
        target (Variant)
        template (dict): indel template
        pileup (list): list of dictized reads (dict)
    Return:
        pileup (list): annotated pileup
    """

    pos, indel_type, indel_seq = (target.pos, target.variant_type, target.indel_seq)

    pileup = [
        find_candidate_softclips(read, pos, indel_type, indel_seq) for read in pileup
    ]
    pileup = [
        is_target_by_sftclp_split(read, pos, indel_type, indel_seq, contig)
        for read in pileup
    ]

    return pileup


def find_candidate_softclips(read, pos, indel_type, indel_seq):
    """Find and annotate softclipped indels to be realigned

    Args:
        read (dict): dictized read with "is_target" key
        pos (int): 1-based
        indel_type (str): "I" for insertion "D" for deletion
        indel_seq (str): inserted or deleted sequence
    Returns:
        read (dict): dictized read with softclip info annotated.
    """

    # already identified as target or seq identical to reference
    if read["is_target"] or read["is_reference_seq"]:
        read["softclip_pattern"] = None
        return read

    # no softclipped bases
    if not "S" in read["cigar_string"]:
        read["softclip_pattern"] = None
        return read

    # softclipped reads covering the locus
    if read["is_covering"]:
        read["softclip_pattern"] = classify_softclip_patterns(read, pos)
        return read

    # reads with large deletion may not cover the locus
    if read["read_end"] < pos:
        if indel_type == "D" and pos < read["read_end"] + len(indel_seq):
            read["softclip_pattern"] = "trailing_deletion"
        else:
            read["softclip_pattern"] = None
    else:
        if indel_type == "D" and read["read_start"] - len(indel_seq) < pos:
            read["softclip_pattern"] = "leading_deletion"
        else:
            read["softclip_pattern"] = None

    return read


def classify_softclip_patterns(read, pos):
    """Annotate if softclip starts before or after indel postion

    Args:
        read (dict): dictized read
        pos (int): 1-based
    Returns:
        softclip pattern (str)
    """
    # event_pos = read["read_start"]  # 1-based genomic pos
    event_pos = read["covering_subread"][0]  # 1-based genomic pos

    last_event = "O"
    for i, c in enumerate(read["cigar_list"]):
        event, event_len = c[-1], int(c[:-1])
        event_pos += event_len

        if pos <= event_pos:
            last_event = event
            is_leading = i == 0
            break

    if last_event == "M":
        return "off_clipping"
    elif last_event == "S" and is_leading:
        return "leading"
    elif last_event == "S" and not is_leading:
        return "trailing"
    else:
        return "other"


def is_target_by_sftclp_split(read, pos, indel_type, indel_seq, contig, slided=False):
    """Find softclipped target indels

    Args:
        read (dict): dictized read
        pos (int): 1-based
        indel_type (str): "I" for insertion and "D" for deletion
        indel_seq (str): inserted or deleted sequence
        template (dict): indel template
        slided (bool): True if "slide_insertion(...)" applied
                       Used internally for recurrence
    Returns:
        read (dict):
    """
    if read["is_target"] or not read["softclip_pattern"]:
        return read

    read = split_softclipped_read(read, pos, indel_type, indel_seq)

    read["is_target"] = is_compatible(read, contig, indel_type)

    # slide insertions
    if not read["is_target"] and not slided and indel_type == "I":
        return is_target_by_sftclp_split(
            slide_insertion(read, contig),
            pos,
            indel_type,
            indel_seq,
            contig,
            slided=True,
        )

    if slided:
        read["read_start"] = read["orig_start"]
        read["read_end"] = read["orig_end"]

        del read["orig_start"], read["orig_end"]

    return read


def split_softclipped_read(read, pos, indel_type, indel_seq):

    indel_len = len(indel_seq)
    cigar_string = read["cigar_string"]
    reverse = True if read["softclip_pattern"] == "leading" else False
    string_pos = read["read_end"] if reverse else read["read_start"]

    if indel_type == "D" and reverse:
        pos += indel_len

    lt_flank, rt_flank = split(
        read["read_seq"],
        cigar_string,
        pos,
        string_pos,
        is_for_ref=False,
        reverse=reverse,
    )
    mid_seq = ""
    lt_qual, rt_qual = split(
        read["read_qual"],
        cigar_string,
        pos,
        string_pos,
        is_for_ref=False,
        reverse=reverse,
    )

    if indel_type == "I":
        mid_seq, rt_flank = rt_flank[:indel_len], rt_flank[indel_len:]
        read["del_seq"] = ""
    else:
        read["del_seq"] = indel_seq
        #del_pos = get_end_pos(read["read_start"], lt_flank, cigar_string)
        #del_pos = pos
        #aln_pos = read["aln_end"] if reverse else read["aln_start"]

        #lt_ref, rt_ref = split(
        #    read["ref_seq"],
        #    cigar_string,
        #    del_pos,
        #    aln_pos,
        #    is_for_ref=True,
        #    reverse=reverse,
        #)

        #read["del_seq"] = "GAATTAAGAGAAGCA"

    read["lt_flank"] = lt_flank
    read["lt_qual"] = lt_qual
    read["indel_seq"] = mid_seq
    read["rt_flank"] = rt_flank
    read["rt_qual"] = rt_qual

    return read


def slide_insertion(read, contig):
    """Move read by the length of insertion
       (this is for BWA aligned reads)

    Args:
        read (dict): dictized read
        template (dict): indel template
    Returns:
        read (dict): dictized read with slided position
    """
    total_slide = sum([int(c[:-1]) for c in contig.gaps if "I" in c])
    read["orig_start"] = read["read_start"]
    read["orig_end"] = read["read_end"]

    if read["softclip_pattern"] == "leading":
        read["read_start"] += total_slide
        read["softclip_pattern"] = "other"
    else:
        read["read_end"] -= total_slide
        read["softclip_pattern"] = "leading"

    return read
