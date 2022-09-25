#cython: profile=False

import re
import random
cimport cython

from cpython cimport array
import array

from  functools import partial
from difflib import get_close_matches, SequenceMatcher


from .variant cimport Variant
from .local_reference cimport UnsplicedLocalReference

from .consensus import consensus_refseq
from .gappedaln import find_by_normalization

from .utilities cimport (
    split,
    locate_indels,
    get_spliced_subreads,
    count_lowqual_non_ref_bases,
)
from .utilities import (
    split_cigar,
    most_common,
    to_flat_list,
    get_mapped_subreads,
    get_local_reference,
    make_insertion_first,
)
from .localn import (
    align,
    make_aligner,
    is_worth_realn,
    findall_indels,
    findall_mismatches,
)

from pysam.libcfaidx cimport FastaFile
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libcalignmentfile cimport AlignmentFile

random.seed(123)

cigar_ptrn = re.compile(r"[0-9]+[MIDNSHPX=]")


cdef tuple make_pileup(
    Variant target,
    AlignmentFile bam,
    UnsplicedLocalReference unspl_loc_ref,
    bint exclude_duplicates,
    int window,
    int downsamplethresh,
    int basequalthresh,
):
    cdef str chrom
    cdef int pos
    cdef FastaFile reference
    cdef AlignedSegment seg
    cdef dict read

    chrom, pos, reference = target.chrom, target.pos, target.reference
    rpos = max(v.pos for v in target.generate_equivalents())

    ref_len = reference.get_reference_length(chrom)

    chroms = bam.references
    if chrom not in chroms:
        if chrom.startswith("chr"):
            _chrom = chrom.replace("chr","")
        else:
            _chrom = "chr" + chrom
    else:
        _chrom = chrom

    pileup = fetch_reads(_chrom, pos, bam, ref_len, window, exclude_duplicates)
    call_back = "all" if exclude_duplicates else "nofilter"
    orig_depth = bam.count(_chrom, pos - 1, pos, read_callback=call_back)
    orig_read_num = len(pileup)

    # downsampling
    if orig_depth > downsamplethresh:
        random.seed(123)
        n_sample = int(orig_read_num * (downsamplethresh / orig_depth))

        ### lower bounded by downsample thresh / 2
        ### to prevent over-downsampling
        ### example: orig_depth = 1197119 (exclude_duplicates=False)
        ###          orig_read_num = 1962
        ###          downsamplethresh = 1000
        ###          n_sample -> 1
        ###          over downsampled.
        if n_sample >= downsamplethresh / 2 > 0:
            pileup = random.sample(pileup,  n_sample)
            sample_factor = orig_read_num / len(pileup)
        else:
            sample_factor = 1.0
    else:
        sample_factor = 1.0

    pileup = [
        dictize_read(seg, chrom, pos, rpos, reference, unspl_loc_ref, basequalthresh) for seg in pileup
    ]

    pileup = [read for read in pileup if not is_within_intron(read, pos, window)]

    return pileup, sample_factor


def is_within_intron(read, pos, window):
    intron = read["intron_pattern"]
    if intron == (0, 0):
        return False
    else:
        intron_start, intron_end = intron[0], intron[1]
        if intron_start < pos - window and pos + window < intron_end:
            return True
        else:
            return False


cdef list fetch_reads(str chrom, int pos, AlignmentFile bam, int ref_len, int window, bint exclude_duplicates):

    cdef AlignedSegment read

    pos = pos - 1  # convert to 0-based

    all_reads = bam.fetch(
        chrom, max(0, pos - window), min(pos + 1 + window, ref_len), until_eof=True
    )

    if exclude_duplicates:
        reads = [
            read
            for read in all_reads
            if not read.is_duplicate
            and not read.is_secondary
            and read.cigarstring
            and read.reference_start
        ]
    else:
        reads = [
            read
            for read in all_reads
            if not read.is_secondary
            and read.cigarstring
        ]

    return reads


cdef dict dictize_read(
    AlignedSegment read,
    str chrom,
    int pos,
    int rpos,
    FastaFile reference,
    UnsplicedLocalReference unspl_loc_ref,
    int basequalthresh
):

    cdef tuple ins, deln

    cigar_string = read.cigarstring
    cigar_list = cigar_ptrn.findall(cigar_string)

    # adjust read start/end if starts or ends with softclipping
    aln_start = read.reference_start + 1
    start_offset = int(cigar_list[0][:-1]) if cigar_list[0].endswith("S") else 0
    read_start = aln_start - start_offset

    aln_end = read.reference_end  # do not add 1
    if aln_end is None:
        aln_end = aln_start + sum(int(c[:-1]) for c in cigar_list if c[-1] in ("M", "N", "D", "=", "X"))

    end_offset = int(cigar_list[-1][:-1]) if cigar_list[-1].endswith("S") else 0

    read_end = aln_end + end_offset

    read_seq = read.query_sequence
    read_qual = read.query_qualities
    ref_seq = get_ref_seq(chrom, aln_start, aln_end, cigar_string, cigar_list, reference, unspl_loc_ref)

    read_dict = {
        "read": read,
        "read_seq": read_seq,
        "read_qual": read_qual,
        "ref_seq": ref_seq,
        "is_reverse": read.is_reverse,
        "read_name": read.query_name,
        "mapq": read.mapping_quality,
        "start_offset": start_offset,
        "aln_start": aln_start,
        "read_start": read_start,
        "end_offset": end_offset,
        "aln_end": aln_end,
        "read_end": read_end,
        "cigar_string": cigar_string,
        "cigar_list": cigar_list,
        "is_reference_seq": (read_seq == ref_seq),
        "I": [],
        "D": [],
    }

    # base qual check
    read_dict["low_qual_base_num"] = count_lowqual_non_ref_bases(read_seq, ref_seq, read_qual, cigar_list, basequalthresh)
    read_dict["is_end_dirty"] = is_end_dirty(read_qual, basequalthresh, pos, read_start, read_end, cigar_string)
    read_dict["is_dirty"] = sum(q <= basequalthresh for q in read_qual) / len(read_seq) > 0.15

    insertions, deletions = locate_indels(cigar_string, read_start)

    for ins in insertions:
        read_dict["I"].append(
            leftalign_indel_read(
                chrom,
                ins[0],
                ins[1],
                "I",
                cigar_string,
                read_start,
                aln_start,
                read_seq,
                ref_seq,
                read_qual,
                reference,
            )
        )

    for deln in deletions:
        read_dict["D"].append(
            leftalign_indel_read(
                chrom,
                deln[0],
                deln[1],
                "D",
                cigar_string,
                read_start,
                aln_start,
                read_seq,
                ref_seq,
                read_qual,
                reference,
            )
        )

    indels = [ins[-1] for ins in read_dict["I"]] + [deln[-1] for deln in read_dict["D"]]

    (
        is_covering,
        covering_subread,
        is_spliced,
        splice_ptrn,
        intron_ptrn,
    ) = parse_spliced_read(cigar_string, read_start, read_end, pos, rpos)

    read_dict["is_covering"] = is_covering
    read_dict["covering_subread"] = covering_subread
    read_dict["is_spliced"] = is_spliced
    read_dict["splice_pattern"] = splice_ptrn
    read_dict["intron_pattern"] = intron_ptrn

    return read_dict


cdef str get_ref_seq(
    str chrom,
    int aln_start,
    int aln_end,
    str cigar_string,
    list cigar_list,
    FastaFile reference,
    UnsplicedLocalReference unspl_loc_ref,
):
    cdef int current_pos = aln_start - 1

    if not "N" in cigar_string:
        #reference.fetch(chrom, current_pos, aln_end)
        return unspl_loc_ref.get_ref_seq(current_pos, aln_end)

    cdef str ref_seq = ""
    cdef str event, cigar
    cdef int event_len

    for cigar in cigar_list:
        event, event_len = cigar[-1], int(cigar[:-1])
        if event == "M" or event == "D":
            ref_seq += reference.fetch(chrom, current_pos, current_pos + event_len)
            current_pos += event_len
        elif event in ("I", "S", "H", "P"):
            pass
        else:
            current_pos += event_len

    return ref_seq


cdef tuple leftalign_indel_read(
    str chrom,
    int pos,
    int indel_len,
    str indel_type,
    str cigar_string,
    int read_start,
    int aln_start,
    str read_seq,
    str ref_seq,
    object read_qual,
    FastaFile reference,
):
    lt_flank, rt_flank = split(
        read_seq, cigar_string, pos, read_start, is_for_ref=False, reverse=False
    )
    lt_ref, rt_ref = split(
        ref_seq, cigar_string, pos, aln_start, is_for_ref=True, reverse=False
    )
    lt_qual, rt_qual = split(
        read_qual, cigar_string, pos, read_start, is_for_ref=False, reverse=False
    )

    padding_base = reference.fetch(chrom, pos - 1, pos) if "N" in cigar_string or not lt_ref else lt_ref[-1]
    if indel_type == "I":
        indel_seq = rt_flank[:indel_len]
        rt_flank = rt_flank[indel_len:]
        rt_qual = rt_qual[indel_len:]
        var = Variant(chrom, pos, padding_base, padding_base + indel_seq, reference, skip_validation=True)
    else:
        indel_seq = rt_ref[:indel_len]
        rt_ref = rt_ref[indel_len:]
        var = Variant(chrom, pos, padding_base + indel_seq, padding_base, reference, skip_validation=True)

    return pos, lt_flank, indel_seq, rt_flank, lt_ref, rt_ref, lt_qual, rt_qual, var


cdef bint is_end_dirty(array.array read_qual, int basequalthresh, int pos, int read_start, int read_end, str cigar_string):
    cdef int dist_to_left_end = pos - read_start
    cdef int dist_to_right_end = read_end - pos
    cdef bint is_lefty

    if dist_to_left_end < 0:
        is_lefty = True
    elif dist_to_right_end < 0:
        is_lefty = False
    else:
        is_lefty = (dist_to_left_end <= dist_to_right_end)

    if cigar_string.count("N") > 1:
        return False
    else:
        if is_lefty:
            return min(read_qual[: 3]) < basequalthresh
        else:
            return min(read_qual[-3 :]) < basequalthresh


cdef str leftalign_cigar(str cigarstring, Variant target, int read_start):
    target.normalize(inplace=True)

    pos = target.pos

    lt_cigar_lst, rt_cigar_lst = split_cigar(cigarstring, pos, read_start)

    if len(rt_cigar_lst) < 3:
        return cigarstring

    tmp0, tmp1, tmp2 = rt_cigar_lst[0], rt_cigar_lst[1], rt_cigar_lst[2]
    if "M" in tmp0 and "M" in tmp2:
        tmp0, tmp2 = int(tmp0[:-1]), int(tmp2[:-1])
    else:
        return cigarstring

    new_cigar = tmp1 + str(tmp0 + tmp2) + "M" + "".join(rt_cigar_lst[3:])

    return "".join(lt_cigar_lst) + new_cigar


cdef tuple parse_spliced_read(str cigar_string, int read_start, int read_end, int pos, int rpos):
    spliced_subreads = get_spliced_subreads(cigar_string, read_start, read_end)


    is_covering = False
    covering_subread = None
    for subread in spliced_subreads:
        if subread[0] <= pos <= subread[1]:
            is_covering = True
            covering_subread = tuple(subread)
        elif subread[0] <= rpos <= subread[1]:
            is_covering = True
            covering_subread = tuple(subread)
            pos = rpos


    intron_ptrn = (0, 0)
    if len(spliced_subreads) > 1:
        is_spliced = True

        lt_ptrn, rt_ptrn = "", ""
        positions = to_flat_list(spliced_subreads)
        positions = positions[1:-1]

        i = 0
        while i < len(positions):
            start = positions[i] + 1
            end = positions[i + 1] - 1

            if end < pos:
                if not lt_ptrn:
                    lt_ptrn += str(start) + "-" + str(end)
                else:
                    lt_ptrn += ":" + str(start) + "-" + str(end)

            elif pos < start - 1:
                if not rt_ptrn:
                    rt_ptrn += str(start) + "-" + str(end)
                else:
                    rt_ptrn += ":" + str(start) + "-" + str(end)

            # overhang reads
            if start - 4 <= pos <= end:
                intron_ptrn = (start, end)

            i += 2

        splice_ptrn = (lt_ptrn, rt_ptrn)
    else:
        is_spliced = False
        splice_ptrn = ("", "")

    return is_covering, covering_subread, is_spliced, splice_ptrn, intron_ptrn


def check_overhangs(pileup, splice_rate=0.2):

    intron_ptrns = [read["intron_pattern"] for read in pileup if is_junctional(read)]
    introns = [ptrn for ptrn in intron_ptrns if ptrn != (0, 0)]
    if not introns:
        return None
    else:
        intron = most_common(introns)
        if intron_ptrns.count(intron) / len(intron_ptrns) < splice_rate:
            return None

    intron_start, intron_end = intron[0], intron[1]
    overhangs = [read for read in pileup if is_overhang(read, intron_start, intron_end)]
    if overhangs:
        return intron, overhangs
    else:
        return None


def is_junctional(read):
    if read["intron_pattern"] == (0, 0):
        return read["is_covering"]
    else:
        return True


def is_overhang(read, intron_start, intron_end):
    covering_subread = read["covering_subread"]
    if not covering_subread:
        return False

    lt_read_lim = max(covering_subread[0], read["aln_start"])
    rt_read_lim = min(covering_subread[1], read["aln_end"])

    if lt_read_lim < intron_start and rt_read_lim < intron_end:
        return True
    elif intron_start < lt_read_lim and intron_end < rt_read_lim:
        return True
    else:
        return False


def overhang_aligners(target, intron, match_score, mismatch_penalty):
    genome_ref = target.reference.fetch(
        target.chrom, target.pos - 100, target.pos + 100
    )
    genome_aligner = make_aligner(genome_ref, match_score, mismatch_penalty)

    lt_exon_end, rt_exon_start = intron[0] - 1, intron[1]

    junction_ref = target.reference.fetch(
        target.chrom, lt_exon_end - 100, lt_exon_end
    ) + target.reference.fetch(target.chrom, rt_exon_start, rt_exon_start + 100)

    junction_aligner = make_aligner(junction_ref, match_score, mismatch_penalty)

    return genome_aligner, junction_aligner


def filter_spurious_overhangs(
    target,
    intron,
    overhangs,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
):

    genome_aligner, junctional_aligner = overhang_aligners(
        target, intron, match_score, mismatch_penalty
    )
    non_spurious_overhangs = [
        read
        for read in overhangs
        if not read["is_reference_seq"]
        and is_non_spurious_overhang(
            read,
            target,
            intron,
            genome_aligner,
            junctional_aligner,
            match_score,
            mismatch_penalty,
            gap_open_penalty,
            gap_extension_penalty,
        )
    ]

    return non_spurious_overhangs


def is_non_spurious_overhang(
    read,
    target,
    intron,
    genome_aligner,
    junction_aligner,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
):
    read_seq = read["read_seq"]

    genome_aln = align(
        genome_aligner, read_seq, gap_open_penalty, gap_extension_penalty
    )
    junction_aln = align(
        junction_aligner, read_seq, gap_open_penalty, gap_extension_penalty
    )

    genome_score = genome_aln.optimal_score
    junction_score = junction_aln.optimal_score
    if genome_score <= junction_score:
        return False

    genome_cigar = make_insertion_first(genome_aln.CIGAR)
    gap_cnt = genome_cigar.count("I") + genome_cigar.count("D")

    if gap_cnt > 3:
        return False
    elif 1 < gap_cnt <= 3:
        if genome_score / junction_score < 1.2 or genome_score < match_score * 50:
            return False
    elif gap_cnt == 0:
        aln_len = genome_aln.read_end - genome_aln.read_start + 1
        if aln_len / len(read_seq) > 0.98:
            return False

    lt_exon_end, rt_exon_start = intron[0] - 1, intron[1]
    indels_within_intron = [
        lt_exon_end < var[-1].pos < rt_exon_start for var in read["D"] and read["I"]
    ]

    if indels_within_intron:
        return True

    read = findall_mismatches(read)
    return is_worth_realn(read, target)


def retarget(
    target,
    pileup,
    window,
    mapq4retarget,
    within,
    retargetcutoff,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    unspl_loc_ref,
    require_exact_for_shiftable,
):
    target_seq, target_type, target_pos  = target.indel_seq, target.variant_type, target.pos


    if target.is_ins:
        non_refs = [
            read_dict
            for read_dict in pileup
            if not read_dict["is_reference_seq"]
            and read_dict["is_covering"]
            and read_dict["mapq"] > mapq4retarget
        ]
    else:
        non_refs = [
            read_dict
            for read_dict in pileup
            if not read_dict["is_reference_seq"]
            and read_dict["mapq"] > mapq4retarget
        ]

    if not non_refs:
        return None

    # cutoff auto_adjustment
    if len(target.indel_seq) < 3:
        cutoff = 1.0
    else:
        cutoff = retargetcutoff

    #best_matches = None
    # check by realn if no gapped alignments or no good gapped alignmetns
    #if not candidates or not best_matches:
    #    is_gapped_aln = False

    tmp_non_refs = non_refs.copy()

    non_refs = [
        read_dict
        for read_dict
        in non_refs
        if read_dict["low_qual_base_num"] < 6
        and not read_dict["is_dirty"]
        and not read_dict["is_end_dirty"]
        and read_dict.get("is_worth_realn", True)
    ]

    if not non_refs:
        non_refs = [read_dict for read_dict in tmp_non_refs if not read_dict["is_dirty"]]

    ref_starts, ref_alns, ref_seqs, aligners = [], [], [], []
    for read in non_refs:
        ref_seq, lt_len = get_local_reference(target, [read], window, unspl_loc_ref)
        ref_seqs.append(ref_seq)

        aligner = make_aligner(ref_seq, match_score, mismatch_penalty)
        aligners.append(aligner)

        ref_alns.append(align(aligner, read["read_seq"], gap_open_penalty, gap_extension_penalty))
        ref_starts.append(target.pos + 1 - lt_len)

    complex_flags = []
    candidates, candidate_reads, candidate_ref_seqs, candidate_ref_starts, candidate_aligners = [], [], [], [], []
    for read, aln, ref_seq, ref_start, aligner in zip(non_refs, ref_alns, ref_seqs, ref_starts, aligners):
        genome_aln_pos = ref_start + aln.reference_start

        aligned_read_len = (aln.read_end - aln.read_start)
        window_len = window * 6
        aligned_frac = aligned_read_len / min(len(read["read_seq"]), window_len)

        gap_cnt = aln.CIGAR.count("I") + aln.CIGAR.count("D")

        if 0 < gap_cnt < 6 and aligned_frac > 0.7:
            indels = findall_indels(aln, genome_aln_pos, ref_seq, read["read_seq"])
            positions = [d["pos"] for d in indels]
            complex_positions = set([p for p in positions if positions.count(p) == 2])

            target_type_indels = [
                indel for indel in indels if indel["indel_type"] == target_type
            ]

            # check for possible complex case
            if complex_positions:
                complex_flags.append(1)
            else:
                pass

            for indel in target_type_indels:
                if indel["pos"] in complex_positions:
                    complex_del = [jndel for jndel in indels if jndel["pos"] == indel["pos"] and jndel["indel_type"] == "D"][0]
                    complex_ins = [jndel for jndel in indels if jndel["pos"] == indel["pos"] and jndel["indel_type"] == "I"][0]
                    ref = complex_del["lt_ref"][-1] + complex_del["del_seq"]
                    alt = complex_ins["lt_ref"][-1] + complex_ins["indel_seq"]
                else:
                    if target_type == "I":
                        ref = indel["lt_ref"][-1]
                        alt = ref + indel["indel_seq"]
                    else:
                        alt = indel["lt_ref"][-1]
                        ref = alt + indel["del_seq"]

                var = Variant(target.chrom, indel["pos"], ref, alt, target.reference, skip_validation=True)

                # non-target indel found in read end -> do not consider
                read_end_thresh = max(len(read["read_seq"]) / 30, 3)
                if var.pos - read["read_start"] <= read_end_thresh or read["read_end"] - var.pos <= read_end_thresh:
                    if var == target or (complex_positions and var.pos not in complex_positions):
                        candidates.append(var)
                        candidate_reads.append(read)
                        candidate_ref_seqs.append(ref_seq)
                        candidate_ref_starts.append(ref_start)
                        candidate_aligners.append(aligner)
                else:
                    candidates.append(var)
                    candidate_reads.append(read)
                    candidate_ref_seqs.append(ref_seq)
                    candidate_ref_starts.append(ref_start)
                    candidate_aligners.append(aligner)

    if not candidates:
        # long reference may align ins as del
        if target.is_ins and window > 3:
            window = int(window/3)
            return retarget(
                        target,
                        pileup,
                        window,
                        mapq4retarget,
                        within,
                        retargetcutoff,
                        match_score,
                        mismatch_penalty,
                        gap_open_penalty,
                        gap_extension_penalty,
                        unspl_loc_ref,
                        require_exact_for_shiftable,
                   )
        else:
            return None
    elif len(target.indel_seq) <= 3:
        if not sum(complex_flags) and target not in candidates:
            return None

    u_candidates = to_flat_list([var._generate_equivalents_private() for var in set(candidates)])
    u_candidates.sort(key=lambda x: abs(x.pos - target.pos))
    candidate_seqs = [var._get_indel_seq(how=target_type) for var in u_candidates]

    best_match = get_close_matches(
        target.indel_seq, candidate_seqs, n=1, cutoff=cutoff
    )

    if best_match:
        best_seq = best_match[0]
        match_score = SequenceMatcher(None, target.indel_seq, best_seq).ratio()

        idx = candidate_seqs.index(best_seq)
        hit = u_candidates[idx]

        # require exact match for shiftable indels
        if require_exact_for_shiftable:
            if len(hit.generate_equivalents()) > 1 or len(target.generate_equivalents()) > 1:
                if hit != target:
                    return None

        if abs(target.pos - hit.pos) < within:

            try:
                idx2 = candidates.index(hit)  #get original representation, do not normalize here.
            except:
                hit.pos = hit.pos - len(hit.ref)
                idx2 = candidates.index(hit)

            candidate = candidates[idx2]
            idx = [i for i, var in enumerate(candidates) if var == candidate]

            # check if the candidate can be del/ins component of a complex indel
            if candidate.is_non_complex_indel():
                complex_candidates = [var for var in set(candidates) if not var.is_non_complex_indel()]
                if complex_candidates:
                    for cplx_candidate in complex_candidates:
                        reduced = cplx_candidate._reduce_complex_indel(to=target_type)
                        if candidate == reduced:
                            idx = [i for i, var in enumerate(candidates) if var == cplx_candidate]
                            candidate = reduced
                            break
            else:
                candidate = candidate._reduce_complex_indel(to=target_type)

            candidate_reads = [candidate_reads[i] for i in idx]
            candidate_ref_seqs = [candidate_ref_seqs[i] for i in idx]
            candidate_ref_starts = [candidate_ref_starts[i] for i in idx]
            candidate_aligners = [candidate_aligners[i] for i in idx]

            return candidate, candidate_reads, match_score, candidate_ref_seqs, candidate_ref_starts, candidate_aligners
        else:
            return None
    else:
        if target.is_ins and window > 3:
            window = int(window/3)
            return retarget(
                        target,
                        pileup,
                        window,
                        mapq4retarget,
                        within,
                        retargetcutoff,
                        match_score,
                        mismatch_penalty,
                        gap_open_penalty,
                        gap_extension_penalty,
                        unspl_loc_ref,
                        require_exact_for_shiftable
                    )
        else:
            return None


def update_read_info(
    read,
    candidate,
    is_gapped_aln=True,
    gap_open_penalty=3,
    gap_extension_penalty=1,
    aligner=None,
    ref_seq=None,
    ref_start=None,
):
    if is_gapped_aln:
        parsed = leftalign_indel_read(
            candidate.chrom,
            candidate.pos,
            len(candidate.indel_seq),
            candidate.variant_type,
            read["cigar_string"],
            read["read_start"],
            read["aln_start"],
            read["read_seq"],
            read["ref_seq"],
            read["read_qual"],
            candidate.reference,
        )
        read["lt_flank"] = parsed[1]
        read["indel_seq"] = parsed[2] if candidate.is_ins else ""
        read["rt_flank"] = parsed[3]
        read["lt_ref"] = parsed[4]
        read["rt_ref"] = parsed[5]
        read["lt_qual"] = parsed[6]
        read["rt_qual"] = parsed[7]

        read["lt_cigar"], read["rt_cigar"] = split_cigar(
            read["cigar_string"], candidate.pos, read["read_start"]
        )

        read["is_target"] = True
    else:
        aln = align(aligner, read["read_seq"], gap_open_penalty, gap_extension_penalty)
        genome_aln_pos = ref_start + aln.reference_start

        indels = findall_indels(
            aln, genome_aln_pos, ref_seq, read["read_seq"], basequals=read["read_qual"]
        )

        is_found= False
        for indel in indels:
            if not indel.get("del_seq", False):
                ref = indel["lt_ref"][-1]
                alt = ref + indel["indel_seq"]
            else:
                alt = indel["lt_ref"][-1]
                ref = alt + indel["del_seq"]

            obj = Variant(candidate.chrom, indel["pos"], ref, alt, candidate.reference, skip_validation=True)
            if candidate == obj:
                is_found = True
                indel_pos_in_this_read = indel["pos"]
                break

        if not is_found:
            read["cigar_updated"] = False
            return read

        read["lt_flank"] = indel["lt_flank"]
        read["indel_seq"] = candidate.indel_seq if candidate.is_ins else ""
        read["rt_flank"] = indel["rt_flank"]
        read["lt_qual"] = indel["lt_qual"]
        read["rt_qual"] = indel["rt_qual"]

        realn_lt_cigar, realn_rt_cigar = split_cigar(
            make_insertion_first(aln.CIGAR), indel["pos"], genome_aln_pos
        )

        read["lt_ref"] = trim_ref_flank(indel["lt_ref"], realn_lt_cigar, left=True)
        read["rt_ref"] = trim_ref_flank(indel["rt_ref"], realn_rt_cigar, left=False)

        read["lt_cigar"] = update_cigar(
            read["cigar_string"],
            realn_lt_cigar,
            read["read_start"],
            read["splice_pattern"],
            indel["lt_clipped"],
            left=True,
        )
        read["rt_cigar"] = update_cigar(
            read["cigar_string"],
            realn_rt_cigar,
            candidate.pos,
            read["splice_pattern"],
            indel["rt_clipped"],
            left=False,
        )

        read["cigar_list"] = read["lt_cigar"] + read["rt_cigar"]
        read["cigar_string"] = "".join(read["cigar_list"])
        read["cigar_updated"] = True

        update_read_positions(read, indel_pos_in_this_read)

        read["is_target"] = True

    return read


def trim_ref_flank(ref_flank, flank_cigar, left):

    cum = sum([int(c[:-1]) for c in flank_cigar if c[-1] != "I"])
    if left:
        ref_flank = ref_flank[-cum:]
    else:
        ref_flank = ref_flank[:cum]

    return ref_flank


def update_cigar(
    orig_cigar_string, realn_cigar, start_pos, splice_prtn, clipped_bases, left
):

    splice_ptrn = splice_prtn[0] if left else splice_prtn[1]

    if splice_ptrn:
        spl_spans = [numeric_span(spl_span) for spl_span in splice_ptrn.split(":")]
    else:
        spl_spans = []

    clip_len = len(clipped_bases)

    if left:
        new_cigar = [str(clip_len) + "S"] if clip_len else []
        current_pos = start_pos + clip_len
    else:
        new_cigar = []
        target_event = realn_cigar[0]
        target_type, target_len = target_event[-1], int(target_event[:-1])
        current_pos = (
            start_pos + 1 if target_type == "I" else start_pos + target_len + 1
        )
        trailing_clip = [str(clip_len) + "S"] if clip_len else []
        realn_cigar = realn_cigar[1:]

    for c in realn_cigar:
        event, event_len = c[-1], int(c[:-1])
        if event == "M":
            if spl_spans:
                last = len(spl_spans) - 1
                tmp = spl_spans.copy()
                for i, span in enumerate(tmp):
                    n = span[1] - span[0] + 1
                    if span[0] <= current_pos + event_len:
                        if i != last:
                            m = span[0] - current_pos

                            if m:
                                new_cigar += [str(m) + "M", str(n) + "N"]
                            else:
                                new_cigar += [str(n) + "N"]

                            current_pos += m + n
                            event_len -= m
                        else:
                            m1 = span[0] - current_pos
                            m2 = event_len - m1

                            if m2:
                                if m1:
                                    new_cigar += [str(m1) + "M", str(n) + "N", str(m2) + "M"]
                                else:
                                    new_cigar += [str(n) + "N", str(m2) + "M"]
                            else:
                                new_cigar += [str(event_len) + "M", str(n) + "N"]

                            current_pos += n + event_len
                        spl_spans = spl_spans[1:]
                    else:
                        new_cigar.append(str(event_len) + "M")
                        current_pos += event_len - 1
                        # hotfix -1 should be checked! Mar 12
                        break

            else:
                new_cigar.append(str(event_len) + "M")
                current_pos += event_len

        elif event == "I":

            if spl_spans:
                span = spl_spans[0]
                spl_start, spl_end = span[0], span[1]
                n = spl_end - spl_start + 1

                if spl_start == current_pos:
                    this_event = [str(event_len) + "I", str(n) + "N"]

                    new_cigar += this_event
                    current_pos += n

                    spl_spans = spl_spans[1:]
                else:
                    new_cigar.append(str(event_len) + "I")
                    current_pos += 1
            else:
                new_cigar.append(str(event_len) + "I")
                current_pos += 1

        elif event == "D":
            new_cigar.append(str(event_len) + "D")
            current_pos += event_len

    if left:
        return new_cigar
    else:
        return [target_event] + new_cigar + trailing_clip

    #return updated_cigar


def numeric_span(spl_span):
    spl_span_lst = spl_span.split("-")
    return [int(i) for i in spl_span_lst]


def update_read_positions(read, target_pos):

    left_adjust = sum([-int(c[:-1]) if c[-1] != "I" else 0 for c in read["lt_cigar"]])
    right_adjust = sum([int(c[:-1]) if c[-1] != "I" else 0 for c in read["rt_cigar"]])

    read["read_start"] = target_pos + left_adjust + 1
    read["read_end"] = target_pos + right_adjust

    lt_most_cigar = read["lt_cigar"][0]
    read["start_offset"] = int(lt_most_cigar[:-1]) if "S" in lt_most_cigar else 0

    rt_most_cigar = read["rt_cigar"][-1]
    read["end_offset"] = int(rt_most_cigar[:-1]) if "S" in rt_most_cigar else 0

    read["aln_start"] = read["read_start"] + read["start_offset"]
    read["aln_end"] = read["read_end"] - read["end_offset"]


def update_pileup(
    pileup,
    new_target,
    window,
    match_score,
    mismatch_penalty,
    gap_open_penalty,
    gap_extension_penalty,
    basequalthresh,
    bypass_search=False
):

    rpos = max(v.pos for v in new_target.generate_equivalents())
    for read in pileup:

        (
            is_covering,
            covering_subread,
            is_spliced,
            splice_ptrn,
            intron_ptrn,
        ) = parse_spliced_read(
            read["cigar_string"], read["read_start"], read["read_end"], new_target.pos, rpos
        )

        read["is_covering"] = is_covering
        read["covering_subread"] = covering_subread
        read["is_spliced"] = is_spliced
        read["splice_pattern"] = splice_ptrn
        read["intron_pattern"] = intron_ptrn

    if bypass_search:
        return new_target, pileup
    else:
        return find_by_normalization(
            new_target,
            pileup,
            window,
            match_score,
            mismatch_penalty,
            gap_open_penalty,
            gap_extension_penalty,
            basequalthresh,
        )
