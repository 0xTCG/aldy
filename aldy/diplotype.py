# 786
# Aldy source: diplotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Tuple, Dict

import collections
import re

from natsort import natsorted
from .common import allele_sort_key, td
from .gene import Gene, Mutation
from .coverage import Coverage
from .solutions import MinorSolution
from .version import __version__ as version


OUTPUT_COLS = [
    "Sample",
    "Gene",
    "SolutionID",
    "Major",
    "Minor",
    "Copy",
    "Allele",
    "Location",
    "Type",
    "Coverage",
    "Effect",
    "dbSNP",
    "Code",
    "Status",
]
"""list[str]: Output column descriptions"""


def write_decomposition(
    sample: str, gene: Gene, sol_id: int, minor: MinorSolution, f
) -> None:
    """
    Write an allelic decomposition to the file `f`.

    Args:
        sample (str): Sample name.
        gene (:obj:`aldy.gene.Gene`): Gene instance.
        sol_id (int): Solution ID (each solution must have a different ID).
        minor (:obj:`aldy.solutions.MinorSolution`):
            Minor star-allele solution to be written.
        f (file): Output file.
    """

    for copy, a in enumerate(minor.solution):
        mutations = set(gene.alleles[a.major].func_muts) | set(
            gene.alleles[a.major].minors[a.minor].neutral_muts
        )
        mutations |= set(a.added)
        mutations -= set(a.missing)
        items = []
        if len(mutations) > 0:
            for m in sorted(mutations):
                items.append(
                    [
                        sample,
                        gene.name,
                        sol_id,
                        minor.diplotype,
                        ";".join(ay.minor for ay in minor.solution),
                        copy,
                        a.minor,
                        m.pos,
                        m.op,
                        -1,
                        ["NEUTRAL", "DISRUPTING"][gene.is_functional(m)],
                        gene.get_dbsnp(m),
                        "",
                    ]
                )
        else:
            items.append(
                [
                    sample,
                    gene.name,
                    sol_id,
                    minor.diplotype,
                    ";".join(ay.minor for ay in minor.solution),
                    copy,
                    a.minor,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                ]
            )
        for it in items:
            print("\t".join(str(i) for i in it), file=f)


def write_vcf(
    sample: str, gene: Gene, coverage: Coverage, minors: List[MinorSolution], f
):
    header = f"""
    ##fileformat=VCFv4.2
    ##source=aldy-v{version}
    ##INFO=<ID=ANN,Number=1,Type=String,Description="Location within {gene.name}">
    ##INFO=<ID=TYPE,Number=1,Type=String,Description="Mutation kind">
    ##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=MA,Number=1,Type=String,Description="Major genotype star-allele calls">
    ##FORMAT=<ID=MI,Number=1,Type=String,Description="Minor genotype star-allele calls">
    """
    all_mutations: Dict[Mutation, List[dict]] = {
        m: [collections.defaultdict(int)] * len(minors)
        for minor in minors
        for a in minor.solution
        for m in set(gene.alleles[a.major].func_muts)
        | set(gene.alleles[a.major].minors[a.minor].neutral_muts)
        | set(a.added)
    }
    for mi, minor in enumerate(minors):
        for ai, a in enumerate(minor.solution):
            mutations = set(gene.alleles[a.major].func_muts) | set(
                gene.alleles[a.major].minors[a.minor].neutral_muts
            )
            mutations |= set(a.added)
            for m in mutations:
                all_mutations[m][mi][ai] = 1
    print(
        td(header).strip()
        + "\n"
        + "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + [f"{sample}:{mi}:{m.diplotype}" for mi, m in enumerate(minors)]
        ),
        file=f,
    )
    for m in sorted(mutations):
        ref = gene.seq[m.pos - gene.region.start]
        if m.op[:3] == "SNP":
            alt = m.op[5]
        elif m.op[:3] == "INS":
            alt = ref + m.op[4:]
        else:
            ref = gene.seq[m.pos - 1 - gene.region.start]
            alt = ref + m.op[4:], ref

        info = [
            "TYPE={}".format(["NEUTRAL", "DISRUPTING"][gene.is_functional(m)]),
            "GENE=" + gene.name,
        ]
        data = []
        for mi, minor in enumerate(minors):
            nall = len(minor.solution)
            data.append(
                {
                    "GT": "|".join(str(all_mutations[m][mi][i]) for i in range(nall)),
                    "DP": str(coverage[m]),
                    "MA": ",".join(
                        f"*{minor.solution[i].major}"
                        if all_mutations[m][mi][i] > 0
                        else "-"
                        for i in range(nall)
                    ),
                    "MI": ",".join(
                        f"*{minor.solution[i].minor}"
                        if all_mutations[m][mi][i] > 0
                        else "-"
                        for i in range(nall)
                    ),
                }
            )
        pattern = (
            "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t"
            + "{format}\t{data}"
        )
        print(
            pattern.format(
                chrom=gene.region.chr,
                pos=m.pos + 1,
                id=gene.get_dbsnp(m),
                ref=ref,
                alt=alt,
                qual=0,
                filter="PASS",
                info=";".join(info),
                format=":".join(data[0].keys()),
                data="\t".join(":".join(d.values()) for d in data),
            ),
            file=f,
        )


def estimate_diplotype(gene: Gene, solution: MinorSolution) -> str:
    """
    Calculate the diplotype assignment for a minor solution.
    Assigns a ``diplotype`` attribute of the :obj:`aldy.solutions.MinorSolution`.
    Relies on the diplotype assignment heuristics to assign correct diplotypes.
    This heuristics has no biological validity whatsoever--- it is purely used
    to pretty-print the final solutions.

    Returns:
        str: Diplotype assignment.
    """

    del_allele = gene.deletion_allele()

    # solution is the array of (major, minor) tuples
    major_dict: Dict[str, List[str]] = collections.defaultdict(list)
    for a in solution.solution:
        n = str(a.major).split("#")[0:1]  # chop off fusion suffix
        for m in a.added:
            if gene.is_functional(m, infer=False):
                n.append(gene.get_dbsnp(m))
        # get "real" name
        components = re.split(r"(\d+)", n[0])
        real = components[0] if components[0] != "" else components[1]
        major_dict[real].append("+".join(n))
    if len(solution.solution) == 1 and del_allele:
        major_dict[del_allele].append(del_allele)

    diplotype: Tuple[List[str], List[str]] = ([], [])
    dc = 0

    # Handle tandems (heuristic that groups common tandems together,
    #                 e.g. 1, 2, 13 -> 13+1/2 if [13,1] is a common tandem)
    if len(solution.solution) > 2:
        for ta, tb in gene.common_tandems:
            while major_dict[ta] and major_dict[tb]:
                diplotype[dc % 2].append(f"{major_dict[ta][0]}+{major_dict[tb][0]}")
                dc += 1
                del major_dict[ta][0], major_dict[ta][1]

    # Handle duplicates (heuristics that groups duplicate alleles together,
    #                    e.g. 1, 1, 2 -> 1+1/2)
    # First check should we split them (e.g. 1, 1, 1, 1 -> 1+1/1+1)?
    if len(major_dict) == 1:
        items = next(iter(major_dict.values()))
        if len(items) % 2 == 0:
            diplotype = items[: len(items) // 2], items[len(items) // 2 :]
            major_dict.clear()
    for allele, items in major_dict.items():
        if len(items) > 1:
            if len(diplotype[dc % 2]) > len(diplotype[(dc + 1) % 2]):
                dc += 1
            diplotype[dc % 2].extend(items)
            items.clear()
            dc += 1

    # Handle the rest
    for allele, items in major_dict.items():
        if items:
            if len(diplotype[dc % 2]) > len(diplotype[(dc + 1) % 2]):
                dc += 1
            assert len(items) == 1
            diplotype[dc % 2].extend(items)
            dc += 1

    # Each diplotype should have at least one item
    # e.g. 1, 1 -> becomes 1+1/_ due to duplicate heuristic -> fixed to 1/1
    if len(diplotype[1]) == 0:
        diplotype = diplotype[0][:-1], [diplotype[0][-1]]

    # Make sure that the elements are sorted and that the tandems are grouped together
    result = natsorted([natsorted(diplotype[0]), natsorted(diplotype[1])])
    res = " / ".join(" + ".join("*{}".format(y) for y in x) for x in result)
    solution.diplotype = res
    return res
