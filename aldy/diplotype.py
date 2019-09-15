# 786
# Aldy source: diplotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Tuple

import collections

from .common import allele_sort_key, td
from .gene import Gene
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
                        ["NEUTRAL", "DISRUPTING"][bool(m.is_functional)],
                        m.aux.get("dbsnp", ""),
                        m.aux.get("old", ""),
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
    sample: str,
    gene: Gene,
    coverage: Coverage,
    minors: List[MinorSolution],
    f,
):
    header = f"""
    ##fileformat=VCFv4.2
    ##source=aldy-v{version}
    ##INFO=<ID=ANN,Number=1,Type=String,Description="Location within {gene.name}">
    ##INFO=<ID=TYPE,Number=1,Type=String,Description="Mutation kind">
    ##INFO=<ID=TYPE,Number=1,Type=String,Description="Gene">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=MA,Number=1,Type=String,Description="Major genotype star-allele calls">
    ##FORMAT=<ID=MI,Number=1,Type=String,Description="Minor genotype star-allele calls">
    """
    minors=[minors[0], minors[0]]
    all_mutations = {
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
        td(header).strip() + "\n"
        + "\t".join(
            [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
            ]
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
            "ANN={}".format(m.aux.get("old", ".")),
            "TYPE={}".format(["NEUTRAL", "DISRUPTING"][bool(m.is_functional)]),
            "GENE=" + gene.name
        ]
        data = []
        for mi, minor in enumerate(minors):
            nall = len(minor.solution)
            data.append(
                {
                    "GT": "|".join(str(all_mutations[m][mi][i]) for i in range(nall)),
                    "DP": str(coverage[m]),
                    "MA": "|".join(
                        f'*{minor.solution[i].major}' if all_mutations[m][mi][i] > 0 else "-"
                        for i in range(nall)
                    ),
                    "MI": "|".join(
                        f'*{minor.solution[i].minor}' if all_mutations[m][mi][i] > 0 else "-"
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
                id=m.aux.get("dbsnp", "."),
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
    majors = [
        str(allele_sort_key(a.major)[0])
        + ("-like" if sum(1 for m in a.added if m.is_functional) > 0 else "")
        for a in solution.solution
    ]
    diplotype: Tuple[List[str], List[str]] = ([], [])

    if len(majors) == 1 and del_allele:
        majors.append(del_allele)

    major_dict = collections.Counter(majors)
    dc = 0

    # Handle tandems (heuristic that groups common tandems together,
    #                 e.g. 1, 2, 13 -> 1+13/2 if [1,13] is a common tandem)
    for ta, tb in gene.common_tandems:
        while major_dict[ta] > 0 and major_dict[tb] > 0:
            diplotype[dc % 2].extend([ta, tb])
            dc += 1
            major_dict[ta] -= 1
            major_dict[tb] -= 1

    # Handle duplicates (heuristics that groups duplicate alleles together,
    #                    e.g. 1, 1, 2 -> 1+1/2)
    # First check should we split them (e.g. 1, 1, 1, 1 -> 1+1/1+1)?
    el = list(major_dict.elements())
    if len(major_dict) == 1 and len(el) % 2 == 0:
        p, a = len(el) // 2, el[0]
        diplotype = (p * [a], p * [a])
        major_dict = collections.Counter()
    for allele, count in major_dict.items():
        if count > 1:
            if len(diplotype[dc % 2]) > len(diplotype[(dc + 1) % 2]):
                dc += 1
            diplotype[dc % 2].extend(count * [str(allele)])
            major_dict[allele] -= count
            dc += 1

    # Handle the rest
    for allele, count in major_dict.items():
        if count > 0:
            if len(diplotype[dc % 2]) > len(diplotype[(dc + 1) % 2]):
                dc += 1
            assert count == 1
            diplotype[dc % 2].append(str(allele))
            dc += 1

    # Each diplotype should have at least one item
    # e.g. 1, 1 -> becomes 1+1/_ due to duplicate heuristic -> fixed to 1/1
    if len(diplotype[1]) == 0:
        diplotype = (diplotype[0][:-1], [diplotype[0][-1]])

    # Make sure that the elements are sorted and that the tandems are grouped together
    result = sorted(diplotype)
    for i in range(2):
        nd: List[tuple] = []
        for ta, tb in gene.common_tandems:
            if ta in result[i] and tb in result[i]:
                nd.append((ta, tb))
                result[i].remove(ta)
                result[i].remove(tb)
        nd += [(x,) for x in result[i]]
        result[i] = [f for e in sorted(nd) for f in e]

    res = "/".join("+".join("*{}".format(y) for y in x) for x in result)
    solution.diplotype = res
    return res
