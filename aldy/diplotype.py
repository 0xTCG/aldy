# 786
# Aldy source: diplotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict, Any
from natsort import natsorted
import collections
import re

from .common import td
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
"""Output column descriptions"""


def write_decomposition(
    sample: str, gene: Gene, coverage: Coverage, sol_id: int, minor: MinorSolution, f
):
    """Write an allelic decomposition to the given file.

    :param sample: Sample name.
    :param gene: Gene instance.
    :param sol_id: Solution ID (each solution must have a different ID).
    :param minor: Minor star-allele solution to be written.
    :param f: Output file.
    """

    for copy, a in enumerate(minor.solution):
        assert a.minor
        mutations = set(gene.alleles[a.major].func_muts) | set(
            gene.alleles[a.major].minors[a.minor].neutral_muts
        )
        mutations |= set(a.added)
        mutations -= set(a.missing)
        items = []
        if len(mutations) > 0:
            for m in sorted(mutations):
                fn = gene.get_functional(m, False)
                items.append(
                    [
                        sample,
                        gene.name,
                        sol_id,
                        minor.get_major_diplotype().replace(" ", ""),
                        ";".join(ay.minor for ay in minor.solution if ay.minor),
                        copy,
                        a.minor,
                        m.pos,
                        m.op,
                        coverage[m],
                        fn if fn else "none",
                        gene.get_rsid(m, default=False),
                        "",
                    ]
                )
        else:
            items.append(
                [
                    sample,
                    gene.name,
                    sol_id,
                    minor.get_major_diplotype().replace(" ", ""),
                    ";".join(ay.minor for ay in minor.solution if ay.minor),
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
    """
    Write an allelic decomposition in the VCF format to the given file.

    :param sample: Sample name.
    :param gene: Gene instance.
    :param sol_id: Solution ID (each solution must have a different ID).
    :param minor: Minor star-allele solution to be written.
    :param f: Output file.
    """

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
        if a.minor
        for m in set(gene.alleles[a.major].func_muts)
        | set(gene.alleles[a.major].minors[a.minor].neutral_muts)
        | set(a.added)
    }
    for mi, minor in enumerate(minors):
        for ai, a in enumerate(minor.solution):
            assert a.minor
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
            + [
                f"{sample}:{mi}:{m.get_major_diplotype().replace(' ', '')}"
                for mi, m in enumerate(minors)
            ]
        ),
        file=f,
    )
    for m in sorted(all_mutations):
        ref = m.op[0]  # TODO: should be genome nucleotide, not the RefSeq nucleotide?!
        if m.op[1] == ">":
            alt = m.op[2]
        elif m.op[:3] == "ins":
            alt = ref + m.op[3:]
        else:
            ref, alt = ".", f"{m.op[3:]}, ."  # TODO: this is wrong?!

        fm = gene.get_functional(m)
        fm = fm.replace(" ", "_").replace("\t", "_").replace(";", "_") if fm else "none"
        info = [
            f"EFFECT={fm}",
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
                chrom=gene.chr,
                pos=m.pos + 1,
                id=gene.get_rsid(m, default=False),
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
    """Calculate the diplotype assignment for a minor solution.
    Set the `diplotype` attribute of the :py:class:`aldy.solutions.MinorSolution`.

    Uses the diplotype assignment heuristics to assign correct diplotypes.
    This heuristics has no biological validity whatsoever---it is purely used
    for pretty-printing the final solutions.

    :returns: Diplotype assignment.
    """

    del_allele = gene.deletion_allele()

    # solution is the array of (major, minor) tuples
    major_dict: Dict[str, List[int]] = collections.defaultdict(list)
    for ai, a in enumerate(solution.solution):
        n = str(a.major).split("#")[0]  # chop off fusion suffix
        components = re.split(r"(\d+)", n)
        real = components[0] if components[0] != "" else components[1]
        major_dict[real].append(ai)
    if del_allele:
        if len(solution.solution) == 0:
            major_dict[del_allele].append(-1)
            major_dict[del_allele].append(-1)
        elif len(solution.solution) == 1:
            major_dict[del_allele].append(-1)

    diplotype: Any = [[], []]
    dc = 0

    # Handle tandems (heuristic that groups common tandems together,
    #                 e.g. 1, 2, 13 -> 13+1/2 if [13,1] is a common tandem)
    if len(solution.solution) > 2:
        for ta, tb in gene.common_tandems:
            while major_dict[ta] and major_dict[tb]:
                diplotype[dc % 2].append((major_dict[ta][0], major_dict[tb][0]))
                dc += 1
                del major_dict[ta][0], major_dict[tb][0]

    # Handle duplicates (heuristics that groups duplicate alleles together,
    #                    e.g. 1, 1, 2 -> 1+1/2)
    # First check should we split them (e.g. 1, 1, 1, 1 -> 1+1/1+1)?
    xlen = lambda d: sum(2 if isinstance(n, tuple) else 1 for n in d)  # noqa

    if len(major_dict) == 1:
        items = next(iter(major_dict.values()))
        if len(items) % 2 == 0:
            diplotype[dc % 2] += items[: len(items) // 2]
            dc += 1
            diplotype[dc % 2] += items[len(items) // 2 :]
            dc += 1
            major_dict.clear()
    for _, items in major_dict.items():
        if len(items) > 1:
            if xlen(diplotype[dc % 2]) > xlen(diplotype[(dc + 1) % 2]):
                dc += 1
            diplotype[dc % 2] += items
            items.clear()
            dc += 1

    # Handle the rest
    for _, items in major_dict.items():
        if items:
            if xlen(diplotype[dc % 2]) > xlen(diplotype[(dc + 1) % 2]):
                dc += 1
            assert len(items) == 1
            diplotype[dc % 2] += items
            dc += 1

    # Each diplotype should have at least one item
    # e.g. 1, 1 -> becomes 1+1/_ due to duplicate heuristic -> fixed to 1/1
    if len(diplotype[1]) == 0:
        if len(diplotype[0]) > 1:
            diplotype = diplotype[0][:-1], [diplotype[0][-1]]
        elif len(diplotype[0]) == 1 and isinstance(diplotype[0][0], tuple):
            diplotype = diplotype[0][0]

    def flatten(d):
        def key(x):
            if isinstance(x, tuple):
                return solution.get_major_name(x[0])
            return solution.get_major_name(x)

        for i in natsorted(d, key=key):
            if isinstance(i, tuple):
                yield from i
            else:
                yield i

    # Make sure that the elements are sorted and that the tandems are grouped together
    diplotype = natsorted(
        [list(flatten(diplotype[0])), list(flatten(diplotype[1]))],
        key=lambda x: [solution.get_major_name(y) for y in x],
    )
    solution.set_diplotype(diplotype)
    return diplotype
