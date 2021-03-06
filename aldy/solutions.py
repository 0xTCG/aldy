# 786
# Aldy source: solutions.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import collections
from dataclasses import dataclass, field
from natsort import natsorted
from typing import List, Dict, Optional
from .gene import Gene, Mutation


@dataclass
class CNSolution:
    r"""
    Valid copy-number configuration assignment.

    :param gene: Gene instance.
    :param score: ILP model objective score (0 for user-provided solutions).
    :param solution:
        Copy-number configurations mapped to the corresponding copy number
        (e.g. ``{1: 2}`` means that there are two copies of \*1 configuration).
    :param region_cn: Gene region copy numbers inferred by this solution.

    .. note:: Has custom printer (``__str__``).
    """

    gene: Gene
    score: float
    solution: Dict[str, int]
    region_cn: List[Dict[str, int]]

    def __init__(self, gene: Gene, score: float, solution: List[str]):
        self.region_cn = [
            {r: 0 for r in gene.regions[0]} for _ in gene.cn_configs["1"].cn
        ]
        for conf in solution:
            for gi, g in enumerate(gene.cn_configs[conf].cn):
                for r in g:
                    self.region_cn[gi][r] += g[r]
        self.gene = gene
        self.score = score
        self.solution = collections.Counter(solution)

    def __hash__(self):
        return str(self).__hash__()

    def position_cn(self, pos: int) -> float:
        """:return: Copy number at the locus ``pos``."""
        r = self.gene.region_at(pos)
        return self.region_cn[r[0]][r[1]] if r else 0

    def _solution_nice(self):
        return ",".join(f"{v}x*{k}" for k, v in natsorted(self.solution.items()))

    def __str__(self):
        return "CNSol[{:.2f}; sol=({}); cn={}]".format(
            self.score,
            self._solution_nice(),
            "|".join(
                "".join(str(self.region_cn[g][r]) for r in self.gene.regions[0])
                for g, _ in enumerate(self.region_cn)
            ),
        )

    def max_cn(self):
        """Maximum copy-number in this solution"""
        return sum(self.solution.values())


@dataclass
class SolvedAllele:
    """
    Solved star-allele assignment.

    :param gene: Gene instance.
    :param major: Major star-allele identifier.
    :param minor: Minor star-allele identifier. Can be None.
    :param added:
        Mutations that are added to the star-allele
        (present in the allele database definition).
    :param missing:
        Mutations that are omitted from the star-allele
        (present in the allele database definition but not here).

    .. note:: Has custom printer (``__str__``).
    """

    gene: Gene
    major: str
    minor: Optional[str] = None
    added: List[Mutation] = field(default_factory=list)
    missing: List[Mutation] = field(default_factory=list)

    def mutations(self):
        m = self.gene.alleles[self.major].func_muts
        m |= self.gene.alleles[self.major].minors[self.minor].neutral_muts
        m |= set(self.added)
        m -= set(self.missing)
        return m

    def major_repr(self):
        """
        Pretty-formats major star-allele name.
        """
        return "*{}{}".format(
            self.major,
            "".join(
                " +" + self.gene.get_rsid(m)
                for m in sorted(
                    m for m in self.added if self.gene.is_functional(m, False)
                )
            ),
        )

    def __str__(self):
        """
        Pretty-formats minor star-allele name.
        """
        return "*{}{}{}".format(
            self.minor if self.minor else self.major,
            "".join(" +" + self.gene.get_rsid(m) for m in sorted(self.added)),
            "".join(" -" + self.gene.get_rsid(m) for m in sorted(self.missing)),
        )

    def __hash__(self):
        return hash(
            (self.gene.name, self.major, self.minor, *self.added, *self.missing)
        )


@dataclass
class MajorSolution:
    r"""
    Valid major star-allele assignment.

    :param score: ILP model objective score.
    :param solution:
        Major star-alleles and the corresponding copy numbers
        (e.g. ``{1: 2}`` means that we have two copies of \*1).
    :param cn_solution: Copy-number solution that was used to assign major star-alleles.
    :param added:
        List of added mutations. Will be assigned to :obj:`SolvedAllele` in the
        minor star-allele calling step if phasing is enabled.

    .. note:: Has custom printer (``__str__``).
    """

    score: float
    solution: Dict[SolvedAllele, int]
    cn_solution: CNSolution
    added: List[Mutation]

    def _solution_nice(self):
        x = ", ".join(
            f"{v}x{s}"
            for s, v in natsorted(self.solution.items(), key=lambda x: x[0].major)
        )
        y = ", ".join(self.cn_solution.gene.get_rsid(m) for m in sorted(self.added))
        return " & ".join([x, y] if y else [x])

    def __str__(self):
        return (
            f"MajorSol[{self.score:.2f}; "
            + f"sol=({self._solution_nice()}); "
            + f"cn={self.cn_solution}"
        )

    def __hash__(self):
        return str(self).__hash__()


@dataclass
class MinorSolution:
    """
    Valid minor star-allele assignment.

    :param score: ILP model objective score.
    :param solution:
        List of solved minor star-alleles.
        Modifications to minor alleles are represented in :obj:`SolvedAllele` format.
    :param major_solution:
        Major star-allele solution used for calculating minor star-allele assignments.

    .. note:: Has custom printer (``__str__``).
    """

    score: float
    solution: List[SolvedAllele]
    major_solution: MajorSolution

    def _solution_nice(self):
        return ", ".join(
            str(s) for s in natsorted(self.solution, key=lambda x: x.minor)
        )

    def __str__(self):
        return (
            f"MinorSol[{self.score:.2f}; "
            + f"sol=({self._solution_nice()}); "
            + f"major={self.major_solution}"
        )

    def get_diplotype(self):
        try:
            return self.diplotype
        except AttributeError:
            self.diplotype = ([], [])
            return self.diplotype

    def set_diplotype(self, d):
        self.diplotype = d

    def get_major_name(self, i):
        gene = self.major_solution.cn_solution.gene
        if i == -1:
            return gene.deletion_allele()
        n = str(self.solution[i].major).split("#")[0:1]
        for m in sorted(self.solution[i].added):
            if gene.is_functional(m, infer=False):
                n.append(gene.get_rsid(m))
        return "+".join(n)

    def get_minor_name(self, i, legacy=False):
        gene = self.major_solution.cn_solution.gene
        if i == -1:
            return gene.deletion_allele()

        m = gene.alleles[self.solution[i].major]
        t = [mi for mi in m.minors if mi == self.solution[i].minor]
        assert len(t) == 1
        n = [m.minors[t[0]].alt_name if legacy and m.minors[t[0]].alt_name else t[0]]
        for m in sorted(self.solution[i].added):
            n.append("+" + gene.get_rsid(m))
        for m in sorted(self.solution[i].missing):
            n.append("-" + gene.get_rsid(m))
        return " ".join(n)

    def get_major_diplotype(self):
        return " / ".join(
            " + ".join(f"*{self.get_major_name(i)}" for i in d) for d in self.diplotype
        )

    def get_minor_diplotype(self, legacy=False):
        return " / ".join(
            " + ".join(f"[*{self.get_minor_name(i, legacy)}]" for i in d)
            for d in self.diplotype
        )
