# 786
# Aldy source: solutions.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import collections
from natsort import natsorted
from typing import List
from .gene import Gene


class CNSolution(
    collections.namedtuple("CNSolution", ["gene", "score", "solution", "region_cn"])
):
    r"""
    Valid copy-number configuration assignment.
    Immutable.

    Attributes:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        score (float):
            ILP model objective score (0 for user-provided solutions).
        solution (dict[str, int]):
            Copy-number configurations mapped to the corresponding copy number
            (e.g. ``{1: 2}`` means that there are two copies of \*1 configuration).
        region_cn (list[dict[str, int]]):
            Gene region copy numbers inferred by this solution.

    Notes:
        Has custom printer (``__str__``).
    """

    def __new__(self, gene: Gene, score: float, solution: List[str]):
        vec = [{r: 0 for r in gene.regions[0]} for _ in gene.cn_configs["1"].cn]
        for conf in solution:
            for gi, g in enumerate(gene.cn_configs[conf].cn):
                for r in g:
                    vec[gi][r] += g[r]
        return super(CNSolution, self).__new__(  # type:ignore
            self, gene, score, collections.Counter(solution), vec,
        )

    def __hash__(self):
        return str(self).__hash__()

    def position_cn(self, pos: int) -> float:
        """
        Returns:
            float: Copy number at the locus ``pos``.
        """
        try:
            g, region = self.gene.region_at(pos)
            return self.region_cn[g][region]
        except KeyError:
            return 0

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


class SolvedAllele(
    collections.namedtuple(
        "SolvedAllele", ["gene", "major", "minor", "added", "missing"]
    )
):
    """
    Solved star-allele assignment.
    Immutable.

    Attributes:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        major (str):
            Major star-allele identifier.
        minor (str, optional):
            Minor star-allele identifier. Can be None.
        added (tuple[:obj:`aldy.gene.Mutation`]):
            Mutations that are added to the star-allele
            (i.e. these mutations are not present in the allele database definition).
        missing (tuple[:obj:`aldy.gene.Mutation`]):
            Mutations that are omitted from the star-allele
            (i.e. these mutations are present in the allele database definition
            but not here).

    Notes:
        Has custom printer (``__str__``).
    """

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
                " +" + self.gene.get_dbsnp(m)
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
            "".join(" +" + self.gene.get_dbsnp(m) for m in sorted(self.added)),
            "".join(" -" + self.gene.get_dbsnp(m) for m in sorted(self.missing)),
        )


class MajorSolution(
    collections.namedtuple(
        "MajorSolution", ["score", "solution", "cn_solution", "added"]
    )
):
    r"""
    Valid major star-allele assignment.
    Immutable.

    Attributes:
        score (float):
            ILP model objective score.
        solution (dict[:obj:`SolvedAllele`, int]):
            Major star-alleles and the corresponding copy numbers
            (e.g. ``{1: 2}`` means that we have two copies of \*1).
        cn_solution (:obj:`aldy.solutions.CNSolution`):
            Copy-number solution that was used to assign major star-alleles.
        added (list[:obj:`Mutation`]):
            List of added mutations. Will be assigned to :obj:`SolvedAllele` in the
            minor star-allele calling step if phasing is enabled.

    Notes:
        Has custom printer (``__str__``).
    """

    def _solution_nice(self):
        x = ", ".join(
            f"{v}x{s}"
            for s, v in natsorted(self.solution.items(), key=lambda x: x[0].major)
        )
        y = ", ".join(self.cn_solution.gene.get_dbsnp(m) for m in self.added)
        return " & ".join([x, y] if y else [x])

    def __str__(self):
        return (
            f"MajorSol[{self.score:.2f}; "
            + f"sol=({self._solution_nice()}); "
            + f"cn={self.cn_solution}"
        )

    def __hash__(self):
        return str(self).__hash__()


class MinorSolution(
    collections.namedtuple("MinorSolution", ["score", "solution", "major_solution"])
):
    """
    Valid minor star-allele assignment.
    Immutable.

    Attributes:
        score (float):
            ILP model objective score.
        solution (list[:obj:`SolvedAllele`]):
            List of solved minor star-alleles.
            Modifications to the minor alleles are represented in
            :obj:`SolvedAllele` format.
        major_solution (:obj:`aldy.solutions.MajorSolution`):
            Major star-allele solution used for calculating the
            minor star-allele assignment.
        diplotype (tuple[list[int], list[int]]):
            Diplotype assignment

    Notes:
        Has custom printer (``__str__``).
    """

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
        for m in self.solution[i].added:
            if gene.is_functional(m, infer=False):
                n.append(gene.get_dbsnp(m))
        return "+".join(n)

    def get_minor_name(self, i, legacy=False):
        gene = self.major_solution.cn_solution.gene
        if i == -1:
            return gene.deletion_allele()

        m = gene.alleles[self.solution[i].major]
        t = [mi for mi in m.minors if mi == self.solution[i].minor]
        assert len(t) == 1
        n = [m.minors[t[0]].alt_name if legacy and m.minors[t[0]].alt_name else t[0]]
        for m in self.solution[i].added:
            n.append("+" + gene.get_dbsnp(m))
        for m in self.solution[i].missing:
            n.append("-" + gene.get_dbsnp(m))
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
