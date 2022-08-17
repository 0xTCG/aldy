# 786
# Aldy source: solutions.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict
from dataclasses import dataclass, field
from natsort import natsorted
from collections import Counter, defaultdict

from .gene import Gene, Mutation


@dataclass
class CNSolution:
    """Valid copy-number configuration assignment."""

    gene: Gene
    """Gene instance."""
    score: float
    """ILP model objective score (0 for user-provided solutions)."""
    solution: Dict[str, int]
    """
    Copy-number configurations mapped to the corresponding copy number
    (e.g. `{1: 2}` means that there are two copies of `*1` configuration).
    """
    region_cn: List[Dict[str, int]]
    """Gene region copy numbers inferred by this solution."""

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
        self.solution = Counter(solution)

    def __hash__(self):
        return str(self).__hash__()

    def position_cn(self, pos: int) -> float:
        """:returns: Copy number at the given location."""
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
        """Maximum copy-number in the solution."""
        return sum(self.solution.values())


@dataclass
class SolvedAllele:
    """Valid star-allele assignment."""

    gene: Gene
    """Gene instance."""
    major: str
    """Major star-allele identifier."""
    minor: str = ""
    """Minor star-allele identifier. Empty string when not assigned."""
    added: List[Mutation] = field(default_factory=list)
    """Mutations that are added to the star-allele (from the allele database)."""
    missing: List[Mutation] = field(default_factory=list)
    """
    Mutations that are omitted from the star-allele (but present in the star-allele
    definition).
    """

    def mutations(self):
        """Set of allele mutations."""
        m = self.gene.alleles[self.major].func_muts
        if self.minor:
            m |= self.gene.alleles[self.major].minors[self.minor].neutral_muts
        m |= set(self.added)
        m -= set(self.missing)
        return m

    def major_repr(self):
        """Pretty-formats major star-allele name."""
        extra = "".join(
            " +" + self.gene.get_rsid(m)
            for m in sorted(m for m in self.added if self.gene.is_functional(m, False))
        )
        s = f"{self.major}{extra}"
        return f"*({s})" if extra else f"*{s}"

    def __str__(self):
        """Pretty-formats minor star-allele name."""
        extra = "".join(" +" + self.gene.get_rsid(m) for m in sorted(self.added))
        miss = "".join(" -" + self.gene.get_rsid(m) for m in sorted(self.missing))
        s = f"{self.minor if self.minor else self.major}{extra}{miss}"
        return f"*({s})" if extra or miss else f"*{s}"

    def __hash__(self):
        return hash(
            (self.gene.name, self.major, self.minor, *self.added, *self.missing)
        )


@dataclass
class MajorSolution:
    """Valid major star-allele assignment."""

    score: float
    """ILP model objective score."""
    solution: Dict[SolvedAllele, int]
    """
    Major star-alleles and the corresponding copy numbers
    (e.g. `{1: 2}` means that we have two copies of `*1`).
    """
    cn_solution: CNSolution
    """Copy-number solution that was used to assign major star-alleles."""
    added: List[Mutation]
    """
    List of added mutations. Will be assigned to :py:class:`SolvedAllele` in the
    minor star-allele calling step if phasing is enabled.
    """

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
    """Valid minor star-allele assignment."""

    score: float
    """ILP model objective score."""
    solution: List[SolvedAllele]
    """
    List of solved minor star-alleles.
    Modifications to minor alleles are represented in :py:class:`SolvedAllele`.
    """
    major_solution: MajorSolution
    """Major star-allele solution used for calculating minor star-allele assignments."""

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
        n = [m.minors[t[0]].alt_name or t[0] if legacy else t[0]]
        for m in sorted(self.solution[i].added):
            n.append("+" + gene.get_rsid(m))
        for m in sorted(self.solution[i].missing):
            n.append("-" + gene.get_rsid(m))
        return " ".join(n)

    def get_major_diplotype(self):
        return " / ".join(
            " + ".join(f"*{self.get_major_name(i)}" for i in d)
            for d in self.diplotype
            if d
        )

    def get_minor_diplotype(self, legacy=False):
        return " / ".join(
            " + ".join(f"[*{self.get_minor_name(i, legacy)}]" for i in d)
            for d in self.diplotype
            if d
        )

    def get_mutation_coverages(self, coverage):
        muts = defaultdict(int)
        for sa in self.solution:
            for m in (
                sa.gene.alleles[sa.major].func_muts
                | sa.gene.alleles[sa.major].minors[sa.minor].neutral_muts
            ):
                if m not in sa.missing:
                    muts[m] += 1
            for m in sa.added:
                muts[m] += 1
        covs = {}
        for m in muts:
            covs[m] = coverage.single_copy(m.pos, self.major_solution.cn_solution)
            covs[m] = coverage[m] / covs[m] if covs[m] > 0 else 0
        return [(m, muts[m], covs[m]) for m in muts]
