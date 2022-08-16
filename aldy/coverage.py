# 786
# Aldy source: coverage.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, Tuple, Callable, List, Any
import copy

from .profile import Profile
from .common import log, AldyException
from .gene import Mutation, Gene
from .solutions import CNSolution


class Coverage:
    """Data structure that maintains the coverage information for a given sample."""

    def __init__(
        self,
        gene: Gene,
        profile: Profile,
        sam,
        coverage: Dict[int, Dict[str, List]],
        cnv_coverage: Dict[int, int],
    ) -> None:
        """
        :param gene: Gene instance.
        :param profile: Profile instance.
        :param sam: Sample instance.
        :param coverage: Coverage for each sample location. Each location is represented
            as a dictionary that maps a mutation (or a reference position indicated by
            `_`) to the list of read quality scores that cover it. For example,
            `coverage[10]['A>G'] = [(10, 20), (10, 10)]` indicates that 2 reads have G
            (instead of A) at the location 10.
        :param cnv_coverage: Coverage of the copy-number neutral region within the
            sample. Each location is represented by the total corresponsing read
            coverage. Used for coverage rescaling.
        """
        self.gene = gene
        self.profile = profile
        self.sam = sam
        self._coverage = coverage
        self._cnv_coverage = cnv_coverage
        self._region_coverage: Dict[Tuple[int, str], float] = {}

    def __getitem__(self, mut: Mutation) -> float:
        """:returns: Mutation coverage."""
        return self.coverage(mut)

    def coverage(self, mut: Mutation) -> float:
        """:returns: Mutation coverage."""
        if mut.pos in self._coverage and mut.op in self._coverage[mut.pos]:
            return len(self._coverage[mut.pos][mut.op])
        else:
            return 0

    def total(self, pos: int) -> float:
        """:returns: Location coverage."""
        if pos not in self._coverage:
            return 0
        return float(
            sum(len(v) for p, v in self._coverage[pos].items() if p[:3] != "ins")
        )

    def percentage(self, m: Mutation) -> float:
        """:returns: Mutation coverage expressed as percentage (0-100%)."""
        total = self.total(m.pos)
        if total == 0:
            return 0
        return 100.0 * self.coverage(m) / total

    def single_copy(self, pos: int, cn_solution: CNSolution) -> float:
        """
        :param pos: Genomic locus.
        :param cn_solution: Copy-number solution.
        :returns: Coverage of a *single* gene copy at the given location.
        """
        if cn_solution.position_cn(pos) == 0:
            return 0
        return max(1, self.total(pos)) / cn_solution.position_cn(pos)

    def region_coverage(self, gene: int, region: str) -> float:
        """:returns: Average coverage of a gene region."""
        return self._region_coverage[gene, region]

    def average_coverage(self) -> float:
        """:returns: Average coverage of the gene."""
        return sum(self.total(pos) for pos in self._coverage) / float(
            len(self._coverage) + 0.1
        )

    def dump(self, out=None):
        """Pretty-print the coverage data."""
        for pos, pos_mut in sorted(self._coverage.items()):
            if len(pos_mut) == 1 and "_" in pos_mut:
                continue
            for _, (op, _) in enumerate(sorted(pos_mut.items(), reverse=True)):
                p = self.percentage(Mutation(pos, op))
                if pos in self.gene:
                    x = self.gene.get_functional((pos, op))
                    if op == "_":
                        x = ""
                    if x and (pos, op) not in self.gene.mutations:
                        x += "**"
                    r = self.gene.region_at(pos)
                    t = f"{x if x else ''}\t{r[1] if r else ''}"
                    t += "\t" + self.gene.get_rsid((pos, op))
                else:
                    t = "\t\t"
                if out and self.sam:
                    out(
                        f"[dump] {self.sam.name}\t{self.gene.name}\t"
                        f"{self.gene.chr_to_ref.get(pos, -1) + 1}\t{op}\t{p:.1f}\t{t}"
                    )

    def filtered(self, filter_fn: Callable[[Any, Mutation], List]):
        """
        :param filter_fn: Function that performs mutation filtering with the following
            arguments:

                1. mut (:py:class:`aldy.gene.Mutation`): mutation to be filtered
                2. cov (float): coverage of the mutation
                3. total (float): total coverage of the mutation locus
                4. thres (float): filtering threshold

            `filter_fn` returns `False` if a mutation should be filtered out.

        :returns: Filtered coverage.
        """

        new_cov = copy.copy(self)
        new_cov._coverage = {}
        for pos, pos_mut in self._coverage.items():
            new_cov._coverage[pos] = {}
            for o in pos_mut:
                f = filter_fn(self, Mutation(pos, o))
                if f:
                    new_cov._coverage[pos][o] = f
        return new_cov

    def diploid_avg_coverage(self) -> float:
        """:returns: Average coverage of the copy-number neutral region."""
        assert self.profile.cn_region, "CN region not set"
        return float(sum(self._cnv_coverage.values())) / abs(
            self.profile.cn_region.end - self.profile.cn_region.start
        )

    def _normalize_coverage(self) -> None:
        """Normalize the sample coverage with the profile coverage."""

        assert self.profile.cn_region and self.profile.data, "CN region not set"
        sam_ref = sum(
            self._cnv_coverage[i]
            for i in range(self.profile.cn_region.start, self.profile.cn_region.end)
        )
        if sam_ref == 0:
            raise AldyException(
                f"CN-neutral region {self.profile.cn_region} has no reads. "
                + "Double check your input file for CYP2D8 (are you using hg19?), "
                + "or pass an alternative CN-neutral region via -n parameter."
            )
        ratio = self.profile.neutral_value / sam_ref
        if ratio == 0:
            raise AldyException("Invalid CN-neutral region in the provided profile.")
        log.debug("[coverage] scale_ratio: {:.1f}", 1 / ratio)

        self._region_coverage = {}
        for gene, gr in enumerate(self.gene.regions):
            for region, rng in gr.items():
                s = sum(self.total(i) for i in range(rng.start, rng.end))
                p = self.profile.data[self.gene.name][region][gene]
                p /= 2  # profile has 2 copies, so divide it with 2 for normalization
                self._region_coverage[gene, region] = (ratio * s / p) if p != 0 else 0.0

    def basic_filter(self, mut: Mutation, cn=None, thres=None) -> List:
        """Basic threshold-based filter."""
        thres = (thres or self.profile.threshold) / (cn or 1)
        quals = self._coverage[mut.pos][mut.op]
        min_cov = max(self.profile.min_coverage, self.total(mut.pos) * thres)
        return quals if len(quals) >= min_cov else []

    def quality_filter(self, mut: Mutation) -> List:
        """Basic quality filter."""
        quals = self._coverage[mut.pos].get(mut.op, [])
        return [
            (m, q)
            for m, q in quals
            if q >= self.profile.min_quality
            if m >= self.profile.min_mapq
        ]
