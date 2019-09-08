# 786
# Aldy source: coverage.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, Tuple, Callable

import collections

from .common import log
from .gene import Mutation, GeneRegion, GRange


class Coverage:
    """
    Data structure that maintains the coverage information for a given sample.
    """

    def __init__(
        self,
        coverage: Dict[int, Dict[str, int]],
        threshold: float,
        cnv_coverage: Dict[int, int],
    ) -> None:
        """
        Coverage initialization.

        Args:
            coverage (dict[int, dict[str, int]]):
                Coverage for each locus within a sample represented as a dictionary
                that maps mutation (or a reference position indicated by `_`)
                to the number of reads supporting that mutation.
                For example, ``coverage[10]['SNP.AG'] = 2`` means that there are 2
                reads that have G (instead of A) at the genomic locus 10.
            threshold (float):
                Threshold `t` used for filtering out low-quality mutations.
                Any mutation with the coverage less than `t`% is filtered out.
                Ranges from 0 to 1 (normalized percentage).
            cnv_coverage (dict[int, int]):
                Coverage of the copy-number neutral region of the sample:
                each genomic locus within that region points to the corresponsing read
                coverage.
                Used for coverage rescaling.
        """
        self._coverage = coverage
        self._threshold = threshold
        self._cnv_coverage = cnv_coverage

        self._rescaled: Dict[int, float] = {}
        self._region_coverage: Dict[Tuple[int, GeneRegion], float] = {}

    def __getitem__(self, mut: Mutation) -> float:
        """
        Returns:
            float: Coverage of the mutation ``mut``.
        """
        return self.coverage(mut)

    def coverage(self, mut: Mutation) -> float:
        """
        Returns:
            float: Coverage of the mutation ``mut``.
        """
        if mut.op in self._coverage[mut.pos]:
            return self._coverage[mut.pos][mut.op]
        else:
            return 0

    def total(self, pos: int) -> float:
        """
        Returns:
            float: Total coverage at the locus ``pos``.
        """
        if pos not in self._coverage:
            return 0
        return float(sum(v for p, v in self._coverage[pos].items() if p[:3] != "INS"))

    def percentage(self, m: Mutation) -> float:
        """
        Returns:
            float: Coverage of the mutation ``mut`` expressed as a percentage
            in the range 0-100.
        """
        total = self.total(m.pos)
        if total == 0:
            return 0
        return 100.0 * self.coverage(m) / total

    def loci_cn(self, pos: int) -> float:
        """
        Returns:
            float: Copy number of the locus ``pos``.
        """
        if self._rescaled[pos] == 0:
            return 0
        return self.total(pos) * (1 / self._rescaled[pos])

    def single_copy(self, pos: int, cn_solution) -> float:
        """
        Args:
            pos (int): genomic locus
            cn_solution (:obj:`aldy.solutions.CNSolution`): copy-number solution
        Returns:
            float: Coverage of a single gene copy at the locus ``pos``.
        """
        if cn_solution.position_cn(pos) == 0:
            return 0
        return max(1, self.total(pos)) / cn_solution.position_cn(pos)

    def region_coverage(self, gene: int, region: GeneRegion) -> float:
        """
        Returns:
            float: Average coverage of the region ``region`` in ``gene``.
        """
        return self._region_coverage[gene, region]

    def average_coverage(self) -> float:
        """
        Returns:
            float: Average coverage of the sample.
        """
        return sum(self.total(pos) for pos in self._coverage) / float(
            len(self._coverage) + 0.1
        )

    def filtered(
        self, filter_fn: Callable[[Mutation, float, float, float], bool]
    ):  # -> Coverage
        """
        Args:
            filter_fn (callable):
                Function that performs mutation filtering with the following arguments:

                    1. mut (:obj:`aldy.gene.Mutation`): mutation to be filtered
                    2. cov (float): coverage of the mutation
                    3. total (float): total coverage of the mutation locus
                    4. thres (float): filtering threshold

                ``filter_fn`` returns ``False`` if a mutation is filtered out.

        Returns:
            :obj:`Coverage`: Filtered coverage.
        """

        cov = collections.defaultdict(
            lambda: collections.defaultdict(int),
            {
                pos: collections.defaultdict(
                    int,
                    {
                        o: c
                        for o, c in pos_mut.items()
                        if filter_fn(
                            Mutation(pos, o),  # type: ignore
                            c,
                            self.total(pos),
                            self._threshold,
                        )
                    },
                )
                for pos, pos_mut in self._coverage.items()
            },
        )

        new_cov = Coverage(cov, self._threshold, self._cnv_coverage)  # type: ignore
        new_cov._rescaled = self._rescaled
        new_cov._region_coverage = self._region_coverage
        return new_cov

    def diploid_avg_coverage(self) -> float:
        """
        Returns:
            float: Average coverage of the copy-number neutral region.
        """
        return float(sum(self._cnv_coverage.values())) / abs(
            self._cn_region.end - self._cn_region.start
        )

    def _normalize_coverage(
        self,
        profile: Dict[str, Dict[int, float]],
        gene_regions: Dict[int, Dict[GeneRegion, GRange]],
        cn_region: GRange,
    ) -> None:
        """
        Normalize the sample coverage to match the profile coverage.

        Args:
            profile (dict[str, dict[int, float]]):
                Profile coverage in the form `chromosome: (position -> coverage)`.
            gene_regions
            (dict[int, dict[:obj:`aldy.common.GeneRegion`, :obj:`aldy.common.GRange`]]):
                List of genic regions for each gene.
            cn_region (:obj:`aldy.common.GRange`):
                Copy-number neutral region.
        """

        #: GRange: store the CN-neutral region
        self._cn_region: GRange = cn_region
        sam_ref = sum(
            self._cnv_coverage[i] for i in range(cn_region.start, cn_region.end)
        )
        cnv_ref = sum(
            profile[cn_region.chr][i] for i in range(cn_region.start, cn_region.end)
        )

        cn_ratio = float(cnv_ref) / sam_ref
        log.debug("CNV factor: {} ({})", cn_ratio, 1.0 / cn_ratio)

        self._rescaled: Dict[int, float] = {}
        self._region_coverage: Dict[Tuple[int, GeneRegion], float] = {}
        for gene, gr in gene_regions.items():
            for region, rng in gr.items():
                s = sum(self.total(i) for i in range(rng.start, rng.end))  # !IMPORTANT
                p = sum(profile[rng.chr][i] for i in range(rng.start, rng.end))
                self._rescaled.update(
                    {
                        i: profile[rng.chr][i] / cn_ratio
                        for i in range(rng.start, rng.end)
                    }
                )
                self._region_coverage[gene, region] = (
                    (cn_ratio * float(s) / p) if p != 0 else 0.0
                )

    @staticmethod
    def basic_filter(mut: Mutation, cov: float, total: float, thres: float) -> bool:
        """
        Basic filtering function.
        """
        return cov >= max(1, total * thres)

    @staticmethod
    def cn_filter(
        mut: Mutation, cov: float, total: float, thres: float, cn_solution
    ) -> bool:
        """
        Filtering function that takes into the account the copy number of the mutation.
        """
        cn = cn_solution.position_cn(mut.pos)
        total = total / cn if cn > 0 else total
        return mut.op == "_" or cov > max(1, total * thres)
