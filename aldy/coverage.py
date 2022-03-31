# 786
# Aldy source: coverage.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, Tuple, Callable, List

import collections

from .common import log, AldyException, script_path
from .gene import Mutation, Gene


class Coverage:
    """
    Data structure that maintains the coverage information for a given sample.
    """

    def __init__(
        self,
        gene: Gene,
        profile,
        sam,
        coverage: Dict[int, Dict[str, int]],
        cnv_coverage: Dict[int, int],
        threshold: float,
        min_cov: float = 1.0,
    ) -> None:
        """
        Coverage initialization.

        :param coverage:
            Coverage for each locus within a sample represented as a dictionary that
            maps mutation (or a reference position indicated by `_`) to the number of
            reads supporting that mutation.
            For example, ``coverage[10]['A>G'] = 2`` means that there are 2 reads that
            have G (instead of A) at the genomic locus 10.
        :param threshold:
            Threshold `t` used for filtering out low-quality mutations.
            Any mutation with the coverage less than `t`% is filtered out.
            Ranges from 0 to 1 (normalized percentage).
        :param cnv_coverage:
            Coverage of the copy-number neutral region of the sample: each genomic locus
            within that region points to the corresponsing read coverage.
            Used for coverage rescaling.
        :param sample: Sample name. Default: 'sample'.
        """
        self.gene = gene
        self.profile = profile
        self.sam = sam

        self._coverage = coverage
        self._cnv_coverage = cnv_coverage
        self._threshold = threshold
        self.min_cov = min_cov

        self._region_coverage: Dict[Tuple[int, str], float] = {}

    def __getitem__(self, mut: Mutation) -> float:
        """:return: Coverage of the mutation ``mut``."""
        return self.coverage(mut)

    def coverage(self, mut: Mutation) -> float:
        """:return: Coverage of the mutation ``mut``."""
        if mut.op in self._coverage[mut.pos]:
            return self._coverage[mut.pos][mut.op]
        else:
            return 0

    def total(self, pos: int) -> float:
        """:return: Total coverage at the locus ``pos``."""
        if pos not in self._coverage:
            return 0
        return float(sum(v for p, v in self._coverage[pos].items() if p[:3] != "ins"))

    def percentage(self, m: Mutation) -> float:
        """:return: Coverage of the mutation ``mut`` as a percentage (0-100%)."""
        total = self.total(m.pos)
        if total == 0:
            return 0
        return 100.0 * self.coverage(m) / total

    def single_copy(self, pos: int, cn_solution) -> float:
        """
        :param pos: Genomic locus,
        :param cn_solution: Copy-number solution.
        :return: Coverage of a single gene copy at the locus ``pos``.
        """
        if cn_solution.position_cn(pos) == 0:
            return 0
        return max(1, self.total(pos)) / cn_solution.position_cn(pos)

    def region_coverage(self, gene: int, region: str) -> float:
        """:return: Average coverage of the region ``region`` in ``gene``."""
        return self._region_coverage[gene, region]

    def average_coverage(self) -> float:
        """:return: Average coverage of the sample."""
        return sum(self.total(pos) for pos in self._coverage) / float(
            len(self._coverage) + 0.1
        )

    def dump(self, out):
        for pos, pos_mut in sorted(self._coverage.items()):
            if not (len(pos_mut) == 1 and "_" in pos_mut):
                for i, (op, cov) in enumerate(sorted(pos_mut.items(), reverse=True)):
                    p = self.percentage(Mutation(pos, op))
                    if pos in self.gene:
                        x = self.gene.get_functional((pos, op))
                        if op == "_":
                            x = ""
                        if x and (pos, op) not in self.gene.mutations:
                            x += "**"
                        t = f"{x if x else ''}\t{self.gene.region_at(pos)[1]}"
                        t += "\t" + self.gene.get_rsid((pos, op))
                    else:
                        t = "\t\t"
                    # out(
                    #     f"[dump] {self.sample}\t{gene.name}\t{gene.chr_to_ref.get(pos, '-')+1}\t{op}\t{p:.1f}\t{t}"
                    # )

    def filtered(
        self, filter_fn: Callable[[Mutation, float, float, float], bool]
    ):  # -> Coverage
        """
        :param filter_fn:
            Function that performs mutation filtering with the following arguments:

                1. mut (:obj:`aldy.gene.Mutation`): mutation to be filtered
                2. cov (float): coverage of the mutation
                3. total (float): total coverage of the mutation locus
                4. thres (float): filtering threshold

            ``filter_fn`` returns ``False`` if a mutation is filtered out.

        :return: Filtered coverage.
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

        new_cov = Coverage(
            self.gene,
            self.profile,
            self.sam,
            cov,
            self._cnv_coverage,
            self._threshold,
            self.min_cov,
        )
        new_cov._region_coverage = self._region_coverage
        return new_cov

    def diploid_avg_coverage(self) -> float:
        """:return: Average coverage of the copy-number neutral region."""
        return float(sum(self._cnv_coverage.values())) / abs(
            self.profile.cn_region.end - self.profile.cn_region.start
        )

    def _normalize_coverage(self) -> None:
        """
        Normalize the sample coverage to match the profile coverage.
        """

        sam_ref = sum(
            self._cnv_coverage[i]
            for i in range(self.profile.cn_region.start, self.profile.cn_region.end)
        )
        # print(self._cnv_coverage, cn_region, sam_ref)

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
                s = sum(self.total(i) for i in range(rng.start, rng.end))  # !IMPORTANT
                p = self.profile.data[self.gene.name][region][gene]
                self._region_coverage[gene, region] = ((ratio * s / p) if p != 0 else 0.0)

    @staticmethod
    def basic_filter(
        mut: Mutation, cov: float, total: float, thres: float, min_cov: float
    ) -> bool:
        """
        Basic filtering function.
        """
        return cov >= max(min_cov, total * thres)

    @staticmethod
    def cn_filter(
        mut: Mutation,
        cov: float,
        total: float,
        thres: float,
        cn_solution,
        min_cov: float,
    ) -> bool:
        """
        Filtering function that takes into the account the copy number of the mutation.
        """
        cn = cn_solution.position_cn(mut.pos)
        total = total / cn if cn > 0 else total

        return mut.op == "_" or cov >= max(min_cov, total * thres)

    # @deprecated
    def _load_phase(self, gene, ploidy=2):
        """Loads a HapTree-X/HapCUT2 phase file."""

        if len(self.fragments) == 0:
            return

        import subprocess, os, tempfile

        phases: List[List[Mutation]] = []
        with tempfile.TemporaryDirectory() as tmp:
            snps = {}
            for read in self.fragments:
                for pos, op in read:
                    snps.setdefault(pos, set()).add(op)
            snps = {pos: ops for pos, ops in snps.items() if len(ops) > 1}
            frags = {}
            for read in self.fragments:
                key = tuple(sorted([(pos, op) for pos, op in read if pos in snps]))
                if key not in frags:
                    frags[key] = 0
                frags[key] += 1

            frag_file = f"{tmp}/fragments.txt"
            with open(frag_file, "w") as fo:
                for s, c in frags.items():
                    if not s:
                        continue
                    fo.write(f"{c}")
                    for p, a in s:
                        fo.write(f";{p},{a}")
                    fo.write("\n")
            phase_file = f"{tmp}/phase.txt"

            htx = script_path("aldy.resources/haptreex")
            # subprocess.run(f"cp {frag_file} _y", shell=True)
            cmd = [htx, "--aldy", frag_file, "-o", phase_file, "-p", str(ploidy)]
            ret = subprocess.run(cmd, env={**os.environ, "OMP_NUM_THREADS": "1"})
            if ret.returncode != 0:
                raise AldyException(
                    "HapTree-X failed"
                    + (": " + ret.stderr.decode("ascii") if ret.stderr else ""),
                )
            log.debug("HapTree-X done")

            haplotypes = [[] for _ in range(ploidy)]
            _, g_s, g_e = gene.get_wide_region()
            with open(phase_file) as hap:
                for li, line in enumerate(hap):
                    if line[:5] == "BLOCK":
                        if haplotypes[0]:
                            phases += list(haplotypes)
                        haplotypes = [], []
                    elif line[:5] == "*****":
                        continue
                    else:
                        ls = line.strip().split("\t")
                        if len(ls) < 7:
                            raise AldyException(f"Invalid phasing line {li + 1}")
                        _, gt0_, gt1_, _, pos_, al0, al1, *_ = ls
                        gt0, gt1, pos = int(gt0_), int(gt1_), int(pos_) - 1
                        if gt0 + gt1 != 1:
                            continue
                        if pos < g_s or pos >= g_e:
                            continue

                        if (pos, f"{al0}>{al1}") not in gene.mutations:
                            if (pos, f"{al1}>{al0}") in gene.mutations:
                                log.warn(f"reorienting {pos} {al0} {al1}")
                                al1, al0 = al0, al1
                                gt0, gt1 = gt1, gt0
                        haplotypes[0].append(
                            Mutation(pos, f"{al0}>{al1}" if gt0 else "_")
                        )
                        haplotypes[1].append(
                            Mutation(pos, f"{al0}>{al1}" if gt1 else "_")
                        )
            if haplotypes[0]:
                phases += [haplotypes[0], haplotypes[1]]
        self.phases = phases
