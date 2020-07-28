# 786
# Aldy source: cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, List, Tuple, Optional

import copy
from natsort import natsorted

from . import lpinterface

from .common import log, json_print, sorted_tuple, AldyException
from .gene import CNConfig, Gene
from .coverage import Coverage
from .solutions import CNSolution
from .lpinterface import escape_name


MAX_CN = 20.0
"""float: Maximum allowed copy number."""

# Model parameters

LEFT_FUSION_PENALTY = 0.1
"""float: Extra penalty applied to left fusions to account for their rarity
   (0 for no penalty)."""

PCE_PENALTY_COEFF = 1.5
"""float: Error penalty applied to the PCE region (1 for no penalty)."""

MAX_CN_ERROR = MAX_CN
"""float: Upper bound for absolute error."""

PARSIMONY_PENALTY = 0.5
"""float: Extra penalty applied to each present gene copy to account for
   solution parsimony (0 for no penalty)."""


def estimate_cn(
    gene: Gene,
    coverage: Coverage,
    solver: str,
    gap: float = 0,
    fusion_penalty: float = LEFT_FUSION_PENALTY,
    user_solution=None,
    debug: Optional[str] = None,
) -> List[CNSolution]:
    """
    Estimate optimal copy number configurations given a gene and read data.

    Args:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        coverage (:obj:`aldy.coverage.Coverage`):
            Read data coverage instance.
        solver (str):
            ILP solver. Check :obj:`aldy.lpinterface` for the list of supported solvers.
        gap (float):
            Relative optimality gap. Use non-zero values to allow non-optimal solutions.
            Default is 0 (reports only optimal solutions).
        fusion_penalty (float):
            Fusion penalty. Use higher values to avoid fusions.
            Default is 0.1.
        user_solution (list[str], optional):
            User-specified list of copy number configurations.
            ILP solver will not run if this parameter is provided.
            Default is ``None``.
        debug (str, optional):
            If set, create a "`debug`.cn.lp" file for debug purposes.
            Default is ``None``.

    Returns:
        list[:obj:`aldu.solutions.CNSolution`]: Optimal copy number configurations.
    """

    log.debug("\n" + "*" * 80)

    if user_solution is not None:
        return [_parse_user_solution(gene, user_solution)]
    elif not gene.do_copy_number:
        basic = list(gene.cn_configs.keys())[0]
        return [_parse_user_solution(gene, [basic, basic])]
    else:
        # TODO: filter CN configs with non-present alleles
        max_observed_cn = 1 + max(
            int(round(coverage.region_coverage(g, r)))
            for g in gene.regions
            for r in gene.regions[g]
        )
        region_cov = _region_coverage(gene, coverage)
        total_cov = sum(r0 + r1 for r0, r1 in region_cov.values())
        min_cov = min(
            sum(sum(v.values()) for v in gene.cn_configs[c].cn) for c in gene.cn_configs
        )
        if total_cov < min_cov / 2.0:
            raise AldyException(
                f"Coverage for {gene.name} too low for copy number calling."
            )

        configs = _filter_configs(gene, coverage)
        log.debug("[cn] candidates= {}", ", ".join(natsorted(configs)))
        log.debug("[cn] max_cn= {}", max_observed_cn)
        _print_coverage(gene, coverage)
        sol = solve_cn_model(
            gene,
            configs,
            max_observed_cn,
            region_cov,
            solver,
            gap,
            fusion_penalty,
            debug,
        )

        return sol


def solve_cn_model(
    gene: Gene,
    cn_configs: Dict[str, CNConfig],
    max_cn: int,
    region_coverage: Dict[str, Tuple[float, float]],
    solver: str,
    gap: float = 0,
    fusion_penalty: float = LEFT_FUSION_PENALTY,
    debug: Optional[str] = None,
) -> List[CNSolution]:
    """
    Solve the copy number estimation problem (an instance of closest vector problem).

    Args:
        cn_configs (dict[str, :obj:`aldy.gene.CNConfig`]):
            Available copy number configurations (vectors).
        max_cn (int):
            Maximum allowed copy number.
        region_coverage (dict[str, tuple[float, float]]):
            Observed copy number of the main gene and the pseudogene
            for each genic region.
        solver (str):
            ILP solver to use.
            Check :obj:`aldy.lpinterface` for the list of supported solvers.
        gap (float):
            Relative optimality gap. Use non-zero values to allow non-optimal solutions.
            Default is 0 (reports only optimal solutions).
        fusion_penalty (float):
            Fusion penalty. Use higher values to avoid fusions.
            Default is 0.1.
        debug (str, optional):
            If set, create a "`debug`.cn.lp" file for debug purposes.
            Default is ``None``.

    Returns:
        list[:obj:`aldy.solutions.CNSolution`]: Optimal copy-number solutions.

    Notes:
        Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_
        (section Methods/Copy number and structural variation estimation)
        for the detailed description of ILP model.
    """

    model = lpinterface.model("AldyCN", solver)

    # List of CN configurations (a.k.a. structures):
    # dict of (`name`, `number`): structure, where multiple copies of the same `name`
    # get different numbers (binary LP transformation).
    # Note that `number` = {0, -1} represents full configuration with *all* pseudogenes.
    # Thus a diploid genome must have exactly *2* full configurations.
    # Any `number` > 0 describes only the main gene configuration
    # (a.k.a. "weak" configurations), and there can be many such configurations.
    structures = {(name, 0): structure for name, structure in cn_configs.items()}
    for a, ai in list(structures.keys()):
        structures[a, -1] = copy.deepcopy(structures[a, 0])
        if cn_configs[a].kind != CNConfig.CNConfigType.DEFAULT:
            continue
        for i in range(1, max_cn):
            structures[a, i] = copy.deepcopy(structures[a, 0])
            for g in range(1, len(structures[a, i].cn)):
                # If this is a pseudogene, remove a copy (create a "weak configuration")
                structures[a, i].cn[g] = {
                    r: v - 1 for r, v in structures[a, i].cn[g].items()
                }

    # Add a binary variable for each CN structure.
    VCN = {
        (a, ai): model.addVar(vtype="B", name=f"CN_{a}_{ai}") for a, ai in structures
    }

    # We assume diploid genome, so the number of haplotype-inducing configurations
    # must be 2.
    diplo_inducing = model.quicksum(VCN[a] for a in VCN if a[1] <= 0)
    model.addConstr(diplo_inducing == 2, name="CDIPLO")

    # Ensure that we cannot link any allele to the whole-gene deletion.
    del_allele = gene.deletion_allele()
    for (a, ai), v in VCN.items():
        if del_allele and a != del_allele:
            model.addConstr(v + VCN[del_allele, -1] <= 1, name=f"CDEL_{a}_{ai}")

    # Ensure that binary transformation is properly formed (i.e. A_i <= A_{i-1}).
    for a, ai in structures:
        # The second haplotype (-1) is only present if the first one (0) is there.
        if ai == -1:
            model.addConstr(VCN[a, ai] <= VCN[a, 0], name=f"CORD_{a}_{ai}")
        # Ignore 1, because A[1] = 1 && A[0] = 0 is valid (e.g. *13/*13+*1)
        elif ai > 1:
            model.addConstr(VCN[a, ai] <= VCN[a, ai - 1], name=f"CORD_{a}_{ai}")

    # Add error variables
    VERR = {}
    json_print(debug, '    "data": {', end="")
    for r, (exp_cov0, exp_cov1) in region_coverage.items():
        json_print(debug, f"'{r}': ({exp_cov0}, {exp_cov1}), ", end="")
        expr = 0
        for s, structure in structures.items():
            if r in structure.cn[0]:
                expr += structure.cn[0][r] * VCN[s]
            if len(structure.cn) > 1 and r in structure.cn[1]:
                expr -= structure.cn[1][r] * VCN[s]
        VERR[r] = model.addVar(name=f"E_{r}", lb=-MAX_CN_ERROR, ub=MAX_CN_ERROR)
        model.addConstr(expr + VERR[r] == exp_cov0 - exp_cov1, name=f"CCOV_{r}")
    json_print(debug, "},")
    # Objective: minimize the sum of absolute errors.
    # PCE_REGION (in CYP2D7) is penalized with an extra score as it is important
    # fusion marker.
    objective = model.abssum(VERR.values(), coeffs={"E_pce": PCE_PENALTY_COEFF})
    # Objective: also minimize the total number of present alleles (maximum parsimony)
    # Also penalize left fusions as they are not likely to occur.
    objective += model.quicksum(
        # MPICL requires only one term for each variable in the objective function,
        # thus this ugly expression:
        (fusion_penalty + PARSIMONY_PENALTY) * v
        if cn_configs[s].kind == CNConfig.CNConfigType.LEFT_FUSION
        else PARSIMONY_PENALTY * v
        for (s, k), v in VCN.items()
    )
    model.setObjective(objective)
    if debug:
        model.dump(f"{debug}.cn.lp")

    # Solve the model
    try:
        lookup = {model.varName(v): a for (a, ai), v in VCN.items()}
        result: dict = {}
        for status, opt, sol in model.solutions(gap):
            log.debug(f"[cn] status= {status}; opt= {opt:.2f}")
            sol_tuple = sorted_tuple(lookup[v] for v in sol)
            # Because A[1] can be 1 while A[0] is 0, we can have biologically
            # duplicate solutions
            if sol_tuple not in result:
                result[sol_tuple] = CNSolution(  # type: ignore
                    opt, solution=list(sol_tuple), gene=gene
                )
                log.debug("[cn] solution= {}", result[sol_tuple])
        json_print(
            debug, '    "sol": ' + str([dict(r.solution) for r in result.values()])
        )
        return list(result.values())
    except lpinterface.NoSolutionsError:
        log.debug("[cn] solution= []")
        return []


def _filter_configs(gene: Gene, coverage: Coverage) -> Dict[str, CNConfig]:
    """
    Filter out low-quality mutations and copy number configurations
    that are not supported by the remaining mutations.

    Returns:
        dict[str, :obj:`aldy.gene.CNConfig`]
    """
    cov = coverage.filtered(
        lambda mut, cov, total, thres: Coverage.basic_filter(
            mut, cov, total, thres / MAX_CN
        )
    )
    configs = copy.deepcopy(gene.cn_configs)
    for an in natsorted(gene.cn_configs):
        if an not in gene.alleles:
            continue  # This is just a CN configuration w/o any mutations
        bad_alleles = []
        for a in gene.cn_configs[an].alleles:
            if any(cov[m] <= 0 for m in gene.alleles[a].func_muts):
                bad_alleles.append(a)
        if len(bad_alleles) == len(gene.cn_configs[an].alleles):
            log.trace(f"Removing CN {a} due to low support")
            del configs[an]
    return configs


def _region_coverage(gene: Gene, coverage: Coverage) -> Dict[str, Tuple[float, float]]:
    """
    Calculate the coverage  of the main gene and the pseudogene in each genic region.
    Returns dictionary where

    Returns:
        dict[str, tuple[float, float]]: Region coverage
        that links genic regions (e.g. exon 1) to the coverage of the main gene
        and the pseudogene (expressed as a tuple).
    """

    return {
        r: (
            coverage.region_coverage(0, r),
            coverage.region_coverage(1, r) if 1 in gene.regions else 0,
        )
        for r in gene.unique_regions
    }


def _print_coverage(gene: Gene, coverage: Coverage) -> None:
    """
    Pretty-print the region coverage.
    """
    for ri in range(0, len(gene.regions[0]), 13):
        regions = list(gene.regions[0])[ri : ri + 13]
        log.debug("[cn]          {}", " ".join(f"{r:5}" for r in regions))
        g = " ".join(f"{coverage.region_coverage(0, r):5.2f}" for r in regions)
        log.debug(f"[cn] {gene.name:7} {g}")
        if gene.pseudogenes:
            g = " ".join(f"{coverage.region_coverage(1, r):5.2f}" for r in regions)
            log.debug(f"[cn] {gene.pseudogenes[0]:7} {g}")
            diffs = {
                r: coverage.region_coverage(0, r) - coverage.region_coverage(1, r)
                for r in regions
            }
            log.debug(
                "[cn]      =  {}",
                " ".join(
                    f"{d:5.2f}" if r in gene.unique_regions else "     "
                    for r, d in diffs.items()
                ),
            )


def _parse_user_solution(gene: Gene, sols: List[str]) -> CNSolution:
    """
    Parse user-provided copy number solutions.

    Args:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        sols (list[str]):
            List of valid copy number configurations.

    Returns:
        :obj:`aldy.solutions.CNSolution`: User-provided copy number solution.

    Raises:
        :obj:`aldy.common.AldyException` if a user-provided solution does not match
        the gene database.
    """
    for sol in sols:
        if sol not in gene.cn_configs:
            raise AldyException(
                f"Given copy number solution contains unknown copy number configuration"
                + f"{sol}. Please run 'aldy show --gene {gene.name}' for the list the"
                + f"valid configurations"
            )
    s = CNSolution(0, solution=sols, gene=gene)  # type: ignore
    log.debug("[cn] result= {} (provided)", s)
    return s
