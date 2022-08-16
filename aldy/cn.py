# 786
# Aldy source: cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, List, Tuple, Optional
from math import ceil
from natsort import natsorted
from functools import partial
import copy
import re

from . import lpinterface
from .common import log, json, sorted_tuple, AldyException
from .gene import CNConfig, CNConfigType, Gene
from .coverage import Coverage
from .solutions import CNSolution
from .profile import Profile


def estimate_cn(
    gene: Gene,
    profile: Profile,
    coverage: Optional[Coverage],
    solver: str,
    debug: Optional[str] = None,
) -> List[CNSolution]:
    """
    Estimate the optimal copy number configuration for a sample given a gene and
    coverage information.

    :param gene: Gene instance.
    :param profile: Profile instance.
    :param coverage: Read coverage instance.
    :param solver: ILP solver (see :py:mod:`aldy.lpinterface` for supported solvers).
    :param debug: When set, create a `{debug}.cn.lp` model description file for
        debug purposes.
        Default: `None` (no debug dumps).

    :returns: List of optimal copy number configurations.
    """

    log.debug("\n" + "*" * 80)

    if profile.cn_solution:
        return [_parse_user_solution(gene, profile.cn_solution)]
    elif not gene.do_copy_number:
        basic = list(gene.cn_configs.keys())[0]
        cn = 2
        if profile.male and gene.chr in ["X", "Y"]:
            cn = 1
        return [_parse_user_solution(gene, [basic] * cn)]
    else:
        assert coverage, "Coverage not provided"

        # TODO: filter CN configs with non-present alleles
        max_observed_cn = 1 + max(
            ceil(coverage.region_coverage(gi, r))
            for gi, g in enumerate(gene.regions)
            for r in g
        )
        region_cov = {
            r: (
                coverage.region_coverage(0, r),
                coverage.region_coverage(1, r) if len(gene.regions) > 1 else 0.0,
            )
            for r in gene.unique_regions
        }
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

        fusion_support = None
        if coverage.sam._fusion_counter:
            fusion_support = {
                fn: (a / b) if b else 0.0
                for fn, [a, b] in coverage.sam._fusion_counter.items()
            }
        sol = solve_cn_model(
            gene,
            profile,
            configs,
            max_observed_cn,
            region_cov,
            solver,
            debug,
            fusion_support,
        )

        return sol


def solve_cn_model(
    gene: Gene,
    profile: Profile,
    cn_configs: Dict[str, CNConfig],
    max_cn: int,
    region_coverage: Dict[str, Tuple[float, float]],
    solver: str,
    debug: Optional[str] = None,
    fusion_support: Optional[Dict[str, float]] = None,
) -> List[CNSolution]:
    """
    Solve the copy number estimation problem (instance of the closest vector problem).

    :param gene: Gene instance.
    :param profile: Profile instance.
    :param cn_configs: Available copy number configurations (vectors).
    :param max_cn: Maximum allowed copy number.
    :param region_coverage: Observed copy number of the main gene and the pseudogene
        for each genic region.
    :param solver: ILP solver (see :py:mod:`aldy.lpinterface` for supported solvers).
    :param gap: Optimality gap. Non-zero values enable non-optimal solutions.
        Default: 0 (only optimal solutions).
    :param debug: When set, create a `{debug}.cn.lp` model description file for
        debug purposes.
        Default: `None` (no debug dumps).
    :param fusion_support: Dictionary that contains read support of each available
        fusion. Used only for the long-read fusion calling.
        Default is `None` (all fusions are treated equally).
    :returns: List of optimal copy-number solutions.

    .. note::
        Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_
        (section Methods/Copy number and structural variation estimation)
        for the detailed description of ILP model.
    """

    model = lpinterface.model("AldyCN", solver)
    debug_info = json[gene.name]["cn"]

    del_allele = gene.deletion_allele()

    # List of CN configurations (a.k.a. structures). Each configuration is a binary
    # variable. Each structure is defined by `('structure_name', number)`.
    # When `number` ∈ {0, -1}, the configuration is complete with *all* pseudogenes
    # included (a diploid genome must contain exactly *2* complete configurations).
    # When `number` > 0, the configuration describes only the main gene and not the
    # pseudogene (there can be many such configurations).
    # N.B. (3/2022) this step filters "weak" fusions without long-read support
    structures: Dict[Tuple[str, int], CNConfig] = {
        (name, 0): structure
        for name, structure in cn_configs.items()
        if not fusion_support
        or name == "1"
        or (del_allele and name == del_allele)
        or (name in fusion_support and fusion_support[name] >= 1 / (2 * max_cn))
    }
    for a, ai in list(structures.keys()):
        structures[a, -1] = copy.deepcopy(structures[a, 0])
        if cn_configs[a].kind != CNConfigType.DEFAULT:
            continue
        for i in range(1, max_cn):
            structures[a, i] = copy.deepcopy(structures[a, 0])
            for g in range(1, len(structures[a, i].cn)):
                # If this is a pseudogene, remove a copy (create a "weak configuration")
                structures[a, i].cn[g] = {
                    r: v - 1 for r, v in structures[a, i].cn[g].items()
                }
    # Add "fake" pseudogenes to handle pseudogene CN changes or CN noise
    if len(gene.regions) > 1 and del_allele:
        for i in range(max_cn):
            structures["PSEUDO", i + 1] = CNConfig(
                copy.deepcopy(structures[del_allele, 0].cn),
                CNConfigType.DELETION,
                set(),
                "pseudogene",
            )

    # Add a binary variable for each CN structure
    VCN = {
        (a, ai): model.addVar(vtype="B", name=f"CN_{a}_{ai}") for a, ai in structures
    }

    # We assume diploid genome, so the number of the complete configurations must be 2
    diplo_inducing = model.quicksum(VCN[a] for a in VCN if a[1] <= 0)
    model.addConstr(diplo_inducing <= 2, name="CDIPLO")
    model.addConstr(diplo_inducing >= 2, name="CDIPLO")

    # Ensure that we cannot link any allele to the whole-gene deletion
    if del_allele:
        for (a, ai), v in VCN.items():
            if a != del_allele:
                model.addConstr(v + VCN[del_allele, -1] <= 1, name=f"CDEL_{a}_{ai}")

    # Ensure that binary transformation is properly formed (i.e. A_i <= A_{i-1})
    for a, ai in structures:
        # The second haplotype (-1) is only present if the first one (0) is there
        if ai == -1:
            model.addConstr(VCN[a, ai] <= VCN[a, 0], name=f"CORD_{a}_{ai}")
        # Ignore 1, because A[1] = 1 && A[0] = 0 is valid (e.g. *13/*13+*1)
        elif ai > 1:
            model.addConstr(VCN[a, ai] <= VCN[a, ai - 1], name=f"CORD_{a}_{ai}")

    # Add error variables
    VERR, VERR_GENE = {}, {}
    for r, (exp_cov0, exp_cov1) in region_coverage.items():
        debug_info["data"][r] = (exp_cov0, exp_cov1)
        expr, expr_gene = 0, 0
        for s, structure in structures.items():
            if r in structure.cn[0]:
                expr += structure.cn[0][r] * VCN[s]
                expr_gene += structure.cn[0][r] * VCN[s]
            if len(structure.cn) > 1 and r in structure.cn[1]:
                expr -= structure.cn[1][r] * VCN[s]

        if r not in gene.unique_regions:
            continue

        VERR_GENE[r] = model.addVar(
            name=f"EG_{r}", lb=-profile.cn_max, ub=profile.cn_max
        )
        model.addConstr(expr_gene + VERR_GENE[r] <= exp_cov0, name=f"CG_COV_{r}")
        model.addConstr(expr_gene + VERR_GENE[r] >= exp_cov0, name=f"CG_COV_{r}")

        scale = max(exp_cov0, exp_cov1) + 1
        VERR[r] = model.addVar(name=f"E_{r}", lb=-profile.cn_max, ub=profile.cn_max)
        model.addConstr(
            expr / scale + VERR[r] <= (exp_cov0 - exp_cov1) / scale, name=f"C_COV_{r}"
        )
        model.addConstr(
            expr / scale + VERR[r] >= (exp_cov0 - exp_cov1) / scale, name=f"C_COV_{r}"
        )

    # Objective component 1: minimize the sum of absolute errors.
    # PCE_REGION (in CYP2D7) is penalized with an extra score (important fusion marker).
    DIFF_COEFF = profile.cn_diff / len(gene.unique_regions)
    o_diff = DIFF_COEFF * model.abssum(
        VERR.values(), coeffs={"E_pce": profile.cn_pce_penalty}
    )

    # Objective component 2: minimize the main gene fitness.
    FIT_COEFF = profile.cn_fit / len(gene.unique_regions)
    o_fit = FIT_COEFF * model.abssum(VERR_GENE.values())

    # Objective component 3: minimize the total number of present alleles
    #                        (maximum parsimony)
    PARSIMONY_PENALTY = 10.0 / len(gene.unique_regions)
    PARSIMONY_PENALTY *= 0.75
    penalty = {s: PARSIMONY_PENALTY for s, _ in VCN}
    # Also penalize the left fusions as they are less likely
    for n, sv in gene.cn_configs.items():
        if n in penalty and sv.kind == CNConfigType.RIGHT_FUSION:
            penalty[n] += PARSIMONY_PENALTY * profile.cn_fusion_right
        if n in penalty and sv.kind == CNConfigType.LEFT_FUSION:
            penalty[n] += PARSIMONY_PENALTY * profile.cn_fusion_left
    o_pars = profile.cn_parsimony * model.quicksum(
        penalty[s] * v for (s, _), v in VCN.items()
    )

    model.setObjective(o_diff + o_fit + o_pars)
    if debug:
        model.dump(f"{debug}.{gene.name}.cn.lp")

    # Solve the model
    lookup = {model.varName(v): a for (a, _), v in VCN.items()}
    result: dict = {}
    for status, opt, sol in model.solutions(profile.gap):
        sol_tuple = sorted_tuple(
            lookup[v] for v in sol if lookup[v] not in [del_allele, "PSEUDO"]
        )
        if sol_tuple not in result:
            result[sol_tuple] = CNSolution(gene, opt, list(sol_tuple))
            log.debug(
                f"[cn] status= {status}; opt= {opt:.2f} "
                + "(diff= {:.2f}, fit= {:.2f}, pars= {:.2f}) "
                + f"solution= {result[sol_tuple]}",
                model.getValue(o_diff),
                model.getValue(o_fit),
                model.getValue(o_pars),
            )
    if not result:
        log.debug("[cn] solution= []")

    debug_info["sol"] = [dict(r.solution) for r in result.values()]
    return list(result.values())


def _filter_configs(gene: Gene, coverage: Coverage) -> Dict[str, CNConfig]:
    """
    Filter out low-quality mutations and copy number configurations that are not
    supported by the remaining mutations.
    """

    cov = coverage.filtered(
        partial(
            Coverage.basic_filter,
            thres=coverage.profile.threshold / coverage.profile.cn_max,
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
            log.trace(f"[cn] removing *{an} due to low support")
            del configs[an]
    return configs


def _print_coverage(gene: Gene, coverage: Coverage) -> None:
    """Pretty-print the region coverage."""

    log.debug("[cn] coverage=")
    gname_ = re.split(r"(\d.+)", gene.name)
    gname = gname_[1] if len(gname_) > 1 else gene.name
    pname = ""
    if gene.pseudogenes:
        pname_ = re.split(r"(\d.+)", gene.pseudogenes[0])
        pname = pname_[1] if len(pname_) > 1 else gene.pseudogenes[0]
    log.debug(f"  {'':5} {gname:4} {pname:4}")
    for r in gene.regions[0]:
        g = coverage.region_coverage(0, r)
        if pname:
            p = coverage.region_coverage(1, r)
            if r in gene.unique_regions:
                log.debug(f"  {r:5} {g:4.1f} {p:4.1f} = {g - p:4.1f}")
            else:
                log.debug(f"  {r:5} {g:4.1f} {p:4.1f}")
        else:
            log.debug(f"  {r:5}: {g:4.1f}")


def _parse_user_solution(gene: Gene, sols: List[str]) -> CNSolution:
    """
    Parse user-provided copy number solutions.

    :param gene: Gene instance.
    :param sols: List of valid copy number configurations.
    :returns: User-provided copy number solution.
    :raise: :py:class:`aldy.common.AldyException` if a user-provided solution does
        not match the gene database.
    """

    for sol in sols:
        if sol not in gene.cn_configs:
            raise AldyException(
                "The copy number solution contains unknown copy number configuration "
                f"{sol}. Please run 'aldy query {gene.name}' for the list of "
                "the valid configurations"
            )
    s = CNSolution(gene, 0, sols)  # type: ignore
    log.debug("[cn] result= {} (provided)", s)
    return s
