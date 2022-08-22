# 786
# Aldy source: major.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict, Tuple, Any, Optional
from natsort import natsorted
import collections
import copy

from . import lpinterface
from .common import log, json, sorted_tuple
from .gene import MajorAllele, Mutation, Gene
from .coverage import Coverage
from .solutions import CNSolution, MajorSolution, SolvedAllele


def estimate_major(
    gene: Gene,
    coverage: Coverage,
    cn_solution: CNSolution,
    solver: str,
    identifier: int = 0,
    debug: Optional[str] = None,
) -> List[MajorSolution]:
    """
    Estimate optimal major star-alleles.

    :param gene: Gene instance.
    :param profile: Profile instance.
    :param coverage: Read coverage instance.
    :param cn_solution: Copy-number solution for major star-allele calling.
    :param solver: ILP solver (see :py:mod:`aldy.lpinterface` for supported solvers).
    :param identifier: Unique solution identifier. Used for debug purposes.
        Default: 0.
    :param debug: When set, create a `{debug}.major{identifier}.lp` model description
        file for debug purposes.
        Default: `None` (no debug dumps).
    """

    log.debug("*" * 80)
    log.debug("[major] struct= {}", cn_solution._solution_nice())

    alleles, coverage = _filter_alleles(gene, coverage, cn_solution)
    coverage.dump(log.trace)

    # Check if some CN solution has no matching allele
    if set(cn_solution.solution) - set(a.cn_config for a in alleles.values()):
        results: List[MajorSolution] = []
    else:
        results = solve_major_model(
            gene, coverage, cn_solution, alleles, solver, identifier, debug
        )
    # TODO: Check for novel functional mutations and do something with them

    return results


def solve_major_model(
    gene: Gene,
    coverage: Coverage,
    cn_solution: CNSolution,
    allele_dict: Dict[str, MajorAllele],
    solver: str,
    identifier: int = 0,
    debug: Optional[str] = None,
) -> List[MajorSolution]:
    """
    Solves the major star-allele detection problem via integer linear programming.

    :param gene: Gene instance.
    :param coverage: Read coverage instance.
    :param cn_solution: Copy-number solution for major star-allele calling.
    :param allele_dict: Candidate major star-alleles.
    :param solver: ILP solver (see :py:mod:`aldy.lpinterface` for supported solvers).
    :param identifier: Unique solution identifier. Used for debug purposes.
        Default: 0.
    :param debug: When set, create a `{debug}.major{identifier}.lp` model description
        file for debug purposes.
        Default: `None` (no debug dumps).

    .. note::
        Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_
        (section Methods/Major star-allele identification) for the model explanation.
    """

    model = lpinterface.model("AldyMajor", solver)
    debug_info = json[gene.name]["major"][len(json[gene.name]["major"])]

    # Get the list of _all_ functional mutations present in the sample
    # and the database (intersection)
    func_muts = {
        Mutation(*m)
        for m in gene.mutations
        if gene.is_functional(m) and coverage[Mutation(*m)] > 0
    }
    _print_candidates(gene, allele_dict, coverage, cn_solution, func_muts)

    a: Any = 0

    # Create a binary variable for every possible allele copy
    alleles = {(a, int(0)): allele_dict[a] for a in allele_dict}
    for (an, _), a in list(alleles.items()):
        max_cn = cn_solution.solution[a.cn_config]
        for i in range(1, max_cn):
            alleles[an, i] = alleles[an, 0]
    VA = {a: model.addVar(vtype="B", name=f"A_{a[0]}_{a[1]}") for a in alleles}

    # Make sure that VA[i+1] <= VA[i] (to avoid equivalent solutions)
    for a, ai in alleles.keys():
        if ai > 0:
            model.addConstr(VA[a, ai] <= VA[a, ai - 1], name=f"CORD_{a}_{ai}")

    # Add an error variable for every mutation
    VERR = {
        m: model.addVar(lb=-model.INF, ub=model.INF, name=f"E_{m.pos}_{m.op}")
        for m in func_muts
    }
    constraints = {e: 0 for e in VERR}

    # Add a binary variable for each mutation copy
    # For each present functional mutation, add a binary variable VNEW s.t.
    #   VNEW[m] = 0 <=> sum(VA[a] if m in a) >= 1
    VNEW = {m: model.addVar(vtype="B", name=f"N_{m}") for m in func_muts}

    # Populate the constraints
    for a in alleles:
        for m in alleles[a].func_muts:
            constraints[m] += VA[a]
    for m, v in VNEW.items():
        constraints[m] += v

    # Populate constraints of non-variations (i.e. matches with the reference genome)
    for pos in set(m.pos for m in constraints):
        ref_m = Mutation(pos, "_")  # type: ignore
        constraints[ref_m] = 0
        VERR[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=f"E_{pos}_REF")
        for a in alleles:
            if not gene.has_coverage(a[0], pos):
                continue
            # An insertion still contributes to the total coverage at this loci
            if any(ma[0] == pos and ma[1][:3] != "ins" for ma in alleles[a].func_muts):
                continue
            constraints[ref_m] += VA[a]

        # We ensure that only one additional mutation can be selected here (HACK)
        z = model.quicksum(
            v for m, v in VNEW.items() if m[0] == pos and m[1][:3] != "ins"
        )
        model.addConstr(z <= 1, name=f"CONE_{pos}")

    # Each allele must express all of its functional mutations
    debug_info["id"] = identifier
    debug_info["cn"] = str(dict(cn_solution.solution))
    debug_info["data"] = []
    for m, expr in sorted(constraints.items()):
        if coverage.single_copy(m.pos, cn_solution) == 0:
            cov = 0.0
        else:
            cov = coverage[m] / coverage.single_copy(m, cn_solution)
        model.addConstr(expr + VERR[m] <= cov, name=f"CFUNC_{m.pos}_{m.op}")
        model.addConstr(expr + VERR[m] >= cov, name=f"CFUNC_{m.pos}_{m.op}")
        debug_info["data"].append((m[0], m[1], cov))

    # Each CN configuration must be satisfied by corresponding alleles
    for cnf, cnt in cn_solution.solution.items():
        expr = sum(VA[a] for a in VA if alleles[a].cn_config == cnf)
        model.addConstr(expr <= cnt, name=f"CSAT_{cnf}")
        model.addConstr(expr >= cnt, name=f"CSAT_{cnf}")

    # Each functional mutation must be either chosen by some allele or marked as novel:
    #   1 == (m chosen by allele) XOR (m novel) == OR(VA[a] if m in a) XOR VNEW[m]
    for m in func_muts:
        VOR = model.addVar(vtype="B", name=f"OR_{m}")
        m_all = {a for a in alleles if m in alleles[a].func_muts}
        model.addConstr(VOR <= model.quicksum(VA[a] for a in m_all), name="COR")
        for a in m_all:
            model.addConstr(VOR >= VA[a], name="COR")
        VXOR = model.addVar(vtype="B", name=f"XOR_{m}")
        model.addConstr(VXOR <= VNEW[m] + VOR, name="CXOR")
        model.addConstr(VXOR <= 2 - VNEW[m] - VOR, name="CXOR")
        model.addConstr(VXOR >= VNEW[m] - VOR, name="CXOR")
        model.addConstr(VXOR >= VOR - VNEW[m], name="CXOR")
        model.addConstr(VXOR >= 1, name="CXOR")

    # Objective: minimize the absolute sum of errors and the number of novel mutations
    objective = 0
    objective += model.abssum(e for e in VERR.values())

    z = model.addVar(vtype="B", name="NOVEL")
    for m in VNEW:
        model.addConstr(z >= VNEW[m], name=f"NOVEL_UB_{VNEW[m]}")
    model.addConstr(z <= model.quicksum(VNEW[m] for m in VNEW), name="NOVEL_LB")
    objective += coverage.profile.major_novel * z
    objective += 0.1 * model.quicksum(VNEW[m] for m in VNEW)
    model.setObjective(objective)
    if debug:
        model.dump(f"{debug}.{gene.name}.major{identifier}.lp")

    # Solve the model
    lookup = {
        **{model.varName(v): a for a, v in VA.items()},
        **{model.varName(v): m for m, v in VNEW.items()},
    }
    result: Dict[Any, MajorSolution] = {}
    debug_info["sol"] = []
    for status, opt, sol in model.solutions(coverage.profile.gap):
        solved_alleles = sorted_tuple(
            [lookup[s][0] for s in sol if s in lookup and s.startswith("A_")]
        )
        novel_muts = sorted_tuple(
            [lookup[s] for s in sol if s in lookup and s.startswith("N_")]
        )
        if (solved_alleles, novel_muts) not in result:
            solution = collections.Counter(
                SolvedAllele(gene, major=a) for a in solved_alleles
            )
            debug_info["sol"].append(
                dict(collections.Counter(a for a in solved_alleles))
            )
            sol = MajorSolution(
                score=opt,
                solution=solution,
                cn_solution=cn_solution,
                added=list(novel_muts),
            )
            log.debug(
                f"[major] status= {status}; opt= {opt:.2f}; "
                + f"solution= {sol._solution_nice()}"
            )
            result[solved_alleles, novel_muts] = sol
    if not result:
        log.debug("[major] solution= []")
    return list(result.values())


def _filter_alleles(
    gene: Gene, coverage: Coverage, cn_solution: CNSolution
) -> Tuple[Dict[str, MajorAllele], Coverage]:
    """
    Filter out low-quality mutations and alleles that are not expressed.

    :returns: Expressed alleles and the coverage of high-quality variants.
    """

    def filter_fns(cov, mut):
        cond = cov.basic_filter(mut, cn=coverage.profile.cn_max)
        if mut.op != "_":
            cond = cond and cov.basic_filter(
                mut, cn=cn_solution.position_cn(mut.pos) + 0.5
            )
        return cond

    cov = coverage.filtered(Coverage.quality_filter)
    cov = cov.filtered(filter_fns)

    alleles = copy.deepcopy(gene.alleles)
    for an, a in natsorted(gene.alleles.items()):
        if a.cn_config not in cn_solution.solution:
            del alleles[an]
        elif any(cov[m] <= 0 for m in a.func_muts):
            s = (f"{gene.get_rsid(m)}" for m in a.func_muts if cov[m] <= 0)
            log.trace("[major] removing {} because of {}", an, " and ".join(s))
            del alleles[an]

    return alleles, cov


def _print_candidates(
    gene,
    alleles: Dict[str, Any],
    coverage: Coverage,
    cn_solution: CNSolution,
    muts: set,
):
    """Pretty-prints the list of allele candidates and their functional mutations."""

    def print_mut(m):
        copies = (
            coverage[m] / (coverage.total(m) / cn_solution.position_cn(m.pos))
            if cn_solution.position_cn(m.pos) and coverage.total(m)
            else 0
        )
        g = gene.region_at(m.pos)
        return (
            f"  {gene.get_rsid(m):12} {str(m):15} "
            + f"{gene.get_refseq(m, from_atg=True):10} "
            + f"(cov={coverage[m]:4}, cn= {copies:3.1f}; "
            + f"region={g[1] if g else '?'}; "
            + f"impact={gene.get_functional(m)}; "
            + ")"
        )

    log.debug("[major] candidate mutations=")
    muts = set(m for a in alleles for m in alleles[a].func_muts) | set(muts)
    for m in sorted(muts):
        log.debug(print_mut(m))
        als = natsorted(
            f"*{a:8}"
            for a, al in gene.alleles.items()
            if m in al.func_muts and "#" not in a
        )
        log.debug(
            "    {}",
            "\n    ".join(" ".join(als[i : i + 6]) for i in range(0, len(als), 6)),
        )
    log.debug("[major] candidate alleles=")
    muts = muts.copy()
    for a in natsorted(alleles):
        log.debug("  *{} (struct= *{})", a, alleles[a].cn_config)
        for m in sorted(alleles[a].func_muts):
            if m in muts:
                muts.remove(m)
            log.debug("  {}", print_mut(m))
    if len(muts) > 0:
        log.debug("  Other mutations:")
        for m in sorted(muts):
            a = (f"*{a}" for a, b in gene.alleles.items() if m in b.func_muts)
            log.debug("  {}, alleles={})", print_mut(m)[:-1], ", ".join(a))
