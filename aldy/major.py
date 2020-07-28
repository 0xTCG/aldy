# 786
# Aldy source: major.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict, Tuple, Any, Optional

import collections
import copy

from natsort import natsorted

from . import lpinterface
from .common import log, AldyException, json_print, sorted_tuple, allele_sort_key
from .cn import MAX_CN
from .gene import MajorAllele, Mutation, Gene
from .coverage import Coverage
from .solutions import CNSolution, MajorSolution, SolvedAllele


# Model parameters
NOVEL_MUTATION_PENAL = MAX_CN + 1
"""float: Penalty for each novel functional mutation (0 for no penalty).
   Should be large enough to prevent novel mutations unless really necessary."""


def estimate_major(
    gene: Gene,
    coverage: Coverage,
    cn_solution: CNSolution,
    solver: str,
    gap: float = 0,
    identifier: int = 0,
    debug: Optional[str] = None,
) -> List[MajorSolution]:
    """
    Estimate optimal major star-alleles.

    Args:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        coverage (:obj:`aldy.coverage.Coverage`):
            Read coverage data.
        cn_solution (:obj:`aldy.solutions.CNSolution`):
            Copy-number solution that will be used for major star-allele calling.
        solver (str):
            ILP solver. Check :obj:`aldy.lpinterface` for the list of supported solvers.
        gap (float):
            Relative optimality gap. Use non-zero values to allow non-optimal solutions.
            Default is 0 (reports only optimal solutions).
        identifier (int):
            Unique solution identifier. Used for generating debug information.
            Default is 0.
        debug (str, optional):
            If set, create a "`debug`.major`identifier`.lp" file for debug purposes.
            Default is ``None``.

    Returns:
        list[:obj:`aldy.solutions.MajorSolution`]
    """

    log.debug("*" * 80)
    log.debug("[major] cn= {}", cn_solution._solution_nice())

    if sum(cn_solution.solution.values()) < 2:
        raise AldyException(
            "estimate_major requires at least two valid gene configurations"
        )

    alleles, coverage = _filter_alleles(gene, coverage, cn_solution)
    # Check if some CN solution has no matching allele
    if set(cn_solution.solution) - set(a.cn_config for a in alleles.values()):
        results: List[MajorSolution] = []
    else:
        results = solve_major_model(
            gene, alleles, coverage, cn_solution, solver, gap, identifier, debug
        )

    # TODO: re-implement phasing step from Aldy 1.2
    # TODO: Check for novel functional mutations and do something with them

    return results


def solve_major_model(
    gene: Gene,
    allele_dict: Dict[str, MajorAllele],
    coverage: Coverage,
    cn_solution: CNSolution,
    solver: str,
    gap: float = 0,
    identifier: int = 0,
    debug: Optional[str] = None,
) -> List[MajorSolution]:
    """
    Solves the major star-allele detection problem via integer linear programming.

    Args:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        allele_dict (dict[str, :obj:`aldy.gene.MajorAllele`]):
            Dictionary of the candidate major star-alleles.
        coverage (:obj:`aldy.coverage.Coverage`):
            Mutation coverage data.
        cn_solution (:obj:`aldy.solutions.CNSolution`):
            Copy-number solution that will be used for calling major star-alleles
            (:obj:`aldy.solutions.CNSolution`).
        solver (str):
            ILP solver. Check :obj:`aldy.lpinterface` for the list of supported solvers.
        gap (float):
            Relative optimality gap. Use non-zero values to allow non-optimal solutions.
            Default is 0 (reports only optimal solutions).
        identifier (int):
            Unique solution identifier. Used for generating debug information.
            Default is 0.
        debug (str, optional):
            If set, create a "`debug`.major`identifier`.lp" file for debug purposes.
            Default is ``None``.

    Returns:
        list[:obj:`aldy.solutions.MajorSolution`]

    Notes:
        Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_
        (section Methods/Major star-allele identification) for the model explanation.
    """

    model = lpinterface.model("AldyMajor", solver)

    # Get the list of _all_ functional mutations present in the sample
    # and the database (intersection)
    func_muts = {
        M
        for m, M in gene.mutations.items()
        if gene.is_functional(M) and coverage[M] > 0
    }
    _print_candidates(gene, allele_dict, coverage, cn_solution, func_muts)

    # HACK: silence type checker
    a: Any = 0

    # Create a binary variable for every possible allele copy
    alleles = {(a, 0): allele_dict[a] for a in allele_dict}
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
    # Add a binary variable for each allele/novel mutation pair
    # For each binary variable, add a multiplication variable
    # MUL = VNEW[a, m] * VA[a] (binary product linearization)
    del_allele = gene.deletion_allele()
    VNEW = {
        a: {
            m: (
                model.addVar(vtype="B", name=f"N_{m.pos}_{m.op}_{a[0]}_{a[1]}"),
                model.addVar(vtype="B", name=f"MUL_N_{m.pos}_{m.op}_{a[0]}_{a[1]}"),
            )
            if a[0] != del_allele
            else (0, 0)  # deletion alleles should not be assigned any mutations
            for m in constraints
            if m[0] not in set(mm[0] for mm in alleles[a].func_muts)
            and gene.has_coverage(a[0], m.pos)
        }
        for a in alleles
    }
    # Populate the constraints
    for a in alleles:
        for m in alleles[a].func_muts:
            constraints[m] += VA[a]
    # Add novel mutation constraints
    # Also ensure that MUL = 1 <=> VA * VNEW = 1
    for a in VNEW:
        for m, (vn, vm) in VNEW[a].items():
            constraints[m] += model.prod(vm, [VA[a], vn])

    # Populate constraints of non-variations (i.e. matches with the reference genome)
    for pos in set(m.pos for m in constraints):
        ref_m = Mutation(pos, "_")  # type: ignore
        constraints[ref_m] = 0
        VERR[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=f"E_{pos}_REF")
        for a in alleles:
            if not gene.has_coverage(a[0], pos):
                continue
            # An insertion still contributes to the total coverage at this loci
            if any(
                ma[0] == pos and not ma[1][:3] == "INS" for ma in alleles[a].func_muts
            ):
                continue
            # Make sure that the allele has no novel mutations at this position
            constraints[ref_m] += VA[a]
            muts = [m for m in VNEW[a] if m[0] == pos and m[1][:3] != "INS"]
            for m in muts:
                constraints[ref_m] -= VNEW[a][m][1]
            # We ensure that only one additional mutation can be selected here
            model.addConstr(
                model.quicksum(VNEW[a][m][0] for m in muts) <= 1,
                name=f"CONE_{pos}_{a[0]}_{a[1]}",
            )

    # Each allele must express all of its functional mutations
    json_print(debug, "  { # {}", identifier)
    json_print(debug, f'    "cn": {str(dict(cn_solution.solution))}, ')
    json_print(debug, '    "data": {', end="")
    prev = 0
    for m, expr in sorted(constraints.items()):
        if coverage.single_copy(m.pos, cn_solution) == 0:
            cov = 0.0
        else:
            cov = coverage[m] / coverage.single_copy(m.pos, cn_solution)
        model.addConstr(expr + VERR[m] == cov, name=f"CFUNC_{m.pos}_{m.op}")
        if m.pos != prev and prev != 0:
            json_print(debug, "\n             ", end="")
        prev = m.pos
        json_print(debug, f"({m[0]}, '{m[1]}'): {cov:.4f}, ", end="")
    json_print(debug, "}, ")

    # Each CN configuration must be satisfied by corresponding alleles
    for cnf, cnt in cn_solution.solution.items():
        expr = sum(VA[a] for a in VA if alleles[a].cn_config == cnf)
        model.addConstr(expr == cnt, name=f"CSAT_{cnf}")

    # Each functional mutation must be chosen by some allele
    for m in func_muts:
        expr = model.quicksum(VA[a] for a in alleles if m in alleles[a].func_muts)
        expr += model.quicksum(
            VNEW[a][m][1]
            for a in alleles
            if m in VNEW[a] and m not in alleles[a].func_muts
        )
        model.addConstr(expr >= 1, name=f"CEXP_{m.pos}_{m.op}")

    # Objective: minimize the absolute sum of errors and the number of novel mutations
    objective = model.abssum(e for e in VERR.values())
    objective += NOVEL_MUTATION_PENAL * model.quicksum(
        VNEW[a][m][0] for a in VNEW for m in VNEW[a]
    )
    model.setObjective(objective)
    if debug:
        model.dump(f"{debug}.major{identifier}.lp")

    # Solve the model
    try:
        lookup = {
            **{model.varName(v): a for a, v in VA.items()},
            **{model.varName(v): (a, m) for a in VNEW for m, (v, _) in VNEW[a].items()},
        }
        result: Dict[Any, MajorSolution] = {}
        json_print(debug, '    "sol": [', end="")
        for status, opt, sol in model.solutions(gap):
            log.debug(f"[major] status= {status}; opt= {opt:.2f}")
            # MajorAllele: novel mutations
            solved_alleles: Any = collections.defaultdict(lambda: [])
            for s in sol:
                if s not in lookup:
                    continue
                elif isinstance(lookup[s][0], tuple):  # Novel mutation
                    solved_alleles[lookup[s][0]].append(lookup[s][1])
                elif s not in solved_alleles:
                    solved_alleles[lookup[s]] = []
            sol_tuple = sorted_tuple(
                (a, sorted_tuple(nm)) for (a, _), nm in solved_alleles.items()
            )
            if sol_tuple not in result:
                solution = collections.Counter(
                    SolvedAllele(major=a, minor=None, added=tuple(mut), missing=tuple())
                    for (a, _), mut in solved_alleles.items()
                )
                json_print(
                    debug,
                    dict(
                        collections.Counter(
                            tuple([a] + [(m[0], m[1]) for m in mut])
                            if len(mut) > 0
                            else a
                            for (a, _), mut in solved_alleles.items()
                        )
                    ),
                    end=", ",
                )
                sol = MajorSolution(
                    score=opt, solution=solution, cn_solution=cn_solution
                )
                log.debug("[major] solution= {}", sol._solution_nice())
                result[sol_tuple] = sol
        json_print(debug, "]\n  }, ", end="")
        return list(result.values())
    except lpinterface.NoSolutionsError:
        log.debug("[major] solution= []")
        return []


def _filter_alleles(
    gene: Gene, coverage: Coverage, cn_solution: CNSolution
) -> Tuple[Dict[str, MajorAllele], Coverage]:
    """
    Filter out low-quality mutations and alleles that are not expressed.

    Returns:
        tuple[dict[str, :obj:`aldy.gene.MajorAllele`], :obj:`aldy.coverage.Coverage`]:
        Tuple of allele dictionary describing the feasible alleles,
        and the coverage description of high-confidence variants.
    """

    def filter_fns(mut, cov, total, thres):
        z = Coverage.basic_filter(
            mut, cov, total, thres / MAX_CN
        ) and Coverage.cn_filter(mut, cov, total, thres, cn_solution)
        return z

    cov = coverage.filtered(filter_fns)
    alleles = copy.deepcopy(gene.alleles)
    for an, a in sorted(gene.alleles.items()):
        if a.cn_config not in cn_solution.solution:
            del alleles[an]
        elif any(cov[m] <= 0 for m in a.func_muts):
            s = (
                "{} in {}".format(m, gene.region_at(m.pos))
                for m in a.func_muts
                if cov[m] <= 0
            )
            log.trace("Major solver: removing {} because of {}", an, " and ".join(s))
            del alleles[an]

    return alleles, cov


def _print_candidates(
    gene,
    alleles: Dict[str, Any],
    coverage: Coverage,
    cn_solution: CNSolution,
    muts: set,
) -> None:
    """
    Pretty-prints the list of allele candidates and their functional mutations.
    """
    log.debug("[major] candidates=")
    muts = muts.copy()
    for a in natsorted(alleles):
        log.debug("  *{} (struct= *{})", a, alleles[a].cn_config)
        for m in sorted(alleles[a].func_muts):
            if m in muts:
                muts.remove(m)
            copies = (
                coverage[m] / (coverage.total(m.pos) / cn_solution.position_cn(m.pos))
                if cn_solution.position_cn(m.pos) and coverage.total(m.pos)
                else 0
            )
            log.debug(
                f"    {coverage[m]:4} (cn= {copies:3.1f}) {str(m):20}  "
                + f"{gene.get_dbsnp(m):10} {gene.get_refseq(m, from_atg=True)}",
            )
    if len(muts) > 0:
        log.debug("  Other mutations:")
        for m in sorted(muts):
            copies = (
                coverage[m] / (coverage.total(m.pos) / cn_solution.position_cn(m.pos))
                if cn_solution.position_cn(m.pos) and coverage.total(m.pos)
                else 0
            )
            a = (f"*{a}" for a, b in gene.alleles.items() if m in b.func_muts)
            log.debug(
                f"    {coverage[m]:4} (cn= {copies:3.1f}) {str(m):20}  "
                + f"{gene.get_dbsnp(m):10} {gene.get_refseq(m, from_atg=True):12}"
                + f" ({', '.join(a)})",
            )
