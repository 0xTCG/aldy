# 786
# Aldy source: minor.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.

from typing import List, Set, Callable, Optional, Tuple, Dict

from . import lpinterface
from .common import log, json_print, allele_sort_key
from .gene import Mutation, Gene
from .cn import MAX_CN
from .major import NOVEL_MUTATION_PENAL
from .coverage import Coverage
from .solutions import MajorSolution, SolvedAllele, MinorSolution
from .diplotype import estimate_diplotype


# Model parameters
MISS_PENALTY_FACTOR = 1.5
"""float: Penalty for each missed minor mutation (0 for no penalty).
          Ideally larger than `ADD_PENALTY_FACTOR` as additions should be cheaper."""

ADD_PENALTY_FACTOR = 1.0
"""float: Penalty for each novel minor mutation (0 for no penalty).
          Zero penalty always prefers mutation additions over coverage errors
          if the normalized SNP slack coverage is >= 50%.
          Penalty of 1.0 prefers additions if the SNP slack coverage is >= 75%."""


def estimate_minor(
    gene: Gene,
    coverage: Coverage,
    major_sols: List[MajorSolution],
    solver: str,
    filter_fn: Optional[Callable] = None,
    max_solutions: int = 1,
    debug: Optional[str] = None,
) -> List[MinorSolution]:
    """
    Estimate the optimal minor star-allele.

    Args:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        coverage (:obj:`aldy.coverage.Coverage`):
            Read coverage data.
        major_sol (:obj:`aldy.solutions.MajorSolution`):
            Major allele solution that will be used for minor star-allele calling.
        solver (str):
            ILP solver. Check :obj:`aldy.lpinterface` for the list of supported solvers.
        filter_fn (callable):
            Custom filtering function.
            Default is ``None``.
        max_solutions (int):
            Maximum number of solutions to report.
            Default is 1.
        debug (str, optional):
            If set, create a "`debug`.minor`identifier`.lp" file for debug purposes.
            Default is ``None``.

    Returns:
        list[:obj:`aldy.solutions.MinorSolution`]
    """

    # Get the list of potential alleles and mutations
    alleles: List[SolvedAllele] = list()
    # Consider all major and minor mutations
    # *from all available major solutions* together
    mutations: Set[Mutation] = set()
    for major_sol in major_sols:
        for (ma, _, added, _) in major_sol.solution:
            alleles += [
                SolvedAllele(ma, mi, added, tuple()) for mi in gene.alleles[ma].minors
            ]
            mutations |= set(gene.alleles[ma].func_muts)
            mutations |= set(added)
            for sa in gene.alleles[ma].minors.values():
                mutations |= set(sa.neutral_muts)

    # Filter out low quality mutations
    def default_filter_fn(mut, cov, total, thres):
        # TODO: is this necessary?
        if mut.op != "_" and mut not in mutations:
            return False
        return Coverage.basic_filter(
            mut, cov, total, thres / MAX_CN
        ) and Coverage.cn_filter(mut, cov, total, thres, major_sol.cn_solution)

    if filter_fn:
        cov = coverage.filtered(filter_fn)
    else:
        cov = coverage.filtered(default_filter_fn)

    minor_sols: List[MinorSolution] = []
    for i, major_sol in enumerate(
        sorted(major_sols, key=lambda s: list(s.solution.items()))
    ):
        minor_sols += solve_minor_model(
            gene, alleles, cov, major_sol, mutations, solver, i, max_solutions, debug
        )
    return minor_sols


def solve_minor_model(
    gene: Gene,
    alleles_list: List[SolvedAllele],
    coverage: Coverage,
    major_sol: MajorSolution,
    mutations: Set[Mutation],
    solver: str,
    identifier: int = 0,
    max_solutions: int = 1,
    debug: Optional[str] = None,
) -> List[MinorSolution]:
    """
    Solves the minor star-allele detection problem via integer linear programming.

    Args:
        gene (:obj:`aldy.gene.Gene`):
            Gene instance.
        alleles_list (list[:obj:`aldy.major.SolvedAllele`]):
            List of candidate minor star-alleles.
        coverage (:obj:`aldy.coverage.Coverage`):
            Sample coverage used to find out the coverage of each mutation.
        major_sol (:obj:`aldy.solutions.MajorSolution`):
            Major star-allele solution to be used for detecting minor star-alleles
            (check :obj:`aldy.solutions.MajorSolution`).
        mutations (set[:obj:`aldy.gene.Mutation`]):
            List of mutations to be considered during the solution build-up
            (all other mutations are ignored).
        solver (str):
            ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.
        identifier (int):
            Unique solution identifier. Used for generating debug information.
            Default is 0.
        max_solutions (int):
            Maximum number of solutions to report.
            Default is 1.
        debug (str, optional):
            If set, Aldy will create "<debug>.minor.lp" file for debug purposes.
            Default is ``None``.

    Returns:
        list[:obj:`aldy.solutions.MinorSolution`]

    Notes:
        Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_
        (section Methods/Genotype refining) for the model explanation.
        Currently returns only the first optimal solution.
    """

    log.debug("*" * 80)
    log.debug("Minor solver: major = {}", major_sol)
    model = lpinterface.model("AldyMinor", solver)

    # Establish minor alleles and their mutations
    alleles: Dict[Tuple[SolvedAllele, int], Set[Mutation]] = {
        (a, 0): set(gene.alleles[a.major].func_muts)
        | set(gene.alleles[a.major].minors[a.minor].neutral_muts)
        | set(a.added)
        for a in alleles_list
    }

    _print_candidates(gene, alleles, major_sol, coverage)

    for a, _ in list(alleles):
        max_cn = major_sol.solution[SolvedAllele(a.major, None, a.added, a.missing)]
        for cnt in range(1, max_cn):
            alleles[a, cnt] = alleles[a, 0]

    VA = {
        a: model.addVar(vtype="B", name=f"A_{a[0].major}_{a[0].minor}_{a[1]}")
        for a in alleles
    }
    for a, cnt in alleles:
        if cnt > 0:
            model.addConstr(VA[a, cnt] <= VA[a, cnt - 1], name=f"CORD_{a.minor}_{cnt}")

    # Make sure that the sum of all subaleles matches the count
    # of the corresponding major alleles
    for sa, cnt in major_sol.solution.items():
        expr = model.quicksum(
            v
            for ((ma, _, ad, mi), _), v in VA.items()
            if SolvedAllele(ma, None, ad, mi) == sa
        )
        model.addConstr(expr == cnt, name=f"CCNT_{sa.major}")

    # Add a binary variable for each allele/mutation pair, where mutation belongs
    # to that allele, that will indicate whether such mutation will be kept or not.
    # Second variable stands for the corresponding VA * VKEEP binary product
    # (binary product transformation).
    VKEEP: dict = {
        a: {
            m: (
                model.addVar(
                    vtype="B", name=f"K_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}"
                ),
                model.addVar(
                    vtype="B",
                    name=f"MUL_K_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}",
                ),
            )
            for m in alleles[a]
        }
        for a in alleles
    }
    # Add a binary variable for each allele/mutation pair, where mutation DOES NOT
    # belong to that allele,that will indicate whether such mutation will be assigned
    # to that allele or not.
    # Second variable stands for the corresponding VA * VNEW product
    # (binary product transformation).
    VNEW: dict = {
        a: {
            m: (
                model.addVar(
                    vtype="B", name=f"N_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}"
                ),
                model.addVar(
                    vtype="B",
                    name=f"MUL_N_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}",
                ),
            )
            for m in mutations
            if gene.has_coverage(a[0].major, m.pos) and m not in alleles[a]
        }
        for a in alleles
    }
    # Add an error variable for each mutation, and populate error constraints
    VERR = {
        m: model.addVar(lb=-model.INF, ub=model.INF, name=f"E_{m.pos}_{m.op}")
        for m in mutations
    }
    constraints = {m: 0 for m in mutations}
    for m in mutations:
        for a in alleles:
            if m in alleles[a]:
                constraints[m] += model.prod(VKEEP[a][m][1], [VA[a], VKEEP[a][m][0]])
            # Add this *only* if CN of this region in a given allele is positive
            elif gene.has_coverage(a[0].major, m.pos):
                constraints[m] += model.prod(VNEW[a][m][1], [VA[a], VNEW[a][m][0]])

    # Fill the constraints for non-variations (i.e. reference genome matches)
    for pos in set(m.pos for m in constraints):
        ref_m = Mutation(pos, "_")  # type: ignore
        VERR[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=f"E_{pos}_REF")
        constraints[ref_m] = 0
        for a in alleles:
            if not gene.has_coverage(a[0].major, pos):
                continue
            # Does this allele contain any mutation at the position `pos`?
            # Insertions are not counted as they always contribute to `_`.
            present_muts = [m for m in alleles[a] if m.pos == pos and m[1][:3] != "INS"]
            assert len(present_muts) < 2
            if len(present_muts) == 1:
                constraints[ref_m] += VA[a] - VKEEP[a][present_muts[0]][1]
            else:
                constraints[ref_m] += VA[a]
                muts = [m for m in VNEW[a] if m.pos == pos and m[1][:3] != "INS"]
                for m in muts:
                    constraints[ref_m] -= VNEW[a][m][1]
                # We ensure that only one additional mutation can be selected here
                model.addConstr(
                    model.quicksum(VNEW[a][m][0] for m in muts) <= 1,
                    name=f"CONE_{pos}_{a[0].major}_{a[0].minor}_{a[1]}",
                )

    # Ensure that each constraint matches the observed coverage
    json_print(debug, "  {")
    json_print(debug, f'    "cn": {str(dict(major_sol.cn_solution.solution))}, ')
    json_print(
        debug,
        '    "major": {}'.format(
            {
                tuple([s.major] + [(m[0], m[1]) for m in s.added])
                if len(s.added) > 0
                else s.major: v
                for s, v in major_sol.solution.items()
            }
        ),
        end=", ",
    )
    json_print(debug, '    "data": {', end="")
    prev = 0
    for m, expr in sorted(constraints.items()):
        scov = coverage.single_copy(m.pos, major_sol.cn_solution)
        # If scov = 0, no mutations should be selected at that locus
        # (enforced by other constraints)
        cov = coverage[m] / scov if scov > 0 else 0
        model.addConstr(expr + VERR[m] == cov, name=f"CCOV_{m.pos}_{m.op}")
        if m.pos != prev and prev != 0:
            json_print(debug, "\n             ", end="")
        prev = m.pos
        json_print(debug, f"({m.pos}, '{m.op}'): {coverage[m]:4}, ", end="")
    json_print(debug, "}, ")

    # Enforce the following rules:
    # 1) Each mutation is assigned only to alleles that are present in the solution
    for a, mv in VKEEP.items():
        for m, v in mv.items():
            model.addConstr(
                v[0] <= VA[a],
                name=f"CVK_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}",
            )
    for a, mv in VNEW.items():
        for m, v in mv.items():
            model.addConstr(
                v[0] <= VA[a],
                name=f"CVN_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}",
            )
    # 2) Each allele must express ALL of its functional mutations
    for a in alleles:
        for m in alleles[a]:
            if m.is_functional:
                assert m in VKEEP[a]
                model.addConstr(
                    VKEEP[a][m][0] >= VA[a],
                    name=f"CFUNC_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}",
                )
    # 3) No allele can include a mutation with coverage 0
    for a in VKEEP:
        for m, v in VKEEP[a].items():
            if not gene.has_coverage(a[0].major, m.pos):
                model.addConstr(
                    v[0] <= 0,
                    name=f"CZERO_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}",
                )
    # 4) No allele can include an extra functional mutation
    #    (that should be done in the major model)
    for a in VNEW:
        for m, v in VNEW[a].items():
            if m.is_functional:
                model.addConstr(
                    v[0] <= 0,
                    name=f"CNEWFUNC_{m.pos}_{m.op}_{a[0].major}_{a[0].minor}_{a[1]}",
                )
    # 5) Avoid extra mutations if there is an existing mutation at the corresponding
    #    locus (either from the definition or added via VNEW)
    for pos in set(m.pos for m in constraints):
        for a in alleles:
            mp = [VKEEP[a][m][1] for m in VKEEP[a] if m.pos == pos]
            ma = [VNEW[a][m][1] for m in VNEW[a] if m.pos == pos]
            # TODO: add support for extra insertions!
            if len(ma) > 1:  # useless if == 1 since variables are binary
                model.addConstr(
                    model.quicksum(ma) <= 1,
                    name=f"CSINGLE_{pos}_{a[0].major}_{a[0].minor}_{a[1]}",
                )
            if len(ma) + len(mp) > 1:
                model.addConstr(
                    model.quicksum(mp + ma) <= 1,
                    name=f"CSINGLEFULL_{pos}_{a[0].major}_{a[0].minor}_{a[1]}",
                )
    # 5) Make sure that each mutation copy has at least one read supporting it
    for m in mutations:
        expr = model.quicksum(VKEEP[a][m][1] for a in alleles if m in VKEEP[a])
        expr += model.quicksum(VNEW[a][m][1] for a in alleles if m in VNEW[a])
        if major_sol.cn_solution.position_cn(m.pos) == 0 or coverage[m] == 0:
            model.addConstr(expr <= 0, name=f"CNOCOV_{m.pos}_{m.op}")
        else:
            model.addConstr(expr <= coverage[m], name=f"CMAXCOV_{m.pos}_{m.op}")
            # Ensure that at least one allele picks an existing non-filtered mutation
            model.addConstr(expr >= 1, name=f"CMINONE_{m.pos}_{m.op}")
    # 6) Do the same for non-mutations
    for pos in set(m.pos for m in constraints):
        expr = 0
        for a in alleles:
            e = [VKEEP[a][m][1] for m in VKEEP[a] if m.pos == pos]
            e += [VNEW[a][m][1] for m in VNEW[a] if m.pos == pos]
            expr += len(e) * VA[a] - model.quicksum(e)
        if isinstance(expr, int) and expr == 0:
            continue
        m = Mutation(pos, "_")
        if major_sol.cn_solution.position_cn(m.pos) == 0:
            model.addConstr(expr <= 0, name=f"CNOCOV_{m.pos}_{m.op}")
        else:
            # If there is no coverage at `pos` at all despite CN being > 0 (or if there
            # is coverage only on a functional mutation that cannot be selected),
            # allow the allele to select a non-existent non-mutation to prevent
            # infeasibility. Should not happen with "sane" datasets...
            model.addConstr(
                expr <= max(major_sol.cn_solution.position_cn(m.pos), coverage[m]),
                name=f"CMAXCOV_{m.pos}_{m.op}",
            )

    # Objective: minimize the absolute sum of errors ...
    objective = model.abssum(v for _, v in VERR.items())
    # ... and penalize the misses ...
    objective += MISS_PENALTY_FACTOR * model.quicksum(
        len(VKEEP[a]) * VA[a] for a in VKEEP
    )
    objective -= MISS_PENALTY_FACTOR * model.quicksum(
        v[1] for a in VKEEP for _, v in VKEEP[a].items()
    )
    # ... and additions ...
    objective += ADD_PENALTY_FACTOR * model.quicksum(
        v[0] for a in VNEW for _, v in VNEW[a].items()
    )
    # ... and novel functional mutations from the major model!
    objective += NOVEL_MUTATION_PENAL * model.quicksum(
        v[0]
        for a in VKEEP
        for m, v in VKEEP[a].items()
        if m.is_functional and m not in gene.alleles[a[0].major].func_muts
    )

    model.setObjective(objective)
    if debug:
        model.dump(f"{debug}.minor{identifier}.lp")

    # Solve the model
    try:
        results = {}
        for status, opt, sol in model.solutions():
            log.debug(f"Minor solver: {status}, opt: {opt:.2f}")

            solution = []
            for allele, value in VA.items():
                if model.getValue(value) <= 0:
                    continue
                added: List[Mutation] = []
                missing: List[Mutation] = []
                for m, mv in VKEEP[allele].items():
                    if not model.getValue(mv[0]):
                        missing.append(m)
                for m, mv in VNEW[allele].items():
                    if model.getValue(mv[0]):
                        added.append(m)
                solution.append(
                    SolvedAllele(
                        allele[0].major,
                        allele[0].minor,
                        allele[0].added + tuple(added),
                        tuple(missing),
                    )
                )
            json_print(
                debug,
                '    "sol": {}, '.format(
                    [
                        (
                            s.minor,
                            [(m[0], m[1]) for m in s.added],
                            [(m[0], m[1]) for m in s.missing],
                        )
                        for s in solution
                    ]
                ),
            )
            sol = MinorSolution(score=opt, solution=solution, major_solution=major_sol)
            _ = estimate_diplotype(gene, sol)
            json_print(debug, f'    "diplotype": "{sol.diplotype}"')
            json_print(debug, "  },")
            log.debug(f"Minor solver: {sol}")
            if str(sol) not in results:
                results[str(sol)] = sol
            if len(results) >= max_solutions:
                break
        return results.values()
    except lpinterface.NoSolutionsError:
        log.debug("Minor solver: no solutions")
        json_print(debug, '    "sol": []')
        json_print(debug, "  },")
        if debug and False:  # Enable to debug infeasible models
            model.model.computeIIS()
            model.dump(f"{debug}.iis.ilp")
        return []


def _print_candidates(gene, alleles, major_sol, coverage):
    """
    Pretty-prints the list of allele candidates and their mutations.
    """

    log.debug("Minor solver: candidates:")
    for a in sorted(alleles, key=lambda x: allele_sort_key(x[0].minor)):
        (ma, mi, _, _), _ = a
        log.debug("  *{} (cn=*{})", mi, gene.alleles[ma].cn_config)
        for m in sorted(alleles[a], key=lambda m: m.pos):
            m_gene, m_region = gene.region_at(m.pos)
            scopy = coverage.single_copy(m.pos, major_sol.cn_solution)
            log.debug(
                "    {:26}  {:.2f} ({:4} / {} * {:4.0f}) {}:{:10} {:>20}",
                str(m),
                coverage[m] / scopy if scopy > 0 else 0,
                coverage[m],
                major_sol.cn_solution.position_cn(m.pos),
                coverage.single_copy(m.pos, major_sol.cn_solution),
                m_gene,
                str(m_region),
                m.aux.get("old", ""),
            )
