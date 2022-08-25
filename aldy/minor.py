# 786
# Aldy source: minor.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.

from typing import List, Set, Tuple, Dict
from natsort import natsorted
import collections
import os

from . import lpinterface
from .common import log, json, Timing
from .gene import Mutation, Gene
from .coverage import Coverage
from .solutions import MajorSolution, SolvedAllele, MinorSolution
from .diplotype import estimate_diplotype


def estimate_minor(
    gene: Gene,
    coverage: Coverage,
    major_sols: List[MajorSolution],
    solver: str,
    max_solutions: int = 1,
    novel: bool = False,
) -> List[MinorSolution]:
    """
    Estimate the optimal minor star-allele.

    :param gene: Gene instance.
    :param coverage: Read coverage instance.
    :param major_sol: Major allele solution for minor star-allele calling.
    :param solver: ILP solver (see :py:mod:`aldy.lpinterface` for supported solvers).
    :param max_solutions: Maximum number of solutions to report.
        Default: 1.
    :param novel: Include non-database functional and silent mutations in the model.
        Default: `False`.
    """

    # Get the list of potential alleles and mutations
    alleles: List[SolvedAllele] = list()
    # Consider all major and minor mutations
    # *from all available major solutions* together
    mutations: Set[Mutation] = set()
    for major_sol in major_sols:
        for sa in major_sol.solution:
            alleles += [
                SolvedAllele(gene, sa.major, mi) for mi in gene.alleles[sa.major].minors
            ]
            mutations |= set(gene.alleles[sa.major].func_muts)
            for minor in gene.alleles[sa.major].minors.values():
                mutations |= set(minor.neutral_muts)
        mutations |= set(major_sol.added)
    mutations |= gene.random_mutations

    # Filter out low quality mutations
    def default_filter_fn(cov, mut):
        # TODO: is this necessary?
        r = gene.region_at(mut.pos)
        if mut.op != "_" and not (
            mut in mutations
            or (r and r[1][0] == "e")
            or (r and r[1] in ["utr3", "utr5", "up"])
        ):
            return False
        cond = cov.basic_filter(mut, cn=coverage.profile.cn_max)
        if mut.op != "_":
            cond = cond and cov.basic_filter(
                mut, cn=major_sol.cn_solution.position_cn(mut.pos) + 0.5
            )
        return cond

    cov = coverage.filtered(Coverage.quality_filter)
    cov = cov.filtered(default_filter_fn)

    if novel:
        for pos, c in cov._coverage.items():
            for m in c:
                if m != "_" and Mutation(pos, m) not in mutations:
                    r = gene.region_at(pos)
                    log.info(
                        "[minor] novel {} ({}; coverage= {:.0f}; func= {})",
                        r[1] if r else "-",
                        cov.percentage(Mutation(pos, m)),
                        gene.is_functional((pos, m)),
                    )
                    mutations.add(Mutation(pos, m))

    # Group by CN solutions
    minor_sols: List[MinorSolution] = []
    cn_sols = {m.cn_solution for m in major_sols}
    min_score = min(m.score for m in major_sols)
    for c in sorted(cn_sols, key=lambda x: x._solution_nice()):
        log.debug("*" * 80)
        majors = [m for m in major_sols if m.cn_solution == c]
        _print_candidates(gene, alleles, c, cov, mutations)
        for major_sol in natsorted(majors, key=lambda s: str(s.solution)):
            sols = solve_minor_model(
                gene,
                cov,
                major_sol,
                alleles,
                mutations,
                solver,
                max_solutions,
            )
            for s in sols:
                s.score += major_sol.score - min_score
            minor_sols += sols
    return minor_sols


def solve_minor_model(
    gene: Gene,
    coverage: Coverage,
    major_sol: MajorSolution,
    alleles_list: List[SolvedAllele],
    mutations: Set[Mutation],
    solver: str,
    max_solutions: int = 1,
) -> List[MinorSolution]:
    """
    Solves the minor star-allele detection problem via integer linear programming.

    :param gene: Gene instance.
    :param coverage: Read coverage instance.
    :param major_sol: Major allele solution for minor star-allele calling.
    :param alleles_list: Candidate minor star-alleles.
    :param mutations: Mutations to consider for model building (all other mutations
        are ignored).
    :param solver: ILP solver (see :py:mod:`aldy.lpinterface` for supported solvers).
    :param max_solutions: Maximum number of solutions to report.
        Default: 1.

    .. note::
        Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_
        (section Methods/Genotype refining) for the model explanation.
        Currently returns only the first optimal solution.
    """

    log.debug("[minor] major= {}", major_sol._solution_nice())
    model = lpinterface.model("AldyMinor", solver)
    debug_info = json[gene.name]["minor"][len(json[gene.name]["minor"])]

    # Establish minor alleles and their mutations
    alleles: Dict[Tuple[SolvedAllele, int], Set[Mutation]] = {
        (a, 0): set(gene.alleles[a.major].func_muts)
        | set(gene.alleles[a.major].minors[a.minor].neutral_muts)
        for a in alleles_list
    }

    for a, _ in list(alleles):
        max_cn = major_sol.solution[SolvedAllele(gene, a.major, "", a.added, a.missing)]
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
            for (vs, _), v in VA.items()
            if (vs.major, vs.added, vs.missing) == (sa.major, sa.added, sa.missing)
        )
        model.addConstr(expr <= cnt, name=f"CCNT_{sa.major}_1")
        model.addConstr(expr >= cnt, name=f"CCNT_{sa.major}_2")
    # Disable other alleles
    model.addConstr(
        model.quicksum(VA.values()) <= sum(major_sol.solution.values()),
        name="CCNT_OTHER",
    )

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
            present_muts = [m for m in alleles[a] if m.pos == pos and m[1][:3] != "ins"]
            assert len(present_muts) < 2
            if len(present_muts) == 1:
                constraints[ref_m] += VA[a] - VKEEP[a][present_muts[0]][1]
            else:
                constraints[ref_m] += VA[a]
                muts = [m for m in VNEW[a] if m.pos == pos and m[1][:3] != "ins"]
                for m in muts:
                    constraints[ref_m] -= VNEW[a][m][1]
                # We ensure that only one additional mutation can be selected here
                model.addConstr(
                    model.quicksum(VNEW[a][m][0] for m in muts) <= 1,
                    name=f"CONE_{pos}_{a[0].major}_{a[0].minor}_{a[1]}",
                )

    # Ensure that each constraint matches the observed coverage
    debug_info["cn"] = dict(major_sol.cn_solution.solution)
    debug_info["major"] = {s.major: v for s, v in major_sol.solution.items()}
    debug_info["data"] = []
    for m, expr in sorted(constraints.items()):
        scov = coverage.single_copy(m, major_sol.cn_solution)
        # If scov = 0, no mutations should be selected at that locus
        # (enforced by other constraints)
        cov = coverage[m] / scov if scov > 0 else 0
        model.addConstr(expr + VERR[m] >= cov, name=f"CCOV_{m.pos}_{m.op}")
        model.addConstr(expr + VERR[m] <= cov, name=f"CCOV_{m.pos}_{m.op}")
        debug_info["data"].append((m.pos, m.op, coverage[m]))

    # TODO: downscale ambiguous mutations (right now score == 1 for all mutations)
    score = {}
    for m, v in VERR.items():
        s = 1
        # ops = coverage._coverage[m.pos].get(m.op, [])
        # if ops:
        #     mx = max(mq for mq, _ in ops)
        #     sx = sum(math.log(1 + mq) / math.log(1 + mx) for mq, _ in ops) / len(ops)
        #     s += 1 - sx
        score[model.varName(v)] = s

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
            if gene.is_functional(m):
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
    # 4) Avoid extra mutations if there is an existing mutation at the corresponding
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

    # 7) Respect phasing
    VPHASEERR = []
    VPHASE = {}
    modes = collections.defaultdict(int)
    if coverage.profile.phase and coverage.sam:
        mut_pos = {m.pos for m in mutations}
        for rr, rv in coverage.sam.phases.items():
            c = []
            for k, v in rv.items():
                if k in mut_pos:
                    c.append((k, v))
            c = sorted(c)
            if len(c) > 1:
                modes[tuple(c)] += 1
        log.debug("[minor] number of phases= {}", len(modes))
        log.debug("[minor] number of alleles= {}", len(alleles))
        if len(modes) * len(alleles) > coverage.profile.minor_phase_vars:
            max_sample = len(modes) * (
                coverage.profile.minor_phase_vars / (len(modes) * len(alleles))
            )
            log.debug("[minor] downsampling to= {}", max_sample)
            skip = len(modes) / max_sample
            if max_sample < len(modes):
                mi = list(modes.items())
                modes = dict(mi[i] for i in range(0, len(mi), int(skip)))
        with Timing("[minor] Model setup"):
            addVar = lambda n: model.addVar(vtype="B", name=n, update=False)  # noqa
            vars = 0
            for ri, (rr, cnt) in enumerate(modes.items()):
                r = dict(rr)
                for ai, a in enumerate(alleles):
                    pos, neg, e = [], [], 0
                    for m in mutations:
                        if m.pos not in r or not gene.has_coverage(a[0].major, m.pos):
                            continue
                        if m in VKEEP[a]:
                            (pos if m.op == r[m.pos] else neg).append(VKEEP[a][m][0])
                        elif m in VNEW[a]:
                            (pos if m.op == r[m.pos] else neg).append(VNEW[a][m][0])
                    if len(pos) + len(neg) > 1:
                        v = VPHASE[ai, ri] = addVar(f"PH_{ai}_{ri}")
                        vars += 1
                        model.addConstr(VPHASE[ai, ri] <= VA[a], name=f"PH_{ai}_{ri}")
                        e = 0
                        for i, vm in enumerate(pos):
                            v_vp = model.prod(addVar(f"PHASE2_{ai}_{ri}_{i}"), [v, vm])
                            vars += 1
                            e += v - v_vp
                        for i, vm in enumerate(neg):
                            v_vp = model.prod(addVar(f"PHASE3_{ai}_{ri}_{i}"), [v, vm])
                            vars += 1
                            e += v_vp
                        VPHASEERR.append(cnt * e)
            for ri, _ in enumerate(modes):
                v = [
                    VPHASE[ai, ri] for ai, _ in enumerate(alleles) if (ai, ri) in VPHASE
                ]
                if v:
                    e = model.quicksum(v)
                    model.addConstr(e <= 1, name=f"PHASE4_{ri}_1")
                    model.addConstr(e >= 1, name=f"PHASE4_{ri}_2")
            model.update()
            log.debug("[minor] phasing setup done, vars= {}", vars)

    # Objective: minimize the absolute sum of errors ...
    objective = 0
    o_error = model.abssum((v for v in VERR.values()), coeffs=score)
    objective += o_error
    # ... and penalize the misses ...
    o_penal = coverage.profile.minor_miss * model.quicksum(
        len(VKEEP[a]) * VA[a] for a in VKEEP
    )
    o_penal -= coverage.profile.minor_miss * model.quicksum(
        v[1] for a in VKEEP for _, v in VKEEP[a].items()
    )
    # ... and additions ...
    cnt = 0
    for a in VNEW:
        for _, v in VNEW[a].items():
            # HACK:
            # Add cnt/10000 to select the smallest allele if there is a tie
            o_penal += coverage.profile.minor_add * (1 + cnt / 1000000) * v[0]
            cnt += 1
    # ... and novel functional mutations from the major model...
    for m in {m for a in VNEW for m in VNEW[a]}:
        vars = [
            VNEW[a][m][0]
            for a in VNEW
            if m in VNEW[a]
            if gene.is_functional(m)
            if m not in gene.alleles[a[0].major].func_muts
        ]
        if vars:
            vo = model.addVar(vtype="B", name=f"VNEWOR_{m.pos}_{m.op}")
            model.addConstr(vo <= model.quicksum(vars), name=f"VNEWOR1_{m.pos}_{m.op}")
            for vi, v in enumerate(vars):
                model.addConstr(vo >= v, name=f"VNEWOR2_{m.pos}_{m.op}_{vi}")
            o_penal += coverage.profile.minor_add / 2 * vo
    objective += o_penal
    # ... and phasing error!
    o_phase = coverage.profile.minor_phase * model.quicksum(VPHASEERR)
    objective += o_phase

    model.setObjective(objective)

    # Solve the model
    results = {}
    with Timing("[minor] Model solution"):
        for status, opt, sol in model.solutions():
            solution = []
            for allele, value in VA.items():
                if model.getValue(value) <= 0:
                    continue

                if coverage.profile.phase and coverage.sam:
                    ai = next(i for i, a in enumerate(alleles) if allele == a)
                    _print_phase(
                        model, gene, mutations, modes, VPHASE, VKEEP, VNEW, allele, ai
                    )

                added: List[Mutation] = []
                missing: List[Mutation] = []
                for m, mv in VKEEP[allele].items():
                    if not model.getValue(mv[0]):
                        missing.append(m)
                for m, mv in VNEW[allele].items():
                    if model.getValue(mv[0]):
                        added.append(m)
                    else:
                        copies = coverage.single_copy(m, major_sol.cn_solution)
                        copies = coverage[m] / copies if copies > 0 else 0
                        if abs(copies - major_sol.cn_solution.max_cn()) > 1e-5:
                            continue
                        # HACK: add homozygous mutation to _all_ alleles
                        # Works only if a mutation is unambiguiusly homozygous;
                        # otherwise, it will fall back to the model defaults
                        if m not in alleles[allele]:
                            added.append(m)
                solution.append(
                    SolvedAllele(
                        gene,
                        allele[0].major,
                        allele[0].minor,
                        allele[0].added + added,
                        missing,
                    )
                )

            sol = MinorSolution(score=opt, solution=solution, major_solution=major_sol)
            _ = estimate_diplotype(gene, sol)
            debug_info["sol"] = [(s.minor, s.added, s.missing) for s in solution]
            debug_info["diplotype"] = sol.diplotype
            log.debug(
                f"[minor] status= {status}; opt= {opt:.2f} "
                + "(error= {:.2f}, penal={:.2f}, phase= {:.2f}) "
                + f"solution= {sol._solution_nice()}",
                model.getValue(o_error),
                model.getValue(o_penal),
                model.getValue(o_phase),
            )
            if str(sol) not in results:
                results[str(sol)] = sol
            if len(results) >= max_solutions:
                break
    if not results:
        log.debug("[minor] solution= []")
        debug_info["sol"] = []
        if "ALDY_IIS" in os.environ:  # Set to debug infeasible models
            model.model.computeIIS()
            model.dump("minor.iis.ilp")
    return sorted(results.values(), key=lambda x: str(x.get_minor_diplotype()))


def _print_candidates(gene, alleles, cn_sol, coverage, muts):
    """Pretty-print the list of allele candidates and their mutations."""

    def print_mut(m):
        copies = coverage.single_copy(m, cn_sol)
        copies = coverage[m] / copies if copies > 0 else 0
        impact = gene.get_functional(m)
        return (
            f"  {gene.get_rsid(m):12} {str(m):15} "
            + f"{gene.get_refseq(m, from_atg=True):10} "
            + f"(cov={coverage[m]:4}, cn= {copies:3.1f}; "
            + (f"impact={impact})" if impact else "")
        )

    log.debug("[minor] candidate mutations=")
    for m in sorted(muts):
        if coverage[m]:
            log.debug(print_mut(m))

    log.trace("[minor] candidate alleles=")
    muts = muts.copy()
    for sa in natsorted(alleles, key=lambda x: (x.major, x.minor)):
        am = (
            set(gene.alleles[sa.major].func_muts)
            | set(gene.alleles[sa.major].minors[sa.minor].neutral_muts)
            | set(sa.added)
        )
        miss = sum(1 for a in am if coverage[a] == 0)
        log.trace(f"  *{sa.minor} (major= *{sa.major}; miss= {miss})")
        for m in sorted(am):
            if m in muts:
                muts.remove(m)
            log.trace("  {}", print_mut(m))
    if len(muts) > 0:
        log.trace("  Other mutations:")
        for m in sorted(muts):
            log.trace("  {}", print_mut(m))


def _print_phase(model, gene, mutations, modes, VPHASE, VKEEP, VNEW, allele, ai):
    """Pretty-print phasing results."""
    rds = [
        r
        for ri, r in enumerate(modes.items())
        if (ai, ri) in VPHASE
        if model.getValue(VPHASE[ai, ri])
    ]
    log.trace("[phase] {} (reads= {})", allele[0].minor, len(rds))
    for m in sorted(mutations):
        sm = f"  {m} "
        if m in gene.alleles[allele[0].major].func_muts:
            sm += "*"
        elif m in gene.alleles[allele[0].major].minors[allele[0].minor].neutral_muts:
            sm += "#"
        else:
            sm += " "

        if m in VKEEP[allele]:
            sm += "K" if model.getValue(VKEEP[allele][m][0]) else "-"
        elif m in VNEW[allele]:
            sm += "N" if model.getValue(VNEW[allele][m][0]) else "-"
        else:
            sm += " "

        ph = collections.defaultdict(int)
        for r in rds:
            dr = dict(r[0])
            if m.pos in dr:
                c = dr[m.pos]
                if ">" in c:
                    c = c[2]
                elif c.startswith("del"):
                    c = f"-{len(c)-3}"
                ph[c] += r[1]
        if len(ph):
            p = sorted(ph.items())
            log.trace(
                "[phase]    {}\t{}",
                sm,
                "; ".join(f"{o} ({n})" for o, n in p),
            )
