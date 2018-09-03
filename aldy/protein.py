# 786

# Aldy source: protein.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict, Tuple, Set

import collections
import itertools
import functools
import multiprocessing
import copy

from . import lpinterface
from .common import *
from .cn import CNSolution, MAX_CN
from .gene import Allele, Mutation, Gene
from .sam import Sample, Coverage


"""
Describes a potential (possibly optimal) major star-allele configuration.
Members:
   score (float):
      ILP model error score (0 for user-provided solutions)
   solution (dict of str: int):
      dict of major star-alleles where a value denotes the copy-number
      of each star-allele (e.g. {1: 2} means that we have two copies
      of *1).
   cn_solution (cn.CNSolution):
      associated copy-number solution used for calculating the major 
      star-allele assignment
   novel (dict of str: list of Mutation):
      dict of major-star alleles where a value describes novel functional
      mutations assigned to that major star-allele.
"""
MajorSolution = nt('MajorSolution', 'score solution cn_solution novel'.split())


def estimate_major(gene: Gene, sam: Sample, cn_solution: CNSolution, 
                   solver: str) -> List[MajorSolution]:
   """
   Entry-point to the star-allele detection solver.

   Params:
      gene (gene.Gene): gene description
      sam (sam.Sample): sample description
      cn_solution (cn.CNSolution): copy-number profile of the sample
      solver (str): ILP solver to use. Check lpinterface.py for available solvers
   Returns:
      list of MajorSolution objects describing the optimal solutions to the ILP.
   """

   log.debug('Solving major alleles for cn={}', cn_solution)
   
   # Case of two deletions
   if len(cn_solution.solution) == 0: 
      del_allele = next(a for a in gene.cn_configs if a.kind == gene.CNConfig.CNConfigType.DELETION)
      sol = MajorSolution(score=0, solution={del_allele: 2}, cn_solution=cn_solution)
      return [sol]
   
   alleles, coverage = _filter_alleles(gene, sam, cn_solution)
   # Check if some CN solution has no matching allele
   if set(cn_solution.solution) - set(a.cn_config for a in alleles.values()):
      return MajorSolution(score=float('inf'), solution=[], cn_solution=cn_solution) 

   # TODO: re-implement phasing step from Aldy 1.4   
   # TODO: Check for novel functional mutations and do something with them
   # novel_functional_mutations = _get_novel_mutations(gene, coverage, cn_solution)

   results = solve_major_model(alleles, coverage, cn_solution, solver)
   return results


def solve_major_model(alleles: Dict[str, Allele], coverage: Coverage, 
                      cn_solution: CNSolution, solver: str) -> List[MajorSolution]:
   """
   Solves the major star-allele detection problem via ILP.

   Params:
      alleles (dict of str: gene.Allele):
         dictionary describing candidate alleles. 
      coverage (coverage.Coverage):
         sample coverage that is used to find out the coverage 
         of each major mutation
      cn_solution (cn.CNSolution):
         copy-number profile of the sample
      solver (str): 
         ILP solver to use. Check lpinterface.py for available solvers
   Returns:
      list of MajorSolution objects describing the optimal solutions to the ILP.
   """

   # Model parameters
   # Make sure that each novel mutation gets heavily penalized
   NOVEL_MUTATION_PENAL = 100000

   # Make sure that coverage defaults to 0 on empty values
   model = lpinterface.model('aldy_major_allele', solver)
   _print_candidates(alleles, coverage, cn_solution)
   
   # Create a binary variable for all possible allele copies
   alleles = {(a, 0): alleles[a] for a in alleles}
   for (an, _), a in list(alleles.items()):
      max_cn = cn_solution.solution[a.cn_config]
      log.trace('Maximum CN for {}: {}', an, max_cn)
      for i in range(1, max_cn):
         alleles[(an, i)] = alleles[(an, 0)]
   A = {a: model.addVar(vtype='B', name='A_{}_{}'.format(*a)) for a in alleles}

   # Make sure that A[i+1] <= A[i] (to avoid equivalent solutions)
   for a, ai in alleles:
      if ai > 0:
         log.trace('LP contraint: A_{}_{} <= A_{}_{}', a, ai, a, ai - 1)
         model.addConstr(A[(a, ai)] <= A[(a, ai - 1)])
   
   # Add an error variable to the ILP for any mutation
   error_vars = {m: model.addVar(lb=-model.INF, ub=model.INF, name='MA_{}_{}_{}'.format(m, *a))
      for a in alleles
      for m in alleles[a].func_muts}
   constraints = {e: 0 for e in error_vars}
   # Add a binary variable for any allele/novel mutation pair
   M = {a: {m: model.addVar(vtype='B', name='EXTRA_{}_{}_{}'.format(m, *a))
         for m in constraints
         if m not in alleles[a].func_muts} 
      for a in alleles}
   # Populate constraints
   for a in alleles:
      for m in alleles[a].func_muts:
         cov = max(1, coverage.total(m.pos)) / cn_solution.position_cn(m.pos) \
            if cn_solution.position_cn(m.pos) > 0 else 0
         constraints[m] += cov * A[a]
   # Add novel mutation constraints
   for a in M:
      for m in M[a]:
         if cn_solution.position_cn(m.pos) == 0:
            continue
         cov = max(1, coverage.total(m.pos)) / cn_solution.position_cn(m.pos) \
            if cn_solution.position_cn(m.pos) > 0 else 0
         constraints[m] += cov * A[a] * M[a][m]
   
   # Populate constraints of non-variations (i.e. matches with the reference genome)
   for m in list(constraints):
      if m.op[:3] == 'INS':
         continue
      
      ref_m = Mutation(m.pos, op='REF')
      if ref_m in constraints:
         ref_m = Mutation(m.pos, op=ref_m.op + '#') # TODO: check if needed
      if ref_m not in constraints:
         constraints[ref_m] = 0
         error_vars[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=str(ref_m))

      cov = max(1, coverage.total(m.pos)) / cn_solution.position_cn(m.pos) \
         if cn_solution.position_cn(m.pos) > 0 else 0
      for a in alleles:
         constraints[ref_m] += cov * A[a]

   # Make sure that each constraint equals the observed coverage with some error
   for m, expr in constraints.items():
      log.trace('LP contraint: {} == {} + err for {} with cn={}', coverage[m], expr, m, cn_solution.position_cn(m.pos))
      model.addConstr(expr + error_vars[m] == coverage[m])

   # Make sure that the CN for each CN config matches the CN solution
   for cnf, cnt in cn_solution.solution.items():
      expr = sum(A[a] for a in A if alleles[a].cn_config == cnf)
      log.trace('LP contraint: {} == {} for {}', cnt, expr, cnf)
      model.addConstr(expr == cnt)

   # Each allele must express all of its functional mutations
   func_muts = (m for a in alleles for m in alleles[a].func_muts if coverage[m] > 0)
   for m in func_muts:
      expr = model.quicksum(A[a] for a in alleles if m in alleles[a].func_muts)
      expr += model.quicksum(A[a] * M[a][m] for a in alleles if m not in alleles[a].func_muts)
      log.trace('LP contraint: {} >= 1 for {}', expr, m)
      model.addConstr(expr >= 1)

   # Set objective: minimize the absolute sum of errors   
   objective = \
      model.abssum(e for e in error_vars.values()) + \
      NOVEL_MUTATION_PENAL * model.quicksum(M[a][m] for a in M for m in M[a])
   log.trace('LP objective: {}', objective)

   # Solve the ILP
   try:
      status, opt, solutions = model.solveAll(objective, 
         dict(list(A.items()) + [((a, m), M[a][m]) for a in M for m in M[a]]))
      log.debug('CN Solver status: {}, opt: {}', status, opt)
   except lpinterface.NoSolutionsError:
      return MajorSolution(score=float('inf'), solution=[], cn_solution=cn_solution)

   result = []
   print(pr(solutions))
   for sol in solutions:
      sol_alleles = collections.Counter(ma for ma, _ in sol)
      result.append(MajorSolution(score=opt, solution=sol_alleles, cn_solution=cn_solution, novel={}))
      # TODO: for (a, _), extra in sd.items():
         # if len(extra) > 0:
            # rai
            # log.warn('Novel major star-allele (*{}-like) found!', a)
            # extra = [Mutation(m.pos, m.op, m.functional, dict(list(m.aux.items()) + [('novel', True)])) for m in extra]
            # an = Allele(
            #    gene.alleles[a].name + '+{}'.format(unique_key[gene.alleles[a].name]),
            #    gene.alleles[a].cnv_configuration,
            #    gene.alleles[a].functional_mutations | set(extra),
            #    gene.alleles[a].suballeles
            # )
            # unique_key[gene.alleles[a].name] = 'a' if unique_key[gene.alleles[a].name] == '' else chr(ord(unique_key[gene.alleles[a].name]) + 1)
            # gene.alleles[an.name] = an
            # solutions[si].append(an.name)
   
   return result


def _filter_alleles(gene: Gene, sam: Sample, 
                    cn_solution: CNSolution) -> Tuple[Dict[str, Allele], Coverage]:
   """
   Filters out all low-quality mutations and impossible alleles.
   Returns:
      tuple of allele dictionary describing feasible alleles and 
      Coverage object describing high-confidence variants
   """
   
   def filter_fns(mut, cov, total, thres):
      return Coverage.basic_filter(mut, cov, total, thres / MAX_CN) and \
             Coverage.cn_filter(mut, cov, total, thres, cn_solution)
   cov = sam.coverage.filtered(filter_fns)
   alleles = copy.deepcopy(gene.alleles)
   for an, a in sorted(gene.alleles.items()):
      if a.cn_config not in cn_solution.solution:
         del alleles[an]
      elif any(cov[m] <= 0 for m in a.func_muts):
         s = ('{} in {}'.format(m, gene.region_at(m.pos))
            for m in a.func_muts
            if cov[m] <= 0)
         log.trace('Removing {} because of {}', an, ' and '.join(s))
         del alleles[an]
   
   return alleles, cov


def _get_novel_mutations(gene: Gene, coverage: Coverage, 
                         cn_solution: CNSolution) -> Set[Mutation]:
   """
   Calculates a set of expressed major functional mutations that are not present in the database.
   Returns the set of Mutation objects.
   TODO: integrate it into the model.
   """

   # Require AT LEAST 80% coverage per copy for a nover mutation
   MIN_COVERAGE_PER_COPY = 80.0

   result = set()
   for pos, muts in coverage._coverage.items():
      for op, cov in muts.items():
         if op == '_' or (pos, op) in gene.mutations:
            continue
         try:
            _, region = gene.region_at(pos)
         except KeyError:
            continue
         # TODO: handle non-unique regions as well (remapping)
         if region not in gene.unique_regions:
            continue
         cn = cn_solution.position_cn(pos)
         if cn == 0 or coverage.percentage(Mutation(pos, op)) < MIN_COVERAGE_PER_COPY / cn:
            continue
         if gene.check_functional(Mutation(pos, op)):
            log.debug('Novel mutation: {} {} {} ({} or {}%)', gene.region_at(pos), pos, op, cov, coverage.percentage(Mutation(pos, op)))
            result.add(Mutation(pos, op, is_functional=True))
   return result


def _print_candidates(alleles: Dict[str, Allele], coverage: Coverage, 
                      cn_solution: CNSolution) -> None:
   """Prints the list of allele candidates and their functional mutations"""

   log.debug('Possible candidates:')
   for a in sorted(alleles, key=allele_sort_key):
      log.debug('  *{} (cn=*{})', a, alleles[a].cn_config)
      for m in sorted(alleles[a].func_muts, key=lambda m: m.pos):
         log.debug('    {} {:4} ({:.1f} copies) {} {}',
            #coverage.region_at(m.pos),
            m, coverage[m], 
            coverage[m] / (coverage.total(m.pos) / cn_solution.position_cn(m.pos)),
            'F', m.aux.get('old', ''))
