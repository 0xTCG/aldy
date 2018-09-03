# 786

# Aldy source: refiner.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict, Tuple, Set
from functools import reduce, partial

import collections
import multiprocessing

from . import lpinterface
from .common import *
from .gene import Mutation, Gene, Allele, Suballele
from .sam import Sample
from .cn import MAX_CN
from .coverage import Coverage
from .protein import MajorSolution

"""
Describes a potential (possibly optimal) minor star-allele configuration.
Members:
   score (float):
      ILP model error score (0 for user-provided solutions)
   solution (dict of str: int):
      dict of minor star-alleles where a value denotes the copy-number
      of each star-allele (e.g. {1: 2} means that we have two copies
      of *1).
   major_solution (major.MajorSolution):
      associated major star-allele solution used for calculating the minor 
      star-allele assignment
"""
MinorSolution = nt('MinorSolution', 'score solution major_solution missing added'.split())


def estimate_minor(gene: Gene, sam: Sample, major_sol: MajorSolution, solver: str) -> MinorSolution:
   """
   """

   # Get list of alleles and mutations to consider   
   alleles = list()
   mutations = set()
   for ma in major_sol.solution:
      alleles += [(ma, mi) for mi in gene.alleles[ma].minors]
      mutations |= set(gene.alleles[ma].func_muts)
      for sa in gene.alleles[ma].minors.values():
         mutations |= set(sa.neutral_muts)
   
   # Filter out "bad" mutations
   def filter_fns(mut, cov, total, thres):
      # TODO: is this necessary?
      if mut.op != '_' and mut not in mutations: 
         return False
      return Coverage.basic_filter(mut, cov, total, thres / MAX_CN) and \
             Coverage.cn_filter(mut, cov, total, thres, major_sol.cn_solution) 
   cov = sam.coverage.filtered(filter_fns)

   return solve_minor_model(gene, alleles, cov, major_sol, mutations, solver)


def solve_minor_model(gene: Gene,
                      alleles: Set[Tuple[str, str]], 
                      coverage: Coverage, 
                      major_sol: MajorSolution, 
                      mutations: Set[Mutation], 
                      solver: str) -> MinorSolution:
   """
   """

   MISS_PENALTY_FACTOR = 2
   ADD_PENALTY_FACTOR = 2

   log.debug('Refining {}', major_sol)
   model = lpinterface.model('aldy_refine', solver)
   
   # Establish minor allele binary variables
   alleles = {(ma, mi, 0): set(gene.alleles[ma].func_muts) | \
                           set(gene.alleles[ma].minors[mi].neutral_muts)
              for ma, mi in alleles}
   for ma, mi, _ in list(alleles):
      for cnt in range(1, major_sol.solution[ma]):
         alleles[(ma, mi, cnt)] = alleles[(ma, mi, 0)]

   A = {a: model.addVar(vtype='B', name='{}-{}_{}'.format(*a)) 
        for a in alleles}
   for ma, mi, cnt in alleles:
      if cnt == 0:
         continue
      log.trace('LP constraint: A_{} <= A_{} for {}-{}', cnt, cnt - 1, ma, mi)
      model.addConstr(A[(ma, mi, cnt)] <= A[(ma, mi, cnt - 1)])

   # Make sure that sum of all subaleles is exacly as the count of their major alleles
   for a, cnt in major_sol.solution.items():
      expr = model.quicksum(v for (ma, _, _), v in A.items() if ma == a)
      log.trace('LP constraint: {} == {} for {}', cnt, expr, a)
      model.addConstr(expr == cnt)

   # Add a binary variable for each allele/mutation pair where mutation belongs to that allele
   # that will indicate whether such mutation will be assigned to that allele or will be missing
   MMISS = {a: {m: model.addVar(vtype='B', name='MISS-{}-{}-{}_{}'.format(m, *a)) 
                for m in alleles[a]}
            for a in alleles}
   # Add a binary variable for each allele/mutation pair where mutation DOES NOT belongs to that allele
   # that will indicate whether such mutation will be assigned to that allele or not
   MADD = {a: {m: model.addVar(vtype='B', name='ADD-{}-{}-{}_{}'.format(m, *a)) 
               for m in mutations 
               if m not in alleles[a]}
           for a in alleles}
   # Add an error variable for each mutation and populate the error constraints
   error_vars = {m: model.addVar(lb=-model.INF, ub=model.INF, name=str(m)) 
                 for m in mutations}
   constraints = {m: 0 for m in mutations} 
   for m in mutations:
      m_gene, m_region = gene.region_at(m.pos)
      m_cn = major_sol.cn_solution.position_cn(m.pos)
      cov = coverage.total(m.pos) / m_cn if m_cn > 0 else 0
      
      for a in alleles:
         ma, mi, _ = a
         if m in alleles[a]:
            constraints[m] += cov * MMISS[a][m] * A[a]
         else:
            if gene.cn_configs[gene.alleles[ma].cn_config].cn[m_gene][m_region] > 0:
            # Add this *only* if CN of this region in a given allele is > 0 (i.e. do not add mutation to allele
            # where such mutation is deleted due to the fusion)
               constraints[m] += cov * MADD[a][m] * A[a]

      # Fill the constraints for non-variations (i.e. where nucleotide matches reference genome)
      ref_m = Mutation(pos=m.pos, op='REF')
      if ref_m not in constraints:
         error_vars[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=str(ref_m))
         constraints[ref_m] = 0

      for a in alleles:
         if gene.cn_configs[gene.alleles[ma].cn_config].cn[m_gene][m_region] == 0:
            continue
         if m in alleles[a]:
            constraints[ref_m] += cov * (1 - MMISS[a][m]) * A[a]
         elif not any(m.pos == fm.pos for fm in alleles[a]):
            # Make sure not to change functionality of an allele (only major step can do that!)
            constraints[ref_m] += cov * (1 - MADD[a][m]) * A[a]

   # Ensure that each constraint matches the observed coverage
   for m, expr in constraints.items():
      log.trace('LP constraint: {} == {} + err for {}', coverage[m], expr, m)
      model.addConstr(expr + error_vars[m] == coverage[m])
   # Ensure that a mutation is not assigned to allele that does not exist 
   for a, mv in MMISS.items():
      for m, v in mv.items():
         log.trace('LP contraint: {} >= {} for MMISS[{}, {}]', A[a], v, a, m)
         model.addConstr(v <= A[a])
   for a, mv in MADD.items():
      for m, v in mv.items():
         log.trace('LP contraint: {} <= {} for MADD[{}, {}]', A[a], v, a, m)
         model.addConstr(v <= A[a])

   # Ensure the following rules for all mutations:
   # 1) a minor allele must express ALL its functional mutations
   for a in alleles:
      p = [MMISS[a][m] for m in alleles[a] if m.is_functional]
      if len(p) == 0: 
         continue
      log.trace('LP constraint: {} = {} * A_{} for {}', sum(p), len(p), A[a], a)
      model.addConstr(sum(p) == len(p) * A[a]) # Either all or none
   # 2) No allele can include a novel mutation with coverage 0
   zero_muts = (m for a in alleles for m in alleles[a] if coverage[m] == 0)
   for m in zero_muts:
      expr = model.quicksum(MMISS[a][m] for a in alleles if m in alleles[a])
      log.trace('LP constraint: 0 >= {} for {}', expr, m)
      model.addConstr(expr <= 0)
   # 3) Make sure that CN of each variation does not exceed total supporting allele CN
   for m in mutations:
      m_cn = major_sol.cn_solution.position_cn(m.pos)
      expressed_cn = coverage.percentage(m) / 100.0

      # Get lower/upper CN bound: [floor(expressed_cn), ceil(expressed_cn)]
      if m_cn == 0:
         expressed_cn = (0, 0)
      elif coverage[m] > 0 and int(expressed_cn * m_cn) == 0:
         expressed_cn = (1, 1) # Force minimal CN to be 1
      else:
         expressed_cn = (int(expressed_cn * m_cn), int(expressed_cn * m_cn) + 1)
      
      expr = model.quicksum(MMISS[a][m] * A[a] for a in alleles if m in alleles[a])
      if not m.is_functional:
         # If this is functional mutation, alleles CANNOT include it via MADD
         # TODO: ensure that mutations agree with partially fused alleles 
         expr += model.quicksum(MADD[a][m] * A[a] for a in alleles if m not in alleles[a])

      lo, hi = expressed_cn
      log.trace('LP constraint: {} <= {} <= {} for {}', m, lo, expr, hi, m)
      model.addConstr(expr >= lo); model.addConstr(expr <= hi)

   # Objective: absolute sum of errors
   objective = model.abssum(error_vars.values())
   if solver == 'scip':
      # HACK: Non-linear objective linearization for SCIP:
      #       min f(x) <==> min w s.t. f(x) <= w
      w = model.addVar(name='W')
      nonlinear_obj = 0
      for a in alleles:
         nonlinear_obj += MISS_PENALTY_FACTOR * model.quicksum(
            A[a] * (1 - v) for m, v in MMISS[a].items())
         objective += ADD_PENALTY_FACTOR * model.quicksum(
            v for m, v in MADD[a].items())
      model.addConstr(nonlinear_obj <= w)
      objective += w
   else:
      objective += MISS_PENALTY_FACTOR * model.quicksum(
         A[a] * (1 - v) for a in alleles for _, v in MMISS[a].items())
      objective += ADD_PENALTY_FACTOR * model.quicksum(
         v for a in alleles for _, v in MADD[a].items())
   log.trace('Objective: {}', objective)

   # Solve the model
   try:
      status, opt = model.solve(objective)
      log.debug('CN Solver status: {}, opt: {}', status, opt)
   except lpinterface.NoSolutionsError:
      return MinorSolution(score=float('inf'), solution=[], major_sol=major_sol)


   # Get final minor solutions
   solution = []
   added = collections.defaultdict(list)
   missing = collections.defaultdict(list)
   for allele, value in A.items():
      if model.getValue(value) <= 0:
         continue

      solution.append(allele)
      for m, mv in MMISS[allele].items():
         if model.getValue(mv) == 0:
            missing[a].append(m)
         # TODO printer
         #if 'novel' in m.aux:
         #   mutations.add(('NOVEL', m, sam[m]))
         #else:
         #   mutations.add(('NORMAL', m, sam[m]))
      for m, mv in MADD[allele].items():
         if model.getValue(mv) > 0:
            added[a].append(m)

   return MinorSolution(score=opt, solution=solution, major_solution=major_sol,
                        added=dict(added), missing=dict(missing))

