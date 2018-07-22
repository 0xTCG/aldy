from __future__ import division
# 786

# Aldy source: refiner.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from builtins import map
from builtins import str
from builtins import range
from pprint import pprint

import collections
import functools
import multiprocessing

from . import lpinterface

from .common import *
from .gene import Mutation
from .filtering import cnv_filter


def get_refined_solution(gene, sam, initial_solutions, solver):
   mutations = set()
   for sol in initial_solutions: # tuple (solution, cn_config)
      for an in sol[0]:
         mutations |= gene.alleles[an].functional_mutations
         for sa in list(gene.alleles[an].suballeles.values()):
            mutations |= sa.neutral_mutations

   # TODO fix this!
   # Use only relevant mutations
   for p in sam.coverage:
      for o in sam.coverage[p]:
         if o != '_' and Mutation(p, o) not in mutations:
            sam.coverage[p][o] = 0

   # pool = multiprocessing.Pool(4)
   results = list(map(
      functools.partial(refine_ilp, gene=gene, sam=sam, mutations=mutations, solver=solver),
      initial_solutions
   ))
   if len(results) > 0:
      score = min(results)[0]
      solutions = [(r[1], r[2]) for r in results if abs(r[0] - score) < 1e-6]
   else:
      score, solutions = float('inf'), []
   return score, solutions


def refine_ilp(solution, gene, sam, mutations, solver):
   log.debug('== Refiner ==')

   solution, region_cn = solution
   log.debug('Refining {}', solution)

   c = lpinterface.model('aldy_refine', solver)
   cnv_filter(sam, gene, region_cn)

   # List all alleles!
   alleles = { # key is (major, minor, copy)
      (a, sa, 0): set(gene.alleles[a].functional_mutations) | set(gene.alleles[a].suballeles[sa].neutral_mutations)
      for a in solution
      for sa in gene.alleles[a].suballeles
   }

   count = collections.Counter(solution)
   for a, sa, _ in list(alleles.keys()):
      for cnt in range(1, count[a]):
         alleles[(a, sa, cnt)] = alleles[(a, sa, 0)]

   # add one binary variable to model for any allele copy
   A = {a: c.addVar(vtype='B', name='{}-{}_{}'.format(*a)) for a in alleles}

   # a[i] >= a[i + 1] for all i
   for a, sa, cnt in alleles:
      if cnt != 0:
         continue
      for i in range(1, count[a]):
         log.trace('M  Sentinel constraint for {}-{}: A_{} <= A_{}', a, sa, i, i - 1)
         c.addConstr(A[(a, sa, i)] <= A[(a, sa, i - 1)])

   for a, cnt in count.items():
      expr = c.quicksum(var for aa, var in A.items() if aa[0] == a)
      log.trace('M  CN constraint for {}: {} == {}', a, cnt, expr)
      c.addConstr(expr == cnt)

   # for any mutation, add error variable to expression
   constraints = {
      m: [0, c.addVar(lb=-c.INF, ub=c.INF, name=str(m))]
      for m in mutations
   }
   # add binary variable for any allele/mutation pair indicating whether it's missing or not
   MM = {
      a: {m: c.addVar(vtype='B', name='MISS-{}-{}-{}_{}'.format(m, *a)) for m in alleles[a]}
      for a in alleles
   }
   # indicating is it added or no
   MA = {
      a: {m: c.addVar(vtype='B', name='ADD-{}-{}-{}_{}'.format(m, *a)) for m in mutations if m not in alleles[a]}
      for a in alleles
   }
   # populate constraints
   for m in mutations:
      coverage = 0
      if region_cn[gene.region_at[m.pos]] != 0:
         coverage = sam.total(m.pos) / region_cn[gene.region_at[m.pos]]

      for a in alleles:
         # print m.pos, gene.region_at[m.pos], gene.cnv_configurations[gene.alleles[a[0]].cnv_configuration].keys()
         if m in alleles[a]:
            constraints[m][0] += coverage * MM[a][m] * A[a]
         elif gene.cnv_configurations[gene.alleles[a[0]].cnv_configuration][gene.region_at[m.pos]] > 0:
            constraints[m][0] += coverage * MA[a][m] * A[a]

      # populate constraints for non-variations
      ref_m = Mutation(pos=m.pos, op='REF')
      if ref_m not in constraints:
         constraints[ref_m] = [0, c.addVar(lb=-c.INF, ub=c.INF, name=str(ref_m))]

      for a in alleles:
         if gene.cnv_configurations[gene.alleles[a[0]].cnv_configuration][gene.region_at[m.pos]] == 0:
            continue

         if m in alleles[a]:
            constraints[ref_m][0] += coverage * (1 - MM[a][m]) * A[a]
         elif not any(m.pos == x.pos for x in alleles[a]):
            constraints[ref_m][0] += coverage * (1 - MA[a][m]) * A[a]

   # make sure that each constraint equals its coverage
   for m, (expr, err) in constraints.items():
      log.trace('M  Contraint for {} {}: {} == {} + err', m, gene.region_at[m.pos], sam[m], expr)
      c.addConstr(expr + err == sam[m])
   # make sure that variation is not assigned to allele if allele does not exist
   for a, mv in MM.items():
      for m, v in mv.items():
         log.trace('M  Contraint for MM[missing] {} {}: {} >= {}', a, m, A[a], v)
         c.addConstr(v <= A[a])
   for a, mv in MA.items():
      for m, v in mv.items():
         log.trace('M  Contraint for MA[extra] {} {}: {} >= {}', a, m, A[a], v)
         c.addConstr(v <= A[a])

   # set objective: minimize absolute sum of errors
   objective = c.abssum([cx[1] for cx in constraints.values()])
   # Non-linear objective linearization for SCIP:
   #            min f(x) <==> min w s.t. f(x) <= w
   w = c.addVar(name='W')
   nonlinear_objective = 0
   for a in alleles:
      nonlinear_objective += 2 * A[a] * sum(1 - v for m, v in MM[a].items())
      objective += sum(v for m, v in MA[a].items())
   c.addConstr(nonlinear_objective <= w)
   objective += w
   log.trace('M  Objective: {}', objective)

   # important rules for mutations:
   # - sub-allele must express ALL its functional mutations!
   for a in alleles:
      p = [MM[a][m] for m in alleles[a] if m.functional]
      if len(p) == 0:
         continue
      log.trace('M  Contraint for {}: {} = {} * A_{}', a, sum(p), len(p), A[a])
      c.addConstr(sum(p) == len(p) * A[a])
   # - if coverage is nil, no allele can select mutation
   for m in (m for a in alleles for m in alleles[a] if sam[m] == 0):
      expr = c.quicksum([MM[a][m] for a in alleles if m in alleles[a]])
      log.trace('M  Contraint for {}: 0 >= {}', m, expr)
      c.addConstr(expr <= 0)
   # - make sure that variation is not over-expressed
   for m in mutations:
      region = gene.region_at[m.pos]
      full_cn = region_cn[region]
      mut_cn = sam.percentage(m)
      if full_cn == 0:
         mut_cn = (0, 0) # upper and lower bound
      else:
         mut_cn /= 100.0 / full_cn
         mut_cn = (int(mut_cn), int(mut_cn + 1))
         if sam[m] > 0 and mut_cn[0] == 0:
            mut_cn = (1, 1)
      expr = c.quicksum([
         MM[a][m] * A[a]
         for a in alleles if m in alleles[a]
      ])

      # if it is functional, alleles CANNOT include it via EXTRA
      if not m.functional:
         expr += c.quicksum([
            MA[a][m] * A[a] for a in alleles if m not in alleles[a]
         ])
         log.trace('M  Contraint for {}: {} <= {} <= {}', m, mut_cn[0], expr, mut_cn[1])
         c.addConstr(expr >= mut_cn[0])
         c.addConstr(expr <= mut_cn[1])
      else:
         # TODO think about this constraint harder
         # expr += sum([A[a]
         #  for a in alleles
         #  if a.split('_')[0].split('$')[0] in gene.fusions_left
         #  and gene.fusions[a.split('_')[0].split('$')[0]][gene.region_at[m.pos]] > 0])
         # log.debug('HC {} {} {} == {}', region, m, mut_cn, expr)
         log.trace('M  Contraint for {}: {} <= {} <= {}', m, mut_cn[0], expr, mut_cn[1])
         c.addConstr(expr >= mut_cn[0])
         log.trace('M  Contraint for {}: {} <= {}', m, expr, mut_cn[1])
         c.addConstr(expr <= mut_cn[1])

   # solve ILP
   try:
      status, opt = c.solve(objective)
      log.debug('CN Solver status: {}, opt: {}', status, opt)
   except lpinterface.NoSolutionsError:
      return float('inf'), []

   for allele, value in A.items():
      value = int(round(c.getValue(value)))
      if value <= 0:
         continue
      log.debug('  {}: {}', allele, value)
      for m, mv in MM[allele].items():
         mv = c.getValue(mv)
         if round(mv) <= 0:
            continue
         log.debug(
            '    {} {} {:4} ({:3.0f} cp. {:.0f}) {}',
            gene.region_at[m.pos], m,
            sam[m],
            sam.percentage(m),
            sam.percentage(m) / (100.0 / region_cn[gene.region_at[m.pos]]) if region_cn[gene.region_at[m.pos]] != 0 else 0,
            'F' if m.functional else ''
         )
   log.debug('Missing:')
   for allele, value in A.items():
      value = int(round(c.getValue(value)))
      if value <= 0:
         continue
      log.debug('  {}: {}', allele, value),
      for m, mv in MM[allele].items():
         if round(c.getValue(mv)) > 0:
            continue
         log.debug(
            '    {} {} {:4} ({:3.0f}) {}  {:12} {}',
            gene.region_at[m.pos], m,
            sam[m],
            sam.percentage(m),
            'F' if m.functional else '', m.aux['dbsnp'], m.aux['old']
         )
   log.debug('Extra:')
   for allele, value in A.items():
      value = int(round(c.getValue(value)))
      if value <= 0:
         continue
      log.debug('  {}: {}', allele, value),
      for m, mv in MA[allele].items():
         if round(c.getValue(mv)) <= 0:
            continue
         log.debug(
            '    {} {} {:4} ({:3.0f}) {}  {:12} {}',
            gene.region_at[m.pos], m,
            sam[m],
            sam.percentage(m),
            'F' if m.functional else '', m.aux['dbsnp'], m.aux['old']
         )

   decomposition = []
   for allele, value in A.items():
      if int(round(c.getValue(value))) <= 0:
         continue
      mutations = set()
      for m, mv in MM[allele].items():
         mv = c.getValue(mv)
         if int(round(mv)) > 0:
            if 'novel' in m.aux:
               mutations.add(('NOVEL', m, sam[m]))
            else:
               mutations.add(('NORMAL', m, sam[m]))
         else:
            mutations.add(('MISSING', m, sam[m]))
      for m, mv in MA[allele].items():
         if int(round(c.getValue(mv))) > 0:
            mutations.add(('EXTRA', m, sam[m]))
      decomposition.append((allele, mutations))

   solution = sorted([(a[0], a[1]) for a, v in A.items() if int(round(c.getValue(v))) > 0])
   return opt, tuple(solution), decomposition
