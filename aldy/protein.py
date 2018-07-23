from __future__ import division
# 786

# Aldy source: protein.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from builtins import map
from builtins import chr
from builtins import str
from builtins import range
from functools import reduce

import collections
import itertools
import functools
import multiprocessing

from . import lpinterface

from .common import *
from .cn import CNSolution
from .filtering import cnv_filter
from .gene import Allele, Mutation


#class MajorSolution(CNSolution):
#   def __init__(self, cn):



def get_initial_solution(gene, sam, cn_sols, solver):
   # pool = multiprocessing.Pool(1)
   results = list(map(
      functools.partial(solve_ilp, gene=gene, sam=sam, solver=solver),
      cn_sols.values()
   ))
   if len(results) > 0:
      score, solutions = min(results)[0], results[0][1]
      if len(results) > 1:
         solutions = reduce(lambda x, y: x + y, (r[1] for r in results if abs(r[0] - score) < 1e-6))
   else:
      score, solutions = float('inf'), []

   return score, solutions


def solve_ilp(cn_sol, gene, sam, solver):
   log.debug('== Initial Solver ==')

   structure, region_cn = cn_sol.cn_solution, cn_sol.region_cn
   structure = collections.Counter([y[0] for y in structure])
   log.debug('Solving {}', structure)
   if len(structure) == 0: # only deletions
      return 0, [([gene.deletion_allele] * 2, region_cn)]

   # More filtering
   alleles = []
   cnv_filter(sam, gene, region_cn)
   for an, a in sorted(list(gene.alleles.items()), key=lambda x: sort_key(x[0])):
      remove = any(sam[m] <= 0 for m in a.functional_mutations)
      if remove:
         s = ['Extra Removing {:4}  '.format(an)]
         for m in a.functional_mutations:
            if sam[m] <= 0:
               s.append('{} {}'.format(gene.region_at[m.pos], m))
         log.trace(s[0] + ',  '.join(s[1:]))
      elif a.cnv_configuration in structure:
         alleles.append(an)

   # Check for novel functional mutations
   novel_functional_mutations = set()
   for pos, c in sam.coverage.items():
      for op, cov in c.items():
         if op == '_':
            continue
         if (pos, op) in gene.mutations:
            continue
         # TODO: handle these regions as well
         if gene.region_at[pos][2:] not in gene.unique_regions:
            continue
         if region_cn[gene.region_at[pos]] == 0:
            continue
         # Require AT LEAST 80% coverage per copy
         if sam.percentage(Mutation(pos, op)) < 80.0 / region_cn[gene.region_at[pos]]:
            continue
         if check_functional(gene, Mutation(pos, op)):
            log.debug('Novel mutation: {} {} {} ({} or {}%)', gene.region_at[pos][2:], pos, op, cov, sam.percentage(Mutation(pos, op)))
            novel_functional_mutations.add(Mutation(pos, op, True))
   # log.debug('Novel mutations: {}', list(novel_functional_mutations))

   if set(structure.keys()) - set(gene.alleles[a].cnv_configuration for a in alleles):
      return float('inf'), [] # Some structure has no matching allele; no solutions
   
   rejected = []
   if sam.PHASE:
      for a in alleles:
         muts = sorted((m.pos, m.op) for m in gene.alleles[a].functional_mutations)
         for mi, m in enumerate(muts):
            for mj in range(mi + 1, len(muts)):
               if abs(muts[mi][0] - muts[mj][0]) > 50:
                  continue
               if sam.links[(muts[mi], muts[mj])] <= 1:
                  # Make sure that we have backup mutations to fall back upon
                  # TODO: find better way
                  other_mathcing = [s for s in sam.links if s[0] == muts[mi]]
                  if len(other_mathcing) > 1:
                     rejected.append(a)
                     break
            if len(rejected) > 0 and rejected[-1] == a: 
               break
   log.debug('Rejected alleles: {}', rejected)
   alleles = {(a, 0): gene.alleles[a] for a in alleles if a not in rejected}

   c = lpinterface.model('aldy_major_allele', solver)

   log.debug('Possible candidates:')
   for a, _ in sorted(list(alleles.keys()), key=lambda x: sort_key(x[0])):
      log.debug('  {}*{}', gene.name, a)
      for m in sorted(alleles[(a, 0)].functional_mutations, key=lambda m: m.pos):
         log.debug(
            '    {} {} {:4} ({:3.0f}) {} {}',
            gene.region_at[m.pos],
            m, sam[m], sam.percentage(m),
            'F', m.aux['old']
         )
   # create special allele for any possible copy
   for a, _ in list(alleles.keys()):
      max_cn = sum(structure.values())
      log.trace('M  Max. CN for {} = {}', a, max_cn)
      for i in range(1, max_cn):
         alleles[(a, i)] = alleles[(a, 0)]

   # add one binary variable to model for any allele copy
   A = {a: c.addVar(vtype='B', name='A_{}_{}'.format(*a)) for a in alleles}
   # for any mutation, add error variable to expression
   constraints = {
      # m -> constraint, error_var
      m: [0, c.addVar(lb=-c.INF, ub=c.INF, name='MA_{}_{}_{}'.format(m, *a))]
      for a in alleles
      for m in alleles[a].functional_mutations
   }
   # add binary variable for any allele/mutation pair
   M = {
      a: {
         m: c.addVar(vtype='B', name='EXTRA_{}_{}_{}'.format(m, *a))
         for m in constraints
         if m not in alleles[a].functional_mutations
      } for a in alleles
   }
   # populate constraints
   for a in alleles:
      for m in alleles[a].functional_mutations:
         if region_cn[gene.region_at[m.pos]] == 0:
            coverage = 0
         else:
            coverage = max(1, sam.total(m.pos)) / region_cn[gene.region_at[m.pos]]
         constraints[m][0] += coverage * A[a]
   for a in M:
      for m in M[a]:
         if region_cn[gene.region_at[m.pos]] != 0:
            constraints[m][0] += coverage * A[a] * M[a][m]
   # populate constraints for non-variations
   for m in list(constraints.keys()):
      if m.op[:3] == 'INS':
         continue # no need for insertions...
      ref_m = Mutation(pos=m.pos, op='REF')
      if ref_m in constraints:
         ref_m = Mutation(pos=m.pos, op=ref_m.op + '#')
      if ref_m not in constraints:
         constraints[ref_m] = [0, c.addVar(lb=-c.INF, ub=c.INF, name=str(ref_m))]

      if region_cn[gene.region_at[m.pos]] == 0:
         coverage = 0
      else:
         coverage = max(1, sam.total(m.pos)) / region_cn[gene.region_at[m.pos]]
      for a in alleles:
         constraints[ref_m][0] += coverage * A[a]

   # make sure that each constraint equals its coverage
   for m, (expr, err) in constraints.items():
      log.trace('M  Contraint for {} {}: {} == {} + err', m, gene.region_at[m.pos], sam[m], expr)
      c.addConstr(expr + err == sam[m])

   # force alleles to be as much as needed
   for a, cnt in structure.items():
      expr = sum(A[y] for y in A if alleles[y].cnv_configuration == a)
      log.trace('M  CN contraint for {}: {} == {}', a, cnt, expr)
      c.addConstr(expr == cnt)

   # Make sure that A[i+1] <= A[i] (to avoid equivalent solutions)
   for a, ai in alleles:
      if ai != 0:
         continue
      for i in itertools.count(1):
         if (a, i) not in A:
            break
         log.trace('M  Sentinel contraint for {}: A_{} <= A_{}', a, i, i - 1)
         c.addConstr(A[(a, i)] <= A[(a, i - 1)])

   # allele must express all of its functional functional_mutations
   for m in (m for a in alleles for m in alleles[a].functional_mutations if sam[m] > 0):
      expr = c.quicksum([A[a] for a in alleles if m in alleles[a].functional_mutations])
      expr += c.quicksum([A[a] * M[a][m] for a in alleles if m not in alleles[a].functional_mutations])
      log.trace('M  Contraint for {}: 1 <= {}', m, expr)
      c.addConstr(expr >= 1)

   # set objective: minimize absolute sum of errors
   # also make sure that every non-expressed snp gets heavily penalized
   penal = 100000
   objective = c.abssum([cx[1] for cx in constraints.values()]) + penal * c.quicksum([M[a][m] for a in M for m in M[a]])
   log.trace('M  Objective: {}', objective)

   # solve ILP
   try:
      status, opt, solutions = c.solveAll(objective, dict(list(A.items()) + [((a, m), M[a][m]) for a in M for m in M[a]]))
      log.debug('CN Solver status: {}, opt: {}', status, opt)
   except lpinterface.NoSolutionsError:
      return float('inf'), []

   unique_key = collections.defaultdict(str)
   for si, s in enumerate(solutions):
      sd = {x: [y[1] for y in s if isinstance(y[0], tuple) and y[0][0] == x[0]] for x in s if not isinstance(x[0], tuple)}
      solutions[si] = []
      for (a, _), extra in sd.items():
         if len(extra) > 0:
            log.warn('Novel major star-allele (*{}-like) found!', a)
            extra = [Mutation(m.pos, m.op, m.functional, dict(list(m.aux.items()) + [('novel', True)])) for m in extra]
            an = Allele(
               gene.alleles[a].name + '+{}'.format(unique_key[gene.alleles[a].name]),
               gene.alleles[a].cnv_configuration,
               gene.alleles[a].functional_mutations | set(extra),
               gene.alleles[a].suballeles
            )
            unique_key[gene.alleles[a].name] = 'a' if unique_key[gene.alleles[a].name] == '' else chr(ord(unique_key[gene.alleles[a].name]) + 1)
            gene.alleles[an.name] = an
            solutions[si].append(an.name)
         else:
            solutions[si].append(a)
   log.debug('  solution:')
   for s in solutions:
      log.debug('    {}', s)

   return opt, [(s, region_cn) for s in solutions]
