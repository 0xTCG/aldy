# 786

# Aldy source: minor.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict, Tuple, Set, Callable, Optional
from functools import reduce, partial

import math
import collections
import multiprocessing

from . import lpinterface
from .common import *
from .gene import Mutation, Gene, Allele, Suballele
from .sam import Sample
from .cn import MAX_CN
from .coverage import Coverage
from .major import MajorSolution, SolvedAllele


class MinorSolution(collections.namedtuple('MinorSolution', ['score', 'solution', 'major_solution'])):
   """
   Describes a potential (possibly optimal) minor star-allele configuration.
   Immutable class.

   Attributes:
      score (float):
         ILP model error score (0 for user-provided solutions).
      solution (list[:obj:`SolvedAllele`]):
         List of minor star-alleles in the solution.
         Modifications to the minor alleles are represented in :obj:`SolvedAllele` format.
      major_solution (:obj:`aldy.major.MajorSolution`):
         Major star-allele solution used for calculating the minor star-allele assignment.
      diplotype (str):
         Assigned diplotype string (e.g. `*1/*2`).

   Notes:
      Has custom printer (``__str__``).
   """
   diplotype = ''
   
   def __str__(self):
      return f'MinorSol[{self.score:.2f}; ' + \
              'sol=({}); '.format(', '.join(str(s) for s in self.solution)) + \
             f'major={self.major_solution}'


def estimate_minor(gene: Gene, 
                   coverage: Coverage, 
                   major_sols: List[MajorSolution], 
                   solver: str,
                   filter_fn: Optional[Callable] = None) -> List[MinorSolution]:
   """
   :obj:`MinorSolution`: Detect the minor star-allele in the sample.

   Args:
      gene (:obj:`aldy.gene.Gene`): 
         A gene instance.
      coverage (:obj:`aldy.coverage.Coverage`): 
         Read alignment data.
      major_sol (:obj:`aldy.major.MajorSolution`): 
         Copy-number solution to be used for major star-allele calling.
      solver (str): 
         ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.
      filter_fn (callable):
         Custom filtering function (used for testing only).
   """

   # Get list of alleles and mutations to consider   
   alleles: List[Tuple[str, str]] = list()
   mutations: Set[Mutation] = set()
   for major_sol in major_sols:
      for (ma, _, added, _) in major_sol.solution:
         alleles += [SolvedAllele(ma, mi, added, tuple()) for mi in gene.alleles[ma].minors]
         mutations |= set(gene.alleles[ma].func_muts)
         mutations |= set(added)
         for sa in gene.alleles[ma].minors.values():
            mutations |= set(sa.neutral_muts)

   # Filter out low quality mutations
   def default_filter_fn(mut, cov, total, thres):
      # TODO: is this necessary?
      if mut.op != '_' and mut not in mutations: 
         return False
      return Coverage.basic_filter(mut, cov, total, thres / MAX_CN) and \
             Coverage.cn_filter(mut, cov, total, thres, major_sol.cn_solution) 
   if filter_fn:
      cov = coverage.filtered(filter_fn)
   else:
      cov = coverage.filtered(default_filter_fn)
   
   minor_sols = []
   for major_sol in major_sols:
      minor_sols.append(solve_minor_model(gene, alleles, cov, major_sol, mutations, solver))
   log.debug(f'>> minor_sols = {minor_sols.__repr__()}')
   return minor_sols


def solve_minor_model(gene: Gene,
                      alleles_list: List[SolvedAllele], 
                      coverage: Coverage, 
                      major_sol: MajorSolution, 
                      mutations: Set[Mutation], 
                      solver: str) -> MinorSolution:
   """
   :obj:`MinorSolution`: Solves the minor star-allele detection problem via integer linear programming.

   Args:
      gene (:obj:`aldy.gene.Gene`):
         Gene instance.
      alleles_list (list[:obj:`aldy.major.SolvedAllele`]):
         List of candidate minor star-alleles. 
      coverage (:obj:`aldy.coverage.Coverage`):
         Sample coverage used to find out the coverage of each mutation.
      major_sol (:obj:`aldy.major.MajorSolution`):
         Major star-allele solution to be used for detecting minor star-alleles (check :obj:`aldy.major.MajorSolution`).
      mutations (set[:obj:`aldy.gene.Mutation`]):
         List of mutations to be considered during the solution build-up 
         (all other mutations are ignored).
      solver (str): 
         ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.

   Notes:
      Please see Aldy paper (section Methods/Genotype refining) for the model explanation.
   """

   MISS_PENALTY_FACTOR = 2
   ADD_PENALTY_FACTOR = 1

   log.debug('\nRefining {}', major_sol)
   model = lpinterface.model('aldy_refine', solver)
   
   # Establish minor allele binary variables
   alleles = {(a, 0): set(gene.alleles[a.major].func_muts) | \
                      set(gene.alleles[a.major].minors[a.minor].neutral_muts) | \
                      set(a.added) 
              for a in alleles_list}

   log.debug('Possible candidates:')
   for a in sorted(alleles, key=lambda x: allele_sort_key(x[0].minor)):
      (ma, mi, _, _), _ = a
      log.debug('  *{} (cn=*{})', mi, gene.alleles[ma].cn_config)
      for m in sorted(alleles[a], key=lambda m: m.pos):
         log.debug('    {:26} {:4} ({:.1f} copies) {}',
            str(m), 
            coverage[m], 
            coverage[m] / (coverage.total(m.pos) / major_sol.cn_solution.position_cn(m.pos)) if coverage.total(m.pos) > 0 else 0,
            m.aux.get('old', ''))
   
   for a, _ in list(alleles):
      max_cn = major_sol.solution[SolvedAllele(a.major, None, a.added, a.missing)]
      for cnt in range(1, max_cn):
         alleles[a, cnt] = alleles[a, 0]
      
   A = {a: model.addVar(vtype='B', name=str(a)) for a in alleles}
   for a, cnt in alleles:
      if cnt == 0:
         continue
      log.trace('LP constraint: A_{} <= A_{} for {}', cnt, cnt - 1, a)
      model.addConstr(A[a, cnt] <= A[a, cnt - 1])

   # Make sure that sum of all subaleles is exactly as the count of their major alleles
   for sa, cnt in major_sol.solution.items():
      expr = model.quicksum(v for ((ma, _, ad, mi), _), v in A.items() 
                            if SolvedAllele(ma, None, ad, mi) == sa)
      log.trace('LP constraint: {} == {} for {}', cnt, expr, a)
      model.addConstr(expr == cnt)

   # Add a binary variable for each allele/mutation pair where mutation belongs to that allele
   # that will indicate whether such mutation will be assigned to that allele or will be missing
   MPRESENT = {a: {m: model.addVar(vtype='B', name='MISS-{}-{}'.format(m, a)) 
                   for m in alleles[a]}
               for a in alleles}
   # Add a binary variable for each allele/mutation pair where mutation DOES NOT belongs to that allele
   # that will indicate whether such mutation will be assigned to that allele or not
   MADD = {a: {m: model.addVar(vtype='B', name='ADD-{}-{}'.format(m, a)) 
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
         if m in alleles[a]:
            constraints[m] += cov * MPRESENT[a][m] * A[a]
         else:
            # Add this *only* if CN of this region in a given allele is positive 
            # (i.e. do not add mutation to allele if a region of mutation is deleted due to the fusion)
            (ma, _, _, _), _ = a
            if gene.cn_configs[gene.alleles[ma].cn_config].cn[m_gene][m_region] > 0:
               constraints[m] += cov * MADD[a][m] * A[a]

      # Fill the constraints for non-variations (i.e. where nucleotide matches reference genome)
      ref_m = Mutation(m.pos, 'REF') # type: ignore
      if ref_m not in constraints:
         error_vars[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=str(ref_m))
         constraints[ref_m] = 0

      for a in alleles:
         (ma, _, _, _), _ = a
         if gene.cn_configs[gene.alleles[ma].cn_config].cn[m_gene][m_region] == 0:
            continue
         if m in alleles[a]:
            # TODO: If a mutation is not present in allele, we assign it to REF mutation. 
            #       This should be extended to other potential mutations at the same loci as well.
            constraints[ref_m] += cov * (1 - MPRESENT[a][m]) * A[a]
         elif not any(m.pos == fm.pos for fm in alleles[a]):
            # ^ Make sure that reference position is not introduced 
            #   if there is already some mutation at this locus
            constraints[ref_m] += cov * (1 - MADD[a][m]) * A[a]

   # Ensure that each constraint matches the observed coverage
   for m, expr in constraints.items():
      log.trace('LP constraint: {} == {} + err for {}', coverage[m], expr, m)
      model.addConstr(expr + error_vars[m] == coverage[m])
   # Ensure that a mutation is not assigned to allele that does not exist 
   for a, mv in MPRESENT.items():
      for m, v in mv.items():
         log.trace('LP contraint: {} >= {} for MPRESENT[{}, {}]', A[a], v, a, m)
         model.addConstr(v <= A[a])
   for a, mv in MADD.items():
      for m, v in mv.items():
         log.trace('LP contraint: {} <= {} for MADD[{}, {}]', A[a], v, a, m)
         model.addConstr(v <= A[a])

   # Ensure the following rules for all mutations:
   # 1) A minor allele must express ALL its functional mutations
   for a in alleles:
      p = [MPRESENT[a][m] for m in alleles[a] if m.is_functional]
      if len(p) == 0: 
         continue
      expr = model.quicksum(p)
      log.trace('LP constraint 1: {} = {} * A_{} for {}', expr, len(p), A[a], a)
      model.addConstr(expr == len(p) * A[a]) # Either all or none
   # 2) No allele can include mutation with coverage 0
   zero_muts = (m for a in alleles for m in alleles[a] if coverage[m] == 0)
   for m in zero_muts:
      expr  = model.quicksum(MPRESENT[a][m] for a in alleles if m in alleles[a])
      expr += model.quicksum(MADD[a][m] for a in alleles if m not in alleles[a])
      log.trace('LP constraint 2: 0 >= {} for {}', expr, m)
      model.addConstr(expr <= 0)
   # 3) No allele can include extra functional mutation (this is resolved at step 2)
   for m in mutations: 
      if m.is_functional:
         expr = model.quicksum(MADD[a][m] for a in alleles if m not in alleles[a])
         log.trace('LP constraint 3: 0 >= {} for {}', expr, m)
         model.addConstr(expr <= 0)
   # 4) CNs at each loci must be respected
   for m in mutations:
      m_gene, m_region = gene.region_at(m.pos) 
      exprs = (A[a] 
               for a in alleles 
               if gene.cn_configs[gene.alleles[a[0].major].cn_config].cn[m_gene][m_region] > 0)
      expr = model.quicksum(exprs)
      total_cn = major_sol.cn_solution.position_cn(m.pos)
      log.trace(f'LP constraint 4: {total_cn} == {expr} for {m}')
      # print(f'{m}:{m_gene}/{m_region} --> {total_cn} @ {list(exprs)}')
      # floor/ceil because of total_cn potentially having 0.5 as a summand
      #model.addConstr(expr >= math.floor(total_cn))
      #model.addConstr(expr <= math.ceil(total_cn))
   # 5) Make sure that CN of each variation does not exceed total supporting allele CN
   for m in mutations: 
      m_cn = major_sol.cn_solution.position_cn(m.pos)
      exp_cn = coverage.percentage(m) / 100.0

      # Get lower/upper CN bound: [floor(expressed_cn), ceil(expressed_cn)]
      if m_cn == 0:
         lo, hi = 0, 0
      elif coverage[m] > 0 and int(exp_cn * m_cn) == 0:
         lo, hi = 1, 1 # Force minimal CN to be 1
      else:
         lo, hi = int(exp_cn * m_cn), min(int(exp_cn * m_cn) + 1, m_cn)
            
      expr  = model.quicksum(MPRESENT[a][m] * A[a] for a in alleles if m in alleles[a])
      expr += model.quicksum(MADD[a][m] * A[a] for a in alleles if m not in alleles[a])

      assert(lo >= 0)
      assert(hi <= m_cn)
      log.trace('LP constraint 5: {} <= {} <= {} for {}', m, lo, expr, hi, m)
      #model.addConstr(expr >= lo); model.addConstr(expr <= hi)

   # Objective: absolute sum of errors
   objective = model.abssum(error_vars.values())
   if solver == 'scip':
      # HACK: Non-linear objective linearization for SCIP:
      #       min f(x) <==> min w s.t. f(x) <= w
      w = model.addVar(name='W')
      nonlinear_obj = 0
      for a in alleles:
         nonlinear_obj += MISS_PENALTY_FACTOR * A[a] * \
                          model.quicksum((1 - v) for m, v in MPRESENT[a].items())
         objective += ADD_PENALTY_FACTOR * \
                      model.quicksum(v for m, v in MADD[a].items())
      model.addConstr(nonlinear_obj <= w)
      objective += w
   else:
      objective += MISS_PENALTY_FACTOR * \
                   model.quicksum(A[a] * (1 - v) for a in alleles for _, v in MPRESENT[a].items())
      objective += ADD_PENALTY_FACTOR * \
                   model.quicksum(v for a in alleles for _, v in MADD[a].items())
   log.trace('Objective: {}', objective)

   # Solve the model
   try:
      status, opt = model.solve(objective)
      log.debug('CN Solver status: {}, opt: {}', status, opt)
   except lpinterface.NoSolutionsError:
      log.debug(f'Minor solution: None')
      return MinorSolution(score=float('inf'), solution=[], major_solution=major_sol, missing={}, added={})

   # Get final minor solutions
   solution = []
   for allele, value in A.items():
      if model.getValue(value) <= 0:
         continue
      added: List[Mutation] = []
      missing: List[Mutation] = []
      for m, mv in MPRESENT[allele].items():
         if not model.getValue(mv):
            missing.append(m)
      for m, mv in MADD[allele].items():
         if model.getValue(mv):
            added.append(m)

      solution.append(SolvedAllele(allele[0].major, 
                                   allele[0].minor, 
                                   allele[0].added + tuple(added), 
                                   tuple(missing)))

   sol = MinorSolution(score=opt, 
                       solution=solution, 
                       major_solution=major_sol)
   log.debug(f'Minor solution: {sol}')

   return sol

