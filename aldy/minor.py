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

# Model parameters
MISS_PENALTY_FACTOR = 1.5
"""float: Penalty for each missed minor mutation (0 for no penalty).
          Ideally larger than `ADD_PENALTY_FACTOR` as additions should be cheaper."""

ADD_PENALTY_FACTOR = 1.0
"""float: Penalty for each novel minor mutation (0 for no penalty). 
          Zero penalty always prefers additions over coverage errors 
          if the normalized SNP slack coverage is >= 50%.
          Penalty of 1.0 prefers additions if the SNP slack coverage is >= 75%."""


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
         Assigned diplotype string (e.g. ``*1/*2``).

   Notes:
      Has custom printer (``__str__``).
   """
   diplotype = ''


   def _solution_nice(self):
      return ', '.join(str(s) 
                       for s in sorted(self.solution, 
                                       key=lambda x: allele_sort_key(x.minor)))

   
   def __str__(self):
      return f'MinorSol[{self.score:.2f}; ' + \
             f'sol=({self._solution_nice()}); ' + \
             f'major={self.major_solution}'


def estimate_minor(gene: Gene, 
                   coverage: Coverage, 
                   major_sols: List[MajorSolution], 
                   solver: str,
                   filter_fn: Optional[Callable] = None) -> List[MinorSolution]:
   """
   Detect the minor star-allele in the sample.

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

   Returns:
      list[:obj:`MinorSolution`]
   """

   # Get list of alleles and mutations to consider
   alleles: List[Tuple[str, str]] = list()
   # Consider all major and minor mutations *from all available major solutions* together
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
   for major_sol in sorted(major_sols, key=lambda s: list(s.solution.items())):
      minor_sols += solve_minor_model(gene, alleles, cov, major_sol, mutations, solver)
   return minor_sols


def solve_minor_model(gene: Gene,
                      alleles_list: List[SolvedAllele], 
                      coverage: Coverage, 
                      major_sol: MajorSolution, 
                      mutations: Set[Mutation], 
                      solver: str) -> List[MinorSolution]:
   """
   Solves the minor star-allele detection problem via integer linear programming.

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

   Returns:
      list[:obj:`MinorSolution`]
      
   Notes:
      Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_ (section Methods/Genotype refining) for the model explanation.
      Currently returns only the first optimal solution.
   """

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
         m_gene, m_region = gene.region_at(m.pos)
         log.debug('    {:26}  {:.2f} ({:4} / {} * {:4.0f}) {}:{:10} {:>20}',
            str(m), 
            coverage[m] / coverage.single_copy(m.pos, major_sol.cn_solution),
            coverage[m], 
            major_sol.cn_solution.position_cn(m.pos),
            coverage.single_copy(m.pos, major_sol.cn_solution),
            m_gene, str(m_region),
            m.aux.get('old', ''))
   
   for a, _ in list(alleles):
      max_cn = major_sol.solution[SolvedAllele(a.major, None, a.added, a.missing)]
      for cnt in range(1, max_cn):
         alleles[a, cnt] = alleles[a, 0]
      
   A = {a: model.addVar(vtype='B', name=a[0].minor) for a in alleles}
   for a, cnt in alleles:
      if cnt == 0:
         continue
      model.addConstr(A[a, cnt] <= A[a, cnt - 1])

   # Make sure that sum of all subaleles is exactly as the count of their major alleles
   for sa, cnt in major_sol.solution.items():
      expr = model.quicksum(v for ((ma, _, ad, mi), _), v in A.items() 
                            if SolvedAllele(ma, None, ad, mi) == sa)
      model.addConstr(expr == cnt)

   # Add a binary variable for each allele/mutation pair where mutation belongs to that allele
   # that will indicate whether such mutation will be assigned to that allele or will be missing
   MPRESENT = {a: {m: model.addVar(vtype='B', name=f'P_{m.pos}_{m.op.replace(".", "")}_{a[0].minor}') 
                   for m in alleles[a]} 
               for a in alleles}
   # Add a binary variable for each allele/mutation pair where mutation DOES NOT belongs to that allele
   # that will indicate whether such mutation will be assigned to that allele or not
   MADD = {a: {} for a in alleles}
   for a in MADD:
      for m in mutations:
         if gene.has_coverage(a[0].major, m.pos) and m not in alleles[a]:
            MADD[a][m] = model.addVar(vtype='B', name=f'A_{m.pos}_{m.op.replace(".", "")}_{a[0].minor}')
   # Add an error variable for each mutation and populate the error constraints
   error_vars = {m: model.addVar(lb=-model.INF, ub=model.INF, name=f'E_{m.pos}_{m.op.replace(".", "")}') 
                 for m in mutations}
   constraints = {m: 0 for m in mutations} 
   for m in mutations:
      for a in alleles:
         if m in alleles[a]:
            constraints[m] += MPRESENT[a][m] * A[a]
         elif gene.has_coverage(a[0].major, m.pos):
            # Add this *only* if CN of this region in a given allele is positive 
            constraints[m] += MADD[a][m] * A[a]

   # Fill the constraints for non-variations (i.e. where nucleotide matches reference genome)
   for pos in set(m.pos for m in constraints):
      ref_m = Mutation(pos, '_') # type: ignore
      error_vars[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=f'E_{pos}_REF')
      constraints[ref_m] = 0
      for a in alleles:
         if not gene.has_coverage(a[0].major, pos):
            continue
         # Does this allele contain any mutation at the position `pos`? 
         # Insertions are not counted as they always contribute to `_`.
         present_muts = [m for m in alleles[a] if m.pos == pos and m[1][:3] != 'INS']
         assert(len(present_muts) < 2)
         if len(present_muts) == 1:
            constraints[ref_m] += (1 - MPRESENT[a][present_muts[0]]) * A[a]
         else:
            N = 1
            for m in MADD[a]:
               if m.pos == pos and m[1][:3] != 'INS':
                  N *= 1 - MADD[a][m]
            constraints[ref_m] += N * A[a]

   # Ensure that each constraint matches the observed coverage
   print('  {')
   print(f'    "cn": {str(dict(major_sol.cn_solution.solution))}, ')
   print( '    "major": ' + str({
         tuple([s.major] + [(m[0], m[1]) for m in s.added]) if len(s.added) > 0 else s.major: v
         for s, v in major_sol.solution.items()}) + ", ")
   print('    "data": {', end='')
   for m, expr in sorted(constraints.items()):
      cov = coverage[m] / coverage.single_copy(m.pos, major_sol.cn_solution) 
      model.addConstr(expr + error_vars[m] == cov)
      print(f"({m[0]}, '{m[1]}'): {cov}, ", end='')
   print('}, ')

   # Ensure the following rules for all mutations:
   # 1) Each mutation is assigned only to alleles that are present in the solution
   for a, mv in MPRESENT.items():
      for m, v in mv.items():
         model.addConstr(v <= A[a])
   for a, mv in MADD.items():
      for m, v in mv.items():
         model.addConstr(v <= A[a])
   # 2) Each allele must express ALL its functional mutations
   for a in alleles:
      for m in alleles[a]: 
         if m.is_functional:
            assert(m in MPRESENT[a])
            model.addConstr(MPRESENT[a][m] >= A[a])
   # 3) No allele can include mutation with coverage 0
   for a in MPRESENT:
      for m, v in MPRESENT[a].items():
         if not gene.has_coverage(a[0].major, m.pos):
            model.addConstr(v <= 0)
   # 4) No allele can include additional functional mutation (this should be done in the major model)
   for a in MADD:
      for m, v in MADD[a].items():
         if m.is_functional:
            model.addConstr(v <= 0)
   # 5) Prevent extra mutations if there is already existing mutation at that loci
   #    (either existing or already added)
   for pos in set(m.pos for m in constraints):
      for a in alleles:
         mp = [A[a] * MPRESENT[a][m] for m in MPRESENT[a] if m.pos == pos]
         ma = [A[a] * MADD[a][m] for m in MADD[a] if m.pos == pos]
         # TODO: add support for extra insertions!
         model.addConstr(model.quicksum(ma) <= 1)
         model.addConstr(model.quicksum(mp + ma) <= 1)
   # 5) Make sure that each copy of a mutation has at least one read supporting it
   for m in mutations:
      expr  = model.quicksum(MPRESENT[a][m] * A[a] for a in alleles if m in MPRESENT[a])
      expr += model.quicksum(MADD[a][m] * A[a] for a in alleles if m in MADD[a])
      if major_sol.cn_solution.position_cn(m.pos) == 0 or coverage[m] == 0:
         model.addConstr(expr <= 0)
      else:
         model.addConstr(expr <= coverage[m])
         # Ensure that at least one allele picks an existing non-filtered mutation
         model.addConstr(expr >= 1) 
   # 6) Do the same for non-mutations
   for pos in set(m.pos for m in constraints):
      expr = []
      for a in alleles:
         e  = [(1 - MPRESENT[a][m]) * A[a] for m in MPRESENT[a] if m.pos == pos]
         e += [(1 - MADD[a][m]) * A[a] for m in MADD[a] if m.pos == pos]
         assert(len(e) > 0)
         expr += e
      expr = model.quicksum(expr)
      m = Mutation(pos, '_')
      if major_sol.cn_solution.position_cn(m.pos) == 0:
         model.addConstr(expr <= 0)
      else:
         # If there is no coverage at pos at all despite CN being > 0 
         # (or if there is coverage only on functional mutations that cannot be selected), 
         # allow alleles to select non-existent non-mutation to prevent infeasible model.
         # Should not happen with "sane" datasets...
         model.addConstr(expr <= max(major_sol.cn_solution.position_cn(m.pos), coverage[m]))
         # Ensure that at least one allele picks an existing non-filtered non-mutation
         if coverage[m] > 0:
            model.addConstr(expr >= 1) 
  
   # Objective: absolute sum of errors
   objective = model.abssum(v for _, v in error_vars.items())
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
                   model.quicksum(A[a] * (1 - v) for a in sorted(alleles) for _, v in sorted(MPRESENT[a].items()))
      objective += ADD_PENALTY_FACTOR * \
                   model.quicksum(v for a in sorted(alleles) for _, v in sorted(MADD[a].items()))

   # Solve the model
   try:
      status, opt = model.solve(objective)
      # model.model.write('minor.lp')
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
      print('    "sol": ' + str([
          (s.minor, [(m[0], m[1]) for m in s.added], [(m[0], m[1]) for m in s.missing])
          for s in solution]))
      sol = MinorSolution(score=opt,
                           solution=solution,
                           major_solution=major_sol)
      log.debug(f'Minor solution: {sol}')
      return [sol]
   except lpinterface.NoSolutionsError:
      log.debug('No minor solutions')
      # Enable to debug infeasible models
      # model.model.computeIIS()
      # model.model.write("minor.ilp")
      return []

