# 786

# Aldy source: major.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Dict, Tuple, Set, Any, Optional

import collections
import itertools
import functools
import multiprocessing
import copy
import operator
from functools import reduce

from . import lpinterface
from .common import *
from .cn import CNSolution, MAX_CN
from .gene import Allele, Mutation, Gene, CNConfig
from .coverage import Coverage


# Model parameters
NOVEL_MUTATION_PENAL = MAX_CN + 1
"""float: Penalty for each novel mutation (0 for no penalty)"""


class SolvedAllele(collections.namedtuple('SolvedAllele', ['major', 'minor', 'added', 'missing'])):
   """
   Describes a candidate star-allele configuration.
   Immutable class.

   Attributes:
      major (str):
         Major star-allele identifier.
      minor (str, optional):
         Minor star-allele identifier. Can be None.
      added (tuple[:obj:`aldy.gene.Mutation`]):
         Tuple of mutations that are added to this copy of a major/minor star-allele
         (e.g. these mutations are not present in the database defition of allele).
      missing (tuple[:obj:`aldy.gene.Mutation`]):
         Tuple of mutations that are ommited from this copy of a major/minor star-allele
         (e.g. these mutations are present in the database defition of allele but not in the sample).

   Notes:
      Has custom printer (``__str__``).
   """


   def major_repr(self):
      return '*{}{}'.format(self.major,
         ''.join(' +' + str(m) for m in sorted(m for m in self.added if m.is_functional)))


   def __str__(self):
      return '*{}{}{}'.format(
         self.minor if self.minor else self.major,
         ''.join(' +' + str(m) for m in sorted(self.added,
                                               key=lambda m: (-m.is_functional, m.pos, m.op))),
         ''.join(' -' + str(m) for m in sorted(self.missing)))


class MajorSolution(collections.namedtuple('MajorSolution', ['score', 'solution', 'cn_solution'])):
   """
   Describes a potential (possibly optimal) major star-allele configuration.
   Immutable class.

   Attributes:
      score (float):
         ILP model error score (0 for user-provided solutions).
      solution (dict[:obj:`SolvedAllele`, int]):
         Dictionary of major star-alleles where each major star-allele is
         associated with its copy number
         (e.g. ``{1: 2}`` means that we have two copies of \*1).
      cn_solution (:obj:`aldy.cn.CNSolution`):
         Associated copy-number solution used for calculating the major
         star-alleles.

   Notes:
      Has custom printer (``__str__``).
   """


   def _solution_nice(self):
      return ', '.join(f'{v}x{s}'
                       for s, v in sorted(self.solution.items(),
                                          key=lambda x: allele_sort_key(x[0].major)))


   def __str__(self):
      return f'MajorSol[{self.score:.2f}; ' + \
             f'sol=({self._solution_nice()}); ' + \
             f'cn={self.cn_solution}'


def estimate_major(gene: Gene,
                   coverage: Coverage,
                   cn_solution: CNSolution,
                   solver: str,
                   gap: float = 0,
                   identifier: int = 0,
                   debug: Optional[str] = None) -> List[MajorSolution]:
   """
   Detect the major star-alleles in the sample.

   Args:
      gene (:obj:`aldy.gene.Gene`):
         A gene instance.
      coverage (:obj:`aldy.coverage.Coverage`):
         Read alignment coverage data.
      cn_solution (:obj:`aldy.cn.CNSolution`):
         Copy-number solution to be used for major star-allele calling.
      solver (str):
         ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.
      gap (float):
         Relative optimality gap (percentage).
         Default is 0 (report only optimal solutions).
      identifier (int):
         Unique solution identifier. Used for generating debug information.
         Default is 0.
      debug (str, optional):
         If set, Aldy will create "<debug>.major.lp" file for debug purposes.
         Default is ``None``.

   Returns:
      list[:obj:`MajorSolution`]
   """

   log.debug('Solving major alleles for cn={}', cn_solution)

   if sum(cn_solution.solution.values()) < 2:
      raise AldyException("estimate_major requires at least two valid gene configurations")

   alleles, coverage = _filter_alleles(gene, coverage, cn_solution)
   # Check if some CN solution has no matching allele
   if set(cn_solution.solution) - set(a.cn_config for a in alleles.values()):
      results = []
   else:
      results = solve_major_model(gene, alleles, coverage, cn_solution, solver, gap, identifier, debug)
   # TODO: re-implement phasing step from Aldy 1.4
   # TODO: Check for novel functional mutations and do something with them
   # novel_functional_mutations = _get_novel_mutations(gene, coverage, cn_solution)

   return results


def solve_major_model(gene: Gene,
                      allele_dict: Dict[str, Allele],
                      coverage: Coverage,
                      cn_solution: CNSolution,
                      solver: str,
                      gap: float = 0,
                      identifier: int = 0,
                      debug: Optional[str] = None) -> List[MajorSolution]:
   """
   Solves the major star-allele detection problem via integer linear programming.

   Args:
      gene (:obj:`aldy.gene.Gene`):
         A gene instance.
      allele_dict (dict[str, :obj:`aldy.gene.Allele`]):
         Dictionary of candidate major star-alleles.
      coverage (:obj:`aldy.coverage.Coverage`):
         Sample coverage used to find out the coverage of each major mutation
      cn_solution (:obj:`aldy.cn.CNSolution`):
         Copy-number solution to be used for detecting major star-alleles (check :obj:`aldy.cn.CNSolution`).
      solver (str):
         ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.
      gap (float):
         Relative optimality gap (percentage).
         Default is 0 (report only optimal solutions).
      identifier (int):
         Unique solution identifier. Used for generating debug information.
         Default is 0.
      debug (str, optional):
         If set, Aldy will create "<debug>.major.lp" file for debug purposes.
         Default is ``None``.

   Returns:
      list[:obj:`MajorSolution`]

   Notes:
      Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_ (section Methods/Major star-allele identification) for the model explanation.
   """

   # Make sure that coverage defaults to 0 on empty values

   model = lpinterface.model('aldy_major_allele', solver)

   func_muts = {M for m, M in gene.mutations.items()
                if M.is_functional and coverage[M] > 0}
   _print_candidates(allele_dict, coverage, cn_solution, func_muts)

   # hack to silence type checker
   a: Any = 0

   # Create a binary variable for all possible allele copies
   alleles = {(a, 0): allele_dict[a] for a in allele_dict}
   for (an, _), a in list(alleles.items()):
      max_cn = cn_solution.solution[a.cn_config]
      for i in range(1, max_cn):
         alleles[an, i] = alleles[an, 0]
   VA = {a: model.addVar(vtype='B', name=var_name(f'A_{a[0]}_{a[1]}')) for a in alleles}

   # Make sure that VA[i+1] <= VA[i] (to avoid equivalent solutions)
   for a, ai in alleles.keys():
      if ai > 0:
         model.addConstr(VA[a, ai] <= VA[a, ai - 1])

   # Add an error variable to the ILP for any mutation
   VERR = {m: model.addVar(lb=-model.INF, ub=model.INF, name=var_name(f'E_{m.pos}_{m.op.replace(".", "")}'))
              for m in func_muts}
   constraints = {e: 0 for e in VERR}
   # Add a binary variable for any allele/novel mutation pair
   del_allele = gene.deletion_allele()
   VNEW = {a: {m: model.addVar(vtype='B', name=var_name(f'N_{m.pos}_{m.op.replace(".", "")}_{a[0]}'))
               if a[0] != del_allele else 0 # deletion alleles should not be assigned any mutations
               for m in constraints 
               if m[0] not in set(mm[0] for mm in alleles[a].func_muts) and gene.has_coverage(a[0], m.pos)}
           for a in alleles}
   # Populate constraints
   for a in alleles:
      for m in alleles[a].func_muts:
         constraints[m] += VA[a]
   # Add novel mutation constraints
   for a in VNEW:
      for m in list(VNEW[a].keys()):
         constraints[m] += VA[a] * VNEW[a][m]

   # Populate constraints of non-variations (i.e. matches with the reference genome)
   for pos in set(m.pos for m in constraints):
      ref_m = Mutation(pos, '_') # type: ignore
      constraints[ref_m] = 0
      VERR[ref_m] = model.addVar(lb=-model.INF, ub=model.INF, name=var_name(f'E_{pos}_REF'))

      for a in alleles:
         if not gene.has_coverage(a[0], pos):
            continue
         # If allele has insertion at this loci, it will still contribute to the coverage at this loci
         if any(ma[0] == pos and not ma[1][:3] == 'INS' for ma in alleles[a].func_muts):
            continue
         # Make sure that this allele has no novel alleles at this position
         Ns = [VNEW[a][m] for m in VNEW[a] if m[0] == pos and m[1][:3] != 'INS']
         if len(Ns) > 0:
            Ns = 1 - reduce(operator.mul, Ns, 1)
         else:
            Ns = 1
         constraints[ref_m] += VA[a] * Ns

   # Each allele must express all of its functional mutations
   json_print(debug, '  { # {}', identifier)
   json_print(debug, f'    "cn": {str(dict(cn_solution.solution))}, ')
   json_print(debug,  '    "data": {', end='')
   prev = 0
   for m, expr in sorted(constraints.items()):
      cov = coverage[m] / coverage.single_copy(m.pos, cn_solution)
      model.addConstr(expr + VERR[m] == cov)
      if m.pos != prev and prev != 0:
         json_print(debug, '\n             ', end='')
      prev = m.pos
      json_print(debug, f"({m[0]}, '{m[1]}'): {cov:.4f}, ", end='')
   json_print(debug, '}, ')

   # Each CN config must be satisfied by matching alleles
   for cnf, cnt in cn_solution.solution.items():
      expr = sum(VA[a] for a in VA if alleles[a].cn_config == cnf)
      model.addConstr(expr == cnt)

   # Each functional mutation must be chosen by some allele and expressed
   for m in func_muts:
      expr  = model.quicksum(VA[a] for a in alleles if m in alleles[a].func_muts)
      expr += model.quicksum(VA[a] * VNEW[a][m] for a in alleles if m in VNEW[a] and m not in alleles[a].func_muts)
      model.addConstr(expr >= 1)

   # Set objective: minimize the absolute sum of errors
   objective = \
      model.abssum(e for e in VERR.values()) + \
      NOVEL_MUTATION_PENAL * model.quicksum(VNEW[a][m] for a in VNEW for m in VNEW[a])
   model.setObjective(objective)
   if debug:
      model.dump(f'{debug}.major{identifier}.lp')

   # Solve the ILP
   try:
      lookup = {**{model.varName(v): a for a, v in VA.items()},
                **{model.varName(v): (a, m) for a in VNEW for m, v in VNEW[a].items()}}
      result = {}
      json_print(debug, '    "sol": [', end='')
      for status, opt, sol in model.solutions(gap):
         print(f'Major solver status: {status}, opt: {opt}')
         solved_alleles = collections.defaultdict(lambda: []) # dict of allele IDs -> novel mutations
         for s in sol:
            if isinstance(lookup[s][0], tuple): # Novel allele
               solved_alleles[lookup[s][0]].append(lookup[s][1])
            elif s not in solved_alleles:
               solved_alleles[lookup[s]] = []
         sol_tuple = sorted_tuple((a, sorted_tuple(nm)) for (a, _), nm in solved_alleles.items())
         if sol_tuple not in result:
            solution = collections.Counter(SolvedAllele(major=a,
                                                         minor=None,
                                                         added=tuple(mut),
                                                         missing=tuple())
                                          for (a, _), mut in solved_alleles.items())
            json_print(debug, str(dict(collections.Counter(
               tuple([a] + [(m[0], m[1]) for m in mut]) if len(mut) > 0 else a
               for (a, _), mut in solved_alleles.items()))) + ', ', end='')
            sol = MajorSolution(score=opt,
                                 solution=solution,
                                 cn_solution=cn_solution)
            log.debug('Major solution: {}'.format(sol))
            result[sol_tuple] = sol
      json_print(debug, ']\n  }, ', end='')
      return result.values()
   except lpinterface.NoSolutionsError:
      return []


def _filter_alleles(gene: Gene,
                    coverage: Coverage,
                    cn_solution: CNSolution) -> Tuple[Dict[str, Allele], Coverage]:
   """
   Filters out all low-quality mutations and impossible alleles.

   Returns:
      tuple[dict[str, :obj:`aldy.gene.Allele`], :obj:`aldy.coverage.Coverage`]: Allele
      dictionary describing feasible alleles and high-confidence variants.
   """

   def filter_fns(mut, cov, total, thres):
      return Coverage.basic_filter(mut, cov, total, thres / MAX_CN) and \
             Coverage.cn_filter(mut, cov, total, thres, cn_solution)
   cov = coverage.filtered(filter_fns)
   alleles = copy.deepcopy(gene.alleles)
   for an, a in sorted(gene.alleles.items()):
      if a.cn_config not in cn_solution.solution:
         del alleles[an]
      elif any(cov[m] <= 0 for m in a.func_muts):
         s = ('{} in {}'.format(m, gene.region_at(m.pos))
              for m in a.func_muts if cov[m] <= 0)
         log.trace('Removing {} because of {}', an, ' and '.join(s))
         del alleles[an]

   return alleles, cov


def _print_candidates(alleles: Dict[str, Allele],
                      coverage: Coverage,
                      cn_solution: CNSolution,
                      func_muts: set) -> None:
   """
   Pretty-prints the list of allele candidates and their functional mutations.
   """
   log.debug('Possible candidates:')
   func_muts = func_muts.copy()
   for a in sorted(alleles, key=allele_sort_key):
      log.debug('  *{} (cn=*{})', a, alleles[a].cn_config)
      for m in sorted(alleles[a].func_muts, key=lambda m: m.pos):
         if m in func_muts: func_muts.remove(m)
         log.debug('    {} {:4} ({:.1f} copies) {} {}',
            #coverage.region_at(m.pos),
            m, coverage[m],
            coverage[m] / (coverage.total(m.pos) / cn_solution.position_cn(m.pos)),
            'F', m.aux.get('old', ''))
   if len(func_muts) > 0:
      log.debug('  Other mutations:')
      for m in sorted(func_muts, key=lambda m: m.pos):
         log.debug('    {} {:4} ({:.1f} copies) {} {}',
            m, coverage[m],
            coverage[m] / (coverage.total(m.pos) / cn_solution.position_cn(m.pos)),
            'F', m.aux.get('old', ''))
