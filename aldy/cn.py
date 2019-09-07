#!/usr/bin/env python
# 786

# Aldy source: cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, List, Tuple, Optional

import collections
import copy
import os

from . import lpinterface
from . import sam
from . import gene

from .common import *
from .gene import GeneRegion, CNConfig, Gene
from .coverage import Coverage
from .solutions import CNSolution


PCE_REGION = GeneRegion(11, 'pce')
""":obj:`aldy.common.GeneRegion`: *CYP2D7* PCE region that requires special handling (as *CYP2D6* does not have matching PCE region)."""

MAX_CN = 20.0
"""float: Maximum supported number of copies for a gene."""

# Model parameters

LEFT_FUSION_PENALTY = 0.1
"""float: Penalty applied to left fusions (0 for no penalty)."""

PCE_PENALTY_COEFF = 1.5
"""float: Error scaling applied to PCE regions (1 for no penalty)."""

MAX_CN_ERROR = 10.0
"""float: Absolute error upper bound."""

PARSIMONY_PENALTY = 0.5
"""float: Penalty applied to each gene copy (0 for no penalty)."""


def estimate_cn(gene: Gene,
                coverage: Coverage,
                solver: str,
                gap: float = 0,
                user_solution=None,
                debug: Optional[str] = None) -> List[CNSolution]:
   """
   Estimate optimal copy number configuration given the gene and read data.

   Args:
      gene (:obj:`aldy.gene.Gene`):
         Gene instance.
      coverage (:obj:`aldy.coverage.Coverage`):
         Read data coverage instance.
      solver (str):
         ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.
      gap (float):
         Relative optimality gap (percentage).
         Default is 0 (report only optimal solutions).
      user_solution (list[str], optional):
         User-specified list of copy number configurations.
         ILP will not run is this parameter is provided.
         Default is ``None``.
      debug (str, optional):
         If set, Aldy will create "<debug>.cn.lp" file for debug purposes.
         Default is ``None``.

   Returns:
      list[:obj:`CNSolution`]: List of copy number solutions.
   """
   if user_solution is not None:
      return [_parse_user_solution(gene, user_solution)]
   else:
      _print_coverage(gene, coverage)

      # TODO: filter CN configs with non-resent alleles
      # Calculate max. possible CN of the gene
      max_observed_cn = 1 + max(int(round(coverage.region_coverage(g, r)))
                                for g in gene.regions
                                for r in gene.regions[g])
      log.debug('Maximum CN = {}', max_observed_cn)

      region_cov = _region_coverage(gene, coverage)
      configs = _filter_configs(gene, coverage)
      sol = solve_cn_model(gene, configs, max_observed_cn, region_cov, solver, gap, debug)

      return sol


def solve_cn_model(gene: Gene,
                   cn_configs: Dict[str, CNConfig],
                   max_cn: int,
                   region_coverage: Dict[GeneRegion, Tuple[float, float]],
                   solver: str,
                   gap: float = 0,
                   debug: Optional[str] = None) -> List[CNSolution]:
   """
   Solves the copy number estimation problem (similar to the closest vector problem).

   Args:
      cn_configs (dict[str, :obj:`aldy.gene.CNConfig`]):
         Dictionary of copy number configurations (vectors) to be considered in the model.
      max_cn (int):
         Maximum copy number of the sample.
      region_coverage (dict[:obj:`aldy.common.GeneRegion`, tuple[float, float]]):
         A dictionary that provides the copy number of the main gene and the pseudogene
         for each region.
      solver (str):
         ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.
      gap (float):
         Relative optimality gap (percentage).
         Default is 0 (report only optimal solutions).
      debug (str, optional):
         If set, Aldy will create "<debug>.cn.lp" file for debug purposes.
         Default is ``None``.

   Returns:
      list[:obj:`CNSolution`]: List of copy-number solutions.

   Notes:
      Please see `Aldy paper <https://www.nature.com/articles/s41467-018-03273-1>`_ (section Methods/Copy number and structural variation estimation) for the model explanation.
   """

   log.debug('Detecting CN among: {}', ', '.join(sorted(cn_configs)))

   # Initialize the LP model
   model = lpinterface.model('AldyCN', solver)

   # List of CN configurations (a.k.a. structures):
   # dict of (name, number): structure, where multiple copies of the same name get different numbers
   # (because of binary LP restrictions).
   # Note that number = {0, -1} represents *full configuration* with *all* pseudogenes.
   # Thus a diploid genome must have *2* full configurations.
   # Any number > 0 describes only the main gene configuration (a.k.a. "weak" configurations), and there
   # can be multiple of such configurations.
   structures = {(name, 0): structure for name, structure in cn_configs.items()}
   for a, ai in list(structures.keys()):
      # Add another complete haplotype (assuming diploid genomes)
      structures[a, -1] = copy.deepcopy(structures[a, 0])
      if cn_configs[a].kind != CNConfig.CNConfigType.DEFAULT_CN:
         continue
      for i in range(1, max_cn):
         structures[a, i] = copy.deepcopy(structures[a, 0])
         for g in structures[a, i].cn:
            if g != 0: # if it is pseudogene, remove one copy (i.e. create a "weak configuration")
               structures[a, i].cn[g] = {r: v - 1 for r, v in structures[a, i].cn[g].items()}

   # Add one binary variable to model for each structure copy.
   # Uppercase variables are LP variables
   VCN = {(a, ai): model.addVar(vtype='B', name=f'CN_{a}_{ai}') for a, ai in structures}

   # We assume diploid genome, so the number of haplotype-inducing configurations is at most 2
   diplo_inducing = model.quicksum(VCN[a] for a in VCN if a[1] <= 0)
   model.addConstr(diplo_inducing == 2, name='CDIPLO')

   # Ensure that we cannot associate any allele to deletion allele
   del_allele = gene.deletion_allele()
   for (a, ai), v in VCN.items():
      if a != del_allele:
         model.addConstr(v + VCN[del_allele, -1] <= 1, name=f"CDEL_{a}_{ai}")

   # Ensure that binary LP is properly formed (i.e. A_i <= A_{i-1})
   for a, ai in structures:
      if ai == -1: # second haplotype (-1) is present only if the first one (0) is there
         model.addConstr(VCN[a, ai] <= VCN[a, 0], name=f"CORD_{a}_{ai}")
      elif ai > 1: # ignore 1, as A[1] can be 1 while A[0] is 0 to support cases such as *13/*13+*1
         model.addConstr(VCN[a, ai] <= VCN[a, ai - 1], name=f"CORD_{a}_{ai}")

   # Form the error variables
   VERR = {}
   json_print(debug, '    "data": {', end='')
   for r, (exp_cov0, exp_cov1) in region_coverage.items():
      print(r)
      json_print(debug, f"'{str(r)[3:-1]}': ({exp_cov0}, {exp_cov1}), ", end='')
      expr = 0
      for s, structure in structures.items():
         if r in structure.cn[0]:
            expr += structure.cn[0][r] * VCN[s]
         if len(structure.cn) > 1 and r in structure.cn[1]:
            expr -= structure.cn[1][r] * VCN[s]
      VERR[r] = model.addVar(name='E_{}{}'.format(*r), lb=-MAX_CN_ERROR, ub=MAX_CN_ERROR)
      model.addConstr(expr + VERR[r] == exp_cov0 - exp_cov1, name="CCOV_{}{}".format(*r))
   json_print(debug, '},')
   # Set objective: minimize absolute errors AND the number of alleles (max. parsimony)
   # PCE_REGION (in CYP2D7) is penalized with extra score
   objective = model.abssum(VERR.values(),
                            coeffs={'E_{}{}'.format(*PCE_REGION): PCE_PENALTY_COEFF})
   # Minimize the number of alleles among equal solutions
   # Also penalize left fusions (to further ensure max. parsimony)
   objective += model.quicksum(
      # MPICL requires only one term for each variable in the objective function
      (LEFT_FUSION_PENALTY + PARSIMONY_PENALTY) * v 
      if cn_configs[s].kind == CNConfig.CNConfigType.LEFT_FUSION
      else PARSIMONY_PENALTY * v
      for (s, k), v in VCN.items())
   model.setObjective(objective)
   if debug:
      model.dump(f'{debug}.cn.lp')

   # Solve the model
   try:
      lookup = {model.varName(v): a for (a, ai), v in VCN.items()}
      result = {}
      for status, opt, sol in model.solutions(gap):
         log.debug(f'LP status: {status}, opt: {opt}')
         sol_tuple = sorted_tuple(lookup[v] for v in sol)
         if sol_tuple not in result: # Because A[1] can be 1 while A[0] is 0, we can have duplicates
            result[sol_tuple] = CNSolution(opt, solution=list(sol_tuple), gene=gene)
      json_print(debug, '    "sol": ' + str([dict(r.solution) for r in result.values()]))
      return list(result.values())
   except lpinterface.NoSolutionsError:
      return []


def _filter_configs(gene: Gene, coverage: Coverage) -> Dict[str, CNConfig]:
   """
   Filters out all low-quality mutations and
   copy number configurations that are not supported by high-quality mutations.

   Returns:
      dict[str, :obj:`aldy.gene.CNConfig`]
   """
   cov = coverage.filtered(lambda mut, cov, total, thres: \
      Coverage.basic_filter(mut, cov, total, thres / MAX_CN))
   configs = copy.deepcopy(gene.cn_configs)
   for an in sorted(gene.cn_configs):
      if an not in gene.alleles:
         continue # this is just a CN description w/o any mutations
      if any(cov[m] <= 0 for m in gene.alleles[an].func_muts):
         s = ('{} in {}'.format(m, gene.region_at(m.pos))
              for m in gene.alleles[an].func_muts
              if cov[m] <= 0)
         log.trace('Removing {} because of {}', an, ' and '.join(s))
         del configs[an]
   return configs


def _region_coverage(gene: Gene, coverage: Coverage) -> Dict[GeneRegion, Tuple[float, float]]:
   """
   Calculate the coverage
   of the main gene and pseudogene for each region. Returns dictionary where the key
   is the region (e.g. exon 1) and the value is the coverage of the main gene (0)
   and the pseudogene (1).

   Returns:
      dict[:obj:`aldy.common.GeneRegion, tuple[float, float]]
   """

   if 1 in gene.regions: # do we have pseudogenes at all?
      cov = {r: (coverage.region_coverage(0, r), coverage.region_coverage(1, r))
             for r in gene.unique_regions if r != PCE_REGION}
      # HACK: by default, CYP2D6 does not have PCE region, so insert dummy region to allow easy calculation below
      if PCE_REGION in gene.regions[1] and PCE_REGION not in gene.regions[0]:
         cov[PCE_REGION] = (0, coverage.region_coverage(1, PCE_REGION))
      return cov
   else:
      return {r: (coverage.region_coverage(0, r), 0) for r in gene.unique_regions}


def _print_coverage(gene: Gene, coverage: Coverage) -> None:
   """
   Pretty-prints region coverage vectors.
   """
   regions = set(r for g in gene.regions for r in gene.regions[g])
   for r in sorted(regions):
      gc = coverage.region_coverage(0, r) if r in gene.regions[0] else .0
      if 1 in gene.regions:
         pc = coverage.region_coverage(1, r) if r in gene.regions[1] else .0
      else:
         pc = -1
      log.debug('Region {:>5} {:2}: {:5.2f} {} {}',
                r.kind, r.number,
                gc, '' if pc == -1 else f'{pc:5.2f}',
                f'* with diff = {gc - pc:5.2f}'
                  if pc != -1 and r in gene.unique_regions
                  else '')


def _parse_user_solution(gene: Gene, sols: List[str]) -> CNSolution:
   """
   Parses the list of user-provided solutions.

   Args:
      gene (:obj:`aldy.gene.Gene`):
         Gene instance.
      sols (list[str]):
         List of copy number configurations.

   Returns:
      :obj:`CNSolution`: User-provided solution.

   Raises:
      :obj:`aldy.common.AldyException` if a user-provided solution does not match the gene database.
   """
   for sol in sols:
      if sol not in gene.cn_configs:
         raise AldyException(
            'Given copy number solution contains unknown copy number configuration {}. '.format(sol) +
            f'Please run "aldy show --gene {gene.name}" for the list the valid configurations')
   s = CNSolution(0, solution=sols, gene=gene)
   log.debug('User-provided solution: {}', s)
   return s
