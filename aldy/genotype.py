# 786

# Aldy source: genotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, List, Tuple, Optional

import os
import re
import sys
import copy
import logbook
import logbook.more
import platform

from . import cn
from . import protein
from . import refiner
from . import sam
from . import diplotype
from . import remap

from .common import *
from .gene import Gene, GRange
from .version import __version__


DEFAULT_CN_NEUTRAL_REGION = GRange('22', 42547463, 42548249)


@timing
def genotype(gene_db: str, sam_path: str, profile: str,
             threshold: float = 0.5, 
             cn_region: GRange = DEFAULT_CN_NEUTRAL_REGION,
             solver: str = 'any',
             user_cn: Optional[List[str]] = None) -> List[refiner.MinorSolution]:
   """
   """

   gene = Gene(gene_db)
   sample = sam.Sample(sam_path=sam_path, gene=gene, 
      threshold=threshold, cn_region=cn_region, profile=profile)
      # reference, cache, phase)

   avg_cov = sample.coverage.average_coverage()
   if avg_cov < 2:
      raise AldyException(td("""
         Average coverage of {0} for gene {1} is too low; skipping gene {1}. 
         Please ensure that {1} is present in the input SAM/BAM.""").format(avg_cov, gene.name))
   elif avg_cov < 20:
      log.warn("Average coverage is {}. We recommend at least 20x coverage for optimal results.", avg_cov)
   
   # Get copy-number solutions
   cn_sols = cn.estimate_cn(gene, sam=sample, solver=solver, user_solution=user_cn)

   # Get major solutions and pick the best one
   major_sols = [sol 
                 for cn_sol in cn_sols
                 for sol in protein.estimate_major(gene, sample, cn_sol, solver)]
   print(major_sols)
   min_score = min(major_sols, key=lambda m: m.score).score
   major_sols = [m for m in major_sols if abs(m.score - min_score) < 1e-6]

   minor_sols = [refiner.estimate_minor(gene, sample, ma_sol, solver)
                 for ma_sol in major_sols]
   min_score = min(minor_sols, key=lambda m: m.score).score
   minor_sols = [m for m in minor_sols if abs(m.score - min_score) < 1e-6]

   sample_name = os.path.splitext(os.path.basename(sam_path))[0]
   diplotype.write_header(sys.stdout)
   for sol_id, sol in enumerate(minor_sols):
      diplotype.estimate_diplotype(gene, sol)
      log.warn('{}', sol.diplotype)
      diplotype.write_decomposition(sample_name, gene, sol_id, sol, sys.stdout)

   # if do_remap != 0:
   #    log.critical('Remapping! Stay tuned...')
   #    cn_sol = list(cn_sol.values())[0] #!! TODO IMPORTANT just use furst CN for now
   #    new_path = remap.remap(sample_path, gene, sample, cn_sol, remap_mode=do_remap)
   #    gene = gene_backup # refactor this somehow...
   #    sam.SAM.CACHE = False
   #    sample = sam.SAM(new_path, gene, threshold)
   #    gene.alleles = filtering.initial_filter(gene, sample)
   #    cn_sol = cn.estimate_cn(gene, sample, profile, cn_solution, solver)

   return minor_sols


def genotype_init(sample_path, output, log_output, gene, profile, threshold, solver, cn_region, cn_solution, do_remap):
   """
   Genotypes the provided sample and returns the genotyping results.
   """

   with open(sample_path): # Check does file exist
      pass
   log.info('Gene: {}', gene.upper())
   
   database_file = script_path('aldy.resources.genes', '{}.yml'.format(gene.lower()))

   if cn_region is not None:
      r = re.match(r'^(.+?):(\d+)-(\d+)$', cn_region) 
      if not r:
         raise AldyException('Parameter --cn-neutral={} is not in the format chr:start-end (where start and end are numbers)'.format(cn_region))
      ch = r.group(1)
      if ch.startswith('chr'): 
         ch = ch[3:]
      cn_region = GRange(ch, int(r.group(2)), int(r.group(3)))
      log.info('Using {} as CN-neutral region', cn_region)
   else: 
      cn_region = DEFAULT_CN_NEUTRAL_REGION
   
   # gene_backup = copy.deepcopy(gene)
   gene = Gene(database_file)
   #sample_name = os.path.basename(sample_path)
   #sample_name = os.path.splitext(sample_name)[0]
   sol = genotype(database_file, sample_path, profile, threshold, cn_region, solver)
   print(sol)
   
   exit(0)
   return sol

