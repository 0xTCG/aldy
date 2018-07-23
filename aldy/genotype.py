# 786

# Aldy source: genotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from builtins import str

import os
import re
import sys
import copy
import logbook
import logbook.more
import platform

from . import cn
from . import filtering
from . import protein
from . import refiner
from . import sam
from . import diplotype
from . import remap

from .common import *
from .gene import Gene
from .version import __version__


@timing
def genotype(sample_path, output, log_output, gene, profile, threshold, solver, cn_region, cn_solution, do_remap=False):
   with open(sample_path): # Check does file exist
      pass

   log.info('Gene: {}', gene.upper())
   
   database_file = script_path('aldy.resources.genes', '{}.yml'.format(gene.lower()))
   gene = Gene(database_file)

   if cn_region is not None:
      r = re.match(r'^(.+?):(\d+)-(\d+)$', cn_region) 
      if not r:
         raise AldyException('Parameter --cn-neutral={} is not in the format chr:start-end (where start and end are numbers)'.format(cn_region))
      ch = r.group(1)
      if ch.startswith('chr'): 
         ch = ch[3:]
      gene.cnv_region = (ch, r.group(2), r.group(3))
      log.warn('Using {} as CN-neutral region', cn_region)
   gene_backup = copy.deepcopy(gene)

   sample_name = os.path.basename(sample_path)
   sample_name = os.path.splitext(sample_name)[0]
   sample = sam.SAM(sample_path, gene, threshold)

   if sample.avg_coverage < 2:
      raise AldyException("Average coverage of {0} for gene {1} is too low; skipping gene {1}. Please ensure that {1} is present in the input SAM/BAM.".format(sample.avg_coverage, gene.name))
   elif sample.avg_coverage < 20:
      log.warn("Average coverage is {}. We recommend at least 20x coverage for optimal results.", sample.avg_coverage)
   
   gene.alleles = filtering.initial_filter(gene, sample)
   cn_sol = cn.estimate_cn(gene, sample, profile, cn_solution, solver)

   if do_remap:
      log.critical('Remapping! Stay tuned...')
      cn_sol = list(cn_sol.values())[0] #!! TODO IMPORTANT just use furst CN for now
      new_path = remap.remap(sample_path, gene, sample, cn_sol)
      gene = gene_backup # refactor this somehow...
      sample = sam.SAM(new_path, gene, threshold)
      gene.alleles = filtering.initial_filter(gene, sample)
      cn_sol = cn.estimate_cn(gene, sample, profile, cn_solution, solver)

   score, init_sol = protein.get_initial_solution(gene, sample, cn_sol, solver)
   score, sol = refiner.get_refined_solution(gene, sample, init_sol, solver)

   sol = diplotype.assign_diplotype(sample_name, gene, sol, output)

   return sol

