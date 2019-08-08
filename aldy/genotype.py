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

from . import sam
from . import cn
from . import major
from . import minor
from . import diplotype

from .common import *
from .gene import Gene, GRange
from .lpinterface import model as lp_model
from .version import __version__


#@timing
def genotype(gene_db: str, 
             sam_path: str,
             profile: Optional[str],
             output_file: Optional = sys.stdout, 
             cn_region: Optional[GRange] = None,
             cn_solution: Optional[List[str]] = None,
             threshold: float = 0.5, 
             solver: str = 'any',
             cache: bool = False,
             phase: bool = False,
             reference: Optional[str] = None,
             dump: bool = False) -> List[minor.MinorSolution]:
   """
   Genotype a sample.

   Returns:
      list[:obj:`aldy.minor.MinorSolution`]: List of solutions.

   Args:
      gene_db (str): 
         Gene name (if it is located in Aldy's gene database)
         or the location of gene YML description.
      sam_path (str): 
         Location of SAM/BAM/CRAM/DeeZ file to be genotyped.
      profile (str, optional):
         Coverage profile (e.g. 'illumina'). 
         Can be ``None`` if ``cn_solution`` is provided.
      cn_region (:obj:`aldy.common.GRange`, optional):
         Copy-number neutral region.
         Can be ``None`` (will use default CYP2D8 region or ``None`` if ``cn_solution`` is provided).
      output_file (file, optional):
         Location of the output decomposition file. 
         Provide ``None`` for no output.
         Default is ``sys.stdout``.
      cn_solution (list[str], optional):
         User-specified list of copy number configurations.
         Copy-number detection ILP will not run is this parameter is provided.
         Default is ``None``.
      threshold (float):
         Filtering threshold.
         Default is 0.5 (for 50%).
      solver (str):
         ILP solver to use. Check :obj:`aldy.lpinterface` for available solvers.
         Default is ``'any'``.
      cache (bool):
         Use Aldy caching for faster loading. Internal-use only. 
         Default is ``False``.
      phase (bool):
         Construct basic rudimentary phasing of the reads to aid the genotyping.
         Not recommended (slows down the pipeline with no tangible benefits).
         Default is ``False``.
      reference (str, optional):
         A path to the reference genome that was used to encode DeeZ or CRAM files.
         Default is ``None``.
      dump (bool):
         If true, Aldy will create "<filename>.aldy.dump" file for debug purposes.
         Default is ``False``.

   Raises:
      :obj:`aldy.common.AldyException` if the average coverage is too low (less than 2).
   """

   # Test the existence of LP solver
   _ = lp_model('init', solver)

   # Load the gene specification
   db_file = script_path('aldy.resources.genes/{}.yml'.format(gene_db.lower()))
   if os.path.exists(db_file):
      gene_db = db_file
   with open(gene_db): # Check does file exist
      pass
   gene = Gene(gene_db)

   with open(sam_path): # Check does file exist
      pass
   if not cn_region:
      cn_region = sam.DEFAULT_CN_NEUTRAL_REGION 
   sample = sam.Sample(sam_path=sam_path, 
                       gene=gene, 
                       threshold=threshold, 
                       profile=profile,
                       cache=False,
                       phase=False,
                       reference=reference,
                       cn_region=cn_region, 
                       dump=dump)

   avg_cov = sample.coverage.average_coverage()
   if avg_cov < 2:
      raise AldyException(td("""
         Average coverage of {0} for gene {1} is too low; skipping gene {1}. 
         Please ensure that {1} is present in the input SAM/BAM.""").format(avg_cov, gene.name))
   elif avg_cov < 20:
      log.warn("Average coverage is {}. We recommend at least 20x coverage for optimal results.", avg_cov)
   
   # Get copy-number solutions
   # print(f'>>NAME>> {os.path.basename(sam_path)}')
   cn_sols = cn.estimate_cn(gene, sample.coverage, solver=solver, user_solution=cn_solution)

   # Get major solutions and pick the best one
   # dmp = sample.coverage._dump()
   # print('>>COV>> {}'.format(' '.join(
   #    '{}={}'.format(
   #       p,  
   #       ','.join("{}:{:.0f}".format(m,dmp[p][m]) for m in sorted(dmp[p]))
   #    )
   #    for p in sorted(dmp)
   # )))
   
   log.info(f'Potential copy number configurations for {gene.name}:')
   major_sols = []
   for i, cn_sol in enumerate(cn_sols):
      sols = major.estimate_major(gene, sample.coverage, cn_sol, solver)
      # print('>>MAJOR>> {} {}'.format(
      #    ','.join(','.join([s] * v) for s, v in cn_sol.solution.items()),
      #    ';'.join(
      #       ','.join(','.join([s.major_repr()] * v) for s, v in m.solution.items())
      #       for m in sols)
      # ))
      log.info('  {:2}: {}', i + 1, cn_sol._solution_nice())
      major_sols += sols

   min_score = min(major_sols, key=lambda m: m.score).score
   major_sols = sorted([m for m in major_sols 
                        if abs(m.score - min_score) < SOLUTION_PRECISION], 
                       key=lambda m: m.score)

   log.info(f'Potential major star-alleles for {gene.name}:')
   for i, major_sol in enumerate(major_sols):
      log.info('  {:2}: {}', i + 1, major_sol._solution_nice())

   minor_sols = minor.estimate_minor(gene, sample.coverage, major_sols, solver)

   min_score = min(minor_sols, key=lambda m: m.score).score
   minor_sols = [m for m in minor_sols if abs(m.score - min_score) < SOLUTION_PRECISION]

   log.info(f'Best minor star-alleles for {gene.name}:')
   for i, minor_sol in enumerate(minor_sols):
      log.info('  {:2}: {}', i + 1, minor_sol._solution_nice())

   # print('>>MINOR>> {}'.format(' '.join(
   #    '{};{}'.format(round(sol.score, 2),
   #                   ','.join(s.major_repr() for s in sol.solution))
   #    for sol in sorted(minor_sols, key=lambda m: m.score)
   #    )))
   # exit(0)

   sample_name = os.path.splitext(os.path.basename(sam_path))[0]
   for sol_id, sol in enumerate(minor_sols):
      _ = diplotype.estimate_diplotype(gene, sol)
      if output_file:
         diplotype.write_decomposition(sample_name, gene, sol_id, sol, output_file)

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
