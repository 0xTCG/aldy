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
             phase: bool = False,
             reference: Optional[str] = None,
             debug: Optional[str] = None) -> List[minor.MinorSolution]:
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
      phase (bool):
         Construct basic rudimentary phasing of the reads to aid the genotyping.
         Not recommended (slows down the pipeline with no tangible benefits).
         Default is ``False``.
      reference (str, optional):
         A path to the reference genome that was used to encode DeeZ or CRAM files.
         Default is ``None``.
      debug (str, optional):
         The debug prefix for the debugging information. ``None`` for no debug information.
         Default is ``None``.

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
                       phase=False,
                       reference=reference,
                       cn_region=cn_region, 
                       debug=debug)

   avg_cov = sample.coverage.average_coverage()
   if avg_cov < 2:
      raise AldyException(td("""
         Average coverage of {0} for gene {1} is too low; skipping gene {1}. 
         Please ensure that {1} is present in the input SAM/BAM.""").format(avg_cov, gene.name))
   elif avg_cov < 20:
      log.warn("Average coverage is {}. We recommend at least 20x coverage for optimal results.", avg_cov)
   
   json_print(debug, f'"{os.path.basename(sam_path).split(".")[0]}": {{')

   # Get copy-number solutions
   json_print(debug, '  "cn": {')
   cn_sols = cn.estimate_cn(gene, sample.coverage, solver=solver, user_solution=cn_solution, debug=debug)
   json_print(debug, '  },')

   # Get major solutions and pick the best one
   log.info(f'Potential copy number configurations for {gene.name}:')
   major_sols = []
   
   json_print(debug, '  "major": [', end='')
   for i, cn_sol in enumerate(cn_sols):
      sols = major.estimate_major(gene, sample.coverage, cn_sol, solver, debug=debug)
      log.info('  {:2}: {}', i + 1, cn_sol._solution_nice())
      major_sols += sols
   json_print(debug, '],')

   min_score = min(major_sols, key=lambda m: m.score).score
   major_sols = sorted([m for m in major_sols 
                        if abs(m.score - min_score) < SOLUTION_PRECISION], 
                       key=lambda m: m.score)

   log.info(f'Potential major star-alleles for {gene.name}:')
   for i, major_sol in enumerate(major_sols):
      log.info('  {:2}: {}', i + 1, major_sol._solution_nice())
   
   json_print(debug, '  "minor": [', end='')
   minor_sols = minor.estimate_minor(gene, sample.coverage, major_sols, solver, debug=debug)
   min_score = min(minor_sols, key=lambda m: m.score).score
   minor_sols = [m for m in minor_sols if abs(m.score - min_score) < SOLUTION_PRECISION]
   json_print(debug, ']')

   log.info(f'Best minor star-alleles for {gene.name}:')
   for i, minor_sol in enumerate(minor_sols):
      log.info('  {:2}: {}', i + 1, minor_sol._solution_nice())
   
   sample_name = os.path.splitext(os.path.basename(sam_path))[0]
   for sol_id, sol in enumerate(minor_sols):
      s = diplotype.estimate_diplotype(gene, sol)
      if output_file:
         diplotype.write_decomposition(sample_name, gene, sol_id, sol, output_file)
   
   json_print(debug, '},')

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
