#!/usr/bin/env python
# 786

# Aldy source: test.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from __future__ import print_function
from builtins import str

import re
import os
import sys
import pandas as pd
import logbook
import logbook.more
import traceback
from multiprocessing import Pool
from functools import partial   

from aldy import genotype
from aldy.sam import SAM
from aldy.common import log, LOG_FORMAT, colorize, allele_key

def sortkey(x):
   return int(allele_key(x))

def test_single(sample, location, expected, profile, gene, solver):
   expected = [r for r in expected if not pd.isnull(r)]

   message = '{}::{}::{}'.format(sample, solver[:2], gene)
   expected = [
      [str(x).strip() for x in re.split(r'[/\+]', r) if str(x).strip() != ''] 
      for r in expected]
   expected = [tuple(sorted(r, key=sortkey)) for r in expected]

   def fix(s):
      return re.split(r'(\d+)', s)[1]
   expected = [tuple(sorted((fix(p) for p in s), key=sortkey)) for s in expected]
   
   expected, expected_new = set(expected), set(expected[1:])
   
   solutions = genotype.genotype(location, 
      'tmp/{}_{}_{}.out'.format(sample, gene, profile),
      'tmp/{}_{}_{}.log'.format(sample, gene, profile),
      gene, profile, 0.5, 
      solver,
      cn_solution=None,
      reference='/data/cb/inumanag/aldy/cram-genome.fa',
      cn_neutral_region=None)

   orig_solutions = '; '.join(','.join(s[1]) for s in solutions)
   orig_expected = '; '.join(','.join(s) for s in expected)
   solutions = set(tuple(
      sorted((fix(p) for p in s[1]), key=sortkey)) for s in solutions)

   if solutions == expected:
      logbook.warn('{:20} {} {:25} == {}', message, colorize('OK   ', 'green'), orig_solutions, orig_expected)
      return 1
   elif solutions <= expected and len(solutions) != 0:
      if solutions == expected_new:
         logbook.warn('{:20} {} {:25} == {}', message, colorize('OK=  ', 'green'), orig_solutions, orig_expected)
      else:
         logbook.warn('{:20} {} {:25} <= {}', message, colorize('OK?  ', 'green'), orig_solutions, orig_expected)
      return 2
   elif len(expected & solutions) > 0:
      logbook.warn('{:20} {} {:25} =~ {}', message, colorize('MULT ', 'yellow'), orig_solutions, orig_expected)
      return 3
   else:
      logbook.error('{:20} {} {:25} != {}', message, colorize('FAIL ', 'red'), orig_solutions, orig_expected)
      return 0

def f(s):
   sample, gene, expected, loc, profile, solver = s
   message = '{}::{}::{}'.format(sample, solver[:2], gene)
   if not os.path.exists(loc): # or sample == 'PGXT122':
      log.error('{:20} ERROR File {} not found', message, loc)
      return (sample, gene, profile, solver), 0
   try:
      res = test_single(sample, loc, expected, profile, gene, solver)
   except Exception as e:
      log.error('{:20} ERROR {}\n{}', message, 
         e,
         traceback.format_exc())
      res = 0
   return (sample, gene, profile, solver), res

if __name__ == "__main__":
   path = os.path.dirname(os.path.realpath(__file__))
   
   sh = logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='WARNING')
   sh.push_application()
   log_file = path + '/test-results.txt'
   fh = logbook.FileHandler(log_file, format_string=LOG_FORMAT, mode='w', bubble=True, level='WARNING')
   fh.push_application()

   SAM.CACHE = False

   np = int(sys.argv[1])
   log.warn('>> Using {} processes', np)

   def get_samples(sheet, gene, samples_path, profile, solver):
      if 'CDC' in sheet:
         data = pd.read_excel(path + '/experiments.xlsx', 
            header=[0, 2], sheet_name=sheet, parse_dates=False)
         # print (data.columns)
         # data = data.set_index(['Sample'])
      else:
         data = pd.read_excel(path + '/experiments.xlsx', 
            header=[0, 1], sheet_name=sheet, parse_dates=False)
      samples = []
      for r in data[gene].iterrows():
         sample_id = r[0]
         if isinstance(sample_id, tuple):
            sample_id = sample_id[0].upper()
         p = [sample_id, r[1].Validated]
         if 'PGRNseq' in r[1] and not pd.isnull(r[1].PGRNseq):
            p += r[1].PGRNseq.split(';')
         elif 'Additional' in r[1] and not pd.isnull(r[1].Additional):
            p += r[1].Additional.split(';')
         samples.append(p)

      return [(s[0], gene, tuple(s[1:]), samples_path.format(s[0]), profile, solver) for s in samples]
   
   samples = []
   for gene in 'CYP2D6 CYP2A6 CYP2C19 CYP2C8 CYP2C9 CYP3A4 CYP3A5 CYP4F2 TPMT DPYD'.split():
      loc = '/data/cb/inumanag/aldy/cdc/pgrnseq-v1/bams/{}.cram'
      samples += get_samples('PGRNseq-v1 (CDC)', gene, loc, 'pgrnseq-v1', 'gurobi')
   for gene in ['CYP2D6']:
      loc = '/data/cb/inumanag/aldy/baylor/pgrnseq-v2/bams/{}.cram'
      samples += get_samples('PGRNseq-v2', gene, loc, 'pgrnseq-v2', 'gurobi')
   
   pool = Pool(processes=np)
   result = pool.map(f, samples) #[y for y in samples if y[0] == 'NA17012' and y[1] == 'CYP2D6'])
   result = list(result)
   logbook.warn(
      'Passed {} out of {} ({} subset, {} multi)',
      sum(1 for _, x in result if x > 0),
      len(result),
      sum(1 for x in result if x == 2),
      sum(1 for x in result if x == 3)
   )

   fails = [':'.join(s) for s, x in result if x == 0]
   if len(fails) > 0:
      logbook.warn('Fail:\n{}', '\n   '.join(fails))

   # get_samples('PGRNseq-v2 (Old)', gene='CYP2D6', samples_path='/../data/pgrnseq-old/{}.dz', profile='pgrnseq-v2', threshold=.5)   
   # get_samples('Illumina',      gene='CYP2D6', samples_path='/../data/illumina/{}.bam',        profile='illumina', threshold=.5)
   # get_samples('Illumina (IU)', gene='CYP2D6', samples_path='/../data/illumina-milan/{}.bam',  profile='illumina', threshold=.5)
