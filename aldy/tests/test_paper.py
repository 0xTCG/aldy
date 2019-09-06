#!/usr/bin/env python
# 786

# Aldy source: test_cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from nose.tools import *
import unittest
import logbook.more
import ast
import os

import aldy.genotype
from aldy.common import script_path


def solve(data, path, profile, reference):
   solver = os.getenv('ALDY_SOLVER', default='gurobi')
   sols = aldy.genotype.genotype(
      'cyp2d6',
      sam_path=path,
      profile=profile,
      output_file=None,
      reference=reference,
      solver=solver)
   assert_equal(data, [s.diplotype for s in sols])


def test_pgx1_cyp2d6():
   orig_path = os.environ.get('ALDY_INTERNAL_SAMPLES', '.')
   path = f'{orig_path}/cdc/pgrnseq-v1/bams'
   with open(script_path('aldy.tests.paper/data-pgx1.json')) as f:
      data = f.read()
      data = ast.literal_eval(data)
   for file in data:
      yield (solve,
             data[file],
             f'{path}/{file}.cram',
             'pgrnseq-v1',
             f'{orig_path}/cram-genome.fa')


def test_pgx2_cyp2d6():
   orig_path = os.environ.get('ALDY_INTERNAL_SAMPLES', '.')
   path = f'{orig_path}/baylor/pgrnseq-v2/bams'
   with open(script_path('aldy.tests.paper/data-pgx2.json')) as f:
      data = f.read()
      data = ast.literal_eval(data)
   for file in data:
      yield (solve,
             data[file],
             f'{path}/{file}.cram',
             'pgrnseq-v2',
             f'{orig_path}/cram-genome.fa')


def test_illumina_cyp2d6():
   orig_path = os.environ.get('ALDY_INTERNAL_SAMPLES', '.')
   path = f'{orig_path}/1000genomes-illumina/bams'
   with open(script_path('aldy.tests.paper/data-illumina.json')) as f:
      data = f.read()
      data = ast.literal_eval(data)
   for file in data:
      yield (solve,
             data[file],
             f'{path}/{file}.bam',
             'illumina',
             None)
