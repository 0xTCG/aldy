#!/usr/bin/env python
# 786

# Aldy source: test_cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from nose.tools import *
from collections import Counter
import unittest

import os
import re
import aldy.cn
import logbook.more
from aldy.cn import CNSolution
from aldy.gene import Gene
from aldy.genotype import genotype
from aldy.common import *


_multiprocess_can_split_ = True

      
# def setup():
#    sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', 
#                                             level=8)
#    sh.push_application()


def solve(expected, sample, profile):
   sample_dir = os.getenv('ALDY_TEST_SAMPLE_DIR', '.')
   ref_path = os.getenv('ALDY_TEST_REF', '.')

   result = genotype(gene_db=     'cyp2d6', 
                     sam_path=    f'{sample_dir}/{sample}.cram',
                     cn_region=   aldy.sam.DEFAULT_CN_NEUTRAL_REGION,
                     profile=     profile,
                     output_file= None,
                     cn_solution= None,
                     threshold=   0.5,
                     solver=      'gurobi',
                     cache=       False, 
                     phase=       None,
                     reference=   ref_path)

   expected = [sorted(e.replace('_', ' ').split(',')) 
               for e in expected.split(';')]

   print(sample)
   print(result)
   result = [sorted([m.major_repr()[1:] for m in r.solution]) 
             for r in result]

   assert_equal(sorted(expected), sorted(result))


# The tests below require nose v1.1.1 (not later) for parallel runs

@nottest
def test_pgrnseq_cyp2d6_v2():
   with open('aldy/tests/paper/pgrn-v2.txt') as f:
      for l in f:
         sam, sol = l.strip().split()   
         yield solve, sol, sam, 'pgrnseq-v2'


@nottest
def test_pgrnseq_cyp2d6_v1():
   with open('aldy/tests/paper/pgrn-v1-cyp2d6.txt') as f:
      for l in f:
         sam, sol = l.strip().split()   
         yield solve, sol, sam, 'pgrnseq-v1'

