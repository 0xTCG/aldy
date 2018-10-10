#!/usr/bin/env python
# 786

# Aldy source: test_minor.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from nose.tools import *
import unittest
import collections

import aldy.cn
import aldy.major
import aldy.minor
import aldy.coverage
from aldy.gene import Gene, Mutation
from aldy.common import *


def M(p, c):
   return Mutation(p, c)


def assert_minor(gene, majors, mutations, solution):
   cns = [gene.alleles[a].cn_config for a in majors]
   cn_sol = aldy.cn.CNSolution(0, cns, gene)
   print(cn_sol)
   major_sol = aldy.major.MajorSolution(0, collections.Counter(majors), cn_sol, {})

   cov_dict = collections.defaultdict(lambda: collections.defaultdict(int))
   for m, cov in mutations.items():
      cov_dict[m.pos][m.op] = cov
   cov = aldy.coverage.Coverage(cov_dict, 0.5, {})

   sol = aldy.minor.estimate_minor(
      gene, 
      cov,
      [major_sol],
      'gurobi',
      lambda mut, cov, total, thres: True)
   print(sol)

   assert_equal(len(sol), 1)
   sol = sol[0]
   ref_solution = sorted(allele_number(s) for s in solution)
   usr_solution = sorted(allele_number(s[1]) for s in sol.solution)
   assert_equal(ref_solution, usr_solution)
   return sol
   
class MajorFullTest(unittest.TestCase):
   _multiprocess_can_split_ = True

   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))


