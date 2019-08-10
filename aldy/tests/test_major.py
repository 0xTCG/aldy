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


class MajorSyntheticTest(unittest.TestCase):
   _multiprocess_can_split_ = True

   def setUp(self):
      self.gene = Gene(None, 'GENE', TOY_GENE)


   def make_coverage(self, lst):
      cov = {}
      i = 0
      for r in sorted(self.gene.unique_regions):
         cov[r] = next(lst)
      return cov


   def test_basic(self):
      # Test two copies of *1
      assert_cn(self.gene, 
                [{'1': 2}],
                self.make_coverage(zip([2,2, 2,2, 2], [2,2, 2,2, 2])))


   def test_deletion(self):
      # Test a single copy of *1 (*6 is deletion allele)
      assert_cn(self.gene, 
                [{'1': 1, '6': 1}],
                self.make_coverage(zip([1,1, 1,1, 1], [2,2, 2,2, 2])))

      # Test whole gene deletion
      assert_cn(self.gene, 
                [{'6': 2}],
                self.make_coverage(zip([0,0, 0,0, 0], [2,2, 2,2, 2])))

   def test_left_fusion(self):
      # Test two fused copies (*4 is defined as 00001|11100)
      assert_cn(self.gene, 
                [{'4': 2}], 
                self.make_coverage(zip([0,0, 0,1, 2], [2,2, 2,1, 0])))
      
      # Test one fused and one normal (*1) allele 
      # Note: each left fusion operates on the whole genic region; 
      #       thus, the maximum number of left fusions is 2
      assert_cn(self.gene, 
                 [{'4': 2, '1': 1}], 
                 self.make_coverage(zip([1,1, 1,2, 3], [2,2, 2,1, 0])))


   def test_right_fusion(self):
      # Test one fused and one normal (*1) allele (*5 is defined as 11000|11222)
      assert_cn(self.gene, 
                [{'1': 1, '5': 1}],
                self.make_coverage(zip([2,2, 1,1, 1], [2,2, 3,3, 3])))


   def test_multiplication(self):
      # Test twelve copies of *1
      assert_cn(self.gene, 
                [{'1': 12}],
                self.make_coverage(zip([12,12, 12,12, 12], [2,2, 2,2, 2])))

      # Test seven copies of *1 and one fused *5 copy
      assert_cn(self.gene, 
                [{'1': 7, '5': 1}],
                self.make_coverage(zip([8,8, 7,7, 7], [2,2, 3,3, 3])))



class MajorFullTest(unittest.TestCase):
   _multiprocess_can_split_ = True

   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))


