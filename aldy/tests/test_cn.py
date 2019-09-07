#!/usr/bin/env python
# 786

# Aldy source: test_cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from nose.tools import *
import unittest
import logbook.more
import os

import aldy.cn
from aldy.cn import PARSIMONY_PENALTY, LEFT_FUSION_PENALTY
from aldy.gene import Gene
from aldy.common import *


def assert_cn(gene, expected, cov, expected_obj=None):
   solver = os.getenv('ALDY_SOLVER', default='gurobi')
   sols = aldy.cn.solve_cn_model(gene,
                                 cn_configs=gene.cn_configs,
                                 max_cn=20,
                                 region_coverage=cov,
                                 solver=solver)
   if expected_obj:
      for s in sols:
         assert_less(abs(s.score - expected_obj), SOLUTION_PRECISION)
   sols = sorted([sorted(s.solution.items()) for s in sols])
   expected = sorted([sorted(s.items()) for s in expected])
   print(sols)
   assert_equal(len(sols), len(expected))
   assert_equal(sols, expected)


sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', level='DEBUG')
# sh.push_application()


class CNSyntheticTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.tests/toy.yml'))


   def make_coverage(self, lst):
      cov = {}
      for r in sorted(self.gene.unique_regions):
         cov[r] = next(lst)
      return cov


   def test_basic(self):
      # Test two copies of *1
      assert_cn(self.gene,
                [{'1': 2}],
                self.make_coverage(zip([2,2, 2,2, 2], [2,2, 2,2, 2])),
                2 * PARSIMONY_PENALTY)
      # Test two copies of *1 with slightly perturbed coverage
      assert_cn(self.gene,
                [{'1': 2}],
                self.make_coverage(zip([1.8,2.2, 2.1,1.9, 1.7], [2.05,1.95, 2,2, 2.7])),
                2 * PARSIMONY_PENALTY + (0.2 * 2 + 0.1 * 2 + 0.3 + 0.05 * 2 + 0.7))


   def test_deletion(self):
      # Test a single copy of *1 (*6 is deletion allele)
      assert_cn(self.gene,
                [{'1': 1, '6': 1}],
                self.make_coverage(zip([1,1, 1,1, 1], [2,2, 2,2, 2])),
                2 * PARSIMONY_PENALTY)
      # Test whole gene deletion
      assert_cn(self.gene,
                [{'6': 2}],
                self.make_coverage(zip([0,0, 0,0, 0], [2,2, 2,2, 2])),
                2 * PARSIMONY_PENALTY)
      # TODO: test 2 deletions with no coverage
      # sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', level='DEBUG')
      # sh.push_application()
      # assert_cn(self.gene,
      #           [{'6': 2}],
      #           self.make_coverage(zip([0,0, 0,0, 0], [0,0, 0,0, 0])),
      #           2 * PARSIMONY_PENALTY)


   def test_left_fusion(self):
      # Test two fused copies (*4 is defined as 00011|11100)
      assert_cn(self.gene,
                [{'4': 2}],
                self.make_coverage(zip([0,0, 0,2, 2], [2,2, 2,0, 0])),
                2 * PARSIMONY_PENALTY + 2 * LEFT_FUSION_PENALTY)
      # Test one fused and one normal (*1) allele
      # Note: each left fusion operates on the whole genic region;
      #       thus, the maximum number of left fusions is 2
      assert_cn(self.gene,
                 [{'4': 2, '1': 1}],
                 self.make_coverage(zip([1,1, 1,3, 3], [2,2, 2,0, 0])),
                 3 * PARSIMONY_PENALTY + 2 * LEFT_FUSION_PENALTY)


   def test_right_fusion(self):
      # Test one fused and one normal (*1) allele (*5 is defined as 11000|11222)
      assert_cn(self.gene,
                [{'1': 1, '5': 1}],
                self.make_coverage(zip([2,2, 1,1, 1], [2,2, 3,3, 3])),
                2 * PARSIMONY_PENALTY)


   def test_multiplication(self):
      # Test twelve copies of *1
      assert_cn(self.gene,
                [{'1': 12}],
                self.make_coverage(zip([12,12, 12,12, 12], [2,2, 2,2, 2])),
                12 * PARSIMONY_PENALTY)
      # Test seven copies of *1 and one fused *5 copy
      assert_cn(self.gene,
                [{'1': 7, '5': 1}],
                self.make_coverage(zip([8,8, 7,7, 7], [2,2, 3,3, 3])),
                8 * PARSIMONY_PENALTY)


class CNRealTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))


   def make_coverage(self, d):
      cov = {}
      for k, v in d.items():
         k = re.split(r'(\d+)', k)
         cov[GeneRegion(int(k[1]), k[2])] = v
      return cov


   def test_many_copies_multiple_solutions(self):
      # HG00465
      assert_cn(self.gene,
                [{'1': 2, '36': 2},
                 {'1': 2, '61': 2},
                 {'1': 2, '63': 2},
                 {'1': 2, '36': 1, '61': 1},
                 {'1': 2, '36': 1, '63': 1},
                 {'1': 2, '61': 1, '63': 1}],
                self.make_coverage({
                   '1e': (3.8, 2.1), '1i': (3.3, 2.5), '2e': (3.7, 2.1), '2i': (4.0, 1.9),
                   '5e': (3.9, 2.0), '5i': (4.1, 2.0), '6e': (4.0, 2.0), '6i': (3.9, 1.8),
                   '3e': (4.5, 1.8), '9e': (2.5, 3.6), '11pce': (0, 4.1)
                }))


   def test_right_fusion(self):
      # HG01190
      assert_cn(self.gene,
                [{'1': 1, '68': 1}],
                self.make_coverage({
                   '1e': (1.8, 2.0), '1i': (1.9, 2.0), '2e': (0.9, 2.8), '2i': (1.0, 3.0),
                   '5e': (1.0, 3.3), '5i': (1.1, 3.2), '6e': (1.0, 3.1), '6i': (1.1, 2.9),
                   '3e': (1.1, 2.7), '9e': (1.3, 2.7), '11pce': (0, 3.1)
                }))


   def test_normal(self):
      # HG02260
      assert_cn(self.gene,
                [{'1': 2}],
                self.make_coverage({
                   '1e': (1.8, 2.0), '1i': (1.9, 2.0), '2e': (2.0, 1.9), '2i': (1.9, 1.8),
                   '5e': (1.8, 2.1), '5i': (1.9, 2.0), '6e': (2.0, 2.1), '6i': (2.0, 1.9),
                   '3e': (2.0, 1.8), '9e': (2.1, 1.8), '11pce': (0, 2.0)
                }))


   def test_deletion(self):
      # NA12336
      assert_cn(self.gene,
                [{'1': 1, '5': 1}],
                self.make_coverage({
                   '1e': (1.0, 1.9), '1i': (1.3, 1.6), '2e': (0.9, 1.6), '2i': (1.0, 1.9),
                   '5e': (0.9, 2.2), '5i': (0.9, 2.0), '6e': (0.9, 1.9), '6i': (0.9, 1.8),
                   '3e': (1.1, 1.7), '9e': (1.0, 1.8), '11pce': (0, 1.9)
                }))


   def test_right_fusion_with_copy(self):
      # NA12878
      assert_cn(self.gene,
                [{'1': 2, '68': 1}],
                self.make_coverage({
                   '1e': (2.6, 2.0), '1i': (2.4, 2.3), '2e': (1.8, 2.9), '2i': (2.0, 2.9),
                   '5e': (1.9, 3.0), '5i': (1.9, 3.0), '6e': (1.8, 2.9), '6i': (1.9, 2.9),
                   '3e': (2.2, 2.5), '9e': (2.1, 2.6), '11pce': (0, 3.0)
                }))


   def test_normal2(self):
      # NA19239
      assert_cn(self.gene,
                [{'1': 2}],
                self.make_coverage({
                   '1e': (1.6, 2.2), '1i': (2.0, 2.0), '2e': (2.0, 1.9), '2i': (2.0, 2.1),
                   '5e': (1.9, 2.2), '5i': (2.0, 2.0), '6e': (1.9, 2.0), '6i': (1.9, 2.0),
                   '3e': (2.1, 1.9), '9e': (2.1, 2.0), '11pce': (0, 2.1)
                }))


   def test_left_fusion(self):
      # NA19790
      assert_cn(self.gene,
                [{'1': 2, '78': 1}, {'1': 2, '67': 1}],
                self.make_coverage({
                   '1e': (1.9, 2.0), '1i': (2.0, 1.9), '2e': (2.0, 2.0), '2i': (2.0, 2.0),
                   '5e': (2.8, 1.2), '5i': (2.9, 1.0), '6e': (2.7, 0.9), '6i': (2.8, 0.9),
                   '3e': (2.1, 1.9), '9e': (3.1, 1.0), '11pce': (0, 1.0)
                }))
