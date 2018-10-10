#!/usr/bin/env python
# 786

# Aldy source: test_cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from nose.tools import *
import unittest
import logbook.more

import aldy.cn
from aldy.gene import Gene
from aldy.tests.test_gene import TOY_GENE
from aldy.common import *


def assert_cn(gene, expected, cov, max_cn=20):
   sols = aldy.cn.solve_cn_model(gene, 
                                 cn_configs=gene.cn_configs, 
                                 max_cn=max_cn, 
                                 region_coverage=cov, 
                                 solver='gurobi')
   sols = sorted([sorted(s.solution.items()) for s in sols])
   expected = sorted([sorted(s.items()) for s in expected])
   # print(sols); print(expected)
   assert_equal(len(sols), len(expected))
   assert_equal(sols, expected)


class CNSyntheticTest(unittest.TestCase):
   _multiprocess_can_split_ = True

   def setUp(self):
      self.gene = Gene(None, 'GENE', TOY_GENE)
      # sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', level=9)
      # sh.push_application()


   def make_coverage(self, lst):
      cov = {}
      i = 0
      for r in sorted(self.gene.unique_regions):
         cov[r] = next(lst)
      return cov


   def test_basic(self):
      assert_cn(self.gene, 
                [{'1': 2}],
                self.make_coverage(zip([2,2, 2,2, 2], [2,2, 2,2, 2])))


   def test_deletion(self):
      assert_cn(self.gene, 
                [{'1': 1, '6': 1}],
                self.make_coverage(zip([1,1, 1,1, 1], [2,2, 2,2, 2])))


   def test_left_fusion(self):
      # 4 is 00001|11100 
      assert_cn(self.gene, 
                [{'4': 2}], 
                self.make_coverage(zip([0,0, 0,1, 2], [2,2, 2,1, 0])))
      
      # TODO: handle this case
      # assert_cn(self.gene, 
      #           [{'4': 2, '1': 1}], 
      #           self.make_coverage(zip([1,1, 1,2, 3], [2,2, 2,1, 0])))


   def test_right_fusion(self):
      # 5 is 11000|11222
      assert_cn(self.gene, 
                [{'1': 1, '5': 1}],
                self.make_coverage(zip([2,2, 1,1, 1], [2,2, 3,3, 3])))


   def test_multiplication(self):
      assert_cn(self.gene, 
                [{'1': 7, '5': 1}],
                self.make_coverage(zip([8,8, 7,7, 7], [2,2, 3,3, 3])))



class CNRealTest(unittest.TestCase):
   _multiprocess_can_split_ = True

   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))


   def test_many_copies_multiple_solutions(self):
      # HG00465
      assert_cn(self.gene,
                [{'1': 2, '36': 2},
                 {'1': 2, '61': 2},
                 {'1': 2, '63': 2}],
                {GeneRegion(1, 'e'): (3.8496, 2.1395), GeneRegion(1, 'i'): (3.3383, 2.5479),
                 GeneRegion(2, 'e'): (3.7298, 2.1429), GeneRegion(2, 'i'): (4.0521, 1.9735),
                 GeneRegion(3, 'e'): (4.5473, 1.8748), 
                 GeneRegion(5, 'e'): (3.9511, 2.0961), GeneRegion(5, 'i'): (4.1373, 2.0973),
                 GeneRegion(6, 'e'): (4.0175, 2.0495), GeneRegion(6, 'i'): (3.9902, 1.8785),
                 GeneRegion(9, 'e'): (2.5436, 3.6071), 
                 GeneRegion(1, 'pce'): (0, 4.1413)})


   def test_right_fusion(self):
      # HG01190
      assert_cn(self.gene,
                [{'1': 1, '68': 1}],
                {GeneRegion(1, 'e'): (1.8930, 2.0846), GeneRegion(1, 'i'): (1.9497, 2.0202),
                 GeneRegion(2, 'e'): (0.9163, 2.8526), GeneRegion(2, 'i'): (1.0604, 3.0672),
                 GeneRegion(3, 'e'): (1.1844, 2.7596), 
                 GeneRegion(5, 'e'): (1.0381, 3.3204), GeneRegion(5, 'i'): (1.1508, 3.2331),
                 GeneRegion(6, 'e'): (1.0662, 3.1833), GeneRegion(6, 'i'): (1.1335, 2.9531),
                 GeneRegion(9, 'e'): (1.3104, 2.7224), 
                 GeneRegion(1, 'pce'): (0, 3.1630)})


   def test_normal(self):
      # HG02260
      assert_cn(self.gene,
                [{'1': 2}],
                {GeneRegion(1, 'e'): (1.8404, 2.0857), GeneRegion(1, 'i'): (1.9031, 2.0085),
                 GeneRegion(2, 'e'): (2.0439, 1.9288), GeneRegion(2, 'i'): (1.9131, 1.8656),
                 GeneRegion(3, 'e'): (2.0366, 1.8861), 
                 GeneRegion(5, 'e'): (1.8515, 2.1882), GeneRegion(5, 'i'): (1.9531, 2.0443),
                 GeneRegion(6, 'e'): (2.0923, 2.1187), GeneRegion(6, 'i'): (2.0568, 1.9060),
                 GeneRegion(9, 'e'): (2.1798, 1.8956), 
                 GeneRegion(1, 'pce'): (0, 2.0132)})


   def test_deletion(self):
      # NA12336
      assert_cn(self.gene,
                [{'1': 1, '5': 1}],
                {GeneRegion(1, 'e'): (1.0308, 1.9988), GeneRegion(1, 'i'): (1.3748, 1.6491),
                 GeneRegion(2, 'e'): (0.9924, 1.6342), GeneRegion(2, 'i'): (1.0290, 1.9198),
                 GeneRegion(3, 'e'): (1.1670, 1.7540), 
                 GeneRegion(5, 'e'): (0.9421, 2.2108), GeneRegion(5, 'i'): (0.9298, 2.0529),
                 GeneRegion(6, 'e'): (0.9427, 1.9694), GeneRegion(6, 'i'): (0.9973, 1.8855),
                 GeneRegion(9, 'e'): (1.0982, 1.8465), 
                 GeneRegion(1, 'pce'): (0, 1.9869)})


   def test_right_fusion_with_copy(self):
      # NA12878
      assert_cn(self.gene,
                [{'1': 2, '68': 1}],
                {GeneRegion(1, 'e'): (2.6537, 2.0015), GeneRegion(1, 'i'): (2.4258, 2.3148),
                 GeneRegion(2, 'e'): (1.8623, 2.9641), GeneRegion(2, 'i'): (2.0050, 2.9440),
                 GeneRegion(3, 'e'): (2.2375, 2.5721), 
                 GeneRegion(5, 'e'): (1.9376, 3.0845), GeneRegion(5, 'i'): (1.9339, 3.0267),
                 GeneRegion(6, 'e'): (1.8525, 2.9231), GeneRegion(6, 'i'): (1.9799, 2.9167),
                 GeneRegion(9, 'e'): (2.1669, 2.6403), 
                 GeneRegion(1, 'pce'): (0, 3.0387)})


   def test_normal2(self):
      # NA19239
      assert_cn(self.gene,
                [{'1': 2}],
                {GeneRegion(1, 'e'): (1.6961, 2.2349), GeneRegion(1, 'i'): (2.0603, 2.0186),
                 GeneRegion(2, 'e'): (2.0000, 1.9220), GeneRegion(2, 'i'): (2.0129, 2.1490),
                 GeneRegion(3, 'e'): (2.1748, 1.9104), 
                 GeneRegion(5, 'e'): (1.9325, 2.2576), GeneRegion(5, 'i'): (2.0313, 2.0658),
                 GeneRegion(6, 'e'): (1.9981, 2.0114), GeneRegion(6, 'i'): (1.9534, 2.0146),
                 GeneRegion(9, 'e'): (2.1122, 2.0378), 
                 GeneRegion(1, 'pce'): (0, 2.1584)})


   def test_left_fusion(self):
      # NA19790
      assert_cn(self.gene,
                [{'1': 2, '78': 1}],
                {GeneRegion(1, 'e'): (1.9410, 2.0690), GeneRegion(1, 'i'): (2.0166, 1.9104),
                 GeneRegion(2, 'e'): (2.0613, 2.0358), GeneRegion(2, 'i'): (2.0889, 2.0979),
                 GeneRegion(3, 'e'): (2.1226, 1.9755), 
                 GeneRegion(5, 'e'): (2.8184, 1.2211), GeneRegion(5, 'i'): (2.9438, 1.0068),
                 GeneRegion(6, 'e'): (2.7618, 0.9867), GeneRegion(6, 'i'): (2.8389, 0.9691),
                 GeneRegion(9, 'e'): (3.1541, 1.0116), 
                 GeneRegion(1, 'pce'): (0, 1.0484)})
