#!/usr/bin/env python
# 786

# Aldy source: test_major.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from nose.tools import *
import unittest
import collections
import pandas as pd
import ast
import logbook.more
from collections import Counter

import aldy.cn
import aldy.major
import aldy.minor
import aldy.coverage
from aldy.gene import Gene, Mutation
from aldy.major import NOVEL_MUTATION_PENAL
from aldy.common import *


def assert_major(gene, major):
   cn_sol = aldy.cn.CNSolution(0, list(Counter(major['cn']).elements()), gene)

   cov = collections.defaultdict(dict)
   for (pos, op), c in major["data"].items():
      cov[pos][op] = c
   cov = aldy.coverage.Coverage(cov, 0.5, {})
   sols = aldy.major.estimate_major(gene, cov, cn_sol, 'gurobi')

   if 'score' in major:
      for s in sols:
         assert_less(abs(major['score'] - s.score), SOLUTION_PRECISION)
   sols_expected = [sorted(Counter(c).elements(), key=str) for c in major['sol']]
   sols_parsed = [sorted(Counter({
                     tuple([k.major] + [(m[0], m[1]) for m in k.added]) if len(k.added) > 0 else k.major: v
                     for k, v in s.solution.items()
                  }).elements(), key=str)
                  for s in sols]
   assert_equal(sorted(sols_expected), sorted(sols_parsed))


sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', level='DEBUG')
# sh.push_application()

class MajorSyntheticTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.tests/toy.yml'))


   def test_basic(self):
      # Test two copies of *1
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {},
         "sol": [{'1': 2}],
         "score": 0
      })


   def test_deletion(self):
      # Test a single copy of *1C
      assert_major(self.gene, {
         "cn": {'1': 1, '6': 1},
         "data": {(105, '_'): 0, (105, 'SNP.TA'): 10},
         "sol": [{'1.a': 1, '6': 1}],
         "score": 0
      })


   def test_two_copies(self):
      # Test two copies (*1/*2)
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(111, '_'): 10, (111, 'DEL.AC'): 10,
                  (119, '_'): 20, (119, 'INS.TT'): 10},
         "sol": [{'1': 1, '2': 1}],
         "score": 0
      })
      # Test two copies (*2/*3)
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(111, '_'): 10, (111, 'DEL.AC'): 10,
                  (119, '_'): 20, (119, 'INS.TT'): 10,
                  (151, '_'): 10, (151, 'SNP.CT'): 10},
         "sol": [{'2': 1, '3': 1}],
         "score": 0
      })
      # Test slightly perturbed two copies (*2/*3)
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(111, '_'): 11, (111, 'DEL.AC'): 9,
                  (119, '_'): 22, (119, 'INS.TT'): 8,
                  (151, '_'): 10.5, (151, 'SNP.CT'): 9.5},
         "sol": [{'2': 1, '3': 1}],
         "score": 2/10 + 3/11 + 1/10
      })


   def test_multiple_copies(self):
      # Test four copies (*1+*1/*2+*2)
      assert_major(self.gene, {
         "cn": {'1': 4},
         "data": {(111, '_'): 20, (111, 'DEL.AC'): 20,
                  (119, '_'): 40, (119, 'INS.TT'): 20},
         "sol": [{'1': 2, '2': 2}],
         "score": 0
      })


   def test_left_fusion(self):
      # Test left fusion that has no SNPs (*4[*1]/*1+*2)
      assert_major(self.gene, {
         "cn": {'1': 2, '4': 1},
         "data": {(111, '_'): 10, (111, 'DEL.AC'): 10,
                  (119, '_'): 20, (119, 'INS.TT'): 10},
         "sol": [{'1': 1, '2': 1, '4/1': 1}],
         "score": 0
      })
      # Test left fusion that has SNPs (*4[*3]/*1+*3)
      assert_major(self.gene, {
         "cn": {'1': 2, '4': 1},
         "data": {(151, '_'): 10, (151, 'SNP.CT'): 20},
         "sol": [{'1': 1, '3': 1, '4/3': 1},
                 {'3': 2, '4/1': 1}],
         "score": 0
      })
      # Test left fusion combination (*4[*3]/*1C+*3)
      assert_major(self.gene, {
         "cn": {'1': 2, '4': 1},
         "data": {(105, '_'): 10, (105, 'SNP.TA'): 10,
                  (151, '_'): 10, (151, 'SNP.CT'): 20},
         "sol": [{'1.a': 1, '3': 1, '4/3': 1}],
         "score": 0
      })
      # Test left fusion combination (*4[*3]/*4[*1]+*3)
      assert_major(self.gene, {
         "cn": {'1': 1, '4': 2},
         "data": {(105, '_'): 0,  (105, 'SNP.TA'): 10,
                  (151, '_'): 20, (151, 'SNP.CT'): 10},
         "sol": [{'1.a': 1, '4/3': 1, '4/1': 1}],
         "score": 0
      })


   def test_right_fusion(self):
      # Test right fusion (*5/*1)
      assert_major(self.gene, {
         "cn": {'1': 1, '5': 1},
         "data": {(111, '_'): 10, (111, 'DEL.AC'): 10},
         "sol": [{'1': 1, '5': 1}],
         "score": 0
      })
      # Test right fusion (*5/*2)
      assert_major(self.gene, {
         "cn": {'1': 1, '5': 1},
         "data": {(111, '_'): 0,  (111, 'DEL.AC'): 20,
                  (119, '_'): 20, (119, 'INS.TT'): 10},
         "sol": [{'2': 1, '5': 1}],
         "score": 0
      })
      # Test right fusion (*5/*3)
      assert_major(self.gene, {
         "cn": {'1': 1, '5': 1},
         "data": {(111, '_'): 10, (111, 'DEL.AC'): 10,
                  (151, '_'): 0,  (151, 'SNP.CT'): 10},
         "sol": [{'3': 1, '5': 1}],
         "score": 0
      })


   def test_novel_mutations(self):
      # Test novel mutations within a single gene (other is deleted)
      assert_major(self.gene, {
        "cn": {'1': 1, '6': 1},
        "data": {(111, '_'): 0, (111, 'DEL.AC'): 20},
        "sol": [{('1', (111, 'DEL.AC')): 1, '6': 1}],
        "score": NOVEL_MUTATION_PENAL
      })
      # Test novel mutations
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(111, '_'): 10, (111, 'DEL.AC'): 10},
         "sol": [{'1': 1, ('1', (111, 'DEL.AC')): 1}],
         "score": NOVEL_MUTATION_PENAL
      })


class MajorRealTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))


   def test_normal(self):
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {},
         "sol": [{'1': 2}]
      })


   def test_different_1(self): # NA17102/v1
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(42522612, '_'): 316, (42522612, 'SNP.CG'): 279,
                  (42523942, '_'): 266, (42523942, 'SNP.GA'): 291,
                  (42524929, '_'): 575, (42524929, 'INS.ggggcgaaaggggcgaaa'): 161,
                  (42525771, '_'): 123, (42525771, 'SNP.GA'): 127},
         "sol": [{'1': 1, '40': 1}]
      })


   def test_many_alleles_1(self):
      # HG00436/v1
      assert_major(self.gene, {
         "cn": {'1': 3},
         "data": {(42522612, '_'): 357, (42522612, 'SNP.CG'): 634,
                  (42523942, '_'): 332, (42523942, 'SNP.GA'): 622,
                  (42526668, '_'): 983, (42526668, 'SNP.CT'): 517,
                  (42528381, '_'): 16,  (42528381, 'SNP.GC'): 19},
         "sol": [{'2': 1, '2.a': 1, '71': 1}]
      })
      # NA24217/v1
      assert_major(self.gene, {
         "cn": {'1': 4},
         "data": {(42522612, '_'): 0,   (42522612, 'SNP.CG'): 1231,
                  (42523804, '_'): 287, (42523804, 'SNP.CT'): 693,
                  (42523942, '_'): 0,   (42523942, 'SNP.GA'): 1138,
                  (42528381, '_'): 29,  (42528381, 'SNP.GC'): 8, },
         "sol": [{'2': 1, '41': 3}]
      })


   def test_right_fusion_1(self): # NA23878/v1
      data = {(42522612, '_'): 0,   (42522612, 'SNP.CG'): 953,
              (42524946, '_'): 760, (42524946, 'SNP.CT'): 1146,
              (42525810, '_'): 163, (42525810, 'SNP.TC'): 190,
              (42525820, '_'): 128, (42525820, 'SNP.GT'): 126,
              (42526693, '_'): 721, (42526693, 'SNP.GA'): 1425}
      assert_major(self.gene, {
         "cn": {'1': 2, '61': 1},
         "data": data,
         "sol": [{'4': 2, '61': 1}]
      })
      assert_major(self.gene, {
         "cn": {'1': 2, '36': 1},
         "data": data,
         "sol": [{'4': 2, '83': 1},
                 {'39': 1, '4': 1, '4.i': 1}]
      })


   def test_left_fusion_1(self): # NA19785/v1
      assert_major(self.gene, {
         "cn": {'1': 2, '79': 1},
         "data": {(42522612, '_'): 526, (42522612, 'SNP.CG'): 905,
                  (42523942, '_'): 446, (42523942, 'SNP.GA'): 822,
                  (42528381, '_'): 14,  (42528381, 'SNP.GC'): 23},
         "sol": [{'2': 1, '2.a': 1, '79/1': 1},
                 {'1': 1, '2': 1, '79/2': 1},
                 {'2': 1, '34': 1, '79/10': 1},
                 {'2': 1, '39': 1, '79/34': 1}]
      })


   def test_novel_allele_1(self): # NA17012/v1
      assert_major(self.gene, {
         "cn": {'1': 1, '5': 1},
         "data": {(42522612, '_'): 0,  (42522612, 'SNP.CG'): 288,
                  (42526693, '_'): 20, (42526693, 'SNP.GA'): 436,
                  (42526716, '_'): 19, (42526716, 'SNP.CT'): 437},
         "sol": [{('10', (42526716, 'SNP.CT')): 1, '5': 1}]
      })


   def test_different_2(self):
      # NA19828/v2
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(42522612, '_'): 0,   (42522612, 'SNP.CG'): 698,
                  (42523942, '_'): 0,   (42523942, 'SNP.GA'): 862,
                  (42525771, '_'): 355, (42525771, 'SNP.GA'): 339},
         "sol": [{'17': 1, '2.a': 1}]
      })
      # NA10861/v2
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(42522612, '_'): 0,   (42522612, 'SNP.CG'): 797,
                  (42523942, '_'): 484, (42523942, 'SNP.GA'): 483,
                  (42524946, '_'): 398, (42524946, 'SNP.CT'): 375,
                  (42525810, '_'): 339, (42525810, 'SNP.TC'): 185,
                  (42525820, '_'): 322, (42525820, 'SNP.GT'): 153,
                  (42526693, '_'): 406, (42526693, 'SNP.GA'): 377,
                  (42526762, '_'): 390, (42526762, 'SNP.CT'): 406,
                  (42528381, '_'): 109, (42528381, 'SNP.GC'): 70},
         "sol": [{'35': 1, '4': 1}]
      })
      # NA19239/v2
      assert_major(self.gene, {
         "cn": {'1': 2},
         "data": {(42522612, '_'): 359, (42522612, 'SNP.CG'): 374,
                  (42523942, '_'): 468, (42523942, 'SNP.GA'): 492,
                  (42525771, '_'): 363, (42525771, 'SNP.GA'): 328,
                  (42526656, '_'): 572, (42526656, 'INS.a'): 207},
         "sol": [{'15': 1, '17': 1}]
      })


   def test_multiple_solutions_2(self): # NA18507/v2
      assert_major(self.gene, {
         "cn": {'1': 3},
         "data": {(42522612, '_'): 0,   (42522612, 'SNP.CG'): 1113,
                  (42523942, '_'): 881, (42523942, 'SNP.GA'): 528,
                  (42524946, '_'): 457, (42524946, 'SNP.CT'): 737,
                  (42526693, '_'): 419, (42526693, 'SNP.GA'): 824},
         "sol": [{'39': 1, '4.b': 1, '4.g': 1},
                 {'2.a': 1, '4.b': 2}]
      })


   def test_many_alleles_2(self): # HG00463/v2
      data = {(42522612, '_'): 0, (42522612, 'SNP.CG'): 1008,
              (42526693, '_'): 0, (42526693, 'SNP.GA'): 1977}
      assert_major(self.gene, {
         "cn": {'1': 2, '36': 2},
         "data": data,
         "sol": [{'10': 2, '36': 2}]
      })
      assert_major(self.gene, {
         "cn": {'1': 2, '36': 1, '61': 1},
         "data": data,
         "sol": [{'10': 2, '36': 1, '61': 1}]
      })


   def test_right_fusion_2(self): # NA12878/v2
      assert_major(self.gene, {
         "cn": {'1': 2, '68': 1},
         "data": {(42522612, '_'): 340, (42522612, 'SNP.CG'): 380,
                  (42524243, '_'): 485, (42524243, 'DEL.T'): 436,
                  (42524946, '_'): 543, (42524946, 'SNP.CT'): 353,
                  (42525810, '_'): 322, (42525810, 'SNP.TC'): 173,
                  (42525820, '_'): 314, (42525820, 'SNP.GT'): 153,
                  (42526693, '_'): 362, (42526693, 'SNP.GA'): 636},
         "sol": [{'3': 1, '4': 1, '68': 1}]
      })


   def test_left_fusion_2(self): # NA19790/v2
      assert_major(self.gene, {
         "cn": {'1': 2, '78': 1},
         "data": {(42522612, '_'): 338, (42522612, 'SNP.CG'): 638,
                  (42523942, '_'): 356, (42523942, 'SNP.GA'): 790,
                  (42528381, '_'): 69, (42528381, 'SNP.GC'): 64},
         "sol": [{'2': 1, '34': 1, '78/4': 1},
                 {'1': 1, '2': 1, '78/2': 1},
                 {'2': 1, '2.a': 1, '78/1': 1},
                 {'2': 1, '39': 1, '78/34': 1}]
      })
