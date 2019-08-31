#!/usr/bin/env python
# 786

# Aldy source: test_minor.py
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
from aldy.major import MajorSolution, SolvedAllele
from aldy.minor import ADD_PENALTY_FACTOR, MISS_PENALTY_FACTOR
from aldy.common import *


def assert_minor(gene, data, shallow=False):
   cn_sol = aldy.cn.CNSolution(0, list(Counter(data['cn']).elements()), gene)

   cov = collections.defaultdict(dict)
   for (pos, op), c in data["data"].items():
      cov[pos][op] = c
   cov = aldy.coverage.Coverage(cov, 0.5, {})

   major_solved = {SolvedAllele(maj, None, tuple(), tuple()) if not isinstance(maj, tuple)
                   else SolvedAllele(maj[0], None, 
                                     tuple(Mutation(mp, mo, True) for mp, mo in maj[1]), tuple()): cnt 
                   for maj, cnt in data['major'].items()}
   major = MajorSolution(0, major_solved, cn_sol)
   sols = aldy.minor.estimate_minor(gene, cov, [major], 'gurobi')

   if 'score' in data:
      for s in sols:
         assert_less(abs(data['score'] - s.score), SOLUTION_PRECISION)
   
   if not shallow:
      sols_expected = [sorted_tuple((i[0], sorted_tuple(i[1]), sorted_tuple(i[2])) for i in data['sol'])]
      sols_parsed = [sorted_tuple((i.minor, 
                                 sorted_tuple((m.pos, m.op) for m in i.missing), 
                                 sorted_tuple((m.pos, m.op) for m in i.added)) 
                                 for i in s.solution)
                     for s in sols]
      assert_equal(sorted(sols_expected), sorted(sols_parsed))
   else: # As assignments can vary between multiple optimal solutions, just test minor allele assignments
      eall, emiss, enew = set(), set(), set()
      for i in data['sol']:
         eall.add(i[0])
         emiss |= set(i[1])
         enew |= set(i[2])
      assert_equal(len(sols), 1)
      pall, pmiss, pnew = set(), set(), set()
      for i in sols[0].solution:
         pall.add(i.minor)
         pmiss |= set((m.pos, m.op) for m in i.missing)
         pnew |= set((m.pos, m.op) for m in i.added) 

      print('>>', eall, emiss, enew)
      print('>>', pall, pmiss, pnew)
      assert_equal(eall, pall)
      assert_equal(emiss, pmiss)
      assert_equal(enew, pnew)
   return sols[0].score if len(sols) > 0 else 0


sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', level='DEBUG')
sh.push_application()

class MinorSyntheticTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.tests/toy.yml'))


   def test_basic(self):
      # Test two copies of *1
      assert_minor(self.gene, {
         "cn": {'1': 2},
         "data": {(115, '_'): 20},
         "major": {'1': 2},
         "sol": [('1', [], []), ('1', [], [])],
         "score": 0
      })


   def test_minor(self):
      assert_minor(self.gene, {
         "cn": {'1': 2},
         "data": {(115, '_'): 9, (115, 'SNP.TA'): 11},
         "major": {'1': 2},
         "sol": [('1', [], []), ('1B', [], [])],
         "score": 0.2
      })


   def test_miss(self):
      assert_minor(self.gene, {
         "cn": {'1': 2},
         "data": {(115, '_'): 10, (115, 'SNP.TA'): 10,
                  (148, '_'): 20,
                  (151, '_'): 10, (151, 'SNP.CT'): 10},
         "major": {'1': 1, '3': 1},
         "sol": [('1B', [], []), ('3', [(148, 'INS.A')], [])],
         "score": MISS_PENALTY_FACTOR
      })


   def test_add(self):
      assert_minor(self.gene, {
         "cn": {'1': 2},
         "data": {(115, '_'): 0, (115, 'SNP.TA'): 20,
                  (148, '_'): 20, (148, 'INS.A'): 10,
                  (151, '_'): 10, (151, 'SNP.CT'): 10},
         "major": {'1': 1, '3': 1},
         "sol": [('1B', [], []), ('3', [], [(115, 'SNP.TA')])],
         "score": ADD_PENALTY_FACTOR
      })


class MinorRealTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))
      
   
   def test_normal(self): # HG00436/v1 
      assert_minor(self.gene, {
         "cn": {'1': 3}, 
         "major": {'2': 1, '2.a': 1, '71': 1}, 
         "data": {(42522311, 'SNP.CT'): 643, (42522612, 'SNP.CG'): 634, (42522677, 'SNP.GA'): 0, (42523002, 'SNP.GA'): 861, (42523208, 'SNP.CT'): 525, (42523408, 'SNP.TG'): 851, (42523942, 'SNP.GA'): 622, (42524217, 'SNP.GT'): 0, (42524312, 'SNP.GA'): 0, (42524322, 'SNP.AG'): 0, (42525035, 'SNP.GA'): 0, (42525068, 'SNP.GA'): 0, (42525131, 'SNP.CG'): 447, (42525279, 'SNP.GA'): 0, (42525298, 'SNP.AG'): 0, (42525755, 'SNP.GA'): 0, (42525797, 'SNP.GC'): 0, (42525951, 'SNP.AC'): 729, (42526048, 'SNP.GC'): 594, (42526483, 'SNP.CA'): 708, (42526548, 'SNP.TC'): 861, (42526560, 'SNP.TG'): 886, (42526561, 'SNP.CG'): 879, (42526561, '_'): 221, (42526566, 'SNP.AG'): 891, (42526570, 'SNP.GC'): 896, (42526572, 'SNP.GT'): 911, (42526668, 'SNP.CT'): 517, (42527470, 'SNP.CT'): 834, (42527532, 'SNP.GA'): 688, (42527541, 'DEL.TC'): 0, (42528027, 'SNP.TC'): 21, (42528095, 'SNP.CT'): 0, (42528381, 'SNP.GC'): 19, }, 
         "sol": [('2MW', [], [(42527541, 'DEL.TC')]), ('2M', [], [(42527541, 'DEL.TC')]), ('71', [(42523002, 'SNP.GA')], [(42525298, 'SNP.AG')])]
      }, shallow=True)      


   def test_fusion(self): # HG01190/v1
      assert_minor(self.gene, {
         "cn": {'1': 1, '68': 1}, 
         "major": {'4': 1, '68': 1}, 
         "data": {(42522391, 'SNP.GA'): 227, (42522612, 'SNP.CG'): 308, (42523210, 'SNP.TC'): 453, (42523408, 'SNP.TG'): 882, (42524695, 'SNP.TC'): 276, (42524946, 'SNP.CT'): 319, (42525131, 'SNP.CG'): 377, (42525797, 'SNP.GC'): 84, (42525810, 'SNP.TC'): 48, (42525810, '_'): 3, (42525820, 'SNP.GT'): 30, (42525951, 'SNP.AC'): 494, (42525951, '_'): 160, (42526048, 'SNP.GC'): 521, (42526483, 'SNP.CA'): 507, (42526693, 'SNP.GA'): 749, (42527792, 'SNP.CT'): 433, (42528027, 'SNP.TC'): 31, (42528223, 'SNP.GA'): 0, }, 
         "sol": [('4AW', [], [(42526483, 'SNP.CA')]), ('68', [], [(42528223, 'SNP.GA')])]
      }, shallow=True)      


   def test_multiple(self): # NA07439/v1 : currently test only one solution
      assert_minor(self.gene, {
         "cn": {'1': 3}, 
         "major": {'4.b': 2, '41': 1}, 
         "data": {(42522391, 'SNP.GA'): 353, (42522612, 'SNP.CG'): 689, (42523210, 'SNP.TC'): 347, (42523408, 'SNP.TG'): 811, (42523804, 'SNP.CT'): 160, (42523942, 'SNP.GA'): 208, (42524695, 'SNP.TC'): 489, (42524946, 'SNP.CT'): 483, (42525131, 'SNP.CG'): 640, (42525755, 'SNP.GA'): 340, (42525797, 'SNP.GC'): 0, (42525951, 'SNP.AC'): 640, (42526048, 'SNP.GC'): 431, (42526483, 'SNP.CA'): 643, (42526483, '_'): 0, (42526548, 'SNP.TC'): 284, (42526560, 'SNP.TG'): 283, (42526560, '_'): 303, (42526561, 'SNP.CG'): 288, (42526566, 'SNP.AG'): 298, (42526570, 'SNP.GC'): 302, (42526572, 'SNP.GT'): 297, (42526693, 'SNP.GA'): 665, (42527470, 'SNP.CT'): 377, (42527532, 'SNP.GA'): 259, (42527792, 'SNP.CT'): 333, (42528027, 'SNP.TC'): 46, }, 
         "sol": [('4DW', [(42526566, 'SNP.AG')], [(42525797, 'SNP.GC')]), ('41', [(42525951, 'SNP.AC'), (42523408, 'SNP.TG'), (42526483, 'SNP.CA')], []), ('4DW', [(42527470, 'SNP.CT')], [(42525797, 'SNP.GC')])]
      }, shallow=True)


   def test_comparison(self): # NA06991/v1
      s1 = assert_minor(self.gene, {
         "cn": {'1': 2}, 
         "major": {'1': 1, '4': 1}, 
         "data": {(42522391, 'SNP.GA'): 342, (42522612, 'SNP.CG'): 428, (42522964, 'SNP.CT'): 0, (42523210, 'SNP.TC'): 459, (42523408, 'SNP.TG'): 1165, (42524217, 'SNP.GT'): 0, (42524695, 'SNP.TC'): 434, (42524814, 'SNP.GA'): 0, (42524923, 'SNP.AG'): 0, (42524946, 'SNP.CT'): 493, (42525131, 'SNP.CG'): 346, (42525755, 'SNP.GA'): 0, (42525797, 'SNP.GC'): 92, (42525810, 'SNP.TC'): 42, (42525820, 'SNP.GT'): 15, (42525951, 'SNP.AC'): 528, (42526048, 'SNP.GC'): 594, (42526483, 'SNP.CA'): 369, (42526548, 'SNP.TC'): 0, (42526560, 'SNP.TG'): 0, (42526561, 'SNP.CG'): 0, (42526566, 'SNP.AG'): 0, (42526566, '_'): 579, (42526570, 'SNP.GC'): 0, (42526572, 'SNP.GT'): 0, (42526693, 'SNP.GA'): 688, (42527792, 'SNP.CT'): 360, (42528027, 'SNP.TC'): 11, (42528223, 'SNP.GA'): 2, }, 
         "sol": [('1', [], []), ('4AW', [(42528223, 'SNP.GA')], [])]
      }, shallow=True)
      s2 = assert_minor(self.gene, {
         "cn": {'1': 2}, 
         "major": {'10': 1, '4.h': 1}, 
         "data": {(42522391, 'SNP.GA'): 342, (42522612, 'SNP.CG'): 428, (42522964, 'SNP.CT'): 0, (42523210, 'SNP.TC'): 459, (42523408, 'SNP.TG'): 1165, (42524217, 'SNP.GT'): 0, (42524695, 'SNP.TC'): 434, (42524814, 'SNP.GA'): 0, (42524923, 'SNP.AG'): 0, (42524946, 'SNP.CT'): 493, (42525131, 'SNP.CG'): 346, (42525755, 'SNP.GA'): 0, (42525797, 'SNP.GC'): 92, (42525810, 'SNP.TC'): 42, (42525820, 'SNP.GT'): 15, (42525951, 'SNP.AC'): 528, (42526048, 'SNP.GC'): 594, (42526483, 'SNP.CA'): 369, (42526548, 'SNP.TC'): 0, (42526560, 'SNP.TG'): 0, (42526561, 'SNP.CG'): 0, (42526566, 'SNP.AG'): 0, (42526566, '_'): 579, (42526570, 'SNP.GC'): 0, (42526572, 'SNP.GT'): 0, (42526693, 'SNP.GA'): 688, (42527792, 'SNP.CT'): 360, (42528027, 'SNP.TC'): 11, (42528223, 'SNP.GA'): 2, }, 
         "sol": [('10A', [], []), ('4M', [(42528223, 'SNP.GA'), (42526483, 'SNP.CA'), (42527792, 'SNP.CT')], [(42525131, 'SNP.CG')]), ('1', [], [])]
      }, shallow=True) 
      s3 = assert_minor(self.gene, {
         "cn": {'1': 2}, 
         "major": {'39': 1, '4.f': 1}, 
         "data": {(42522391, 'SNP.GA'): 342, (42522612, 'SNP.CG'): 428, (42522964, 'SNP.CT'): 0, (42523210, 'SNP.TC'): 459, (42523408, 'SNP.TG'): 1165, (42524217, 'SNP.GT'): 0, (42524695, 'SNP.TC'): 434, (42524814, 'SNP.GA'): 0, (42524923, 'SNP.AG'): 0, (42524946, 'SNP.CT'): 493, (42525131, 'SNP.CG'): 346, (42525755, 'SNP.GA'): 0, (42525797, 'SNP.GC'): 92, (42525810, 'SNP.TC'): 42, (42525820, 'SNP.GT'): 15, (42525951, 'SNP.AC'): 528, (42526048, 'SNP.GC'): 594, (42526483, 'SNP.CA'): 369, (42526548, 'SNP.TC'): 0, (42526560, 'SNP.TG'): 0, (42526561, 'SNP.CG'): 0, (42526566, 'SNP.AG'): 0, (42526566, '_'): 579, (42526570, 'SNP.GC'): 0, (42526572, 'SNP.GT'): 0, (42526693, 'SNP.GA'): 688, (42527792, 'SNP.CT'): 360, (42528027, 'SNP.TC'): 11, (42528223, 'SNP.GA'): 2, }, 
         "sol": [('1', [], []), ('39', [(42528223, 'SNP.GA')], []), ('4JW', [], [(42525131, 'SNP.CG')])]
      }, shallow=True) 
      assert_less(s1, s2)
      assert_less(s1, s3)

   