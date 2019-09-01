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
   
   shallow=False
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


   def test_normal(self): # NA07439/v1 : currently test only one solution
      assert_minor(self.gene, {
         "cn": {'1': 3}, 
         "major": {'4.b': 2, '41': 1}, 
         "data": {(42522391, 'SNP.GA'): 1.9981, (42522391, '_'): 1.0019, 
                  (42522612, 'SNP.CG'): 3.0000, (42522612, '_'): 0.0000, 
                  (42523210, 'SNP.TC'): 1.9101, (42523210, '_'): 1.0899, 
                  (42523408, 'SNP.TG'): 3.0000, (42523408, '_'): 0.0000, 
                  (42523804, 'SNP.CT'): 0.8696, (42523804, '_'): 2.1304, 
                  (42523942, 'SNP.GA'): 1.0365, (42523942, '_'): 1.9635, 
                  (42524695, 'SNP.TC'): 1.9905, (42524695, '_'): 1.0095, 
                  (42524946, 'SNP.CT'): 1.8435, (42524946, '_'): 1.1565, 
                  (42525131, 'SNP.CG'): 3.0000, (42525131, '_'): 0.0000, 
                  (42525755, 'SNP.GA'): 1.9844, (42525755, '_'): 1.0156, 
                  (42525797, 'SNP.GC'): 0.0000, (42525797, '_'): 3.0000, 
                  (42525951, 'SNP.AC'): 3.0000, (42525951, '_'): 0.0000, 
                  (42526048, 'SNP.GC'): 1.8008, (42526048, '_'): 1.1992, 
                  (42526483, 'SNP.CA'): 3.0000, (42526483, '_'): 0.0000, 
                  (42526548, 'SNP.TC'): 1.3589, (42526548, '_'): 1.6411, 
                  (42526560, 'SNP.TG'): 1.4488, (42526560, '_'): 1.5512, 
                  (42526561, 'SNP.CG'): 1.4669, (42526561, '_'): 1.5331, 
                  (42526566, 'SNP.AG'): 1.4419, (42526566, '_'): 1.5581, 
                  (42526570, 'SNP.GC'): 1.4732, (42526570, '_'): 1.5268, 
                  (42526572, 'SNP.GT'): 1.4348, (42526572, '_'): 1.5652, 
                  (42526693, 'SNP.GA'): 2.0050, (42526693, '_'): 0.9950, 
                  (42527470, 'SNP.CT'): 1.5536, (42527470, '_'): 1.4464, 
                  (42527532, 'SNP.GA'): 1.2738, (42527532, '_'): 1.7262, 
                  (42527792, 'SNP.CT'): 1.8743, (42527792, '_'): 1.1257, 
                  (42528027, 'SNP.TC'): 3.0000, (42528027, '_'): 0.0000}, 
         "sol": [('4DW', [], [(42525797, 'SNP.GC')]), 
                 ('41',  [(42525951, 'SNP.AC'), (42523408, 'SNP.TG'), (42526483, 'SNP.CA')], []), 
                 ('4DW', [], [(42525797, 'SNP.GC')])]
      }, shallow=True)


   def test_multiple(self): # HG00436/v1 
      assert_minor(self.gene, {
         "cn": {'1': 3}, 
         "major": {'2': 1, '2.a': 1, '71': 1}, 
         "data": {(42522311, 'SNP.CT'): 2.4542, (42522311, '_'): 0.5458, 
                  (42522612, 'SNP.CG'): 1.9193, (42522612, '_'): 1.0807, 
                  (42522677, 'SNP.GA'): 0.0000, (42522677, '_'): 3.0000, 
                  (42523002, 'SNP.GA'): 2.6197, (42523002, '_'): 0.3803, 
                  (42523208, 'SNP.CT'): 1.9737, (42523208, '_'): 1.0263, 
                  (42523408, 'SNP.TG'): 2.1691, (42523408, '_'): 0.8309, 
                  (42523942, 'SNP.GA'): 1.9560, (42523942, '_'): 1.0440, 
                  (42524217, 'SNP.GT'): 0.0000, (42524217, '_'): 3.0000, 
                  (42524312, 'SNP.GA'): 0.0000, (42524312, '_'): 3.0000, 
                  (42524322, 'SNP.AG'): 0.0000, (42524322, '_'): 3.0000, 
                  (42525035, 'SNP.GA'): 0.0000, (42525035, '_'): 3.0000, 
                  (42525068, 'SNP.GA'): 0.0000, (42525068, '_'): 3.0000, 
                  (42525131, 'SNP.CG'): 2.1388, (42525131, '_'): 0.8612, 
                  (42525279, 'SNP.GA'): 0.0000, (42525279, '_'): 3.0000, 
                  (42525298, 'SNP.AG'): 0.0000, (42525298, '_'): 3.0000, 
                  (42525755, 'SNP.GA'): 0.0000, (42525755, '_'): 3.0000, 
                  (42525797, 'SNP.GC'): 0.0000, (42525797, '_'): 3.0000, 
                  (42525951, 'SNP.AC'): 2.3694, (42525951, '_'): 0.6306, 
                  (42526048, 'SNP.GC'): 1.8447, (42526048, '_'): 1.1553, 
                  (42526483, 'SNP.CA'): 2.1261, (42526483, '_'): 0.8739, 
                  (42526548, 'SNP.TC'): 2.3525, (42526548, '_'): 0.6475, 
                  (42526560, 'SNP.TG'): 2.4011, (42526560, '_'): 0.5989, 
                  (42526561, 'SNP.CG'): 2.3973, (42526561, '_'): 0.6027, 
                  (42526566, 'SNP.AG'): 2.3845, (42526566, '_'): 0.6155, 
                  (42526570, 'SNP.GC'): 2.3872, (42526570, '_'): 0.6128, 
                  (42526572, 'SNP.GT'): 2.3683, (42526572, '_'): 0.6317, 
                  (42526668, 'SNP.CT'): 1.0340, (42526668, '_'): 1.9660, 
                  (42527470, 'SNP.CT'): 2.2602, (42527470, '_'): 0.7398, 
                  (42527532, 'SNP.GA'): 2.1841, (42527532, '_'): 0.8159, 
                  (42527541, 'DEL.TC'): 0.0000, (42527541, '_'): 3.0000, 
                  (42528027, 'SNP.TC'): 2.5200, (42528027, '_'): 0.4800, 
                  (42528095, 'SNP.CT'): 0.0000, (42528095, '_'): 3.0000, 
                  (42528381, 'SNP.GC'): 1.6286, (42528381, '_'): 1.3714}, 
         "sol": [('2MW', [], [(42527541, 'DEL.TC')]), 
                 ('2M',  [], [(42527541, 'DEL.TC')]),
                 ('71',  [], [(42525298, 'SNP.AG')])]
      }, shallow=True)      


   def test_fusion(self): # HG01190/v1
      assert_minor(self.gene, {
         "cn": {'1': 1, '68': 1}, 
         "major": {'4': 1, '68': 1}, 
         "data": {(42522391, 'SNP.GA'): 0.9619, (42522391, '_'): 0.0381, 
                  (42522612, 'SNP.CG'): 1.0000, (42522612, '_'): 0.0000, 
                  (42523210, 'SNP.TC'): 1.0000, (42523210, '_'): 0.0000, 
                  (42523408, 'SNP.TG'): 1.0000, (42523408, '_'): 0.0000, 
                  (42524695, 'SNP.TC'): 0.6699, (42524695, '_'): 0.3301, 
                  (42524946, 'SNP.CT'): 0.6787, (42524946, '_'): 0.3213, 
                  (42525131, 'SNP.CG'): 1.0000, (42525131, '_'): 0.0000, 
                  (42525797, 'SNP.GC'): 0.9438, (42525797, '_'): 0.0562, 
                  (42525810, 'SNP.TC'): 0.9412, (42525810, '_'): 0.0588, 
                  (42525820, 'SNP.GT'): 0.9091, (42525820, '_'): 0.0909, 
                  (42525951, 'SNP.AC'): 0.7554, (42525951, '_'): 0.2446, 
                  (42526048, 'SNP.GC'): 0.7464, (42526048, '_'): 0.2536, 
                  (42526483, 'SNP.CA'): 1.0000, (42526483, '_'): 0.0000, 
                  (42526693, 'SNP.GA'): 2.0000, (42526693, '_'): 0.0000, 
                  (42527792, 'SNP.CT'): 2.0000, (42527792, '_'): 0.0000, 
                  (42528027, 'SNP.TC'): 2.0000, (42528027, '_'): 0.0000, 
                  (42528223, 'SNP.GA'): 0.0000, (42528223, '_'): 0.0000}, 
         "sol": [('4AW', [], [(42526483, 'SNP.CA')]), 
                 ('68',  [], [(42528223, 'SNP.GA')])]
      }, shallow=True)      


   def test_deletion(self): # HG00276/v1
      assert_minor(self.gene, {
         "cn": {'1': 1, '5': 1}, 
         "major": {'4': 1, '5': 1}, 
         "data": {(42522391, 'SNP.GA'): 0.9680, (42522391, '_'): 0.0320, 
                  (42522612, 'SNP.CG'): 1.0000, (42522612, '_'): 0.0000, 
                  (42523210, 'SNP.TC'): 0.9688, (42523210, '_'): 0.0312, 
                  (42523408, 'SNP.TG'): 1.0000, (42523408, '_'): 0.0000, 
                  (42524695, 'SNP.TC'): 0.7650, (42524695, '_'): 0.2350, 
                  (42524946, 'SNP.CT'): 0.7719, (42524946, '_'): 0.2281, 
                  (42525131, 'SNP.CG'): 1.0000, (42525131, '_'): 0.0000, 
                  (42525797, 'SNP.GC'): 0.9550, (42525797, '_'): 0.0450, 
                  (42525810, 'SNP.TC'): 0.9605, (42525810, '_'): 0.0395, 
                  (42525820, 'SNP.GT'): 1.0000, (42525820, '_'): 0.0000, 
                  (42525951, 'SNP.AC'): 0.6725, (42525951, '_'): 0.3275, 
                  (42526048, 'SNP.GC'): 0.6659, (42526048, '_'): 0.3341, 
                  (42526483, 'SNP.CA'): 1.0000, (42526483, '_'): 0.0000, 
                  (42526693, 'SNP.GA'): 0.9424, (42526693, '_'): 0.0576, 
                  (42527792, 'SNP.CT'): 0.9738, (42527792, '_'): 0.0262, 
                  (42528027, 'SNP.TC'): 1.0000, (42528027, '_'): 0.0000, }, 
         "sol": []
      }, shallow=True)


   def test_comparison(self): # NA06991/v1
      data =  {(42522391, 'SNP.GA'): 0.9108, (42522391, '_'): 1.0892, 
               (42522612, 'SNP.CG'): 0.9438, (42522612, '_'): 1.0562, 
               (42522964, 'SNP.CT'): 0.0000, (42522964, '_'): 2.0000, 
               (42523210, 'SNP.TC'): 1.1475, (42523210, '_'): 0.8525, 
               (42523408, 'SNP.TG'): 1.4831, (42523408, '_'): 0.5169, 
               (42524217, 'SNP.GT'): 0.0000, (42524217, '_'): 2.0000, 
               (42524695, 'SNP.TC'): 0.8535, (42524695, '_'): 1.1465, 
               (42524814, 'SNP.GA'): 0.0000, (42524814, '_'): 2.0000, 
               (42524923, 'SNP.AG'): 0.0000, (42524923, '_'): 2.0000, 
               (42524946, 'SNP.CT'): 0.7977, (42524946, '_'): 1.2023, 
               (42525131, 'SNP.CG'): 1.2791, (42525131, '_'): 0.7209, 
               (42525755, 'SNP.GA'): 0.0000, (42525755, '_'): 2.0000, 
               (42525797, 'SNP.GC'): 1.0110, (42525797, '_'): 0.9890, 
               (42525810, 'SNP.TC'): 0.9130, (42525810, '_'): 1.0870, 
               (42525820, 'SNP.GT'): 0.5455, (42525820, '_'): 1.4545, 
               (42525951, 'SNP.AC'): 0.9897, (42525951, '_'): 1.0103, 
               (42526048, 'SNP.GC'): 1.0215, (42526048, '_'): 0.9785, 
               (42526483, 'SNP.CA'): 0.9498, (42526483, '_'): 1.0502, 
               (42526548, 'SNP.TC'): 0.0000, (42526548, '_'): 2.0000, 
               (42526560, 'SNP.TG'): 0.0000, (42526560, '_'): 2.0000, 
               (42526561, 'SNP.CG'): 0.0000, (42526561, '_'): 2.0000, 
               (42526566, 'SNP.AG'): 0.0000, (42526566, '_'): 2.0000, 
               (42526570, 'SNP.GC'): 0.0000, (42526570, '_'): 2.0000, 
               (42526572, 'SNP.GT'): 0.0000, (42526572, '_'): 2.0000, 
               (42526693, 'SNP.GA'): 1.0700, (42526693, '_'): 0.9300, 
               (42527792, 'SNP.CT'): 1.0942, (42527792, '_'): 0.9058, 
               (42528027, 'SNP.TC'): 1.2222, (42528027, '_'): 0.7778, 
               (42528223, 'SNP.GA'): 1.3333, (42528223, '_'): 0.6667} 
      s1 = assert_minor(self.gene, {
         "cn": {'1': 2}, 
         "major": {'1': 1, '4': 1}, 
         "data": data,
         "sol": [('1',   [], []), 
                 ('4AW', [(42528223, 'SNP.GA')], [])]
      }, shallow=True)
      s2 = assert_minor(self.gene, {
         "cn": {'1': 2}, 
         "major": {'10': 1, '4.h': 1}, 
         "data": data,
         "sol": [('10A', [(42526483, 'SNP.CA')], []), 
                 ('4M',  [(42527792, 'SNP.CT'), (42528223, 'SNP.GA')], [(42525131, 'SNP.CG')])]
      }, shallow=True) 
      s3 = assert_minor(self.gene, {
         "cn": {'1': 2}, 
         "major": {'39': 1, '4.f': 1}, 
         "data": data, 
         "sol": [('39',  [], []), 
                 ('4JW', [(42528223, 'SNP.GA')], [(42525131, 'SNP.CG')])]
      }, shallow=True) 
      assert_less(s1, s2)
      assert_less(s1, s3)

   