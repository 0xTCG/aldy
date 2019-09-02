#!/usr/bin/env python
# 786

# Aldy source: test_diplotype.py
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
import aldy.diplotype
from aldy.gene import Gene, Mutation
from aldy.major import MajorSolution, SolvedAllele
from aldy.minor import ADD_PENALTY_FACTOR, MISS_PENALTY_FACTOR
from aldy.common import *


def assert_diplotype(gene, test, majors):
   sols = []
   for ma in majors:
      if isinstance(ma, tuple):
         sols.append(SolvedAllele(ma[0], None, [Mutation(m[0], m[1], True) for m in ma[1]], []))
      else:
         sols.append(SolvedAllele(ma, None, [], []))
   minor = aldy.minor.MinorSolution(0, sols, aldy.major.MajorSolution(0, sols, None))
   res = aldy.diplotype.estimate_diplotype(gene, minor)
   assert_equal(test, minor.diplotype)
   assert_equal(test, res)


sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', level='DEBUG')
# sh.push_application()

class DiplotypeSyntheticTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.tests/toy.yml'))


   def test_basic(self):
      assert_diplotype(self.gene, '*1/*1', 
         ['1', '1.a'])


   def test_tandem(self):
      assert_diplotype(self.gene, '*1+*4/*3', 
         ['4', '1.a', '3'])


   def test_multi(self):
      assert_diplotype(self.gene, '*1+*4/*3+*3', 
         ['4', '1', '3', '3'])
      assert_diplotype(self.gene, '*1+*4/*3+*3+*3+*3', 
         ['4', '1', '3', '3', '3', '3'])


   def test_combination(self):
      assert_diplotype(self.gene, '*1+*4/*1+*4+*3', 
         ['4', '1.a', '1', '3', '4'])


   def test_novel(self):
      assert_diplotype(self.gene, '*1+*4+*3/*1-like+*4', 
         ['4', ('1', [(151, 'SNP.CT')]), '1', '3', '4'])


   def test_deletion(self):
      assert_diplotype(self.gene, '*2/*6', 
         ['6', '2'])


class DiplotypeRealTest(unittest.TestCase):
   _multiprocess_can_split_ = True


   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))


   def test_basic(self):
      assert_diplotype(self.gene, '*1/*1', 
         ['1', '1.a'])


   def test_tandem(self):
      assert_diplotype(self.gene, '*2+*2/*71', 
         ['2', '2.a', '71'])
      assert_diplotype(self.gene, '*4+*4/*41', 
         ['4.b', '41', '4.b'])
      assert_diplotype(self.gene, '*4/*4+*4', 
         ['4.b', '4', '4.b'])
      assert_diplotype(self.gene, '*3/*68+*4', 
         ['3', '4', '68'])


   def test_36(self):
      assert_diplotype(self.gene, '*36+*10/*36+*41', 
         ['10', '36', '41', '36'])
      assert_diplotype(self.gene, '*1+*36/*36+*10', 
         ['1', '10', '36', '36'])
      assert_diplotype(self.gene, '*10/*36+*10', 
         ['10', '36', '10'])
      assert_diplotype(self.gene, '*36+*10/*36+*10', 
         ['10', '36', '36', '10'])
      assert_diplotype(self.gene, '*10+*10/*36+*10', 
         ['10', '36', '10', '10'])


   def test_fusion(self):
      assert_diplotype(self.gene, '*2+*2/*68+*4', 
         ['2', '4', '68', '2'])
      assert_diplotype(self.gene, '*1/*79+*2', 
         ['1', '2', '79/2'])
      assert_diplotype(self.gene, '*2/*78+*2', 
         ['2', '78/2', '2'])
      assert_diplotype(self.gene, '*1/*78+*2', 
         ['2', '1', '78/2'])


