#!/usr/bin/env python
# 786

# Aldy source: test_minor.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.

from nose.tools import *

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
   
class MinorFullTest(unittest.TestCase):
   _multiprocess_can_split_ = True

   def setUp(self):
      self.gene = Gene(script_path('aldy.resources.genes/cyp2d6.yml'))


   def test_minor_cyp2d6_pgx2_HG00463():
      majors = ['10', '10', '36', '36']
      mutations = {M(42525131, 'SNP.CG'): 1706, M(42524695, 'SNP.TC'): 2705, M(42523210, 'SNP.TC'): 2513, M(42528223, 'SNP.GA'): 77, M(42525951, 'SNP.AC'): 2067, M(42526561, 'SNP.CG'): 0, M(42526548, 'SNP.TC'): 0, M(42528027, 'SNP.TC'): 332, M(42527792, 'SNP.CT'): 2540, M(42526566, 'SNP.AG'): 0, M(42522612, 'SNP.CG'): 1008, M(42526693, 'SNP.GA'): 1977, M(42525755, 'SNP.GA'): 1911, M(42526570, 'SNP.GC'): 0, M(42523408, 'SNP.TG'): 3678, M(42526572, 'SNP.GT'): 0, M(42526560, 'SNP.TG'): 0, M(42526483, 'SNP.CA'): 2661, M(42525131, 'REF'): 1, M(42524695, 'REF'): 3, M(42523210, 'REF'): 7, M(42528223, 'REF'): 1, M(42525951, 'REF'): 4, M(42526561, 'REF'): 1314, M(42526548, 'REF'): 1530, M(42528027, 'REF'): 7, M(42527792, 'REF'): 6, M(42526566, 'REF'): 1401, M(42522612, 'REF'): 4, M(42526693, 'REF'): 14, M(42525755, 'REF'): 2, M(42526570, 'REF'): 1307, M(42523408, 'REF'): 6, M(42526572, 'REF'): 1327, M(42526560, 'REF'): 1334, M(42526483, 'REF'): 2}
      solution = ('10B', '10B', '36', '36')
      assert_minor(majors, mutations, solution)


   def test_minor_cyp2d6_pgx2_HG00593():
      majors = ['10', '2', '36']
      mutations = {M(42523002, 'SNP.GA'): 452, M(42523210, 'SNP.TC'): 1003, M(42522677, 'SNP.GA'): 0, M(42527541, 'DEL.TC'): 0, M(42528223, 'SNP.GA'): 34, M(42525279, 'SNP.GA'): 0, M(42525755, 'SNP.GA'): 726, M(42528381, 'SNP.GC'): 86, M(42522311, 'SNP.CT'): 476, M(42525951, 'SNP.AC'): 1149, M(42525131, 'SNP.CG'): 1038, M(42526572, 'SNP.GT'): 567, M(42523408, 'SNP.TG'): 1960, M(42526693, 'SNP.GA'): 769, M(42524695, 'SNP.TC'): 1086, M(42525068, 'SNP.GA'): 0, M(42523208, 'SNP.CT'): 449, M(42526560, 'SNP.TG'): 598, M(42523942, 'SNP.GA'): 481, M(42528027, 'SNP.TC'): 208, M(42525035, 'SNP.GA'): 0, M(42526566, 'SNP.AG'): 582, M(42524217, 'SNP.GT'): 0, M(42527532, 'SNP.GA'): 503, M(42526483, 'SNP.CA'): 1676, M(42522612, 'SNP.CG'): 794, M(42527792, 'SNP.CT'): 988, M(42526548, 'SNP.TC'): 614, M(42527470, 'SNP.CT'): 679, M(42524322, 'SNP.AG'): 0, M(42524312, 'SNP.GA'): 0, M(42526561, 'SNP.CG'): 597, M(42526570, 'SNP.GC'): 579, M(42526048, 'SNP.GC'): 683, M(42525797, 'SNP.GC'): 0, M(42528095, 'SNP.CT'): 0, M(42523002, 'REF'): 483, M(42523210, 'REF'): 454, M(42522677, 'REF'): 729, M(42527541, 'REF'): 1226, M(42528223, 'REF'): 11, M(42525279, 'REF'): 1378, M(42525755, 'REF'): 392, M(42528381, 'REF'): 139, M(42522311, 'REF'): 404, M(42525951, 'REF'): 2, M(42525131, 'REF'): 0, M(42526572, 'REF'): 556, M(42523408, 'REF'): 2, M(42526693, 'REF'): 396, M(42524695, 'REF'): 578, M(42525068, 'REF'): 1005, M(42523208, 'REF'): 980, M(42526560, 'REF'): 543, M(42523942, 'REF'): 926, M(42528027, 'REF'): 8, M(42525035, 'REF'): 1021, M(42526566, 'REF'): 570, M(42524217, 'REF'): 1426, M(42527532, 'REF'): 728, M(42526483, 'REF'): 0, M(42522612, 'REF'): 0, M(42527792, 'REF'): 505, M(42526548, 'REF'): 641, M(42527470, 'REF'): 730, M(42524322, 'REF'): 1523, M(42524312, 'REF'): 1508, M(42526561, 'REF'): 542, M(42526570, 'REF'): 550, M(42526048, 'REF'): 529, M(42525797, 'REF'): 908, M(42528095, 'REF'): 55}
      solution = ('10BW', '2MW', '36')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['10', '10', '63']
      solution = ('10B', '10BW', '63')
      sol2 = assert_minor(majors, mutations, solution)
      
      assert_less(sol1.score, sol2.score)


   def test_minor_cyp2d6_pgx2_HG01061():
      majors = ['10', '4.h']
      mutations = {M(42522612, 'SNP.CG'): 470, M(42523210, 'SNP.TC'): 861, M(42526483, 'SNP.CA'): 693, M(42525951, 'SNP.AC'): 525, M(42524814, 'SNP.GA'): 0, M(42523408, 'SNP.TG'): 1680, M(42528027, 'SNP.TC'): 71, M(42524946, 'SNP.CT'): 466, M(42522964, 'SNP.CT'): 0, M(42525131, 'SNP.CG'): 481, M(42528223, 'SNP.GA'): 22, M(42525797, 'SNP.GC'): 276, M(42526570, 'SNP.GC'): 0, M(42526548, 'SNP.TC'): 0, M(42526560, 'SNP.TG'): 0, M(42526561, 'SNP.CG'): 0, M(42526048, 'SNP.GC'): 625, M(42522391, 'SNP.GA'): 623, M(42525820, 'SNP.GT'): 207, M(42524695, 'SNP.TC'): 723, M(42524217, 'SNP.GT'): 0, M(42525755, 'SNP.GA'): 0, M(42526566, 'SNP.AG'): 0, M(42526693, 'SNP.GA'): 475, M(42525810, 'SNP.TC'): 244, M(42526572, 'SNP.GT'): 0, M(42527792, 'SNP.CT'): 632, M(42524923, 'SNP.AG'): 0, M(42522612, 'REF'): 406, M(42523210, 'REF'): 570, M(42526483, 'REF'): 554, M(42525951, 'REF'): 531, M(42524814, 'REF'): 1418, M(42523408, 'REF'): 653, M(42528027, 'REF'): 43, M(42524946, 'REF'): 612, M(42522964, 'REF'): 1215, M(42525131, 'REF'): 422, M(42528223, 'REF'): 24, M(42525797, 'REF'): 369, M(42526570, 'REF'): 622, M(42526548, 'REF'): 740, M(42526560, 'REF'): 612, M(42526561, 'REF'): 612, M(42526048, 'REF'): 627, M(42522391, 'REF'): 616, M(42525820, 'REF'): 341, M(42524695, 'REF'): 1030, M(42524217, 'REF'): 1260, M(42525755, 'REF'): 865, M(42526566, 'REF'): 656, M(42526693, 'REF'): 446, M(42525810, 'REF'): 344, M(42526572, 'REF'): 634, M(42527792, 'REF'): 676, M(42524923, 'REF'): 976}
      solution = ('10A', '4M')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['39', '4.f']
      solution = ('39', '4JW')
      sol2 = assert_minor(majors, mutations, solution)

      majors = ['1', '4']
      solution = ('1', '4AW')
      sol3 = assert_minor(majors, mutations, solution)

      assert_less(sol3.score, min(sol1.score, sol2.score))


   def test_minor_cyp2d6_pgx2_HG01062():
      majors = ['1', '4']
      mutations = {M(42526560, 'SNP.TG'): 0, M(42525797, 'SNP.GC'): 249, M(42526570, 'SNP.GC'): 0, M(42523408, 'SNP.TG'): 1658, M(42526483, 'SNP.CA'): 678, M(42528027, 'SNP.TC'): 83, M(42526572, 'SNP.GT'): 0, M(42522612, 'SNP.CG'): 515, M(42524695, 'SNP.TC'): 690, M(42526048, 'SNP.GC'): 658, M(42526548, 'SNP.TC'): 0, M(42525131, 'SNP.CG'): 485, M(42523210, 'SNP.TC'): 875, M(42524946, 'SNP.CT'): 476, M(42526561, 'SNP.CG'): 0, M(42524814, 'SNP.GA'): 0, M(42528223, 'SNP.GA'): 16, M(42525810, 'SNP.TC'): 209, M(42522391, 'SNP.GA'): 638, M(42526566, 'SNP.AG'): 0, M(42522964, 'SNP.CT'): 0, M(42525755, 'SNP.GA'): 0, M(42526693, 'SNP.GA'): 477, M(42527792, 'SNP.CT'): 668, M(42525820, 'SNP.GT'): 184, M(42524923, 'SNP.AG'): 0, M(42524217, 'SNP.GT'): 0, M(42525951, 'SNP.AC'): 607, M(42526560, 'REF'): 649, M(42525797, 'REF'): 417, M(42526570, 'REF'): 655, M(42523408, 'REF'): 704, M(42526483, 'REF'): 686, M(42528027, 'REF'): 50, M(42526572, 'REF'): 669, M(42522612, 'REF'): 494, M(42524695, 'REF'): 1082, M(42526048, 'REF'): 624, M(42526548, 'REF'): 793, M(42525131, 'REF'): 418, M(42523210, 'REF'): 578, M(42524946, 'REF'): 670, M(42526561, 'REF'): 638, M(42524814, 'REF'): 1505, M(42528223, 'REF'): 29, M(42525810, 'REF'): 398, M(42522391, 'REF'): 682, M(42526566, 'REF'): 684, M(42522964, 'REF'): 1280, M(42525755, 'REF'): 911, M(42526693, 'REF'): 521, M(42527792, 'REF'): 750, M(42525820, 'REF'): 392, M(42524923, 'REF'): 1038, M(42524217, 'REF'): 1311, M(42525951, 'REF'): 585}
      solution = ('1', '4AW')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['39', '4.f']
      solution = ('39', '4JW')
      sol2 = assert_minor(majors, mutations, solution)

      majors = ['10', '4.h']
      solution = ('10A', '4M')
      sol3 = assert_minor(majors, mutations, solution)

      assert_less(sol1.score, min(sol3.score, sol2.score))


   def test_minor_cyp2d6_pgx2_HG01190():
      majors = ['4', '68']
      mutations = {M(42525797, 'SNP.GC'): 181, M(42524695, 'SNP.TC'): 487, M(42524946, 'SNP.CT'): 344, M(42525131, 'SNP.CG'): 331, M(42526693, 'SNP.GA'): 669, M(42525810, 'SNP.TC'): 161, M(42526483, 'SNP.CA'): 998, M(42525951, 'SNP.AC'): 524, M(42528027, 'SNP.TC'): 117, M(42527792, 'SNP.CT'): 854, M(42523408, 'SNP.TG'): 1609, M(42522391, 'SNP.GA'): 395, M(42525820, 'SNP.GT'): 139, M(42526048, 'SNP.GC'): 685, M(42523210, 'SNP.TC'): 927, M(42528223, 'SNP.GA'): 32, M(42522612, 'SNP.CG'): 359, M(42525797, 'REF'): 2, M(42524695, 'REF'): 216, M(42524946, 'REF'): 123, M(42525131, 'REF'): 4, M(42526693, 'REF'): 10, M(42525810, 'REF'): 2, M(42526483, 'REF'): 2, M(42525951, 'REF'): 159, M(42528027, 'REF'): 2, M(42527792, 'REF'): 1, M(42523408, 'REF'): 2, M(42522391, 'REF'): 2, M(42525820, 'REF'): 4, M(42526048, 'REF'): 226, M(42523210, 'REF'): 3, M(42528223, 'REF'): 0, M(42522612, 'REF'): 3}
      solution = ('4AW', '68')
      assert_minor(majors, mutations, solution)


   def test_minor_cyp2d6_pgx2_HG01192():
      majors = ['41', '5']
      mutations = {M(42526548, 'SNP.TC'): 520, M(42526561, 'SNP.CG'): 509, M(42526572, 'SNP.GT'): 484, M(42527470, 'SNP.CT'): 507, M(42523804, 'SNP.CT'): 320, M(42522612, 'SNP.CG'): 289, M(42526560, 'SNP.TG'): 514, M(42526570, 'SNP.GC'): 491, M(42528027, 'SNP.TC'): 53, M(42526566, 'SNP.AG'): 492, M(42523942, 'SNP.GA'): 383, M(42525131, 'SNP.CG'): 264, M(42527532, 'SNP.GA'): 358, M(42526548, 'REF'): 1, M(42526561, 'REF'): 1, M(42526572, 'REF'): 1, M(42527470, 'REF'): 2, M(42523804, 'REF'): 7, M(42522612, 'REF'): 0, M(42526560, 'REF'): 1, M(42526570, 'REF'): 1, M(42528027, 'REF'): 2, M(42526566, 'REF'): 1, M(42523942, 'REF'): 4, M(42525131, 'REF'): 1, M(42527532, 'REF'): 1}
      solution = ('41', '5')
      assert_minor(majors, mutations, solution)


   def test_minor_cyp2d6_pgx2_NA10854():
      majors = ['1', '4']
      mutations = {M(42525810, 'SNP.TC'): 137, M(42526548, 'SNP.TC'): 0, M(42525755, 'SNP.GA'): 0, M(42522612, 'SNP.CG'): 294, M(42523210, 'SNP.TC'): 578, M(42526483, 'SNP.CA'): 433, M(42524814, 'SNP.GA'): 0, M(42522964, 'SNP.CT'): 0, M(42526566, 'SNP.AG'): 0, M(42526561, 'SNP.CG'): 0, M(42526570, 'SNP.GC'): 0, M(42528223, 'SNP.GA'): 10, M(42524217, 'SNP.GT'): 0, M(42528027, 'SNP.TC'): 60, M(42526048, 'SNP.GC'): 424, M(42526693, 'SNP.GA'): 294, M(42524695, 'SNP.TC'): 406, M(42524923, 'SNP.AG'): 0, M(42523408, 'SNP.TG'): 1001, M(42526572, 'SNP.GT'): 0, M(42527792, 'SNP.CT'): 392, M(42525951, 'SNP.AC'): 305, M(42526560, 'SNP.TG'): 0, M(42525820, 'SNP.GT'): 118, M(42525797, 'SNP.GC'): 162, M(42524946, 'SNP.CT'): 289, M(42525131, 'SNP.CG'): 291, M(42522391, 'SNP.GA'): 368, M(42525810, 'REF'): 236, M(42526548, 'REF'): 523, M(42525755, 'REF'): 578, M(42522612, 'REF'): 251, M(42523210, 'REF'): 355, M(42526483, 'REF'): 437, M(42524814, 'REF'): 952, M(42522964, 'REF'): 770, M(42526566, 'REF'): 454, M(42526561, 'REF'): 431, M(42526570, 'REF'): 430, M(42528223, 'REF'): 14, M(42524217, 'REF'): 820, M(42528027, 'REF'): 31, M(42526048, 'REF'): 364, M(42526693, 'REF'): 305, M(42524695, 'REF'): 659, M(42524923, 'REF'): 590, M(42523408, 'REF'): 370, M(42526572, 'REF'): 441, M(42527792, 'REF'): 411, M(42525951, 'REF'): 332, M(42526560, 'REF'): 439, M(42525820, 'REF'): 219, M(42525797, 'REF'): 243, M(42524946, 'REF'): 423, M(42525131, 'REF'): 275, M(42522391, 'REF'): 383}
      solution = ('1', '4AW')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['39', '4.f']
      solution = ('39', '4JW')
      sol2 = assert_minor(majors, mutations, solution)

      majors = ['10', '4.h']
      solution = ('10A', '4M')
      sol3 = assert_minor(majors, mutations, solution)

      assert_less(sol1.score, min(sol3.score, sol2.score))


   def test_minor_cyp2d6_pgx2_NA10860():
      majors = ['4', '4.f', '61']
      mutations = {M(42526483, 'SNP.CA'): 991, M(42526572, 'SNP.GT'): 0, M(42524814, 'SNP.GA'): 0, M(42524695, 'SNP.TC'): 1064, M(42525131, 'SNP.CG'): 644, M(42526570, 'SNP.GC'): 0, M(42524946, 'SNP.CT'): 660, M(42525797, 'SNP.GC'): 388, M(42524217, 'SNP.GT'): 0, M(42524923, 'SNP.AG'): 0, M(42528223, 'SNP.GA'): 31, M(42525820, 'SNP.GT'): 281, M(42526548, 'SNP.TC'): 0, M(42522612, 'SNP.CG'): 357, M(42527792, 'SNP.CT'): 934, M(42522391, 'SNP.GA'): 488, M(42523210, 'SNP.TC'): 957, M(42525951, 'SNP.AC'): 641, M(42526560, 'SNP.TG'): 0, M(42526566, 'SNP.AG'): 0, M(42526048, 'SNP.GC'): 701, M(42525755, 'SNP.GA'): 0, M(42526561, 'SNP.CG'): 0, M(42526693, 'SNP.GA'): 723, M(42523408, 'SNP.TG'): 1782, M(42522964, 'SNP.CT'): 0, M(42528027, 'SNP.TC'): 96, M(42525810, 'SNP.TC'): 327, M(42526483, 'REF'): 485, M(42526572, 'REF'): 713, M(42524814, 'REF'): 1509, M(42524695, 'REF'): 779, M(42525131, 'REF'): 295, M(42526570, 'REF'): 714, M(42524946, 'REF'): 448, M(42525797, 'REF'): 263, M(42524217, 'REF'): 1416, M(42524923, 'REF'): 1039, M(42528223, 'REF'): 6, M(42525820, 'REF'): 259, M(42526548, 'REF'): 848, M(42522612, 'REF'): 286, M(42527792, 'REF'): 511, M(42522391, 'REF'): 481, M(42523210, 'REF'): 427, M(42525951, 'REF'): 443, M(42526560, 'REF'): 692, M(42526566, 'REF'): 735, M(42526048, 'REF'): 449, M(42525755, 'REF'): 970, M(42526561, 'REF'): 675, M(42526693, 'REF'): 365, M(42523408, 'REF'): 420, M(42522964, 'REF'): 1163, M(42528027, 'REF'): 37, M(42525810, 'REF'): 255}
      solution = ('4AW', '4JW', '61')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['4', '4.f', '83']
      solution = ('4AW', '4JW', '83')
      sol2 = assert_minor(majors, mutations, solution)

      majors = ['39', '4.f', '4.i']
      solution = ('39', '4JW', '4N')
      sol3 = assert_minor(majors, mutations, solution)

      majors = ['10', '4.h', '4.i']
      solution = ('10A', '4M', '4N')
      sol4 = assert_minor(majors, mutations, solution)

      majors = ['36', '4', '4.h']
      solution = ('36', '4B', '4M')
      sol5 = assert_minor(majors, mutations, solution)

      majors = ['1', '4', '4.i']
      solution = ('1', '4AW', '4N')
      sol6 = assert_minor(majors, mutations, solution)

      assert_less(sol6.score, min(sol1.score, sol2.score, sol3.score, sol4.score, sol5.score))


   def test_minor_cyp2d6_pgx2_NA12156():
      majors = ['1', '4']
      mutations = {M(42524923, 'SNP.AG'): 0, M(42523408, 'SNP.TG'): 1125, M(42526483, 'SNP.CA'): 461, M(42525797, 'SNP.GC'): 183, M(42528027, 'SNP.TC'): 74, M(42526560, 'SNP.TG'): 0, M(42526561, 'SNP.CG'): 0, M(42526693, 'SNP.GA'): 294, M(42524217, 'SNP.GT'): 0, M(42526048, 'SNP.GC'): 407, M(42522612, 'SNP.CG'): 314, M(42522964, 'SNP.CT'): 0, M(42525755, 'SNP.GA'): 0, M(42526548, 'SNP.TC'): 0, M(42528223, 'SNP.GA'): 4, M(42526566, 'SNP.AG'): 0, M(42523210, 'SNP.TC'): 622, M(42522391, 'SNP.GA'): 380, M(42525820, 'SNP.GT'): 134, M(42525951, 'SNP.AC'): 387, M(42526570, 'SNP.GC'): 0, M(42525810, 'SNP.TC'): 157, M(42524814, 'SNP.GA'): 0, M(42527792, 'SNP.CT'): 445, M(42524695, 'SNP.TC'): 416, M(42525131, 'SNP.CG'): 310, M(42526572, 'SNP.GT'): 0, M(42524946, 'SNP.CT'): 258, M(42524923, 'REF'): 621, M(42523408, 'REF'): 447, M(42526483, 'REF'): 478, M(42525797, 'REF'): 261, M(42528027, 'REF'): 33, M(42526560, 'REF'): 426, M(42526561, 'REF'): 422, M(42526693, 'REF'): 359, M(42524217, 'REF'): 921, M(42526048, 'REF'): 413, M(42522612, 'REF'): 313, M(42522964, 'REF'): 858, M(42525755, 'REF'): 599, M(42526548, 'REF'): 518, M(42528223, 'REF'): 12, M(42526566, 'REF'): 445, M(42523210, 'REF'): 421, M(42522391, 'REF'): 404, M(42525820, 'REF'): 236, M(42525951, 'REF'): 374, M(42526570, 'REF'): 449, M(42525810, 'REF'): 246, M(42524814, 'REF'): 944, M(42527792, 'REF'): 471, M(42524695, 'REF'): 787, M(42525131, 'REF'): 272, M(42526572, 'REF'): 453, M(42524946, 'REF'): 443}
      solution = ('1', '4AW')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['10', '4.h']
      solution = ('10A', '4M')
      sol2 = assert_minor(majors, mutations, solution)

      majors = ['39', '4.f']
      solution = ('39', '4JW')
      sol3 = assert_minor(majors, mutations, solution)

      assert_less(sol1.score, min(sol3.score, sol2.score))


   def test_minor_cyp2d6_pgx2_NA12400():
      majors = ['39', '4.f', '68']
      mutations = {M(42522964, 'SNP.CT'): 0, M(42522391, 'SNP.GA'): 553, M(42526693, 'SNP.GA'): 890, M(42526048, 'SNP.GC'): 774, M(42525755, 'SNP.GA'): 0, M(42524695, 'SNP.TC'): 638, M(42527792, 'SNP.CT'): 1184, M(42522612, 'SNP.CG'): 428, M(42526548, 'SNP.TC'): 0, M(42524814, 'SNP.GA'): 0, M(42526561, 'SNP.CG'): 0, M(42523210, 'SNP.TC'): 1143, M(42528027, 'SNP.TC'): 140, M(42526570, 'SNP.GC'): 0, M(42523408, 'SNP.TG'): 1907, M(42525810, 'SNP.TC'): 188, M(42524946, 'SNP.CT'): 415, M(42524923, 'SNP.AG'): 0, M(42524217, 'SNP.GT'): 0, M(42526566, 'SNP.AG'): 0, M(42526560, 'SNP.TG'): 0, M(42526572, 'SNP.GT'): 0, M(42525951, 'SNP.AC'): 692, M(42526483, 'SNP.CA'): 1148, M(42528223, 'SNP.GA'): 35, M(42525797, 'SNP.GC'): 232, M(42525820, 'SNP.GT'): 179, M(42525131, 'SNP.CG'): 370, M(42522964, 'REF'): 1767, M(42522391, 'REF'): 584, M(42526693, 'REF'): 401, M(42526048, 'REF'): 579, M(42525755, 'REF'): 758, M(42524695, 'REF'): 883, M(42527792, 'REF'): 573, M(42522612, 'REF'): 424, M(42526548, 'REF'): 1032, M(42524814, 'REF'): 1298, M(42526561, 'REF'): 847, M(42523210, 'REF'): 502, M(42528027, 'REF'): 44, M(42526570, 'REF'): 850, M(42523408, 'REF'): 544, M(42525810, 'REF'): 316, M(42524946, 'REF'): 563, M(42524923, 'REF'): 895, M(42524217, 'REF'): 1130, M(42526566, 'REF'): 892, M(42526560, 'REF'): 860, M(42526572, 'REF'): 861, M(42525951, 'REF'): 488, M(42526483, 'REF'): 608, M(42528223, 'REF'): 12, M(42525797, 'REF'): 335, M(42525820, 'REF'): 324, M(42525131, 'REF'): 355}
      solution = ('39', '4JW', '68')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['1', '4', '68']
      solution = ('1', '4AW', '68')
      sol2 = assert_minor(majors, mutations, solution)

      majors = ['10', '4.h', '68']
      solution = ('10A', '4M', '68')
      sol3 = assert_minor(majors, mutations, solution)

      assert_less(sol2.score, min(sol3.score, sol1.score))


   def test_minor_cyp2d6_pgx2_NA18507():
      majors = ['39', '4.b', '4.g']
      mutations = {M(42526566, 'SNP.AG'): 0, M(42524322, 'SNP.AG'): 0, M(42526560, 'SNP.TG'): 0, M(42528027, 'SNP.TC'): 202, M(42522612, 'SNP.CG'): 1113, M(42527792, 'SNP.CT'): 982, M(42524312, 'SNP.GA'): 0, M(42522311, 'SNP.CT'): 0, M(42526561, 'SNP.CG'): 0, M(42527532, 'SNP.GA'): 468, M(42526048, 'SNP.GC'): 666, M(42524695, 'SNP.TC'): 1146, M(42526570, 'SNP.GC'): 0, M(42526572, 'SNP.GT'): 0, M(42526548, 'SNP.TC'): 0, M(42527541, 'DEL.TC'): 0, M(42524946, 'SNP.CT'): 737, M(42525951, 'SNP.AC'): 1232, M(42523208, 'SNP.CT'): 522, M(42526483, 'SNP.CA'): 1698, M(42523210, 'SNP.TC'): 1079, M(42528095, 'SNP.CT'): 21, M(42523408, 'SNP.TG'): 2023, M(42525068, 'SNP.GA'): 0, M(42524217, 'SNP.GT'): 0, M(42525035, 'SNP.GA'): 0, M(42526693, 'SNP.GA'): 824, M(42522391, 'SNP.GA'): 970, M(42523942, 'SNP.GA'): 528, M(42527470, 'SNP.CT'): 0, M(42525755, 'SNP.GA'): 825, M(42522677, 'SNP.GA'): 0, M(42525279, 'SNP.GA'): 0, M(42523002, 'SNP.GA'): 533, M(42525131, 'SNP.CG'): 1170, M(42525797, 'SNP.GC'): 0, M(42526566, 'REF'): 883, M(42524322, 'REF'): 1545, M(42526560, 'REF'): 821, M(42528027, 'REF'): 4, M(42522612, 'REF'): 1, M(42527792, 'REF'): 527, M(42524312, 'REF'): 1513, M(42522311, 'REF'): 1001, M(42526561, 'REF'): 815, M(42527532, 'REF'): 684, M(42526048, 'REF'): 543, M(42524695, 'REF'): 621, M(42526570, 'REF'): 867, M(42526572, 'REF'): 870, M(42526548, 'REF'): 942, M(42527541, 'REF'): 1189, M(42524946, 'REF'): 457, M(42525951, 'REF'): 0, M(42523208, 'REF'): 1058, M(42526483, 'REF'): 0, M(42523210, 'REF'): 527, M(42528095, 'REF'): 71, M(42523408, 'REF'): 1, M(42525068, 'REF'): 1131, M(42524217, 'REF'): 1549, M(42525035, 'REF'): 1085, M(42526693, 'REF'): 419, M(42522391, 'REF'): 1001, M(42523942, 'REF'): 881, M(42527470, 'REF'): 1192, M(42525755, 'REF'): 415, M(42522677, 'REF'): 1121, M(42525279, 'REF'): 1423, M(42523002, 'REF'): 507, M(42525131, 'REF'): 0, M(42525797, 'REF'): 1020}
      solution = ('39', '4DW', '4KW')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['2.a', '4.b', '4.b']
      solution = ('2L', '4DW', '4DW')
      sol2 = assert_minor(majors, mutations, solution)

      assert_less(sol2.score, sol1.score)


   def test_minor_cyp2d6_pgx2_NA19684():
      majors = ['39', '4.f']
      mutations = {M(42524814, 'SNP.GA'): 0, M(42526561, 'SNP.CG'): 0, M(42527792, 'SNP.CT'): 431, M(42524217, 'SNP.GT'): 0, M(42526548, 'SNP.TC'): 0, M(42526560, 'SNP.TG'): 0, M(42526566, 'SNP.AG'): 0, M(42522964, 'SNP.CT'): 0, M(42526483, 'SNP.CA'): 411, M(42524695, 'SNP.TC'): 484, M(42526572, 'SNP.GT'): 0, M(42528223, 'SNP.GA'): 13, M(42526048, 'SNP.GC'): 357, M(42523408, 'SNP.TG'): 1100, M(42525810, 'SNP.TC'): 132, M(42526570, 'SNP.GC'): 0, M(42525951, 'SNP.AC'): 345, M(42524923, 'SNP.AG'): 0, M(42522612, 'SNP.CG'): 316, M(42525797, 'SNP.GC'): 158, M(42528027, 'SNP.TC'): 63, M(42522391, 'SNP.GA'): 367, M(42526693, 'SNP.GA'): 336, M(42525131, 'SNP.CG'): 317, M(42525755, 'SNP.GA'): 0, M(42523210, 'SNP.TC'): 579, M(42524946, 'SNP.CT'): 282, M(42525820, 'SNP.GT'): 108, M(42524814, 'REF'): 963, M(42526561, 'REF'): 428, M(42527792, 'REF'): 466, M(42524217, 'REF'): 822, M(42526548, 'REF'): 516, M(42526560, 'REF'): 435, M(42526566, 'REF'): 443, M(42522964, 'REF'): 813, M(42526483, 'REF'): 463, M(42524695, 'REF'): 632, M(42526572, 'REF'): 428, M(42528223, 'REF'): 6, M(42526048, 'REF'): 386, M(42523408, 'REF'): 404, M(42525810, 'REF'): 260, M(42526570, 'REF'): 420, M(42525951, 'REF'): 407, M(42524923, 'REF'): 615, M(42522612, 'REF'): 293, M(42525797, 'REF'): 283, M(42528027, 'REF'): 38, M(42522391, 'REF'): 408, M(42526693, 'REF'): 340, M(42525131, 'REF'): 240, M(42525755, 'REF'): 623, M(42523210, 'REF'): 391, M(42524946, 'REF'): 384, M(42525820, 'REF'): 264}
      solution = ('39', '4JW')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['10', '4.h']
      solution = ('10A', '4M')
      sol2 = assert_minor(majors, mutations, solution)

      majors = ['1', '4']
      solution = ('1', '4AW')
      sol3 = assert_minor(majors, mutations, solution)

      assert_less(sol3.score, min(sol1.score, sol2.score))


   def test_minor_cyp2d6_pgx2_NA19788():
      majors = ['2', '2', '78/2']
      mutations = {M(42526483, 'SNP.CA'): 1589, M(42528027, 'SNP.TC'): 166, M(42527532, 'SNP.GA'): 1118, M(42525755, 'SNP.GA'): 0, M(42522677, 'SNP.GA'): 0, M(42526048, 'SNP.GC'): 1096, M(42524131, 'SNP.CT'): 0, M(42522391, 'SNP.GA'): 0, M(42525035, 'SNP.GA'): 0, M(42523301, 'SNP.CT'): 0, M(42524695, 'SNP.TC'): 0, M(42524312, 'SNP.GA'): 0, M(42526561, 'SNP.CG'): 1387, M(42522311, 'SNP.CT'): 1760, M(42522612, 'SNP.CG'): 1287, M(42524663, 'SNP.TG'): 0, M(42523208, 'SNP.CT'): 1697, M(42526570, 'SNP.GC'): 1346, M(42524217, 'SNP.GT'): 0, M(42524322, 'SNP.AG'): 0, M(42525951, 'SNP.AC'): 1012, M(42523002, 'SNP.GA'): 1701, M(42527541, 'DEL.TC'): 0, M(42523942, 'SNP.GA'): 1620, M(42526566, 'SNP.AG'): 1373, M(42523210, 'SNP.TC'): 0, M(42523762, 'SNP.CT'): 0, M(42526548, 'SNP.TC'): 1399, M(42523408, 'SNP.TG'): 1616, M(42525068, 'SNP.GA'): 0, M(42528381, 'SNP.GC'): 149, M(42523538, 'SNP.AG'): 0, M(42526572, 'SNP.GT'): 1342, M(42525797, 'SNP.GC'): 0, M(42525131, 'SNP.CG'): 835, M(42522257, 'INS.a'): 0, M(42526560, 'SNP.TG'): 1392, M(42527470, 'SNP.CT'): 1252, M(42525279, 'SNP.GA'): 0, M(42528095, 'SNP.CT'): 0, M(42526483, 'REF'): 3, M(42528027, 'REF'): 4, M(42527532, 'REF'): 4, M(42525755, 'REF'): 867, M(42522677, 'REF'): 1196, M(42526048, 'REF'): 0, M(42524131, 'REF'): 2029, M(42522391, 'REF'): 1768, M(42525035, 'REF'): 862, M(42523301, 'REF'): 1764, M(42524695, 'REF'): 1349, M(42524312, 'REF'): 1624, M(42526561, 'REF'): 2, M(42522311, 'REF'): 3, M(42522612, 'REF'): 1, M(42524663, 'REF'): 1350, M(42523208, 'REF'): 3, M(42526570, 'REF'): 2, M(42524217, 'REF'): 1794, M(42524322, 'REF'): 1594, M(42525951, 'REF'): 7, M(42523002, 'REF'): 1, M(42527541, 'REF'): 1122, M(42523942, 'REF'): 4, M(42526566, 'REF'): 3, M(42523210, 'REF'): 1702, M(42523762, 'REF'): 1562, M(42526548, 'REF'): 2, M(42523408, 'REF'): 2, M(42525068, 'REF'): 875, M(42528381, 'REF'): 2, M(42523538, 'REF'): 1658, M(42526572, 'REF'): 2, M(42525797, 'REF'): 716, M(42525131, 'REF'): 1, M(42522257, 'REF'): 1568, M(42526560, 'REF'): 3, M(42527470, 'REF'): 2, M(42525279, 'REF'): 1058, M(42528095, 'REF'): 64}
      solution = ('2MW', '2MW', '2M')
      assert_minor(majors, mutations, solution)


   def test_minor_cyp2d6_pgx2_NA19790():
      majors = ['2', '34', '78/4']
      mutations = {M(42527470, 'SNP.CT'): 545, M(42526561, 'SNP.CG'): 506, M(42524312, 'SNP.GA'): 0, M(42524814, 'SNP.GA'): 0, M(42528381, 'SNP.GC'): 64, M(42522311, 'SNP.CT'): 864, M(42525035, 'SNP.GA'): 0, M(42524663, 'SNP.TG'): 0, M(42523942, 'SNP.GA'): 790, M(42526572, 'SNP.GT'): 490, M(42526483, 'SNP.CA'): 561, M(42523762, 'SNP.CT'): 0, M(42526560, 'SNP.TG'): 506, M(42526548, 'SNP.TC'): 510, M(42527541, 'DEL.TC'): 0, M(42525068, 'SNP.GA'): 0, M(42524923, 'SNP.AG'): 0, M(42524489, 'SNP.GA'): 0, M(42522257, 'INS.a'): 0, M(42528027, 'SNP.TC'): 60, M(42523210, 'SNP.TC'): 0, M(42525131, 'SNP.CG'): 319, M(42523408, 'SNP.TG'): 1163, M(42525951, 'SNP.AC'): 413, M(42528095, 'SNP.CT'): 0, M(42523002, 'SNP.GA'): 757, M(42527532, 'SNP.GA'): 419, M(42522677, 'SNP.GA'): 0, M(42526048, 'SNP.GC'): 393, M(42523208, 'SNP.CT'): 895, M(42525279, 'SNP.GA'): 0, M(42526566, 'SNP.AG'): 490, M(42523538, 'SNP.AG'): 0, M(42524695, 'SNP.TC'): 0, M(42522612, 'SNP.CG'): 638, M(42523301, 'SNP.CT'): 0, M(42525755, 'SNP.GA'): 0, M(42524131, 'SNP.CT'): 0, M(42523183, 'SNP.CA'): 0, M(42525797, 'SNP.GC'): 0, M(42526570, 'SNP.GC'): 495, M(42522964, 'SNP.CT'): 0, M(42522391, 'SNP.GA'): 0, M(42524217, 'SNP.GT'): 0, M(42524322, 'SNP.AG'): 0, M(42527470, 'REF'): 339, M(42526561, 'REF'): 217, M(42524312, 'REF'): 1186, M(42524814, 'REF'): 928, M(42528381, 'REF'): 69, M(42522311, 'REF'): 213, M(42525035, 'REF'): 584, M(42524663, 'REF'): 1119, M(42523942, 'REF'): 356, M(42526572, 'REF'): 238, M(42526483, 'REF'): 420, M(42523762, 'REF'): 1061, M(42526560, 'REF'): 218, M(42526548, 'REF'): 247, M(42527541, 'REF'): 767, M(42525068, 'REF'): 598, M(42524923, 'REF'): 677, M(42524489, 'REF'): 1126, M(42522257, 'REF'): 906, M(42528027, 'REF'): 40, M(42523210, 'REF'): 1300, M(42525131, 'REF'): 277, M(42523408, 'REF'): 376, M(42525951, 'REF'): 398, M(42528095, 'REF'): 47, M(42523002, 'REF'): 140, M(42527532, 'REF'): 331, M(42522677, 'REF'): 894, M(42526048, 'REF'): 414, M(42523208, 'REF'): 424, M(42525279, 'REF'): 766, M(42526566, 'REF'): 242, M(42523538, 'REF'): 1214, M(42524695, 'REF'): 1149, M(42522612, 'REF'): 338, M(42523301, 'REF'): 1250, M(42525755, 'REF'): 689, M(42524131, 'REF'): 1507, M(42523183, 'REF'): 1248, M(42525797, 'REF'): 554, M(42526570, 'REF'): 226, M(42522964, 'REF'): 1025, M(42522391, 'REF'): 1256, M(42524217, 'REF'): 1247, M(42524322, 'REF'): 1196}
      solution = ('2MW', '34', '88')
      assert_minor(majors, mutations, solution)

      majors = ['2', '39', '78/34']
      solution = ('2MW', '39', '34')
      assert_minor(majors, mutations, solution)

      majors = ['2', '2.a', '78/1']
      solution = ('2MW', '2D', '1')
      assert_minor(majors, mutations, solution)

      majors = ['1', '2', '78/2']
      solution = ('1', '2MW', '2M')
      assert_minor(majors, mutations, solution)


   def test_minor_cyp2d6_pgx2_NA19819():
      majors = ['39', '4.b', '4.g']
      mutations = {M(42523002, 'SNP.GA'): 435, M(42522612, 'SNP.CG'): 982, M(42527532, 'SNP.GA'): 379, M(42522677, 'SNP.GA'): 0, M(42524695, 'SNP.TC'): 996, M(42524217, 'SNP.GT'): 0, M(42526572, 'SNP.GT'): 478, M(42523208, 'SNP.CT'): 419, M(42525279, 'SNP.GA'): 0, M(42526548, 'SNP.TC'): 508, M(42525131, 'SNP.CG'): 878, M(42525755, 'SNP.GA'): 655, M(42525797, 'SNP.GC'): 0, M(42526693, 'SNP.GA'): 690, M(42522391, 'SNP.GA'): 793, M(42522311, 'SNP.CT'): 479, M(42525951, 'SNP.AC'): 1102, M(42527541, 'DEL.TC'): 0, M(42526483, 'SNP.CA'): 1499, M(42526048, 'SNP.GC'): 659, M(42524312, 'SNP.GA'): 0, M(42524322, 'SNP.AG'): 0, M(42523408, 'SNP.TG'): 1656, M(42524946, 'SNP.CT'): 667, M(42525035, 'SNP.GA'): 0, M(42523942, 'SNP.GA'): 391, M(42528095, 'SNP.CT'): 0, M(42528027, 'SNP.TC'): 200, M(42525068, 'SNP.GA'): 0, M(42527470, 'SNP.CT'): 586, M(42527792, 'SNP.CT'): 878, M(42526560, 'SNP.TG'): 489, M(42526561, 'SNP.CG'): 485, M(42526566, 'SNP.AG'): 476, M(42526570, 'SNP.GC'): 477, M(42523210, 'SNP.TC'): 858, M(42523002, 'REF'): 410, M(42522612, 'REF'): 0, M(42527532, 'REF'): 644, M(42522677, 'REF'): 906, M(42524695, 'REF'): 507, M(42524217, 'REF'): 1366, M(42526572, 'REF'): 489, M(42523208, 'REF'): 838, M(42525279, 'REF'): 1179, M(42526548, 'REF'): 536, M(42525131, 'REF'): 0, M(42525755, 'REF'): 350, M(42525797, 'REF'): 785, M(42526693, 'REF'): 356, M(42522391, 'REF'): 465, M(42522311, 'REF'): 398, M(42525951, 'REF'): 0, M(42527541, 'REF'): 1049, M(42526483, 'REF'): 1, M(42526048, 'REF'): 453, M(42524312, 'REF'): 1335, M(42524322, 'REF'): 1333, M(42523408, 'REF'): 1, M(42524946, 'REF'): 343, M(42525035, 'REF'): 920, M(42523942, 'REF'): 759, M(42528095, 'REF'): 59, M(42528027, 'REF'): 8, M(42525068, 'REF'): 900, M(42527470, 'REF'): 713, M(42527792, 'REF'): 479, M(42526560, 'REF'): 448, M(42526561, 'REF'): 451, M(42526566, 'REF'): 512, M(42526570, 'REF'): 487, M(42523210, 'REF'): 427}
      solution = ('39', '4DW', '4KW')
      sol1 = assert_minor(majors, mutations, solution)

      majors = ['2.a', '4.b', '4.b']
      solution = ('2M', '4DW', '4DW')
      sol2 = assert_minor(majors, mutations, solution)

      assert_less(sol2.score, sol1.score)
