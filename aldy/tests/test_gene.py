#!/usr/bin/env python
# 786

# Aldy source: test_gene.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from nose.tools import *
import unittest
from unittest.mock import Mock, patch

from aldy.gene import *
from aldy.common import *


class GeneTest(unittest.TestCase):
   def setUp(self):
      self.gene = Gene(script_path('aldy.tests/toy.yml'))


   def test_gene_basic(self):
      assert_equal(self.gene.name, 'GENE')
      assert_equal(self.gene.seq, 'ACGT' * 50)
      assert_equal(self.gene.region, GRange('1', 0, 200))


   def test_gene_regions(self):
      # Test gene regions
      assert_dict_equal(self.gene.regions, {
         0: { GeneRegion(0, "tmp"): GRange('1', 100, 110),
              GeneRegion(1, "e"):   GRange('1', 110, 120),
              GeneRegion(1, "i"):   GRange('1', 120, 130),
              GeneRegion(2, "e"):   GRange('1', 130, 140),
              GeneRegion(2, "i"):   GRange('1', 140, 150),
              GeneRegion(3, "e"):   GRange('1', 150, 160) },
         1: { GeneRegion(0, "tmp"): GRange('1',   0,  10),
              GeneRegion(1, "e"):   GRange('1',  10,  20),
              GeneRegion(1, "i"):   GRange('1',  20,  30),
              GeneRegion(2, "e"):   GRange('1',  30,  40),
              GeneRegion(2, "i"):   GRange('1',  40,  50),
              GeneRegion(3, "e"):   GRange('1',  50,  60) },
      })
      assert_equal(self.gene.pseudogenes, ['PSEUDO'])

      # Test gene auxiliary functions
      assert_equal(self.gene.region_at(105), (0, GeneRegion(0, "tmp")))
      assert_equal(self.gene.region_at(127), (0, GeneRegion(1, "i")))
      assert_equal(self.gene.region_at(155), (0, GeneRegion(3, "e")))
      assert_raises(KeyError, lambda: self.gene.region_at(165))
      assert_equal(self.gene.region_at(5),  (1, GeneRegion(0, "tmp")))
      assert_equal(self.gene.region_at(27), (1, GeneRegion(1, "i")))
      assert_equal(self.gene.region_at(55), (1, GeneRegion(3, "e")))
      assert_raises(KeyError, lambda: self.gene.region_at(80))

      assert_equal(self.gene.unique_regions, [GeneRegion(1, 'e'),
                                              GeneRegion(1, 'i'),
                                              GeneRegion(2, 'e'),
                                              GeneRegion(2, 'i'),
                                              GeneRegion(3, 'e')])


   def test_gene_alleles(self):
      # Test metadata loading
      assert_equal(self.gene.common_tandems, [('1', '4')])
      assert_equal(self.gene.old_notation(Mutation(115, 'SNP.TA')), '3828:T>A')
      assert_equal(self.gene.old_notation(Mutation(119, 'INS.TT')), '')

      def cnify(s):
         s = s.replace(' ', '').replace('|', '').split(',')
         assert_equal(len(s), len(self.gene.regions[0]) + len(self.gene.regions[1]))
         i = 0
         cn = {}
         for g in range(2):
            cn[g] = {}
            for r in sorted(self.gene.regions[g]):
               cn[g][r] = float(s[i])
               i += 1
         return cn

      # Test copy number configurations
      assert_dict_equal(self.gene.cn_configs, {
         '1': CNConfig(cnify('1,  1,1,  1,1,  1, | 1,  1,1,  1,1,  1'), 
            kind=CNConfig.CNConfigType.DEFAULT_CN, 
            alleles=set(['1', '1.a', '2', '3']), 
            description='Normal allelic configuration: all regions present in both gene and pseudogene'),
         '4': CNConfig(cnify('0,  0,0,  0,0.5,  1, | 1,  1,1,  1,0.5,  0'), 
            kind=CNConfig.CNConfigType.LEFT_FUSION, 
            alleles=set(['4/1', '4/3']), 
            description='Fusion: pseudogene until 2i followed by the gene'),
         '5': CNConfig(cnify('1,  1,1,  0,0,  0, | 1,  1,1,  2,2,  2'), 
            kind=CNConfig.CNConfigType.RIGHT_FUSION, 
            alleles=set(['5']), 
            description='Conservation: Pseudogene retention after 2e within the gene'),
         '6': CNConfig(cnify('0,  0,0,  0,0,  0, | 1,  1,1,  1,1,  1'),
            kind=CNConfig.CNConfigType.DELETION, 
            alleles=set(['6']), 
            description='Gene deletion')})

      # Test allele configurations
      assert_equal(sorted(self.gene.alleles.keys()), ['1', '1.a', '2', '3', '4/1', '4/3', '5', '6'])
      assert_equal(self.gene.alleles['1'], 
                   Allele('1', cn_config='1', 
                          func_muts=set(), 
                          minors={'1': Suballele('1', alt_names=[], neutral_muts=set()),
                                  '1B': Suballele('1B', alt_names=[], neutral_muts=set([Mutation(115, 'SNP.TA', 0, aux={'old': '3828:T>A', 'dbsnp': 'rs28371732 or rs28371741'})]))}))
      assert_equal(self.gene.alleles['1.a'], 
                   Allele('1.a', cn_config='1', 
                          func_muts=set([Mutation(105, 'SNP.TA', 1)]), 
                          minors={'1C': Suballele('1C', alt_names=[], neutral_muts=set())}))
      assert_equal(self.gene.alleles['2'], 
                   Allele('2', cn_config='1', 
                          func_muts=set([Mutation(111, 'DEL.AC', 3), Mutation(119, 'INS.TT', 3)]), 
                          minors={'2': Suballele('2', alt_names=[], neutral_muts=set())}))
      assert_equal(self.gene.alleles['3'], 
                   Allele('3', cn_config='1', 
                          func_muts=set([Mutation(151, 'SNP.CT', 1)]),
                          minors={'3': Suballele('3', alt_names=[], neutral_muts=set([Mutation(148, 'INS.A', 0)]))}))
      assert_equal(self.gene.alleles['4/1'], 
                   Allele('4/1', cn_config='4', 
                          func_muts=set(),
                          minors={'1': Suballele('1', alt_names=['1B', '1C', '2'], neutral_muts=set())}))
      assert_equal(self.gene.alleles['4/3'], 
                   Allele('4/3', cn_config='4', 
                          func_muts=set([Mutation(151, 'SNP.CT', 1)]),
                          minors={'3': Suballele('3', alt_names=[], neutral_muts=set([Mutation(148, 'INS.A', 0)]))}))
      assert_equal(self.gene.alleles['5'], 
                   Allele('5', cn_config='5', 
                          func_muts=set([Mutation(111, 'DEL.AC', 3)]), 
                          minors={'5': Suballele('5', alt_names=[], neutral_muts=set())}))
      assert_equal(self.gene.alleles['6'], 
                   Allele('6', cn_config='6', 
                          func_muts=set(), 
                          minors={'6DEL': Suballele('6DEL', alt_names=[], neutral_muts=set())}))


   def test_structure(self):
      # Test the inferred gene structure
      assert_dict_equal(self.gene.mutations, {
         (105, 'SNP.TA'): Mutation(105, 'SNP.TA', 1),
         (111, 'DEL.AC'): Mutation(111, 'DEL.AC', 3),
         (115, 'SNP.TA'): Mutation(115, 'SNP.TA', 0, aux={'dbsnp': 'rs28371732 or rs28371741', 'old': '3828:T>A'}),
         (119, 'INS.TT'): Mutation(119, 'INS.TT', 3),
         (148, 'INS.A'):  Mutation(148, 'INS.A', 0),
         (151, 'SNP.CT'): Mutation(151, 'SNP.CT', 1)})

      coding_seq = ''.join(c for _, c in sorted(self.gene.coding_region.lookup.items()))
      assert_equal(coding_seq, 'GTACGTACGTGTACGTACGTGTACGTACGT')
      assert_equal(self.gene.coding_region.aminoacid, 'VRTCTYVYVR')

