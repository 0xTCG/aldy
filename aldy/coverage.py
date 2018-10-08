#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 786

# Aldy source: coverage.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, Tuple, Callable, List

from .common import *
from .gene import Mutation


class Coverage:
   """
   Describes the positional and regional coverage of a sample.
   """

   def __init__(self, 
         coverage: Dict[int, Dict[str, int]], 
         threshold: float, 
         cnv_coverage: Dict[int, int]) -> None:
      """
      Coverage initialization.

      Args:
         coverage (dict[int, dict[str, int]]):
            Coverage for each locus in the sample, where such coverage
            is given as a dictionary that maps mutation (or a reference position
            indicated by `_`) to the number of reads supporting such mutation.
            For example, ``coverage[10]['SNP.AG'] = 2`` means that there are 2
            reads that have G (instead of A) at the genomic locus 10.
         threshold (float):
            Threshold :math:`t` used for filtering out low quality mutations.
            Ranging from 0 to 1 (normalized percentage).
            Typically any mutation with the coverage :math:`<t\%` is filtered out.
         cnv_coverage (dict[int, int]):
            Coverage of the copy-number neutral region in the sample used for
            coverage rescaling. Dictionary key is the locus, while the value stands
            for the read coverage.
      """
      self._coverage = coverage
      self._threshold = threshold
      self._cnv_coverage = cnv_coverage

      self._rescaled: Dict[int, float] = {}
      self._region_coverage: Dict[Tuple[int, GeneRegion], float] = {}
      

   def __getitem__(self, mut: Mutation) -> float:
      """
      float: Returns the coverage of the mutation `mut`.
      """
      return self.coverage(mut)
   

   def coverage(self, mut: Mutation) -> float:
      """
      float: Returns the coverage of the mutation `mut`.
      """
      op = '_' if mut.op[:3] == 'REF' else mut.op
      if op in self._coverage[mut.pos]:
         return self._coverage[mut.pos][op]
      else:
         return 0


   def _dump(self):
      return {x: {a: b for a, b in c.items()} for x, c in self._coverage.items() if '_' in c and len(c)>1}


   def total(self, pos: int) -> float:
      """
      float: Returns the total coverage at locus `pos`.
      """
      if pos not in self._coverage: 
         return 0
      return float(sum(v 
         for p, v in self._coverage[pos].items() 
         if p[:3] != 'INS'))


   def percentage(self, m: Mutation) -> float:
      """
      float: Returns the coverage, expressed as the percentage in range 0-100, of the mutation `mut`.
      """
      total = self.total(m.pos)
      if total == 0: return 0
      return 100.0 * self.coverage(m) / total


   def loci_cn(self, pos: int) -> float:
      """
      float: Returns the copy number of the locus `pos`.
      """
      if self._rescaled[pos] == 0: 
         return 0
      return self.total(pos) * (1 / self._rescaled[pos])


   def region_coverage(self, gene: int, region: GeneRegion) -> float:
      """
      float: Returns the average coverage of the region `region` in the gene `gene`.
      """
      return self._region_coverage[gene, region]


   def average_coverage(self) -> float:
      """
      float: Returns the average coverage of the sample.
      """
      return sum(self.total(pos) for pos in self._coverage) / float(len(self._coverage) + 0.1)


   def filtered(self, filter_fn: Callable[[Mutation, float, float, float], bool]): # -> Coverage
      """
      :obj:`Coverage`: Returns the filtered coverage.

      Args:
         filter_fn (callable): 
            Function that performs the filtering.
            Its arguments are:
               - mut (:obj:`aldy.gene.Mutation`): mutation to be filtered
               - cov (float): coverage of the mutation
               - total (float): total locus coverage
               - thres (float): threshold value to be used for filtering
            Returns false if a mutation is to be filtered out.
      """

      cov = collections.defaultdict(lambda: collections.defaultdict(int), {
         pos: collections.defaultdict(int, {
            o: c
            for o, c in pos_mut.items()
            if filter_fn(Mutation(pos, o), c, self.total(pos), self._threshold) # type: ignore
         })
         for pos, pos_mut in self._coverage.items()
      })
      
      new_cov = Coverage(cov, self._threshold, self._cnv_coverage) # type: ignore
      new_cov._rescaled = self._rescaled
      new_cov._region_coverage = self._region_coverage
      return new_cov


   def diploid_avg_coverage(self) -> float:
      """
      float: Returns the average coverage of the copy-number neutral region.
      """
      return float(sum(self._cnv_coverage.values())) / abs(self._cn_region.end - self._cn_region.start)
      

   def _normalize_coverage(self, 
                           profile: Dict[str, Dict[int, float]], 
                           gene_regions: Dict[int, Dict[GeneRegion, GRange]], 
                           cn_region: GRange) -> None:
      """
      Normalizes the coverage of the sample to match the profile coverage.

      Args:
         profile (dict[str, dict[int, float]]):
            Profile coverage in the form `chromosome -> (position -> coverage)`.
         gene_regions (dict[int, dict[:obj:`aldy.common.GeneRegion`, :obj:`aldy.common.GRange`]]):
            List of gene regions for each gene ID.
         cn_region (:obj:`aldy.common.GRange`):
            Copy-number neutral region.
      """

      #: GRange: store the CN-neutral region
      self._cn_region: GRange = cn_region
      sam_ref = sum(self._cnv_coverage[i] for i in range(cn_region.start, cn_region.end))
      cnv_ref = sum(profile[cn_region.chr][i] for i in range(cn_region.start, cn_region.end))

      cn_ratio = float(cnv_ref) / sam_ref
      log.debug('CNV factor: {} ({})', cn_ratio, 1.0 / cn_ratio)

      self._rescaled: Dict[int, float] = {} 
      self._region_coverage: Dict[Tuple[int, GeneRegion], float] = {}
      for gene, gr in gene_regions.items():
         for region, rng in gr.items():
            s = sum(self.total(i) for i in range(rng.start, rng.end)) #!IMPORTANT
            p = sum(profile[rng.chr][i] for i in range(rng.start, rng.end))
            self._rescaled.update({
               i: profile[rng.chr][i] / cn_ratio
               for i in range(rng.start, rng.end)
            })
            self._region_coverage[gene, region] = (cn_ratio * float(s) / p) if p != 0 else 0.0


   @staticmethod
   def basic_filter(mut: Mutation, cov: float, total: float, thres: float) -> bool:
      """
      Basic filtering function.
      """
      return cov >= max(1, total * thres)

   @staticmethod
   def cn_filter(mut: Mutation, cov: float, total: float, thres: float, cn_solution) -> bool: # cn_solution: CNSolution
      """
      Filtering function that takes into the account the copy number of the mutation.
      """
      cn = cn_solution.position_cn(mut.pos)
      total = total / cn if cn > 0 else total
      return mut.op == '_' or cov > max(1, total * thres)



   
