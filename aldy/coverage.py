#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 786

# Aldy source: coverage.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, Tuple, Callable, List

from .common import *
from .gene import Mutation


#: dict of aldy.gene.GeneRegion: float: resclaed (average) CN of the whole region
#: dict of int: float: Rescaled coverage for a loci

# rescaled: Dict[int: float], 
         # region_cov: Dict[Tuple[int, GeneRegion], float]

class Coverage:
   """

   """

   def __init__(self, 
         coverage: Dict[int, Dict[str, float]], 
         threshold: float, 
         cnv_coverage: Dict[int, float]):
      self._coverage = coverage
      self.threshold = threshold
      self._cnv_coverage = cnv_coverage
      

   def __getitem__(self, mut: Mutation) -> float:
      return self.coverage(mut)
   

   def coverage(self, mut: Mutation) -> float:
      op = '_' if mut.op[:3] == 'REF' else mut.op
      if op in self._coverage[mut.pos]:
         return self._coverage[mut.pos][op]
      else:
         return 0


   def total(self, pos: int) -> float:
      if pos not in self._coverage: 
         return 0
      return float(sum(v 
         for p, v in self._coverage[pos].items() 
         if p[:3] != 'INS'))


   def percentage(self, m: Mutation) -> float:
      total = self.total(m.pos)
      if total == 0: return 0
      return 100.0 * self.coverage(m) / total


   def loci_cn(self, pos: int) -> float:
      if self._rescaled[pos] == 0: 
         return 0
      return self.total(pos) * (1 / self._rescaled[pos])


   def region_coverage(self, gene: str, region: GeneRegion) -> float:
      return self._region_coverage[gene, region]


   def average_coverage(self) -> float:
      return sum(self.total(pos) for pos in self._coverage) / float(len(self._coverage) + 0.1)


   def filtered(self, filter_fn: Callable[[Mutation, float, float, float], bool]): # -> Coverage
      """Return filtered coverage"""

      cov = collections.defaultdict(lambda: collections.defaultdict(int), {
         pos: collections.defaultdict(int, {
            o: c
            for o, c in pos_mut.items()
            if filter_fn(mut=Mutation(pos, o), cov=c, total=self.total(pos), thres=self.threshold)
         })
         for pos, pos_mut in self._coverage.items()
      })
      
      new_cov = Coverage(cov, self.threshold, self._cnv_coverage)
      new_cov._rescaled = self._rescaled
      new_cov._region_coverage = self._region_coverage
      return new_cov


   def diploid_avg_coverage(self) -> float:
      return float(sum(self._cnv_coverage.values())) / abs(self._cn_region.end - self._cn_region.start)
      

   def _normalize_coverage(self, 
      profile: Dict[str, Dict[int, float]], 
      gene_regions: List[Tuple[int, GeneRegion]], 
      cn_region: GeneRegion) -> None:
      """
      """

      #: GeneRegion: store the CN-neutral region
      self._cn_region = cn_region
      sam_ref = sum(self._cnv_coverage[i] for i in range(cn_region.start, cn_region.end))
      cnv_ref = sum(profile[cn_region.chr][i] for i in range(cn_region.start, cn_region.end))

      cn_ratio = float(cnv_ref) / sam_ref
      log.debug('CNV factor: {} ({})', cn_ratio, 1.0 / cn_ratio)

      self._rescaled = {} 
      self._region_coverage = {}
      for gene, gr in gene_regions:
         for region, rng in gr.items():
            s = sum(self.total(i) for i in range(rng.start, rng.end)) #!IMPORTANT
            p = sum(profile[rng.chr][i] for i in range(rng.start, rng.end))
            self._rescaled.update({
               i: profile[rng.chr][i] / cn_ratio
               for i in range(rng.start, rng.end)
            })
            self._region_coverage[(gene, region)] = (cn_ratio * float(s) / p) if p != 0 else 0.0


   @staticmethod
   def basic_filter(mut: Mutation, cov: float, total: float, thres: float) -> bool:
      return cov >= max(1, total * thres)

   @staticmethod
   def cn_filter(mut: Mutation, cov: float, total: float, thres: float, cn_solution) -> bool: # cn_solution: CNSolution
      cn = cn_solution.position_cn(mut.pos)
      total = total / cn if cn > 0 else total
      return mut.op == '_' or cov > max(1, total * thres)



   
