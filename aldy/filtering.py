from __future__ import division
# 786

# Aldy source: filtering.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import copy

from .common import *
from .gene import *


def basic_filter(sam):
   for pos in sam.coverage:
      s = sam.total(pos)
      sam.unfiltered_coverage[pos] = s

      indels = [i for i in sam[pos] if i[:3] == 'INS']
      # Make sure not to filter out large indels
      for i in indels:
         if Mutation(pos, i) in sam.indel_sites:
            current_ratio = float(sam[pos][i]) / s
            j = sam.indel_sites[Mutation(pos, i)]

            # calculate j
            j0 = 0
            for orig_start, start, read_len in j[0]:
               ins_len = len(i) - 4
               # This is rather hackish way of detecting *40
               # min 20% sa strane, max ... how much?
               min_margin = .20 * max(1, ins_len / 10.0) * read_len
               if start < pos + max(3 * ins_len, min_margin) or \
                  orig_start > pos - max(2 * ins_len, min_margin):
                  continue
               j0 += j[0][(orig_start, start, read_len)]
            j[0] = j0

            new_ratio = j[1] / float(j[0] + j[1])
            log.debug(
               'Refining indel {}:{}: from {} of {} to {} (indel:total = {}:{})',
               pos, i,
               sam[pos][i], s,
               int(sam[pos][i] * (new_ratio / current_ratio)),
               j[1], j[0] + j[1]
            )
            sam[pos][i] = int(sam[pos][i] * (new_ratio / current_ratio))
            sam.original_coverage[pos][i] = sam[pos][i]

      threshold = sam.threshold / 20.0 # MAX copy number
      sam.coverage[pos] = {
         k: v
         for k, v in sam.coverage[pos].items()
         if v >= s * threshold and v > 1
      }


def cnv_filter(sam, gene, copy_number):
   original_coverage = copy.deepcopy(sam.original_coverage)
   sam.coverage = dict()
   for pos in original_coverage:
      s = sam.unfiltered_coverage[pos]
      if gene.region_at[pos] in copy_number and copy_number[gene.region_at[pos]] > 0: # we don't care about other mutations
         s /= copy_number[gene.region_at[pos]]
      sam.coverage[pos] = {
         k: v
         for k, v in original_coverage[pos].items()
         if k == '_' or (v > 1 and v >= s * sam.threshold)
      }


def initial_filter(gene, sam):
   log.debug('== Allele Filter ==')

   basic_filter(sam)
   sam.dump(gene)

   useless = []
   for an, a in sorted(list(gene.alleles.items()), key=lambda x: sort_key(x[0])):
      remove = any(sam[m] <= 0 for m in a.functional_mutations)
      if remove:
         s = ['Removing {:4}  '.format(an)]
         for m in a.functional_mutations:
            if sam[m] <= 0:
               s.append('{} {}'.format(gene.region_at[m.pos], m))
         log.trace(s[0] + ',  '.join(s[1:]))
         useless.append(an)

   for u in useless:
      del gene.alleles[u]

   log.debug('Initial status')
   for an, a in sorted(list(gene.alleles.items()), key=lambda x: sort_key(x[0])):
      if '/' in an:
         continue
      log.debug('  {}*{}', gene.name, an)
      for m in sorted(a.functional_mutations, key=lambda x: x.pos):
         log.debug(
            '    {:14} {:15} {:4} ({:3.0f}) {} {}',
            gene.region_at[m.pos],
            str(m), sam[m], sam.percentage(m),
            'F', m.aux['old']
         )

   return gene.alleles
