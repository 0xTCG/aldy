#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 786

# Aldy source: gene.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from __future__ import division
from __future__ import print_function
from builtins import next
from builtins import filter
from builtins import chr
from builtins import map
from builtins import range
from builtins import object
from collections import namedtuple as nt
from pprint import pprint

import os
import re
import yaml
import collections

from .common import *


# disrupting_mutations, neutral_mutations
class Suballele(nt('Suballele', ['name', 'alt_names', 'neutral_mutations'])):
   """
      alt_names         = list of other names for this suballele
      neutral_mutations    = unordered set of Mutation
   """
   def __repr__(self):
      return '{}: [{}]'.format('|'.join([self.name] + self.alt_names), ', '.join(map(str, self.neutral_mutations)))

   def __str__(self):
      return self.__repr__()


class Allele(nt('Allele', ['name', 'cnv_configuration', 'functional_mutations', 'suballeles'])):
   """
      cnv_configuration       = key for CNV configuration
      functional_mutations    = unordered set of Mutation
      suballeles           = dict of Suballele
   """
   pass


class Mutation(nt('Mutation', ['pos', 'op', 'functional', 'aux'])):
   """   Class representing the mutation
      pos      = position in the reference genome
      op          = mutation oprtation, one of the following:
         SNP.AB:  SNP from A to B (A, B are in [ACGT])
         INS.xx:  Insertion of xx (xx is [acgt]+)
         DEL.xx:  Deletion of xx (xx is [ACGT]+)
      functional  = does it impact the protein?
      aux      = auxiliary information dict (right now dbSNP rsID and old notation)
   """
   def __new__(self, pos, op, functional=0, aux=None, **kwargs):
      if not aux: 
         aux = dict()
      return super(Mutation, self).__new__(self, pos, op, functional, aux)

   def __repr__(self):
      return '{}:{}'.format(self.pos, self.op)

   def __str__(self):
      return self.__repr__()

   def __eq__(self, other):
      return (self.pos, self.op) == (other.pos, other.op)

   def __cmp__(self, other):
      return (self.pos, self.op).__cmp__((other.pos, other.op))

   def __hash__(self):
      return (self.pos, self.op).__hash__()

   def __ne__(self, other):
      return not self.__eq__(other)


class Gene(object):
   """   Class representing the gene cluster structure
      name     = gene name
      seq      = gene nucleotide sequence
      region      = (chr, start, end)
      cnv_region  = (chr, start, end)
      regions  = { region_name: (start, end) }
         region_name is 'gene.name', where:
         gene is in {6, 7, ...}
         name is {1e, 1i, ..., pcr} etc.
      region_at   = { pos: region_name }
      alleles  = { allele_name: Allele }
      max_cn      = { allele_name: max_cn }
      
      unique_regions = [ region_name ... ]
      
      structures
      deletion_allele = name of deletion allele
      fusions     = colection of fusions! {  }
      old_notation   = lookup { Mutation: karolinska_snp_loci }
         this represents the canonical mutation code
         for the gene (e.g. -77:G>A)
   """
   def __init__(self, path, cnv_region=None):
      self.load(path, cnv_region)

   def load(self, path, cnv_region):
      with open(path) as file:
         j = file.read()
         j = yaml.safe_load(j)

      ## (1) Basics
      self.name = os.path.split(path)[-1].split('.')[0].upper()
      self.seq = j['seq']
      self.region = tuple(j['region'])
      self.cnv_region = tuple(j['cnv_region'])
      if cnv_region:
         r = re.match(r'^(.+?):(\d+)-(\d+)$', cnv_region)
         if not r:
            raise Exception('Parameter --cn-neutral-region={} is not in the format chr:start-end (where start and end are numbers)', cnv_region)
         ch = r.group(1)
         if not ch.startswith('chr'):
            ch = 'chr' + ch
         self.cnv_region = [ch, int(r.group(2)), int(r.group(3))]
         log.info('Using {} as copy-number neutral region', cnv_region)

      ## (2) Set up gene regions
      # Convention: Main gene has prefix 6; pseudogenes are numbered 7, 8
      def calculate_regions(id, exons, special):
         regions = {'{}.{}e'.format(id, e): tuple(ex) for e, ex in exons.items()}
         for i in range(1, len(regions)):
            r1, r2 = regions['{}.{}e'.format(id, i + 1)], regions['{}.{}e'.format(id, i)]
            if r1 > r2:
               r1, r2 = r2, r1
            regions['{}.{}i'.format(id, i)] = (r1[1], r2[0])
         for n, r in special.items():
            regions['{}.{}'.format(id, n)] = tuple(r)
         return regions
      self.regions = calculate_regions(6, j['exons'], j['special_regions'])
      if 'pseudogenes' not in j:
         j['pseudogenes'] = {}
      for gi, g in enumerate(sorted(j['pseudogenes'])):
         g = j['pseudogenes'][g]
         self.regions.update(calculate_regions(7 + gi, g['exons'], g['special_regions']))
      self.region_at = collections.defaultdict(str, {
         i: n
         for n, j in self.regions.items()
         for i in range(j[0], j[1])
      })
      self.max_cn = j['max_cn']
      self.unique_regions = j['unique_regions']

      ## (3) Set up alleles
      # is_left_fusion: 6.0 begin = 0
      self.common_tandems = []
      if 'common_tandems' in j:
         self.common_tandems = list(map(tuple, j['common_tandems']))
      alleles = {}
      fusions_left = {}
      fusions_right = {}
      self.descriptions = collections.defaultdict(str)
      self.descriptions['1'] = 'Normal allelic configuration: all regions present in both gene and pseudogene'
      self.old_notation = collections.defaultdict(str)
      self.deletion_allele = '0'
      for allele_name, allele in j['alleles'].items():
         allele_name = allele_name.split('*')[1]
         mutations = []
         if {'op': 'deletion'} in allele['mutations']:
            self.deletion_allele = allele_name
            self.descriptions[allele_name] = 'Gene deletion (only pseudogene is present)'
            mutations = []
         else:
            for m in allele['mutations']:
               # if 'pos' not in m:
                  # continue
               if m['pos'] == 'pseudogene': # has only one mutation indicator
                  if '-' in m['op']:
                     fusions_left[allele_name] = m['op'][:-1][::-1] # eN -> Ne
                     self.descriptions[allele_name] = "Fusion: pseudogene until {} followed by the gene".format(m['op'][:-1][::-1])
                  else:
                     if m['op'][-1] == '+':
                        m['op'] = m['op'][:-1]
                     fusions_right[allele_name] = m['op'][::-1]
                     self.descriptions[allele_name] = "Conservation: Pseudogene retention after {} within the gene".format(m['op'][::-1])
               else:
                  if 'old' in m:
                     self.old_notation[Mutation(**m)] = m['old']
                  m['aux'] = {'dbsnp': ','.join(m['dbsnp']), 'old': m['old']}
                  mutations.append(Mutation(**m))
         alleles[allele_name] = Suballele(allele_name, [], mutations)

      ## (4) Set up copy number structures
      regions = sorted(set(r[2:] for r in self.regions) - set('cnv'))
      fusions_left = { # pseudogene
         a: tuple(1 for r in regions if r < l) + (.5,) + tuple(0 for r in regions if r > l)
         for a, l in fusions_left.items()
      }
      fusions_left[self.deletion_allele] = tuple([1] * len(regions))
      fusions_left = {
         a: tuple(1 - y for y in l) + l
         for a, l in fusions_left.items()
      }
      fusions_right = { # pseudogene
         a: tuple(0 for r in regions if r < l) + tuple(1 for r in regions if r >= l)
         for a, l in fusions_right.items()
      }
      fusions_right = {
         a: tuple(1 - y for y in l) + tuple(1 + y for y in l)
         for a, l in fusions_right.items()
      }
      self.fusions = dict(list(fusions_left.items()) + list(fusions_right.items()))

      self.cnv_configurations = collections.defaultdict(set)
      for a, c in self.fusions.items():
         self.cnv_configurations[c].add(a)
      reference_cnv = tuple([1] * len(regions) * 2)
      self.cnv_configurations[reference_cnv] |= set(a for a in alleles if a not in self.fusions)
      self.cnv_configurations = {
         ('1' if c == reference_cnv else allele_key(min(a))): (dict(
            [('6.' + regions[yi], y) for yi, y in enumerate(c[:len(c) // 2])] +
            [('7.' + regions[yi], y) for yi, y in enumerate(c[len(c) // 2:])]
         ), a) for c, a in self.cnv_configurations.items()
      }
      self.descriptions = {allele_key(a): d for a, d in self.descriptions.items()}
      if self.name == 'CYP2D6': # This is necessary hack
         self.cnv_configurations['1'][0]['6.pce'] = 0
      self.fusions = { # TODO: check assumption that each fusion has unique allele_key
         allele_key(a): dict(
            [('6.' + regions[yi], y) for yi, y in enumerate(t[:len(t) // 2])] +
            [('7.' + regions[yi], y) for yi, y in enumerate(t[len(t) // 2:])]
         )
         for a, t in self.fusions.items()
      }
      self.fusions_left = {allele_key(a): c for a, c in fusions_left.items()}

      ## (5) Set up major and minor allele structures
      self.alleles = collections.defaultdict(set)
      for a in alleles:
         cnv_configuration = next(cn for cn, (conf, alls) in self.cnv_configurations.items() if a in alls)
         self.alleles[tuple([cnv_configuration] + sorted(m for m in alleles[a].neutral_mutations if m.functional))].add(a)
      self.alleles = {
         min(suballeles): Allele('', key[0], key[1:], suballeles)
         for key, suballeles in self.alleles.items()
      }
      self.cnv_configurations = {cn: config[0] for cn, config in self.cnv_configurations.items()}
      used_names = {}
      for an, a in self.alleles.items():
         name = allele_key(an)
         if name in used_names: # append ' to the proteins with same numbers to distinguish them
            used_names[name] += 1
            name += '.' + chr(used_names[name] + ord('a'))
         else:
            used_names[name] = -1
         self.alleles[an] = Allele(
            name,
            a.cnv_configuration,
            set(a.functional_mutations),
            {sa: Suballele(sa, [], set(alleles[sa].neutral_mutations) - set(a.functional_mutations)) for sa in a.suballeles}
         )
      self.alleles = {a.name: a for a in self.alleles.values()}

      ## (6) Add possible partial fusion products
      # TODO: currently only empty left fusions: extend to all fusions (CHECK does it happen)?
      def filter_missing(f, m):
         def filter_f(m):
            return self.fusions[f][self.region_at[m.pos]] > 0
         return set(filter(filter_f, m))
      for f in self.fusions_left:
         if f == self.deletion_allele:
            continue
         # do not extend left alleles which are already defined by some SNVs
         if len(self.alleles[f].functional_mutations) > 0:
            continue
         add = {}
         for an, a in sorted(list(self.alleles.items()), key=lambda x: sort_key(x[0])):
            if a.cnv_configuration != '1':
               continue
            n = f + '/' + an
            key = sorted_tuple(filter_missing(f, a.functional_mutations))
            if key in add:
               add[key] = Allele(
                  add[key].name, add[key].cnv_configuration,
                  add[key].functional_mutations,
                  dict(
                     list(add[key].suballeles.items()) +
                     [(sa, Suballele(sa, [], filter_missing(f, set(alleles[sa].neutral_mutations) - set(a.functional_mutations)))) for sa in a.suballeles]
                  )
               )
            else:
               add[key] = Allele(
                  n, f,
                  filter_missing(f, a.functional_mutations),
                  {sa: Suballele(sa, [], filter_missing(f, set(alleles[sa].neutral_mutations) - set(a.functional_mutations))) for sa in a.suballeles}
               )
         del self.alleles[f]
         self.alleles.update({a.name: a for a in add.values()})
      for an, a in self.alleles.items():
         suballeles = collections.defaultdict(list)
         for s in a.suballeles:
            key = sorted_tuple(a.suballeles[s].neutral_mutations)
            suballeles[key].append(s)
         self.alleles[an] = Allele(
            self.alleles[an].name,
            self.alleles[an].cnv_configuration,
            self.alleles[an].functional_mutations,
            {min(sa): Suballele(min(sa), list(set(sa) - set([min(sa)])), set(nm)) for nm, sa in suballeles.items()}
         )

      ## (7) Set up coding region structure for aminoacid calculation
      self.coding_region = {}
      self.rev_comp = True if j['rev_comp'] == 1 else False
      for ri, (start, end) in self.regions.items():
         p = re.match(r'6\.(\d+)e', ri)
         if not p:
            continue
         for i in range(start, end):
            self.coding_region[i] = self.seq[i - self.region[1]]

      seq = ''.join(self.coding_region[i] for i in sorted(self.coding_region.keys()))
      self.aminoacid = seq_to_amino(rev_comp(seq) if self.rev_comp else seq)

      self.mutations = {}
      for _, a in self.alleles.items():
         self.mutations.update({(m.pos, m.op): m for m in a.functional_mutations})
         for _, sa in a.suballeles.items():
            self.mutations.update({(m.pos, m.op): m for m in sa.neutral_mutations})

   def print_configurations(self):
      log.error('Gene {}\n', self.name)
      for name, config in sorted(self.cnv_configurations.items()):
         c6 = [c[2:] for c in config if c.startswith('6.')]
         c7 = [c[2:] for c in config if c.startswith('7.')]
         labels = sorted(list(set(c6) | set(c7)))
         max_width = str(max(3, max(len(x) for x in labels)) + 1)

         description = '{}\n'.format(self.descriptions[name]) if name in self.descriptions else ''
         alleles = sorted(set('*' + a.split('/')[0] for a in self.alleles if self.alleles[a].cnv_configuration == name))
         description += 'Alleles: ' 
         for i in range(0, len(alleles), 12):
            if i != 0: description += '\n         ';
            description += ', '.join(alleles[i:i + 12]) 

         log.warn('Name: {}\n{}'.format(name, description))
         print('            ', end='')
         for c in labels:
            print(('{:>' + max_width + '}').format(c), end='')
         print(' ')
         print('Gene:       ', end='')
         for c in labels:
            if '6.' + c in config:
               print(('{:' + max_width + '}').format(config['6.' + c]), end='')
         print(' ')
         print('Pseudogene: ', end='')
         for c in labels:
            if '7.' + c in config:
               print(('{:' + max_width + '}').format(config['7.' + c]), end='')
         print('\n')

