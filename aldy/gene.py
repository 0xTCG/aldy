#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 786

# Aldy source: gene.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.

from typing import Tuple, Dict, List

import os
import enum
import yaml 
import functools
import collections

from .common import *


EXON = 'e'
INTRON = 'i'


class Allele(collections.namedtuple('Allele', ['name', 'cn_config', 'func_muts', 'minors'])):
   """
   Immutable class describing major allele.
   Members:
      name (str):
         Name of the major allele
      cn_config (str):
         CN configuration key
      func_muts (set of `Mutation`):  
         Set of functional mutations that describe this major allele
      minors (dict of str: `Suballele`):
         Dictionary of all minor alleles that inherit this major allele
   """
   pass


class Suballele(collections.namedtuple('Suballele', ['name', 'alt_names', 'neutral_muts'])):
   """
   Immutable class describing minor allele. 
   + Custom printer
   Members:
      name (str):
         Name of the minor allale
      alt_names (list of str):   
         Alternative names for this suballele
      neutral_muts (set of `Mutation`):
         Set of netural mutations defining this suballele
   """
   def __repr__(self):
      return 'Sub({}; [{}])'.format('|'.join([self.name] + self.alt_names), ', '.join(map(str, self.neutral_muts)))

   def __str__(self):
      return self.__repr__()


@functools.total_ordering
class Mutation(collections.namedtuple('Mutation', ['pos', 'op', 'is_functional', 'aux'])):
   """   
   Immutable class representing the mutation.
   + Custom printer
   + Comparable and hashable via ``(pos, op)``
   Members:
      pos (int):
         position in the reference genome
      op (str):
         mutation operation that can be one of the following:
            ``SNP.AB``:  SNP from A to B (A, B are in [ACGT])
            ``INS.xx``:  Insertion of xx (xx is [acgt]+)
            ``DEL.xx``:  Deletion of xx (xx is [ACGT]+)
      functional (bool):
         is this mutation functional?
      aux (dict of str: str):
         auxiliary information dictionary 
         Currently stores dbSNP rsID and Karolinska-style mutation notation
   """
   def __new__(self, pos: int, op: str, is_functional=0, functional=None, aux=dict(), **kwargs):
      if functional is not None:
         is_functional = functional
      return super(Mutation, self).__new__(self, pos, op, is_functional, aux)

   def __repr__(self): return '{}:{}{}'.format(self.pos, self.op, '*' if self.is_functional else '')
   def __str__(self): return self.__repr__()
   def __eq__(self, other): return (self.pos, self.op) == (other.pos, other.op)
   def __lt__(self, other): return (self.pos, self.op) < (other.pos, other.op)
   def __hash__(self): return (self.pos, self.op).__hash__()


class CNConfig(collections.namedtuple('CNConfig', ['cn', 'kind', 'alleles', 'description'])):
   """
   Describes the copy number (CN) configuration of the gene.
   Members:
      cn (dict of int: dict of GeneRegion: int):
         for each gene g, holds the value of expected copy number of each region in the gene g
         e.g. `cn[0][GeneRegion(0, 1, EXON)] == 1` means that the exon 1 of the main gene should be present in this configuration
      kind (`CNConfigType`): 
         Type of the configuration. Can be `LEFT_FUSION`, `RIGHT_FUSION`, `DELETION` and `DEFAULT_CN`.
      alleles (set of str):
         IDs of all alleles that have this CN configuratuon
      description (str):
         human-readable description of the configuration
   """
   CNConfigType = enum.Enum('CNConfigType', 'DEFAULT_CN LEFT_FUSION RIGHT_FUSION DELETION')

   def __repr__(self):
      regions = sorted(set(r for g in self.cn for r in self.cn[g]))
      return 'CNConfig({}; vector={}; alleles=[{}])'.format(
         str(self.kind)[13:],
         '|'.join(
            ''.join(
               (str(self.cn[g][r]) if self.cn[g][r] != .5 else '½')
               if r in self.cn[g] else '_'
               for r in regions)
            for g in sorted(self.cn)),
         ' '.join(sorted(self.alleles)))


class Gene:
   """   
   Provides an absraction over the gene that Aldy targets and associated pseudogenes.
   Holds the list of alleles, mutations, copy number configurations etc.
      """

   def __init__(self, path: str) -> None:
      """
      Initializes the Gene class with the database description in `path` YML file.
      Args:
         path (str): location of YML file
      """
      
      with open(path) as file:
         yml = file.read()
         yml = yaml.safe_load(yml)
         gene_name = os.path.split(path)[-1].split('.')[0].upper()
         self.__parse_yml(gene_name, yml)


   def __parse_yml(self, gene_name: str, yml) -> None:
      """
      Initializes the gene structure from a YML data
      """

      self.__init_basic(gene_name, yml)
      self.__init_regions(yml)
      self.__init_alleles(yml)
      self.__init_partials()
      self.__init_structure(yml)


   def __init_basic(self, gene: str, yml) -> None:
      """
      Reads basic gene properties (`name`, `seq` and `region`)
      """

      #: str: name of the gene (e.g. CYP2D6)
      self.name = gene.upper() 
      #: str: reference sequence that includes *1 allele and all associated suballeles
      #: location-compatible with hg19 (i.e. contains no indels w.r.t. the reference genome)
      self.seq = yml['seq']
      #: `GenomeRange`: region of the ``seq`` in the reference genome
      self.region = GRange(*yml['region']) 
      

   def __init_regions(self, yml) -> None:
      """
      Calculates the genic regions and pseudogenes (`regions`, `unique_regions` and `pseudogenes`).
      Also initializes `region_at()` call.
      """

      assert hasattr(self, 'region')

      def calculate_regions(gene, exons, special):
         """
         Given the list of exons, calculate the intronic regions and return the dictionary of all regions
         """
         regions = {
            GeneRegion(e, EXON): GRange(self.region.chr, *ex) 
            for e, ex in exons.items()
         }
         # Fill introns
         for i in range(1, len(regions)):
            r1, r2 = regions[GeneRegion(i + 1, EXON)], regions[GeneRegion(i, EXON)]
            if r1 > r2:
               r1, r2 = r2, r1
            regions[GeneRegion(i, INTRON)] = GRange(self.region.chr, r1[2], r2[1])
         for n, r in special.items():
            regions[GeneRegion(n, r[0])] = GRange(self.region.chr, r[1], r[2])
         return regions
      
      #: dict of int: (dict of `GeneRegion`: `GenomeRegion`) 
      #: describing the gene exonic, intronic and special (e.g. UTR) regions
      #: Gene 0 is the main gene (key is the gene ID)
      self.regions = {0: calculate_regions(0, yml['exons'], yml['special_regions'])}

      #: list of str: set of pseudogene names
      #: Each pseudogene is associated with index > 0
      self.pseudogenes: List[str] = list()
      if 'pseudogenes' in yml:
         for gi, g in enumerate(sorted(yml['pseudogenes'])):
            self.pseudogenes.append(g)
            self.regions[gi + 1] = calculate_regions(gi + 1, 
               yml['pseudogenes'][g]['exons'], yml['pseudogenes'][g]['special_regions'])
      
      #: dict of int: (int, `GeneRegion`): reverse lookup (gene, region) of the gene regions given a location within reference genome
      self._region_at = {
         i: (g, r) 
         for g, d in self.regions.items()
         for r, rng in d.items()
         for i in range(rng.start, rng.end)
      }

      #: list of (int, str): list of stable regions that will be used for copy number calling
      #: each region is given via tuple (number, type) as in GeneRegion
      def parse_unique(u):
         u = list(filter(None, re.split(r'^(\d+)', u)))
         if len(u) > 1: 
            return GeneRegion(int(u[0]), u[1])
         else: # find matching kind number (legacy support)
            for g in self.regions:
               for r in self.regions[g]:
                  if u[0] == r.kind: return r
         raise KeyError
      self.unique_regions = [parse_unique(y) for y in yml['unique_regions']]

      #print('regions', pr(self.regions))    
      #print('unique_regions', pr(self.unique_regions))   
      #print('pseudogenes', pr(self.pseudogenes))   


   def __init_alleles(self, yml) -> None:
      """
      Initializes gene alleles (`alleles`, `common_tandems`) and copy number configurations (`cn_configs`).
      Initializes old mutation lookup (`old_notation`).
      """

      alleles = {}
      fusions_left = {}
      fusions_right = {}
      
      # Human-readable CN descriptions
      descriptions = {'1': 'Normal allelic configuration: all regions present in both gene and pseudogene'}
      
      #: dict of `Mutation`: string: old Karolinska notation (e.g. 32C>T) for a `Mutation` key
      self._old_notation: Dict[Mutation, str] = collections.defaultdict(str)
      
      # allele ID of the deletion allele (i.e. whole gene is missing). By default it is 0.
      deletion_allele = None

      for allele_name, allele in yml['alleles'].items():
         allele_name = allele_name.split('*')[1]
         mutations: List[Mutation] = []
         if {'op': 'deletion'} in allele['mutations']:
            deletion_allele = allele_name
            descriptions[allele_name] = 'Gene deletion'
            mutations = []
         else:
            for m in allele['mutations']:
               if m['pos'] == 'pseudogene': # has only one mutation indicator
                  if '-' in m['op']:
                     fusions_left[allele_name] = (int(m['op'][1:-1]), (EXON if m['op'][0] == 'e' else INTRON))
                     descriptions[allele_name] = "Fusion: pseudogene until {} followed by the gene".format(m['op'][:-1][::-1])
                  else:
                     if m['op'][-1] == '+':
                        m['op'] = m['op'][:-1]
                     fusions_right[allele_name] = (int(m['op'][1:]), (EXON if m['op'][0] == 'e' else INTRON))
                     descriptions[allele_name] = "Conservation: Pseudogene retention after {} within the gene".format(m['op'][::-1])
               else:
                  if 'old' in m:
                     self._old_notation[Mutation(**m)] = m['old']
                  m['aux'] = {'dbsnp': ' or '.join(m['dbsnp']), 'old': m['old']}
                  mutations.append(Mutation(**m))
         alleles[allele_name] = Suballele(allele_name, [], mutations)

      #: list of tuples of string: list of common major allele pairings used for diplotype calling
      #: e.g. the fact that *13 is always followed by *1 is encoded as `('13', '1')` 
      self.common_tandems: List[tuple] = []
      if 'common_tandems' in yml:
         self.common_tandems = [tuple(y) for y in yml['common_tandems']]

      # Set copy number configurations
      # TODO: currently fusions works only between the main gene and the first pseudogene

      #: dict of str: CNConfig: describes all CN-configurations in the gene with '1' as default CN config
      self.cn_configs: Dict[str, CNConfig] = dict()

      def freezekey(x): # hashing for dictionaries
         return tuple(i[1] for i in sorted(x[0].items())) + tuple(i[1] for i in sorted(x[1].items()))
      inverse_cn: Dict[tuple, str] = dict()

      # Left fusions are fusions of type PSEUDOGENE + GENE (breakpoint has CN of 0.5 for no good reason)
      for a, brk in fusions_left.items():
         cn = dict()
         cn[0] = {r: float(0 if (r.number, r.kind) < brk else 1) for r in self.regions[0]}
         cn[1] = {r: float(1 if (r.number, r.kind) < brk else 0) for r in self.regions[1]}
         cn[0][GeneRegion(*brk)] = cn[1][GeneRegion(*brk)] = 0.5

         key = freezekey(cn)
         if key not in inverse_cn:
            self.cn_configs[a] = CNConfig(cn, CNConfig.CNConfigType.LEFT_FUSION, {a}, descriptions[a])
            inverse_cn[key] = a
         else:
            self.cn_configs[inverse_cn[key]].alleles.add(a)
     
      # Deletion is a special kind of left fusion
      if deletion_allele is not None:
         cn = {0: {r: 0 for r in self.regions[0]}, 1: {r: 1 for r in self.regions[1]}}
         self.cn_configs[deletion_allele] = CNConfig(cn, CNConfig.CNConfigType.DELETION, {deletion_allele}, descriptions[deletion_allele])
      
      # Right fusions are fusions of type GENE + PSEUDOGENE + whole copy of PSEUDOGENE
      for a, brk in fusions_right.items():
         cn = dict()
         cn[0] = {r: (1 if (r.number, r.kind) < brk else 0) for r in self.regions[0]}
         cn[1] = {r: (1 if (r.number, r.kind) < brk else 2) for r in self.regions[1]} 
         key = freezekey(cn)
         if key not in inverse_cn:
            self.cn_configs[a] = CNConfig(cn, CNConfig.CNConfigType.RIGHT_FUSION, {a}, descriptions[a])
            inverse_cn[key] = a
         else:
            self.cn_configs[inverse_cn[key]].alleles.add(a)
      
      # Normal CN case
      used_alleles = {a for _, (_, _, alleles, _) in self.cn_configs.items() for a in alleles}
      default_cn = {g: {r: 1 for r in self.regions[g]} for g in range(2)}
      if self.name == 'CYP2D6': # This is necessary hack for CYP2D6
         default_cn[1][GeneRegion(number=11, kind='pce')] = 0

      self.cn_configs = {allele_number(min(v.alleles)): v for k, v in self.cn_configs.items()}
      self.cn_configs['1'] = CNConfig(default_cn, CNConfig.CNConfigType.DEFAULT_CN, 
         set(alleles.keys()) - used_alleles, descriptions['1'])

      # Set up major and minor allele structures
      alleles_inverse: Dict[tuple, set] = collections.defaultdict(set)
      for a in alleles:
         cn_config = next(cn for cn, conf in self.cn_configs.items() if a in conf.alleles)
         fn_muts = sorted(m for m in alleles[a].neutral_muts if m.is_functional)
         alleles_inverse[(cn_config, tuple(fn_muts))].add(a)

      #: dict of str: Allele: dictionary of major star-alleles associated with this gene
      self.alleles: Dict[str, Allele] = dict()
      # Karolinska DB has a lot of ambiguities, and two alleles with same 'number' can correspond to two different major alleles
      # This step fixes this ambiguity by prepending different letter to each different major allele with the same number
      used_names: Dict[str, int] = {}
      for key, minors in alleles_inverse.items():
         #for an, a in alleles.items():
         an = min(minors)
         a = Allele(name='', cn_config=key[0], func_muts=key[1], minors=minors)
         name = allele_number(an)
         if name in used_names: # Append letter
            used_names[name] += 1
            name += '.' + chr(used_names[name] + ord('a'))
         else:
            used_names[name] = 0
         self.alleles[name] = Allele(name=name, cn_config=a.cn_config, func_muts=set(a.func_muts),
            minors={sa: Suballele(name=sa, alt_names=[], 
               neutral_muts=set(alleles[sa].neutral_muts) - set(a.func_muts)) 
               for sa in a.minors
            })

      # print('alleles', pr(self.alleles))    
      # print('common_tandems', pr(self.common_tandems))   
      # print('cn_configs', pr(self.cn_configs))   
      # print('old_notation', pr(self._old_notation))  
      

   def __init_partials(self) -> None:
      """
      Construct "partial" major alleles:
      i.e. if a major allele is cut in half by a fusion, it will create "new" major allele defined
      with functional mutations that have survided after fusion event
      NOTE: This currently is supported only for left fusions
      """

      def preserved_mutations(f, m):
         def filter_f(m):
            gene, region = self.region_at(m.pos)
            # print(f, gene, region, '\n', self.cn_configs[f])
            # print(sorted(self.cn_configs[f].cn[gene].keys()))
            return self.cn_configs[f].cn[gene][region] > 0
         return set(filter(filter_f, m))
      
      for f in filter(lambda x: self.cn_configs[x].kind == CNConfig.CNConfigType.LEFT_FUSION, self.cn_configs):
         # Do not extend left fusions that are already defined by some functional mutations
         # HACK: This is hacky but works well for now
         if len(self.alleles[f].func_muts) > 0:
            continue
         add: Dict[tuple, Allele] = {}
         for an, a in self.alleles.items():
            if a.cn_config != '1': # We only make partial major alleles from non-fused major alleles
               continue
            new_name = '{}/{}'.format(f, an)
            new_muts = preserved_mutations(f, a.func_muts)
            key = sorted_tuple(new_muts)
            if key in add:
               add[key].minors.update({
                  san: Suballele(san, [], preserved_mutations(f, sa.neutral_muts))
                  for san, sa in a.minors.items()
               })
            else:
               add[key] = Allele(name=new_name, cn_config=f, func_muts=new_muts, 
                  minors={san: Suballele(
                     name=san, alt_names=[], 
                     neutral_muts=preserved_mutations(f, sa.neutral_muts))
                     for san, sa in a.minors.items()})
         
         # Remove fusion (will be replaced at least by allele "1/{f}")
         del self.alleles[f]
         self.alleles.update({a.name: a for a in add.values()})
      
      for an, a in self.alleles.items():
         # Clean up minor alleles (as many might be equivalent after fusions)
         # Put references to the cleaned-up alleles in `alt_names` field
         minors: Dict[tuple, List[str]] = collections.defaultdict(list)
         for s in a.minors:
            key = sorted_tuple(a.minors[s].neutral_muts)
            minors[key].append(s)
         self.alleles[an] = Allele(
            self.alleles[an].name, 
            self.alleles[an].cn_config, 
            self.alleles[an].func_muts,
            {min(sa): Suballele(min(sa), 
                                alt_names=list(set(sa) - {min(sa)}), 
                                neutral_muts=nm) 
             for nm, sa in minors.items()})

      # TODO: prune identical post-fusion alleles (e.g. *2 after *13A and *13B is the same--no need to have 2 items)
      # print('alleles', pr(self.alleles))    


   def __init_structure(self, yml) -> None:
      """
      Initialize `mutations` lookup and protein coding region structure `coding_region`
      """

      # Set up coding region structure for aminoacid calculation
      CodingRegion = nt('CodingRegion', ['rev_comp', 'aminoacid', 'lookup'])

      is_rev_comp = True if yml['rev_comp'] == 1 else False
      lookup = {
         i: self.seq[i - self.region.start]
         for (_, kind), (_, start, end) in self.regions[0].items()
         if kind == EXON
         for i in range(start, end)
      }
      seq = ''.join(lookup[i] for i in sorted(lookup.keys()))
      aminoacid = seq_to_amino(rev_comp(seq) if is_rev_comp else seq)

      #: CodingRegion(is_rev_comp=bool, aminoacid=str, lookup=dict of int: char): 
      #: Describes coding region of the gene with its aminoacid and the reverse complement status in the referenc genome
      self.coding_region = CodingRegion(is_rev_comp, aminoacid, lookup)

      #: dict of (int, str): `Mutation`: get `Mutation` for any query (position, mutation_type)
      self.mutations: Dict[int, str] = {}
      for _, a in self.alleles.items():
         self.mutations.update({(m.pos, m.op): m for m in a.func_muts})
         for _, sa in a.minors.items():
            self.mutations.update({(m.pos, m.op): m for m in sa.neutral_muts})

      # print('coding_region', pr(self.coding_region))    
      # print('mutations', pr(self.mutations))    

   
   ####################################################################################################

   def region_at(self, pos: int) -> Tuple[int, GeneRegion]:
      """
      Returns tuple consisting of a gene ID and a region of the given position in the gene. 
      Args:
         pos (int): position in the reference genome 
      """
      return self._region_at[pos]


   def old_notation(self, mut: Mutation) -> str:
      """
      Returns old Karolinska-style notation (e.g. 32C>T) for a `Mutation` key 
      Args:
         mut (`Mutation`): mutation to be queried
      """
      return self._old_notation[mut]


   def sequence(self, gene: int, pad=0) -> str:
      """
      Returns the sequence of a gene
      Args:
         id (int): gene ID (0 for the main gene, 1+ for pseudogenes)
         pad (int, default: 0): extra padding on both sides in bps
      """
      regs = [(st, ed) for _, (_, st, ed) in self.regions[gene].items()]
      mi, ma = min(r[0] for r in regs), max(r[1] for r in regs)
      off = self.region.start
      return self.seq[max(0, mi - off - pad):min(ma - off + pad, len(self.seq))]


   def check_functional(self, m: Mutation) -> bool:
      """Check is a mutation functional (i.e. does it change gene's aminoacid or not)"""
      if m.pos not in self.coding_region.lookup:
         return False
      if m.op[:3] != 'SNP':
         return True
      seq = ''.join(self.coding_region.lookup[i] if i != m.pos else m.op[5] for i in sorted(self.coding_region.lookup))
      amino = seq_to_amino(rev_comp(seq) if self.coding_region.rev_comp else seq)
      return amino != self.coding_region.aminoacid


   def print_configurations(self):
      raise NotImplementedError

      # log.error('Gene {}\n', self.name)
      # for name, config in sorted(self.cnv_configurations.items()):
      #    c6 = [c[2:] for c in config if c.startswith('6.')]
      #    c7 = [c[2:] for c in config if c.startswith('7.')]
      #    labels = sorted(list(set(c6) | set(c7)))
      #    max_width = str(max(3, max(len(x) for x in labels)) + 1)

      #    description = '{}\n'.format(self.descriptions[name]) if name in self.descriptions else ''
      #    alleles = sorted(set('*' + a.split('/')[0] for a in self.alleles if self.alleles[a].cnv_configuration == name))
      #    description += 'Alleles: ' 
      #    for i in range(0, len(alleles), 12):
      #       if i != 0: description += '\n         ';
      #       description += ', '.join(alleles[i:i + 12]) 

      #    log.warn('Name: {}\n{}'.format(name, description))
      #    print('            ', end='')
      #    for c in labels:
      #       print(('{:>' + max_width + '}').format(c), end='')
      #    print(' ')
      #    print('Gene:       ', end='')
      #    for c in labels:
      #       if '6.' + c in config:
      #          print(('{:' + max_width + '}').format(config['6.' + c]), end='')
      #    print(' ')
      #    print('Pseudogene: ', end='')
      #    for c in labels:
      #       if '7.' + c in config:
      #          print(('{:' + max_width + '}').format(config['7.' + c]), end='')
      #    print('\n')


