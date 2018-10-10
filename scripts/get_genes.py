#!/usr/bin/env python

# 786

# Short notes
# (1) This script might break if anything changes (Karolinska's DB, NCBI etc)
# (2) It has been tested with hg19 only
# (3) Few assumptions:
#     - dupA_B = insA_B
#     - all AxN alleles are ignored (persumed to be equal with A)

import sys, json, re, copy, collections, os
from pprint import pprint
from collections import OrderedDict
from noindent import NoIndent, NoIndentEncoder
from common import log
import karolinska, ncbi, dbsnp, merger

def process_mut_helper(m, dbsnp):
   result = []
   rs = re.findall(r'\s*?([-0-9_]+)\s*?((ins|del|dup|[ACGT])(\S*))', ' ' + m) # Match 100C>T, -100C>T, -100 C>T etc.
   if len(rs) == 0:
      return None
   for r in rs:
      pos, op = int(r[0].split('_')[0]), r[1].split()[0]
      op = op.replace('*', '')
      if len(op) == 1: # e.g. 51A
         op = 'A>{}'.format(op)
      elif '>' in op and len(op) != 3: # e.g. -1601_-1600GA>TT;
         op = op.split('>')
         assert(len(op) == 2)
         if len(op[0]) == len(op[1]):
            result += [[pos + i, '{}>{}'.format(op[0][i], op[1][i]), dbsnp] for i in xrange(len(op[1])) if op[0][i] != op[1][i]]
            continue
         else: # e.g. 3030G>G/A
            if '/' in op[1] and len(op[0]) == 1:
               result += [[pos, '{}>{}'.format(op[0], c), dbsnp] for c in op[1].split('/') if op[0] != c]
            else:
               log.warn('Main: Ignoring {}', m)
            continue
      if op[:3] == 'dup': 
         op = 'ins' + op[3:]
      if op[-2:] == 'x2': # detect ins<something>x2 
         op = op[:-2] + op[3:-2]
      result.append([pos, op, dbsnp])
   return result

def parser_helper(row, dbsnp, process_mut, column, ignored_columns):
   text = row.get_text().encode('ascii', 'ignore').strip()
   items = []
   for i, col in enumerate(row.find_all('td')):
      text = col.get_text().encode('ascii', 'ignore').strip()
      if i == column:
         for a in col.find_all('a'):
            r = re.search(r'snp_ref.cgi\?rs=(.+)$', a['href'])
            if r is not None:
               snp = a.get_text().encode('ascii', 'ignore').strip()
               rs = r.group(1)
               if not rs.startswith('rs'): rs = 'rs' + rs
               if snp in dbsnp and dbsnp[snp] != rs:
                  print >>sys.stderr, 'error {} != {} for {}'.format(rs, dbsnp[rs], snp)
               dbsnp[snp] = rs
         # split mutations
         text = filter(None, re.split(r'[;,:]\s*', text))
         muts = [process_mut(x, dbsnp[x]) for x in text]
         muts = [x if type(x) is list else [x] for x in muts if x is not None]
         items.append([n for m in muts for n in m])
         text = ' '.join(text)
      if i not in ignored_columns:
         items.append(text)
   return items

def get_data(force, gene, pseudogene, reverse_complement, parser, 
   fix_karolinska, genome_range, gene_ids, coordinate, patch, 
   post_process, functional_exceptions, unique_regions, max_cn,
   custom_url = None):
   def sf(x):
      y = re.split(r'(\d+)', x[len(gene):])
      return int(y[1]), y[2]

   # Get Karolinska's data
   cypdata = karolinska.get_karolinska_database(gene, parser, force, custom_url)
   if fix_karolinska is not None:
      fix_karolinska(cypdata)
   
   #pprint (cypdata)

   # Get NCBI data for genes and reference genome
   genes, hg19 = ncbi.get_genomes(
      gene_ids[0], genome_range, gene_ids[1:], force=force, reverse_complement=reverse_complement
   )

   new_seq = genes[gene].seq.tomutable()
   for c, n in patch:
      new_seq[coordinate(c, genes[gene])] = n
   genes[gene] = genes[gene]._replace(seq=new_seq.toseq())

   # Fix Karolinska's coordinates
   result = merger.merge(cypdata, genes[gene], coordinate, functional_exceptions, reverse_complement)
   ## pprint(genes['CYP21'].translation)
   ## pprint(genes['CYP21P'].translation)

   mx = collections.defaultdict(lambda: ['', []])
   for a in result:
      for m in result[a]['mutations']:
         mx[(m['pos'], m['op'])][0] = m
         mx[(m['pos'], m['op'])][1].append(a)
   for m in genes[gene].pseudo_mutations.values():
      m['functional'] = merger.is_functional(genes[gene], m, genes[gene].pseudo_mutations.values(), True)
      # if (m['pos'], m['op']) in mx:
      # 	log.warn('[{}] {} (from {}) originates from {}',
      # 		' F'[mx[(m['pos'], m['op'])][0]['functional']], 
      # 		mx[(m['pos'], m['op'])][0]['old'],
      # 		','.join(set(mx[(m['pos'], m['op'])][1])),
      # 		m['old']
      # 	)

   # Remove mutations not present in hg19 and fix the coordinates
   for a in result:
      for m in result[a]['mutations']:
         if m['pos'] == 'pseudogene': continue
         if m['pos'] not in genes[gene].translation:
            log.warn('Main: Translation not found for {}: {} ({})', a, m['old'], m['pos'])
            m['pos'] = None
         else:
            m['pos'] = genes[gene].translation[m['pos']] 
      result[a]['mutations'] = [m for m in result[a]['mutations'] if not m['pos'] is None]

   # Fetch missing dbSNP links
   result = dbsnp.get_dbsnp(result, genome_range, force)

   # Fix exon and intron coordinates
   for _, g in genes.iteritems():
      g.exons[:] = map(lambda x: (g.translation[int(x.start)], g.translation[int(x.end)]), g.exons)
      g.introns[:] = map(lambda x: (g.translation[int(x.start)], g.translation[int(x.end)]), g.introns)

   # Patch hg19 with reference SNPs
   hg19 = list(hg19)
   for gi, hi in genes[gene].translation.iteritems():
      if hg19[hi - genome_range[1]] != genes[gene].seq[gi]:
         hg19[hi - genome_range[1]] = genes[gene].seq[gi]
   hg19 = ''.join(hg19)

   result.update({
      gene + '*1': { 
         'mutations': [], 
         'phenotype': { 'invivo': 'Normal', 'invitro': 'Normal' }
      }
   })

   # Add missing regions
   post_process(genes, result)

   hoi=collections.OrderedDict()
   for pos, m in genes[gene].pseudo_translation.iteritems():
      hoi[genes[gene].translation[pos]] = NoIndent((
         genes[pseudogene].translation[m['old_pos']],
            m['op'] if 'op' in m else ''))
   return dict(
      #map=hoi,
      seq=hg19, 
      region=NoIndent(genome_range),
      name=gene,
      exons={'{}'.format(ei + 1): NoIndent(e) for ei, e in enumerate(genes[gene].exons)},
      special_regions={g: NoIndent(gg) for g, gg in genes[gene].special_regions.iteritems()},
      pseudogenes={
         g: {
            'exons': {'{}'.format(ei + 1): NoIndent(e) for ei, e in enumerate(genes[g].exons)},
            'special_regions': {g: NoIndent(gg) for g, gg in genes[g].special_regions.iteritems()}
         }
         for g in [pseudogene]
      } if pseudogene is not None else {},
      # Regions used for CNV detection of each gene
      unique_regions=NoIndent(unique_regions),
      # Unique CYP2D8 region used for CNV detection
      # Based on BLAT, that is [5e-4i-4e]
      cnv_region=NoIndent(('chr22', 42547463, 42548249)),
      alleles=OrderedDict([
         (a, {
            'phenotype': NoIndent(result[a]['phenotype']),
            'mutations': [
               NoIndent(OrderedDict([(x, y[x]) for x in sorted(y, reverse=True)]))
               for y in result[a]['mutations']
            ] 
         })
         for a in sorted(result, key=sf) 
      ]),
      max_cn=max_cn
   )

# TODO: "variable number of A's in the region -1258 to -1237a" is ignored
# TODO: "CYP2D7-like 3'-flanking region" in *10B ignored
# TODO: "CYP2D7 conversion upstream of exon 1 (-225 to -431)" in *35B is ignored
def get_cyp2d6(force):
   gene, pseudogene = 'CYP2D6', 'CYP2D7'
   
   def process_mut(m, dbsnp):
      if pseudogene not in m:
         return process_mut_helper(m, dbsnp)
      else:
         r = re.search(r'(\w+) (\d)', m)
         if r is None or 'upstream' in m:
            return None 
         pos = 'pseudogene'
         op = '{}{}'.format(r.group(1).lower()[0], r.group(2))
         if op[0] not in 'ie': return None # only introns/exons
         if 'onwards' in m: op += '+'
         r = re.search(r'\((\d+)-(\d+)\)', m)
         if r is not None: op += ':{}-{}'.format(r.group(1), r.group(2))
         return [[pos, op, dbsnp]]
   
   def parser(row, dbsnp):
      text = row.get_text().encode('ascii', 'ignore').strip()
      if not re.search(r'CYP2D6\*13([A-Z]+\d*)*\W', text) != None:
         return parser_helper(row, dbsnp, process_mut, 2, [])
      else:
         orig, name, region = '', '', ['']
         w = re.search(r'[Oo]riginally\s+called\s+([\*\w]+)', text)
         if w != None: orig = w.group(1)
         w = re.search(r'[Tt]emporarily\s+called\s+([\*\w]+)', text)
         if w != None: name = w.group(1)
         w = re.search(r'switch\s+region in ([\w -]+)', text)
         if w != None: 
            region = w.group(1).split('-')
            region = [r.split() for r in region]
            region =  [r[0][0] + r[1][0] + '-' for r in region]
         if orig == '': orig = name
         # TODO incorporate all possitibilites (if the split region is in the middle)
         return [orig, name, [['pseudogene', region[0], '*']], 'None', 'None', '']
      
   def fix_karolinska(cypdata):
      cypdata[gene + '*82']['mutations'] = [m for m in cypdata[gene + '*82']['mutations'] if m[0] != 'pseudogene']
      del cypdata[gene + '*68B']

      for g in ['21A', '21B']:
         a = cypdata[gene + '*' + g]['mutations']
         for mi, m in enumerate(a):
            if '{}:{}'.format(m[0], m[1]) == '2573:insC':
               a[mi][0] += 6

      extra_alleles = [m for m in cypdata[gene + '*4N']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in 
         '4401:C>T 3384:A>C 3582:A>G 2097:A>G 997:C>G 843:T>G 746:C>G 310:G>T -1000:G>A -1235:A>G'.split()]
      for i in '4A 4B 4C 4D 4E 4F 4G 4H 4J 4K 4L'.split():
         cypdata[gene + '*' + i + 'W'] = copy.deepcopy(cypdata[gene + '*' + i])
         cypdata[gene + '*' + i + 'W']['mutations'] += [f for f in extra_alleles if f not in cypdata[gene + '*' + i + 'W']['mutations']]

      extra_alleles = [m for m in cypdata[gene + '*2A']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in
         '-1584:C>G -1235:A>G -740:C>T -678:G>A'.split()]
      extra_alleles += [m for m in cypdata[gene + '*4N']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in 
         '310:G>T 746:C>G 843:T>G'.split()]
      for i in '2A 2B 2C 2D 2E 2F 2G 2H 2K 2L 2M'.split():
         cypdata[gene + '*' + i + 'W'] = copy.deepcopy(cypdata[gene + '*' + i])
         cypdata[gene + '*' + i + 'W']['mutations'] += [f for f in extra_alleles if f not in cypdata[gene + '*' + i + 'W']['mutations']]

      extra_alleles = [m for m in cypdata[gene + '*2A']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in '-1584:C>G'.split()]
      for i in '11 14A 14B 28'.split():
         cypdata[gene + '*' + i + 'W'] = copy.deepcopy(cypdata[gene + '*' + i])
         cypdata[gene + '*' + i + 'W']['mutations'] += [f for f in extra_alleles if f not in cypdata[gene + '*' + i + 'W']['mutations']]

      extra_alleles = [m for m in cypdata[gene + '*2A']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in 'pseudogene:i1'.split()]
      extra_alleles += [m for m in cypdata[gene + '*36']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in '-1235:A>G -1000:G>A 310:G>T 843:T>G 3582:A>G 3384:A>C'.split()]
      for i in '10A 10B 10D'.split():
         cypdata[gene + '*' + i + 'W'] = copy.deepcopy(cypdata[gene + '*' + i])
         cypdata[gene + '*' + i + 'W']['mutations'] += [f for f in extra_alleles if f not in cypdata[gene + '*' + i + 'W']['mutations']]

      extra_alleles = [m for m in cypdata[gene + '*36']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in '-1235:A>G 310:G>T 843:T>G 3384:A>C'.split()]
      for i in '17'.split():
         cypdata[gene + '*' + i + 'W'] = copy.deepcopy(cypdata[gene + '*' + i])
         cypdata[gene + '*' + i + 'W']['mutations'] += [f for f in extra_alleles if f not in cypdata[gene + '*' + i + 'W']['mutations']]

      # For CYP2D6, all intron 1 conservations are interpreted as i1:214-245
      # All exon 9 conservations include all regions after exon 9 
      for a in cypdata:
         for mi, m in enumerate(cypdata[a]['mutations']):
            if m[0] == 'pseudogene' and m[1] == 'i1':
               m[1] += ':214-245'
               # cypdata[a]['mutations'][mi] = (m[0], m[1] + ':214-245', m[2])
            if m[0] == 'pseudogene' and m[1] == 'e9':
               m[1] += '+'
               # cypdata[a]['mutations'][mi] = (m[0], m[1] + '+', m[2])

   def coordinate(p, g): 
      p = -p - 1619 + len(g.seq)
      return p + (1 if p < 0 else 0)

   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[8][0]
      g.special_regions.update({
         '0i':  (g.exons[0][1], g.translation[len(g.seq) - 1]),
         'ins': (gs - 600,        gs),
         'rep': (gs - 600 - 2800, gs - 600)
      })
      g = genes[pseudogene]
      gs = g.exons[8][0]
      #pprint(g.translation)
      g.special_regions.update({
         '0i':  (g.exons[0][1], g.translation[len(g.seq) - 1]),
         'ins': (gs - 600,               gs),
         'pce': (gs - 600 - 1600,        gs - 600), # CYP2D7 spacer; ins <= pce <= rep
         'rep': (gs - 600 - 1600 - 2800, gs - 600 - 1600)
      })

      # Fix database mistakes and nonsenses
      result[gene + '*77'] = result.pop(gene + '*13G1') # 77 is same as 13; ignoring 77
      result[gene + '*68'] = result.pop(gene + '*68A')
      del result[gene + '*13G2'] # 13G1 = 13G2
      # remove e9 mutation since it is already registered in conversion
      for i in ['36', '57', '4N', '83']:
         i = gene + '*' + i
         result[i]['mutations'] = [m for m in result[i]['mutations'] if m['old'] != '4180:G>C']

      result.update({
         'CYP2D6*5': { 
            'mutations': [{'op': 'deletion'}],  
            'phenotype': { 'invivo': 'None (d, s)', 'invitro': '' }
         },
      })

   hg19_range = ('chr22', 42518900, 42546000)
   return get_data(force, gene, pseudogene, True, parser, fix_karolinska, 
      hg19_range, ['224384747', (gene, 'M33388.1'), (pseudogene, 'M33387.1', {'seq_start': 7000})],
      coordinate, [], post_process,
      functional_exceptions=['2988:G>A', '-1584:C>G', '2939:G>A'],
      unique_regions=[
         '1e', '1i',
         '2e', '2i',
         '3e', 
         '5e', '5i',
         '6e', '6i',
         '9e',
         'pce'
      ], 
      max_cn={
         '1':  20,
         '2':  20,
         '41': 20,
         '17': 20,
         '68': 20,
         'X9': 20,
         '4':  2,
         '9':  2,
         '10': 2,
         '35': 2,
         '36': 2,
         '5':  2
      }
   )

# TODO: gene conversion in the 3' flanking region
# TODO: gene conversion in 3'-UTR
# TODO: gene conversion at 7073-7082
# TODO: -1199_-1198ins316bpAlu
def get_cyp2a6(force):
   gene, pseudogene = 'CYP2A6', 'CYP2A7'
   
   def process_mut(m, dbsnp):
      if gene in m: # exons 3-9 of CYP2A6 origin etc
         return None
      elif pseudogene in m:
         r = re.search(r'[Ee]xons (\d)-(\d)', m)
         pos = 'pseudogene'
         op = 'e{}-'.format(r.group(2))
         return [(pos, op, dbsnp)]
      else:
         return process_mut_helper(m, dbsnp)
   
   def parser(row, dbsnp):
      return parser_helper(row, dbsnp, process_mut, 3, [2])
      
   def fix_karolinska(cypdata):
      cypdata[gene + '*38']['mutations'][0] = (5023, 'T>C', 'rs148166815')
      cypdata[gene + '*27']['mutations'] += [
         [2162, 'delG', 'rs28399445'],
         [2163, 'C>A', 'rs28399445']
      ]
      cypdata.update({
         gene + '*3': { # This one is retention in e3, e6 and e8 (check original paper)
            'mutations': {
               ('pseudogene', 'e3:1-151', ''),
               ('pseudogene', 'e6:1-143', ''),
               ('pseudogene', 'e8:1-143', '')
            }, 
            'phenotype': {'invivo': '', 'invitro': ''},
         },
         gene + '*12A': {
            'mutations': { ('pseudogene', 'e2-', '') },
            'phenotype': {'invivo': 'Decr.', 'invitro': 'Decr.'},
         },
         gene + '*34': {
            'mutations': { ('pseudogene', 'i4-', '') },
            'phenotype': {'invivo': '', 'invitro': ''},
         }
      })

   def coordinate(p, g): 
      p = -p - 5021 + len(g.seq)
      return p + (1 if p < 0 else 0)

   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[8][0]
      g.special_regions.update({
         '0i':  (g.exons[0][1], g.translation[len(g.seq) - 1]),
         '3utr': (gs - 258, gs),
         '3flank': (gs - 2000, gs - 258),
      })
      g = genes[pseudogene]
      gs = g.exons[8][0]
      g.special_regions.update({
         '0i':  (g.exons[0][1], g.translation[len(g.seq) - 1]),
         '3utr': (gs - 542, gs),
         '3flank': (gs - 2000, gs - 542),
      })

      # Manual fixes
      assert(result[gene + '*26']['mutations'][9]['old'] == '1710:C>T')
      result[gene + '*26']['mutations'][9]['functional'] = 0

      result.update({
         'CYP2A6*4': { 
            'mutations': [{'op': 'deletion'}],  
            'phenotype': { 'invivo': 'None', 'invitro': '' }
         },
      })

   hg19_range = ('chr19', 41347500, 41400000)
   return get_data(force, gene, pseudogene, True, parser, fix_karolinska, 
      hg19_range, ['224384750', (gene, 'NG_008377.1'), (pseudogene, 'NG_007960.1')],
      coordinate, [(51, 'C')], post_process,
      functional_exceptions=['-48:T>G', '2163:C>A'],
      unique_regions=[
         '1e', '1i',
         '2e', '2i',
         '3e', '3i',
         '4e', '4i',
         '5e', '5i',
         '6e', '6i',
         '7e', '7i',
         '8e', '8i',
         '9e', 
      ], 
      max_cn={}
   )

def get_cyp2c19(force):
   gene, pseudogene = 'CYP2C19', None
   
   def parser(row, dbsnp):
      return parser_helper(row, dbsnp, process_mut_helper, 3, [2])
      
   def fix_karolinska(cypdata):
      cypdata[gene + '*11']['mutations'] = cypdata[gene + '*11']['mutations'][:3]
      cypdata[gene + '*30']['mutations'] = cypdata[gene + '*30']['mutations'][:1]

      extra_alleles = [m for m in cypdata[gene + '*17']['mutations'] 
         if '{}:{}'.format(m[0], m[1]) in ['80161:A>G']]
      cypdata[gene + '*8X'] = copy.deepcopy(cypdata[gene + '*8'])
      cypdata[gene + '*8X']['mutations'] += extra_alleles

   def coordinate(p, g): 
      p = p + 4999
      return p + (1 if p < 0 else 0)

   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[0][0]
      g.special_regions.update({
         '0i':  (gs - 4000, gs),
      })

   hg19_range = ('chr10', 96445000, 96615000)
   return get_data(force, gene, pseudogene, False, parser, fix_karolinska, 
      hg19_range, ['224589801', (gene, 'NG_008384.2')],
      coordinate, [], post_process,
      functional_exceptions=['19154:G>A', '-806:C>T', '-1041:G>A', '12662:A>G'],
      unique_regions=['1e', '2e', '3e', '4e', '5e', '6e', '7e', '8e', '9e'],
      max_cn={}
   )

def get_cyp2c9(force):
   gene, pseudogene = 'CYP2C9', None
   
   def parser(row, dbsnp):
      return parser_helper(row, dbsnp, process_mut_helper, 3, [2])

   def coordinate(p, g): 
      p = p + 5024
      return p + (1 if p < 0 else 0)

   def fix_karolinska(cypdata):
      cypdata[gene + '*6']['mutations'][0][0] -= 1

   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[0][0]
      g.special_regions.update({
         '0i':  (gs - 5000, gs),
      })

   hg19_range = ('chr10', 96691000, 96754000)
   return get_data(force, gene, pseudogene, False, parser, fix_karolinska, 
      hg19_range, ['224589801', (gene, 'NG_008385.1')],
      coordinate, [], post_process,
      functional_exceptions=[],
      unique_regions=['1e', '2e', '3e', '4e', '5e', '6e', '7e', '8e', '9e'],
      max_cn={}
   )

def get_cyp2c8(force):
   gene, pseudogene = 'CYP2C8', None
   
   def parser(row, dbsnp):
      return parser_helper(row, dbsnp, process_mut_helper, 3, [2])

   def coordinate(p, g): 
      p = -p - 3023 + len(g.seq)
      return p + (1 if p < 0 else 0)
      
   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[8][0]
      g.special_regions.update({
         '0i':  (gs - 1000, gs),
      })

   hg19_range = ('chr10', 96796000, 96830000)
   return get_data(force, gene, pseudogene, True, parser, None, 
      hg19_range, ['224589801', (gene, 'NC_000010.9', {'seq_start': 96822172, 'seq_stop': 96783519, 'strand': 2})],
      coordinate, [], post_process,
      functional_exceptions=[],
      unique_regions=['1e', '2e', '3e', '4e', '5e', '6e', '7e', '8e', '9e'],
      max_cn={}
   )

def get_cyp3a5(force):
   gene, pseudogene = 'CYP3A5', None
   
   def parser(row, dbsnp):
      return parser_helper(row, dbsnp, process_mut_helper, 2, [])

   def fix_karolinska(cypdata):
      del cypdata['Home']
      del cypdata[gene + '*9']['mutations'][1] 
      cypdata[gene + '*9X'] = copy.deepcopy(cypdata[gene + '*9'])
      del cypdata[gene + '*9X']['mutations'][1] 

   def coordinate(p, g): 
      p = -p - 182 + len(g.seq)
      return p + (1 if p < 0 else 0)
      
   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[12][0]
      g.special_regions.update({
         '0i':  (gs - 4500, gs),
      })

   hg19_range = ('chr7', 99245000, 99278000)
   return get_data(force, gene, pseudogene, True, parser, fix_karolinska, 
      hg19_range, ['224589819', (gene, 'NG_000004.3', {'seq_start': 240000, 'seq_stop': 272000, 'strand': 1})],
      coordinate, [(6986, 'T'), (31611, 'G')], post_process,
      functional_exceptions=['6986:A>G', '14690:G>A'],
      unique_regions=['1e', '2e', '3e', '4e', '5e', '6e', '7e', '8e', '9e'],
      max_cn={}
   )

def get_cyp3a4(force):
   gene, pseudogene = 'CYP3A4', None
   
   def parser(row, dbsnp):
      return parser_helper(row, dbsnp, process_mut_helper, 3, [2])

   def fix_karolinska(cypdata):
      #15389C>T
      m = cypdata[gene + '*3']['mutations'][0]
      cypdata[gene + '*3']['mutations'][0] = (m[0] + 1, m[1], m[2])
      m = cypdata[gene + '*22']['mutations'][0]
      cypdata[gene + '*22']['mutations'][0] = (m[0] - 9, m[1], m[2])
      pass

   def coordinate(p, g): 
      p = -p - 2037 + len(g.seq)
      return p + (1 if p < 0 else 0)
      
   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[12][0]
      g.special_regions.update({
         '0i':  (gs - 1500, gs),
      })

   hg19_range = ('chr7', 99354000, 99465000)
   return get_data(force, gene, pseudogene, True, parser, fix_karolinska, 
      hg19_range, ['224589819', (gene, 'AF280107.1', {'seq_start': 60000, 'seq_stop': 90000, 'strand': 1})],
      coordinate, [], post_process,
      functional_exceptions=['15389:C>T'],
      unique_regions=['1e', '2e', '3e', '4e', '5e', '6e', '7e', '8e', '9e'],
      max_cn={}
   )

def get_cyp4f2(force):
   gene, pseudogene = 'CYP4F2', None
   
   def parser(row, dbsnp):
      return parser_helper(row, dbsnp, process_mut_helper, 3, [2])

   def fix_karolinska(cypdata):
      # m = cypdata[gene + '*2']['mutations'][0]
      # cypdata[gene + '*2']['mutations'][0] = (m[0]+1, m[1], m[2])
      pass

   def coordinate(p, g): 
      p = -p - 1663 - 900 - 34 + len(g.seq)
      return p + (1 if p < 0 else 0)
      
   def post_process(genes, result):
      g = genes[gene]
      gs = g.exons[11][0]
      g.special_regions.update({
         '0i':  (gs - 1000, gs),
      })

   hg19_range = ('chr19', 15619000, 16009500)
   return get_data(force, gene, pseudogene, True, parser, fix_karolinska, 
      hg19_range, ['224589810', (gene, 'AF467894.1')],
      coordinate, [], post_process,
      functional_exceptions=[],
      unique_regions=['1e', '2e', '3e', '4e', '5e', '6e', '7e', '8e', '9e'],
      max_cn={}
   )


def main(force=False):
   for n, f in globals().iteritems():
      if not n.startswith('get_') or n == 'get_data': continue
      if n == 'get_cyp2c8': continue
      if n == 'get_cyp4f2': continue

      n = n[4:]
      print n
      result = f(force)
      os.system('mkdir -p output')
      with open('output/' + n + '.yml', 'w') as f:
         print >>f, json.dumps(result, indent=2, cls=NoIndentEncoder)

if __name__ == '__main__':
   main(len(sys.argv) == 2 and sys.argv[1] == '-f')
