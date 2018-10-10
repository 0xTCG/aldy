# 786

import os, sys
from cStringIO import StringIO

from Bio import Entrez, SeqIO, SeqFeature, pairwise2
from Bio.Seq import Seq

from common import log
from collections import namedtuple as nt

Gene = nt('Gene', ['name', 'protein', 'introns', 'exons', 'special_regions', 'seq', 'translation', 'pseudo_mutations', 'pseudo_translation'])

def blat(a, b):
   from Bio import SearchIO

   with open('_a', 'w') as f:
      print >>f, '>A\n' + a 
   with open('_b', 'w') as f:
      print >>f, '>B\n' + b
   # If using OS X, make sure to compile BLAT with GNU g++; clang version does not produce correct results
   ret = os.system('./blat _a _b _out -out=pslx -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -fine >/dev/null')
   if ret != 0:
      raise 'BLAT failed'
   with open('_out') as f:
      ret = f.read()
   os.unlink('_a')
   os.unlink('_b')
   os.unlink('_out')

   ret = SearchIO.read(StringIO(ret), 'blat-psl', pslx=True)
   hsp = max(ret.hsps, key=lambda h: h.score)
   return hsp 

def get_pseudo_mutations(gene, pseudogene, force=False):
   mutations = {}
   translation = {}
   # print gene.exons
   # print pseudogene.exons
   # print gene.introns
   # print pseudogene.introns
   # exit(0)
   for ei, (e6, e7) in enumerate(zip(gene.exons, pseudogene.exons) + zip(gene.introns, pseudogene.introns)):
      y6, y7 = e6.start, e7.start
      e6 = gene.seq[e6.start:e6.end]
      e7 = pseudogene.seq[e7.start:e7.end]
      
      
      def yy(x,s=10):
         return ' '.join(x[i:i + s] for i in xrange(0, len(x), s))
      log.debug('NCBI: ALN {}{} (len {} / {})', 
         'E' if ei < len(gene.exons) else 'I',
         ei + 1 if ei < len(gene.exons) else ei - len(gene.exons) + 1,
         len(e6), len(e7))
      if max(len(e6), len(e7)) < 1000:
         al = pairwise2.align.globalxs(e6, e7, -1, 0)[0]
         s = str(pairwise2.format_alignment(*al)).split('\n')
      else:
         log.debug('BLAT')
         a = blat(e6, e7)

         s = [a[0].query.seq, '', a[0].hit.seq]
         for i, f in enumerate(a):
             if i == 0: continue
             if a[i].query_range[0] == a[i - 1].query_range[1]:
                 s[0] += e6[a[i - 1].hit_range[1]:a[i].hit_range[0]]
                 s[2] += '-' * (a[i].hit_range[0] - a[i - 1].hit_range[1])
             if a[i].hit_range[0] == a[i - 1].hit_range[1]:
                 s[0] += '-' * (a[i].query_range[0] - a[i - 1].query_range[1])
                 s[2] += e7[a[i - 1].query_range[1]:a[i].query_range[0]]
             s[0] += a[i].query.seq
             s[2] += a[i].hit.seq
         s = map(str.upper, map(str, s))
         assert(len(s[0])==len(s[2]))

      log.debug('NCBI: ALN {}', yy(s[0]))
      log.debug('NCBI: ALN {}', yy(''.join(['*' if s[0][y]!=s[2][y] else '-' for y in xrange(len(s[0]))])))
      log.debug('NCBI: ALN {}', yy(s[2]))

      gaps6, gaps7 = 0, 0
      i = 0
      while i < len(al[0]):
         c6 = y6 + i - gaps6
         c7 = y7 + i - gaps7
         if c6 not in gene.translation:
            # print >>sys.stderr, 'CYP2D7 mutation {}:{}{} ignored'.format(y7 + i, a, b)
            i += 1
            continue
         # print al[0][i]
         if al[0][i] == '-': 
            seq = ''
            while i < len(al[0]) and al[0][i] == '-': 
               seq += al[1][i]
               i += 1
               gaps6 += 1
            # (deletion in 6 is actually insertion in mapping)
            c6 = y6 + i - gaps6 # INS je ispred!
            translation[c6] = mutations[c6] = dict(
               pos=c6, 
               op='INS.{}'.format(seq.lower()), 
               dbsnp='*',
               old='{}:{}{}:{}'.format(pseudogene.name, 'e' if ei < len(gene.exons) else 'i', 
                                               ei + 1 if ei < len(gene.exons) else ei - len(gene.exons) + 1,  
                                               c7),
               old_pos=c7
            )
            continue
         if al[1][i] == '-': 
            seq = ''
            while i < len(al[0]) and al[1][i] == '-': 
               seq += al[0][i]
               i += 1
               gaps7 += 1
            # (deletion in 7 is actually deletion in mapping)
            translation[c6] = mutations[c6] = dict(
               pos=c6, 
               op='DEL.{}'.format(seq), 
               dbsnp='*',
               old='{}:{}{}:{}'.format(pseudogene.name, 'e' if ei < len(gene.exons) else 'i', 
                                               ei + 1 if ei < len(gene.exons) else ei - len(gene.exons) + 1,  
                                               c7),
               old_pos=c7
            )
            continue
         if al[0][i] != al[1][i]:
            translation[c6] = mutations[c6] = dict(
               pos=c6, 
               op='SNP.{}{}'.format(al[0][i], al[1][i]), 
               dbsnp='*',
               old='{}:{}{}:{}'.format(pseudogene.name, 'e' if ei < len(gene.exons) else 'i', 
                                               ei + 1 if ei < len(gene.exons) else ei - len(gene.exons) + 1,  
                                               c7),
               old_pos=c7
            )
         elif c6 not in mutations: # do not overwrite insertions
            translation[c6] = dict(old_pos=c7)
         i += 1   

   return mutations, translation

def get_genomes(genome_id, genome_region, gene_ids, reverse_complement=True, entrez_mail='test@test.ca', force=False):
   Entrez.email = entrez_mail
   chromosome, start, end = genome_region 
   # NCBI uses 1 based indexing and closed intervals [a,b]
   handle = Entrez.efetch(db='nucleotide', 
                  id=genome_id,
                  rettype='fasta', 
                  strand=1, 
                  seq_start=start + 1,
                  seq_stop=end + 1 + 1)
   record = SeqIO.read(handle, 'fasta')
   hg19 = record.seq

   genomes = {}
   handle = Entrez.read(Entrez.esearch(db='nucleotide', term=' '.join(g[1] for g in gene_ids), retmode='xml'))
   for gi, gid in enumerate(handle['IdList']):
      params = {}
      if len(gene_ids[gi]) > 2:
         params = gene_ids[gi][2]
      genome = Entrez.efetch(db='nucleotide', id=gid, rettype='gb', retmode='text', **params).read()
      genome = SeqIO.read(StringIO(genome), 'genbank')
      
      if reverse_complement:
         genome.seq = genome.seq.reverse_complement()
      alignment = blat(hg19, genome.seq)
      
      log.trace('NCBI: Gene {} BLAT results: hit {}, query {}', genome.id, alignment.hit_range, alignment.query_range)
      translation = dict((i[0], i[1] + start)
         for f in alignment 
         for i in zip(range(*f.query_range), range(*f.hit_range))
      )
      cds = [c for c in genome.features if c.type == 'CDS']
      if len(cds) == 0:
         cds = [c for c in genome.features if c.type == 'misc_RNA']
      for cd in cds:
         protein = ''
         if 'translation' in cd.qualifiers:
            protein = cd.qualifiers['translation']
         
         if reverse_complement:
            exons = [SeqFeature.FeatureLocation(len(genome.seq) - e.end, len(genome.seq) - e.start, 1) for e in cd.location.parts]
            introns = [SeqFeature.FeatureLocation(e2.end, e1.start, 1) for e1, e2 in zip(exons[:-1], exons[1:])]
         else:
            exons = [SeqFeature.FeatureLocation(e.start, e.end, 1) for e in cd.location.parts]
            introns = [SeqFeature.FeatureLocation(e1.end, e2.start, 1) for e1, e2 in zip(exons[:-1], exons[1:])]
         
         genomes[cd.qualifiers['gene'][0]] = Gene(
            name=cd.qualifiers['gene'][0],
            protein=protein, 
            introns=introns, 
            exons=exons, 
            seq=genome.seq,
            translation=translation,
            pseudo_mutations={},
            pseudo_translation={},
            special_regions={})
   
   if len(gene_ids) > 1:
      g, p = gene_ids[0][0], gene_ids[1][0]
      p, pt = get_pseudo_mutations(genomes[g], genomes[p], force)
      genomes[g].pseudo_mutations.update(p)
      genomes[g].pseudo_translation.update(pt)

   return genomes, hg19
