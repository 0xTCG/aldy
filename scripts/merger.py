# 786

import sys, copy
from Bio.Seq import Seq
from common import log

NO_CHANGE = 0
PROTEIN_CHANGE = 1
SPLICING = 2
FRAMESHIFT = 3

def is_functional(gene, m, exceptions, reverse_complement, neigh=None):
   original_transcript = ''
   new_transcript = ''

   for e in gene.exons:
      seq = gene.seq[e.start:e.end]
      if reverse_complement:
         original_transcript += seq.reverse_complement()
      else: 
         original_transcript += seq

      # exon site
      if e.start <= m['pos'] < e.end:
         if m['op'][:3] in ['DEL', 'INS']:
            if len(m['op'][4:]) % 3 == 0:
               return PROTEIN_CHANGE
            else:
               return FRAMESHIFT 
         if seq[m['pos'] - e.start] != m['op'][4]:
            log.debug('Merge: Karolinska\'s SNP {} does not match reference {}', m, seq[m['pos'] - e.start])
         s = seq.tomutable()
         s[m['pos'] - e.start] = m['op'][5]

         if neigh != None:
            nt = len(new_transcript)
            codon = (nt + e.end - m['pos']) / 3
            for n in neigh:
               if isinstance(n['pos'], int) and (n != m) and (nt + e.end - n['pos']) / 3 == codon:
                  log.debug('Merge: {}:{} and {}:{} are in same codon'.format(m['pos'], m['op'], n['pos'], n['op']))
                  s[n['pos'] - e.start] = n['op'][5]

         if reverse_complement:
            new_transcript += s.toseq().reverse_complement()
         else:
            new_transcript += s.toseq()
      # splicing site 
      elif e.start - 3 <= m['pos'] < e.start \
         or e.end <= m['pos'] < e.end + 3 \
         or m['old'] in exceptions:
         return SPLICING 
      else:
         if reverse_complement:
            new_transcript += seq.reverse_complement()
         else:
            new_transcript += seq

   original_transcript = original_transcript.translate()
   new_transcript = new_transcript.translate()
   #print original_transcript, new_transcript
   return int(original_transcript != new_transcript)

from pprint import pprint
def merge(cypdata, gene, coordinate, exceptions, reverse_complement=True):
   result = copy.deepcopy(cypdata)

   for a in sorted(result):
      mutations = []
      for mi, m in enumerate(result[a]['mutations']):
         pos, op, db = m
         if op == 'A6': #hack
            log.info('Merge: hack for {}', m)
            continue
         elif pos == 'pseudogene':
            if op[-1] in '+-':
               mutations.append(dict(pos=pos, op=op, dbsnp=[], old=''))
            else:
               parts = op.split(':') # do we have ranges?
               if len(parts) == 2:
                  start, end = tuple(map(lambda x: coordinate(int(x), gene), parts[1].split('-'))[::-1])
               else:
                  if op[0] == 'e':
                     r = gene.exons[int(op[1:]) - 1]
                  else:
                     r = gene.introns[int(op[1:]) - 1]
                  start, end = r.start, r.end
               mutations += [copy.deepcopy(gene.pseudo_mutations[i]) 
                  for i in xrange(start, end) if i in gene.pseudo_mutations] 
            continue
         if pos < 0: 
            pos += 1
         if op[:3] == 'ins':  
            #pos += 2
            if reverse_complement:
               op = 'INS.{}'.format(Seq(op[3:]).reverse_complement().lower())
            else:
               op = 'INS.{}'.format(Seq(op[3:]).lower())
         elif op[:3] == 'del': 
            if reverse_complement:
               op = 'DEL.{}'.format(Seq(op[3:]).reverse_complement())
            else:
               op = 'DEL.{}'.format(Seq(op[3:]))
            pos += len(op) - 5
         else:
            if reverse_complement:
               op = 'SNP.{}{}'.format(Seq(op[0]).reverse_complement(), Seq(op[2]).reverse_complement())
            else:
               op = 'SNP.{}{}'.format(Seq(op[0]), Seq(op[2]))
         mutations.append(dict(pos=coordinate(pos, gene), op=op, dbsnp=db, old='{}:{}'.format(m[0],m[1])))
      for m in mutations:
         m['functional'] = is_functional(gene, m, exceptions, reverse_complement, mutations)
      result[a]['mutations'] = mutations
   return result
