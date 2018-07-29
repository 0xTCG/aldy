# 786

# Aldy source: genotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.

from __future__ import print_function
from __future__ import division
from builtins import str

import os
import re
import sys
import logbook
import logbook.more
import collections
import platform
import textwrap
import pysam
import subprocess
import tempfile
import shutil
from Bio import SeqIO
from collections import defaultdict
from collections import OrderedDict
from gurobipy import *

from . import cn
from . import gene
from . import sam
from . import lpinterface
from . import diplotype

from .common import *


def remap(sam_path, gene, sam, cn_sol, tempdir, hg19_path="genome.fa"):
   samtools_path = 'samtools'
   bowtie2_path = ''

   infile = sam_path
   outbamfile = '{}.remap.bam'.format(os.path.basename(sam_path))

   avg_coverage = sum(sam.total(pos) for pos in sam.coverage) / float(len(sam.coverage))
   start = gene.region[1]

   readlength = 0
   samfile = pysam.AlignmentFile(infile, 'r')
   
   chromosome = ''
   chrs = [x['SN'] for x in samfile.header['SQ']]
   if 'chr1' not in chrs and '1'in chrs:
       chromosome = '22'
   elif 'chr1' in chrs and '1' not in chrs:
      chromosome = 'chr22'
   
   for segment in samfile.fetch():
       readlength = len(segment.query_sequence)
       break
   samfile.close()
   assert readlength != 0, 'Read Length is 0'

   def read_genome(filename):
       genome = ''
       with open(filename, 'r') as f:
           for line in f:
               if not line[0] == '>':
                   genome += line.rstrip()
       return genome

   os.system('{} faidx {} {}:{}-{} > {tmp}/seq.fa'.format(samtools_path, hg19_path,
                                                           'chr22',
                                                           int(gene.region[1])+1,
                                                           int(gene.region[2])+2, tmp=tempdir))
   os.system('sed -i.bak "1 s/^.*$/>{}/" {tmp}/seq.fa'.format(chromosome, tmp=tempdir))
   seq = read_genome('{tmp}/seq.fa'.format(tmp=tempdir))

   cnv = cn_sol.region_cn

   regions = {'2D6':[], '2D7':[]}
   cyp2d6 = {}
   cyp2d7 = {}
   for k, v in cnv.iteritems():
      if v == 0: continue
      if '6.' in k:
         cyp2d6[k] = v
      elif '7.' in k:
         cyp2d7[k] = v

   def get_regions_helper(cyp2d, which):
      keys = sorted(cyp2d)
      a = keys[0]
      b = cyp2d[a]
      merges = [gene.regions[a]]
      for i in range(1, len(keys)):
         k = keys[i]
         v = cyp2d[k]
         if b == v:
            merges.append(gene.regions[k])
         else:
            regions[which].insert(0, (merges[-1][0], merges[0][1], b))
            merges = [gene.regions[k]]
            a = k
            b = v
         if i == len(keys)-1:
            if b == v:
               regions[which].insert(0, (merges[-1][0], merges[0][1], b))
            else:
               regions[which].insert(0, (gene.regions[k][0], gene.regions[k][1], v))

   get_regions_helper(cyp2d6, '2D6')
   get_regions_helper(cyp2d7, '2D7')

   for which in regions.keys():
      for i in range(len(regions[which])):
         if i == 0:
            regions[which][i] = (regions[which][i][0], 
                            regions[which][i][1]+int(readlength/2), 
                            regions[which][i][2])
         elif i == len(regions[which])-1:
            regions[which][i] = (regions[which][i][0]-int(readlength/2), 
                            regions[which][i][1], 
                            regions[which][i][2])
         else:
            regions[which][i] = (regions[which][i][0]-int(readlength/2), 
                            regions[which][i][1]+int(readlength/2), 
                            regions[which][i][2])


   def create_reads(total):
      command = 'cat'
      for i in range(len(total)):
         os.system('{} view -h {} {}:{}-{} > {tmp}/temp{}.sam'.format(samtools_path, infile, 
                                                        chromosome,
                                                        total[i][0], 
                                                        total[i][1], i, tmp=tempdir))
         os.system('{} fastq {tmp}/temp{}.sam > {tmp}/read{}.fq'.format(samtools_path, i, i, tmp=tempdir))
         command += ' {tmp}/read{}.fq'.format(i, tmp=tempdir)
      command += ' > {tmp}/read.fq'.format(tmp=tempdir)
      os.system(command)

   def align(regions, infile):
      def write_genome(name, genome):
         breakline = 70
         with open('{tmp}/{}.fa'.format(name, tmp=tempdir), 'w') as f:
            f.write('>'+name+'\n')
            for i in range(len(genome)/breakline+1):
               f.write(genome[i*breakline:i*breakline+breakline]+'\n')

      def merge_fasta(keys, filename):
         with open('{tmp}/{}.fa'.format(filename, tmp=tempdir), 'w') as outfile:
            for fname in keys:
               with open('{tmp}/{}.fa'.format(fname, tmp=tempdir)) as infile:
                  for line in infile:
                     outfile.write(line)

      def merge_sam_bowtie2(keys, filename):
         first = "@HD\tVN:1.0\tSO:unsorted\n"
         with open('{tmp}/{}.sam'.format(filename, tmp=tempdir), 'w') as outfile:
            for fname in keys:
               with open('{tmp}/{}.sam'.format(fname, tmp=tempdir)) as infile:
                  i = 0
                  for line in infile:
                     i += 1
                     if i >= 4:
                        outfile.write(line)
                     elif i == 2:
                        first += line
         with open('{tmp}/{}.sam'.format(filename, tmp=tempdir), 'r+') as outfile:
            file_data = outfile.read()
            outfile.seek(0, 0)
            outfile.write(first + file_data)

      keys = regions.keys()
      write_genome(keys[0], seq[regions[keys[0]][0]-start:regions[keys[0]][1]-start])
      write_genome(keys[1], seq[regions[keys[1]][0]-start:regions[keys[1]][1]-start])

      merge_fasta(keys, 'candidate')

      for k in keys:
         os.system('{}bowtie2-build {tmp}/{}.fa {tmp}/{} > /dev/null 2>&1'.format(bowtie2_path, k, k, tmp=tempdir))
         os.system('{}bowtie2 -x {tmp}/{} {tmp}/read.fq -S {tmp}/{}.sam'.format(bowtie2_path, k, k, tmp=tempdir))
      merge_sam_bowtie2(keys, 'candidate')
      
      os.system('{} sort {tmp}/candidate.sam > {tmp}/candidate.bam'.format(samtools_path, tmp=tempdir))
      os.system('{} index {tmp}/candidate.bam'.format(samtools_path, tmp=tempdir))

   def preprocess():
      class SNP(object):
         def __init__(self):
            self.kind = ''
            self.base_ref = ''
            self.base_read = ''
            self.pos_ref = -1
            self.pos_read = -1

         def __repr__(self):
            return '{}.{}{}-{},{}'.format(self.kind, self.base_ref, self.base_read, self.pos_ref, self.pos_read)

      class Read(object):
         def __init__(self):
            self.ref = ''
            self.name = ''
            self.seq = ''
            self.start = -1
            self.end = -1
            self.snp = []

         def __str__(self):
            return '{}, {}, ({}-{}), {}'.format(self.ref, self.name, self.start, self.end, self.snp)

      def data(file_bam, file_fa):
         all_reads = []
         score_read_ref = {}
         which_read_at = defaultdict(list)
         reads = defaultdict(list)
         samfile = pysam.AlignmentFile(file_bam, 'r')
         fasta_parser = open(file_fa, 'r')
         seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_parser, 'fasta'))
         alleles = {}
         for x in samfile.header['SQ']:
            alleles[x['SN']] = x['LN']
         for segment in samfile.fetch():
            nm = dict(segment.tags)['NM']
            r = Read()
            r.ref = segment.reference_name
            r.name = segment.query_name
            r.seq = segment.query_sequence
            r.start = segment.reference_start
            r.end = segment.reference_end

            start, s_start = segment.reference_start, 0
            chr_start = 0
            seq = segment.query_sequence
            ref = str(seq_dict[segment.reference_name][:].seq)
            for op, size in segment.cigartuples:
               if op == 2:
                  s = SNP()
                  s.base_ref = ref[start - chr_start:start - chr_start + size]
                  s.pos_ref = start
                  s.pos_read = s_start
                  s.kind = 'DEL'
                  r.snp.append(s)
                  start += size
               elif op == 1:
                  s = SNP()
                  s.base_read = seq[s_start:s_start + size]
                  s.pos_ref = start
                  s.pos_read = s_start
                  s.kind = 'INS'
                  r.snp.append(s)
                  s_start += size
               elif op == 4:
                  s_start += size
               elif op in [0, 7, 8]:
                  for i in range(size):
                     if 0 <= start + i - chr_start < len(ref) and ref[start + i - chr_start] != seq[s_start + i]:
                        s = SNP()
                        s.base_ref = ref[start + i - chr_start]
                        s.base_read = seq[s_start + i]
                        s.pos_ref = start + i
                        s.pos_read = i
                        s.kind = 'SNP'
                        r.snp.append(s)
                  start += size
                  s_start += size
            score_read_ref[segment.query_name, segment.reference_name] = int(len(segment.query_sequence)-nm)
            reads[segment.reference_name].append(r)
            all_reads.append(segment.query_name)
            for i in range(r.start, r.end):
               which_read_at[segment.reference_name, i].append(r.name)

         samfile.close()
         fasta_parser.close()
         return alleles, reads, score_read_ref, which_read_at, list(OrderedDict.fromkeys(all_reads))
      
      file_bam = '{tmp}/candidate.bam'.format(tmp=tempdir)
      file_fa = '{tmp}/candidate.fa'.format(tmp=tempdir)
      alleles, reads, score_read_ref, which_read_at, all_reads = data(file_bam, file_fa)
      return all_reads, alleles, score_read_ref, which_read_at


   def optimization(all_reads, alleles, score_read_ref, which_read_at, e, r, n, a):
      def _e(k, g, r):
         if k < r:
            return 1.0*e/(r-k)
         elif g-r < k:
            return 1.0*e/(g-k)
         else:
            return e

      all_alleles = alleles.keys()
      model = Model()

      x = {}
      for i in range(len(all_reads)):
         for j in range(len(all_alleles)):
            x[all_reads[i],all_alleles[j]] = model.addVar(vtype=GRB.BINARY, name='x[{},{}]'.format(all_reads[i], all_alleles[j]))

      d = {}
      for i in range(len(all_reads)):
         d[all_reads[i]] = model.addVar(vtype=GRB.BINARY, name='d[{}]'.format(all_reads[i]))

      y = {}
      z = {}
      for g, v in alleles.items():
         for k in range(1,v):
            z[g, k] = model.addVar(vtype=GRB.CONTINUOUS, name='z[{},{}]'.format(g, k), lb=-GRB.INFINITY)
            model.addConstr(z[g, k] == (quicksum([x[i, g] for i in which_read_at[g, k]])-n[g]*_e(k, v, r)))
            y[g, k] = model.addVar(vtype=GRB.CONTINUOUS, name='y[{},{}]'.format(g, k), lb=0)
            model.addGenConstrAbs(y[g, k], z[g, k])


      #Objective
      model.setObjective(quicksum([score_read_ref.get((all_reads[i], all_alleles[j]), 0)*x[all_reads[i], all_alleles[j]] 
                            for i in range(len(all_reads)) for j in range(len(all_alleles))]) -
                     quicksum([y[g, k] for g, v in alleles.items() for k in range(1,v)]),
                     GRB.MAXIMIZE)

      #Constraints_1
      for i in range(len(all_reads)):
         for j in range(len(all_alleles)):
            model.addConstr(x[all_reads[i], all_alleles[j]] * n[all_alleles[j]] >= x[all_reads[i], all_alleles[j]])

      #Constraints_2
      for i in range(len(all_reads)):
         model.addConstr((d[all_reads[i]] + quicksum([x[all_reads[i], all_alleles[j]] for j in range(len(all_alleles))])) == 1)

      model.params.outputFlag = 0
      model.params.logFile = ''
      model.optimize()

      out = open('{tmp}/output{}.txt'.format(a, tmp=tempdir), 'w')
      all_alleles = alleles.keys()
      for i in range(len(all_reads)):
         if model.getVarByName('d[{}]'.format(all_reads[i])).X > 0:
            out.write('{} is not kept\n'.format(all_reads[i]))
         for j in range(len(all_alleles)):
            temp = model.getVarByName('x[{},{}]'.format(all_reads[i], all_alleles[j]))
            if temp.X > 0:
               name = temp.VarName.replace('x[','').replace(']','').split(',')
               out.write('{},{}\n'.format(all_reads[i],all_alleles[j]))
      out.close()



   def produceBam(i, startdict):
      results = []
      with open('{tmp}/output{}.txt'.format(i, tmp=tempdir), 'r') as f:
         for line in f:
            if 'not kept' not in line:
               results.append((line.split(',')[0], line.rstrip().split(',')[1]))

      samfile = pysam.AlignmentFile(infile, 'r')
      header = samfile.header.copy()
      ref_id = 0
      for h in range(len(header['SQ'])):
         if '22' in header['SQ'][h]['SN']:
            ref_id = h
            break
      samfile = pysam.AlignmentFile('{tmp}/candidate.bam'.format(tmp=tempdir), 'r')
      name_indexed = pysam.IndexedReads(samfile)
      name_indexed.build()
      outfile = pysam.Samfile('{tmp}/output{}.sam'.format(i, tmp=tempdir), 'w', header=header)
      for result in results:
         try:
            iterator = name_indexed.find(result[0])
            for x in iterator:
               if x.reference_name == result[1]:
                  a = pysam.AlignedSegment()
                  a.query_name = x.query_name
                  a.query_sequence = x.query_sequence
                  a.flag = x.flag
                  # a.reference_name = 'chr22'
                  a.reference_id = ref_id
                  # a.next_reference_id = x.next_reference_id
                  a.reference_start = x.reference_start + startdict[result[1]]
                  a.mapping_quality = x.mapping_quality 
                  a.cigar = x.cigar
                  a.next_reference_start = x.next_reference_start
                  a.template_length = x.template_length
                  a.query_qualities = x.query_qualities
                  a.tags = x.tags
                  outfile.write(a)
         except:
            pass
      outfile.close()


   def mergeBam(n):
      command = 'cat '
      for i in range(n):
         os.system('{} sort {tmp}/output{}.sam > {tmp}/output{}.bam'.format(samtools_path, i, i, tmp=tempdir))
         os.system('{} index {tmp}/output{}.bam'.format(samtools_path, i, tmp=tempdir))
         command += '{tmp}/output{}.txt '.format(i, tmp=tempdir)
      os.system(command+'> {tmp}/output.txt'.format(tmp=tempdir))
      
      results = []
      with open('{tmp}/output.txt'.format(tmp=tempdir), 'r') as f:
         for line in f:
            if 'not kept' not in line:
               results.append(line.rstrip().split(',')[0])
      
      multiple = list(set([x for x in results if results.count(x) <= 1]))

      bamfile = pysam.AlignmentFile('{tmp}/output0.bam'.format(tmp=tempdir), 'r')
      header = bamfile.header.copy()
      out = pysam.Samfile('{tmp}/output.sam'.format(tmp=tempdir), 'w', header=header)
      for i in range(n):
         bamfile = pysam.AlignmentFile('./temp/output{}.bam'.format(i, tmp=tempdir), 'r')
         name_indexed = pysam.IndexedReads(bamfile)
         name_indexed.build()
         header = bamfile.header.copy()
         for name in multiple:
            try:
               name_indexed.find(name)
            except KeyError:
               pass
            else:
               iterator = name_indexed.find(name)
               for x in iterator:
                  # out.write(x)
                  a = pysam.AlignedSegment()
                  a.query_name = x.query_name.split('/')[0]
                  a.query_sequence = x.query_sequence
                  a.flag = x.flag
                  a.reference_id = x.reference_id
                  a.reference_start = x.reference_start
                  a.mapping_quality = x.mapping_quality 
                  a.cigar = x.cigar
                  a.next_reference_start = x.next_reference_start
                  a.template_length = x.template_length
                  a.query_qualities = x.query_qualities
                  a.tags = x.tags
                  out.write(a)
      out.close()


   def merge_in_out():
      os.system('{} sort {tmp}/output.sam > {tmp}/output.bam'.format(samtools_path, tmp=tempdir))
      os.system('{} index {tmp}/output.bam'.format(samtools_path, tmp=tempdir))
      reads_out = []
      samfile = pysam.AlignmentFile('{tmp}/output.bam'.format(tmp=tempdir), 'r')
      for segment in samfile.fetch():
         reads_out.append(segment.query_name)

      reads_in = []
      samfile = pysam.AlignmentFile(infile, 'r')
      for segment in samfile.fetch():
         reads_in.append(segment.query_name)

      results = list(set(reads_in)-set(reads_out))

      name_indexed = pysam.IndexedReads(samfile)
      name_indexed.build()
      header = samfile.header.copy()
      out = pysam.Samfile('{tmp}/outputtemp.sam'.format(tmp=tempdir), 'w', header=header)
      for name in results:
         try:
            name_indexed.find(name)
         except KeyError:
            pass
         else:
            iterator = name_indexed.find(name)
            for x in iterator:
               out.write(x)
      out.close()

      os.system('{} sort {tmp}/outputtemp.sam > {tmp}/outputtemp.bam'.format(samtools_path, tmp=tempdir))
      os.system('{} merge -f {} {tmp}/outputtemp.bam {tmp}/output.bam'.format(samtools_path, outbamfile, tmp=tempdir))
      os.system('{} index {}'.format(samtools_path, outbamfile))

   total = [(regions['2D6'][0][0], regions['2D6'][-1][1]), 
          (regions['2D7'][0][0], regions['2D7'][-1][1])]
   create_reads(total)
   for a in range(len(regions['2D6'])):
      test = {
         '2D6': regions['2D6'][a],
         '2D7': regions['2D7'][a]
      }
      align(test, infile)
      
      all_reads, alleles, score_read_ref, which_read_at = preprocess()
      n = {
         '2D6' : regions['2D6'][a][2],
         '2D7' : regions['2D7'][a][2]
      }
      optimization(all_reads, alleles, score_read_ref, which_read_at, int(avg_coverage/2), readlength, n, a)
      
      startdict = {
         '2D6' : regions['2D6'][a][0],
         '2D7' : regions['2D7'][a][0]
      }
      produceBam(a, startdict)
      
   mergeBam(len(regions['2D6']))
   merge_in_out()
   
   return outbamfile
   