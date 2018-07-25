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

from . import cn
from . import gene
from . import sam
from . import lpinterface
from . import diplotype

from .common import *


def cmd(cmd):
   log.debug('Executing ' + cmd)
   return subprocess.check_output(cmd, shell=True)


def optimize(gene, reads, alleles, copy_number):
   model = lpinterface.model('reads', 'gurobi')

   # Binary variables for reads
   x = collections.defaultdict(dict)
   for r in reads:
      for a in reads[r]:
         x[r][a] = model.addVar(vtype='B', name='x_{}_{}'.format(a, r), update=False)

   # Should we include the read or not?
   for r in reads:
      x[r][-1] = model.addVar(vtype='B', name='d_{}'.format(r), update=False)

   # Coverage constraints
   z = collections.defaultdict(int)
   for r in reads:
      for a in reads[r]:
         read = reads[r][a][0]
         for k in range(read.reference_start, read.reference_end):
            z[a, k] += x[r][a]
   for (a, k), lhs in z.items():
      a_start = alleles[a][0][0]
      region = gene.region_at[a_start + k]
      assert(region != '')
      z[a, k] = model.addVar(lb=-model.INF, name='z_{}_{}'.format(a, k), update=False)
      model.addConstr(z[a, k] == lhs - copy_number[region] * copy_number.baseline[a_start + k])
   model.update()

   # Objective
   def score(r):
      if not r[0].has_tag('NM'):
         return 0
      return len(r[0].query_sequence) - r[0].get_tag('NM')
   obj = model.quicksum(score(reads[r][a]) * x[r][a] for r in reads for a in reads[r])
   abss = model.abssum(z.values())
   obj -= abss

   # Existence constraints
   for r in reads:
      for a in reads[r]:
         region = gene.region_at[alleles[a][0][0] + reads[r][a][0].reference_start]
         model.addConstr(x[r][a] * copy_number[region] >= x[r][a])
      model.addConstr(model.quicksum(x[r].values()) == 1)

   # Solve
   def params(m):
      # m.params.timeLimit = 60
      m.params.outputFlag = 1
   model.solve(objective=obj, method='max', init=params)

   # Remove bad reads
   result = {}
   deleted = 0
   for r, xr in x.items():
      if model.getValue(xr[-1]):
         deleted += 1
         continue
      correct_allele = [a for a in xr if model.getValue(xr[a])]
      assert(len(correct_allele) == 1)
      correct_allele = correct_allele[0]
      result[r] = reads[r][correct_allele]
   log.warn('Discarded {} reads'.format(deleted))

   return result


def write_reads(sam_path, gene, alleles, reads, out_path):
   y = list(reads.values())
   y = sorted(y, key=lambda r: (alleles[r[0].reference_name][0][0] + r[0].reference_start, r[0].query_name))

   with pysam.AlignmentFile(sam_path) as sam:
      with pysam.AlignmentFile(out_path, "wb", template=sam) as out:
         for x, _ in y:
            a = pysam.AlignedSegment()
            a.query_name = x.query_name.split('/')[0]
            a.query_sequence = x.query_sequence
            a.flag = x.flag
            a.reference_id = sam.get_tid('chr22')
            a.reference_start = x.reference_start + alleles[x.reference_name][0][0]
            a.mapping_quality = x.mapping_quality 
            a.cigar = x.cigar
            a.next_reference_start = x.next_reference_start
            a.template_length = x.template_length
            a.query_qualities = x.query_qualities
            a.tags = x.tags
            out.write(a)
         cnv_chromosome, cnv_start, cnv_end = gene.cnv_region
         region = 'chr{}:{}-{}'.format(cnv_chromosome, cnv_start - 500, cnv_end + 1)
         for read in sam.fetch(region=region):
            out.write(read)
   cmd('samtools index {}'.format(out_path))


def remap(sam_path, gene, sam, cn_sol, tempdir=None, force=True, cleanup=True):
   if tempdir is None:
      tempdir = tempfile.mkdtemp(suffix='aldy', dir='.')

   def sequence(gene, id, pad=0):
      regs = [b for a, b in gene.regions.items() if a.startswith(str(id) + '.')]
      mi, ma = min(r[0] for r in regs), max(r[1] for r in regs)
      st = gene.region[1]
      return (mi, ma), gene.seq[max(0, mi - st - pad):min(ma - st + pad, len(gene.seq))]

   alleles = {
      '6': sequence(gene, '6'),
      '7': sequence(gene, '7'),
   }

   def makeref(gene):
      for i in '67':
         with open('{}/{}.fa'.format(tempdir, i), 'w') as f:
            print('>{}'.format(i), file=f)
            print('\n'.join(textwrap.wrap(alleles[i][1])), file=f)
         out = cmd('bowtie2-build {0}/{1}.fa {0}/{1}'.format(tempdir, i))
         log.trace(out)
   log.warn('Creating bowtie2 reference')
   makeref(gene)

   def realign(path, regions):
      newpath = '{}/{}.fq'.format(tempdir, os.path.basename(path))
      cmd('samtools view -h {} {} | samtools fastq - > {}'.format(path, ' '.join(regions), newpath))
      for i in '67':
         out = cmd('bowtie2 -x {0}/{1} {2} -S {2}.{1}.sam'.format(tempdir, i, newpath))
         log.debug(out)
   log.warn('Aligning reads via bowtie2')
   regions = ['chr22:{}-{}'.format(*alleles[i][0]) for i in alleles]
   realign(sam_path, regions)

   log.warn('Reading realinged data...')
   reads = collections.defaultdict(dict) # read -> allele
   which_read_at = collections.defaultdict(list)
   for ref_id in '67':
      ref = alleles[ref_id][1]
      with pysam.AlignmentFile('{}/{}.fq.{}.sam'.format(tempdir, os.path.basename(sam_path), ref_id)) as sam:
         for read in sam.fetch():
            start = read.reference_start
            start, s_start = read.reference_start, 0
            if not read.cigartuples or read.reference_name != ref_id:
               continue

            seq = read.query_sequence
            rname = read.query_name
            muts = list()
            
            for op, size in read.cigartuples:
               if op == 2:
                  mut = (start, 'DEL.{}'.format(ref[start:start + size]))
                  muts.append(mut)
                  start += size
               elif op == 1:
                  mut = (start, 'INS.{}'.format(seq[s_start:s_start + size].lower()))
                  muts.append(mut)
                  s_start += size
               elif op == 4:
                  s_start += size
               elif op in [0, 7, 8]:
                  for i in range(size):
                     if 0 <= start + i < len(ref) and ref[start + i] != seq[s_start + i]:
                        mut = (start + i, 'SNP.{}{}'.format(ref[start + i], seq[s_start + i]))
                        muts.append(mut)
                  start += size
                  s_start += size
            reads[rname][ref_id] = (read, muts)
   
   log.warn('Optimizing...')
   reads = optimize(gene, reads, alleles, cn_sol)

   out = '{}.remap.bam'.format(os.path.basename(sam_path)) 
   if os.path.exists(out) and not force:
      raise AldyException('{} already exists--- make sure to use --force!'.format(out))
   write_reads(sam_path, gene, alleles, reads, out)

   if cleanup:
      shutil.rmtree(tempdir)
   return out



   