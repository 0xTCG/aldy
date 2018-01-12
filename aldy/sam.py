#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 786

# Aldy source: sam.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import range # noqa
from builtins import object # noqa
from pprint import pprint # noqa

import pysam # noqa
import os # noqa
import copy # noqa
import subprocess # noqa
import tempfile # noqa
import collections # noqa
import pickle as pickle # noqa

from .common import * # noqa
from .gene import Mutation # noqa


class SAM(object):
	CACHE = False
	PHASE = False
	PROFILE = False

	def get_cache(self, path, gene):
		return path + '-{}.aldycache'.format(gene.name)

	def __init__(self, sam_path, gene, threshold):
		if SAM.PROFILE:
			self.load_profile(sam_path, threshold)
			exit(0)

		self.coverage = dict()
		self.cnv_coverage = collections.defaultdict(int)
		self.region_coverage = dict()
		self.sequence = ''
		self.unfiltered_coverage = collections.defaultdict(int)
		
		if SAM.CACHE and os.path.exists(self.get_cache(sam_path, gene)):
			log.debug('Loading cache')
			d = pickle.load(open(self.get_cache(sam_path, gene)))
			self.__dict__.update(d)
		else:
			self.load(sam_path, gene)
			if SAM.CACHE:
				pickle.dump(self.__dict__, open(self.get_cache(sam_path, gene), 'wb'))

		avg_coverage = sum(self.total(pos) for pos in self.coverage) / float(len(self.coverage))
		log.debug('Coverage is {}', avg_coverage)
		self.threshold = threshold

	def total(self, pos):
		return float(sum(v for p, v in self[pos].items() if p[:3] != 'INS'))

	def percentage(self, m):
		total = self.total(m.pos)
		if total == 0:
			return 0
		return 100.0 * float(self[m]) / total

	def __getitem__(self, index):
		if isinstance(index, Mutation):
			if index.pos not in self.coverage:
				self.coverage[index.pos] = {}
			op = index.op
			if op[:3] == 'REF':
				op = '_'
			if op not in self.coverage[index.pos]:
				self.coverage[index.pos][op] = 0
			return self.coverage[index.pos][op]
		else:
			if index not in self.coverage:
				self.coverage[index] = {}
			return self.coverage[index]

	def __setitem__(self, index, value):
		if isinstance(index, Mutation):
			self.coverage[index.pos][index.op] = value
		else:
			self.coverage[index] = value

	def dump(self, gene, limit=-1, full=False):
		for pos in sorted(self.coverage):
			if limit != -1 and pos > limit:
				break
			for op, cov in self.coverage[pos].items():
				m = Mutation(pos=pos, op=op)
				if self.percentage(m) < 5:
					continue
				if (full or op != '_') and gene.region_at[pos] != '':
					log.debug(
						' {:6} {:5}: {:10} {:10} {:5} ({:3.0f}%)',
						gene.region_at[pos], pos, op,
						gene.old_notation[m],
						cov, self.percentage(m)
					)

	def load(self, sam_path, gene):
		log.debug('SAM file: {}', sam_path)

		muts = collections.defaultdict(int)
		norm = collections.defaultdict(int)

		self.links = collections.defaultdict(int)

		is_deez_file = sam_path[-3:] == '.dz'
		pipe = ''
		if is_deez_file:
			log.debug('Using DeeZ file {}', sam_path)

			def is_exe(path): # http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028
				return os.path.isfile(path) and os.access(path, os.X_OK)
			if not is_exe('deez'):
				found = False
				for path in os.environ["PATH"].split(os.pathsep):
					path = path.strip('"')
					if is_exe(os.path.join(path, 'deez')):
						found = True
						break
				if not found:
					log.error('DeeZ not in PATH')
					exit(1)

			tmpdir = tempfile.mkdtemp()
			pipe = os.path.join(tmpdir, 'dzpipe')
			os.mkfifo(pipe)
			chromosome, chr_start, chr_end = gene.region
			cnv_chromosome, cnv_start, cnv_end = gene.cnv_region
			cnv_region = '{}:{}-{}'.format(cnv_chromosome, cnv_start - 500, cnv_end + 1)[3:]
			region = '{}:{}-{}'.format(chromosome, chr_start - 500, chr_end + 1)[3:]

			command = "deez --stats {} 2>&1 | grep 'Block' | awk '{{print $3}}'".format(sam_path)
			p = subprocess.check_output(command, shell=True)
			if all(x[:3] == 'chr' or x[:2] == '*:' for x in p.split()):
				cnv_region = 'chr' + cnv_region
				region = 'chr' + region

			command = 'deez {} -h -c -Q -r /cs/compbio3/ibrahim/mpeg/ref/hs37d5.fa "{};{}" > {} 2>/dev/null'.format(
				sam_path, region, cnv_region, pipe)
			log.debug(command)
			p = subprocess.Popen(command, shell=True)
			sam_path = pipe

		ref = gene.seq

		# TODO: improve over non-functional mutations
		self.indel_sites = {}
		for a in gene.alleles:
			for m in gene.alleles[a].functional_mutations:
				if m.op[:3] == 'INS':
					self.indel_sites[m] = [collections.defaultdict(int), 0]

		total_reads = 0
		filtered_reads = 0
		with pysam.AlignmentFile(sam_path) as sam:
			try:
				sam.check_index()
			except AttributeError:
				pass # SAM files do not have index
			except ValueError as be:
				raise be

			region = ''
			chromosome, chr_start, chr_end = gene.region
			cnv_chromosome, cnv_start, cnv_end = gene.cnv_region
			if not chromosome.startswith('chr'):
				chromosome = 'chr' + chromosome
				cnv_chromosome = 'chr' + cnv_chromosome

			# TODO: saner way 
			# if sum(x['SN'][:3] == 'chr' for x in sam.header['SQ']) < 2.0/3.0 * len(sam.header['SQ']):
			chrs = [x['SN'] for x in sam.header['SQ']]
			if 'chr1' not in chrs and '1'in chrs:
				chromosome = chromosome[3:]
				cnv_chromosome = cnv_chromosome[3:]
			
			fetched_cnv = False
			if sam.has_index():
				log.debug('File {} has index', sam_path)
				
				region = '{}:{}-{}'.format(cnv_chromosome, cnv_start - 500, cnv_end + 1)
				for read in sam.fetch(region=region):
					start = read.reference_start
					if read.cigartuples is None or (read.flag & 0x40) == 0:
						continue
					for op, size in read.cigartuples:
						if op in [0, 7, 8, 2]:
							for i in range(size):
								self.cnv_coverage[start + i] += 1
							start += size
				fetched_cnv = True

				region = '{}:{}-{}'.format(chromosome, chr_start - 500, chr_end + 1)

			mutation_sites = set()
			for read in sam.fetch(region=region):
				start = read.reference_start
				if not fetched_cnv and read.reference_id != -1 and read.reference_name == cnv_chromosome \
					and cnv_start - 500 <= start <= cnv_end:
					if read.cigartuples is None or (read.flag & 0x40) == 0:
						pass
					else:
						for op, size in read.cigartuples:
							if op in [0, 7, 8, 2]:
								for i in range(size):
									self.cnv_coverage[start + i] += 1
								start += size

				start, s_start = read.reference_start, 0
				if not read.cigartuples or read.reference_name != chromosome or not (chr_start <= start < chr_end):
					continue

				if 'S' in read.cigarstring:
					sk = 0
					for op, size in read.cigartuples:
						if op == 4:
							sk += size
					filtered_reads += 1
				total_reads += 1

				seq = read.query_sequence
				
				has_indel = set()
				orig_start = start

				paired_muts = []
				for op, size in read.cigartuples:
					if op == 2:
						mut = (start, 'DEL.{}'.format(ref[start - chr_start:start - chr_start + size]))
						muts[mut] += 1
						if SAM.PHASE: 
							mutation_sites.add(start)
							paired_muts.append(mut)
						start += size
					elif op == 1:
						mut = (start, 'INS.{}'.format(seq[s_start:s_start + size].lower()))
						muts[mut] += 1
						if SAM.PHASE:
							mutation_sites.add(start)
							paired_muts.append(mut)
						has_indel.add((start, size))
						s_start += size
					elif op == 4:
						s_start += size
					elif op in [0, 7, 8]:
						for i in range(size):
							if 0 <= start + i - chr_start < len(ref) and ref[start + i - chr_start] != seq[s_start + i]:
								mut = (start + i, 'SNP.{}{}'.format(ref[start + i - chr_start], seq[s_start + i]))
								muts[mut] += 1
								if SAM.PHASE: 
									mutation_sites.add(start + i)
									paired_muts.append(mut)
							else:
								norm[start + i] += 1
								if SAM.PHASE and start + i in mutation_sites:
									paired_muts.append((start + i, '_'))
						start += size
						s_start += size

				# TODO: check is real insertion? (does it match reference/read)
				for ins in self.indel_sites:
					ins_len = len(ins.op) - 4
					if (ins.pos, ins_len) in has_indel:
						self.indel_sites[ins][1] += 1
					else:
						self.indel_sites[ins][0][(orig_start, start, len(read.query_sequence))] += 1

				if SAM.PHASE:
					paired_muts = sorted(paired_muts)
					for pi, p in enumerate(paired_muts):
						for pj in range(pi + 1, len(paired_muts)):
							if p[1] == '_' and paired_muts[pj][1] == '_':
								continue
							self.links[(p, paired_muts[pj])] += 1

		for i, c in norm.items():
			if c == 0:
				continue
			if i not in self.coverage:
				self.coverage[i] = {}
			self.coverage[i]['_'] = c
		for k, c in muts.items():
			p, o = k
			if p not in self.coverage:
				self.coverage[p] = {}
			self.coverage[p][o] = c

		for pos in self.coverage:
			self.unfiltered_coverage[pos] = self.total(pos)
		self.original_coverage = copy.deepcopy(self.coverage)

		if pipe != '' and os.path.exists(pipe):
			os.unlink(pipe)

		# pprint(self.links)

	# PGRNseq-v1: PGXT104 for all except CYP2B4; PGXT147 with rescale 2.52444127771 for CYP2B6
	#    cypiripi.py --generate-profile bams/PGXT149.bam --profile-factor 2.52444127771 | grep CYP2B6 > _ha &&
	#    cypiripi.py --generate-profile bams/PGXT104.bam | grep -v CYP2B6 >> _ha
	# PGRNseq-v2: NA19789 for all
	# Illumina: all ones
	def load_profile(self, sam_path, factor):
		log.debug('SAM file: {}', sam_path)

		R = [
			('CYP3A5', '7', 99245000, 99278000),
			('CYP3A4', '7', 99354000, 99465000),
			('CYP2C19', '10', 96445000, 96615000),
			('CYP2C9', '10', 96691000, 96754000),
			('CYP2C8', '10', 96796000, 96830000),
			('CYP4F2', '19', 15619000, 16009500),
			('CYP2A6', '19', 41347500, 41400000),
			('CYP2D6', '22', 42518900, 42553000),
			('TPMT', '6', 18126541, 18157374),
			('DPYD', '1', 97541298, 98388615)
		#	('CYP2B6', '19', 41428000, 41525000),
		#	('SLCO1B1', '12', 21282127, 21394730),
		# 	('CYP21A', '6', 31970000, 32010000),
		# 	('CLN3', '16', 28488000, 28504000),
		]
		for r in R:
			cov = collections.defaultdict(lambda: collections.defaultdict(int))
			rep_cov = collections.defaultdict(lambda: collections.defaultdict(int))

			gene, r = r[0], (r[1], r[2], r[3])
			with pysam.AlignmentFile(sam_path) as sam:
				chrs = [x['SN'] for x in sam.header['SQ']]
				if r[0] not in chrs and 'chr' + r[0] in chrs:
					r = ('chr' + r[0], r[1], r[2])
				log.info('Generating profile for {} ({}:{}-{})', gene, *r)
				try:
					for read in sam.fetch(region='{}:{}-{}'.format(r[0], r[1], r[2])):
						start, s_start = read.reference_start, 0
						if not read.cigartuples:
							continue
						chr = read.reference_name
						if chr.startswith('chr'):
							chr = chr[3:]
						
						for op, size in read.cigartuples:
							if op == 2:
								for i in range(size):
									cov[chr][start + i] += 1
									if read.mapping_quality == 0:
										rep_cov[chr][start + i] += 1
								start += size
							elif op == 1:
								s_start += size
							elif op == 4:
								s_start += size
							elif op in [0, 7, 8]:
								for i in range(size):
									cov[chr][start + i] += 1
									if read.mapping_quality == 0:
										rep_cov[chr][start + i] += 1
								start += size
								s_start += size
				except ValueError as e:
					log.warn('Cannot fetch gene {} ({}:{}-{})', gene, *r)
					
			for c in sorted(cov.keys()):
				for p in sorted(cov[c].keys()):
					if not r[1] - 500 <= p <= r[2] + 500:
						continue
					print(gene, c, p, cov[c][p] * (factor / 2.0))
