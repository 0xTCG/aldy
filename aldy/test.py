#!/usr/bin/env python
# 786

# Aldy source: test.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from __future__ import print_function
from builtins import str

import re
import os
import sys
import pandas as pd
import logbook
import logbook.more

import genotype
import solver

from gene import SAM
from common import log, LOG_FORMAT, colorize, allele_key


def sortkey(x):
	return int(allele_key(x))


def test_single(sample, location, expected, profile, threshold, gene):
	expected = [r for r in expected if not pd.isnull(r)]
	message = '{} - {}'.format(sample, ' or '.join(expected))
	message = colorize('{:30}'.format(message, 'teal'))
	if '' not in expected:
		expected = [[str(x).strip() for x in re.split('[/\+]', r)] for r in expected]
		expected = set(tuple(sorted(r, key=sortkey)) for r in expected)
	else:
		expected = set()

	try:
		solutions = genotype.genotype(location, gene, profile, threshold)
	except:
		logbook.error('{} {}', message, colorize('CRASH ', 'red'))
		exit(1)

	def fix(s):
		return re.split('(\d+)', s)[1]
	orig_solutions = solutions
	solutions = set(tuple(sorted((fix(p) for p in s), key=sortkey)) for s in solutions)
	expected = set(tuple(sorted((fix(p) for p in s), key=sortkey)) for s in expected)

	if solutions == expected:
		logbook.info('{} {} {}', message, colorize('OK   ', 'green'), list(orig_solutions))
		return 1
	elif solutions <= expected and len(solutions) != 0:
		logbook.info('{} {} {}', message, colorize('OK<  ', 'green'), list(orig_solutions))
		return 2
	elif len(expected & solutions) > 0:
		logbook.warn('{} {} {}', message, colorize('MULTI', 'yellow'), list(orig_solutions))
		return 3
	else:
		logbook.error('{} {} {}', message, colorize('FAIL ', 'red'), list(orig_solutions))
		return 0


def test_samples(samples, location, profile, threshold, gene='cyp2d6'):
	from pathos.multiprocessing import ProcessingPool as Pool
	
	pool = Pool(processes=16)
	
	def f(s):
		global counter
		sample, expected = s[0], s[1:]
		# if sample not in T: continue
		loc = location.format(sample)
		if not os.path.exists(loc): # or sample == 'PGXT122':
			return 0
		try:
			res = test_single(sample, loc, expected, profile, threshold, gene)
		except Exception as e:
			print(s, e)
			res = 0
		# sys.stderr.write('-')
		return res
	
	# sys.stderr.write(' ' * len(samples) + '|\r')
	result = pool.map(f, samples)
	logbook.warn(
		'Passed {} out of {} ({} subset, {} multi)\n',
		colorize(str(len([x for x in result if x > 0])), 'green'),
		colorize(str(len(result)), 'blue'),
		colorize(str(len([x for x in result if x == 2])), 'yellow'),
		colorize(str(len([x for x in result if x == 3])), 'yellow')
	)


if __name__ == "__main__":
	path = os.path.dirname(os.path.realpath(__file__))
	
	if '-f' in sys.argv:
		os.system(
			'curl -L "https://docs.google.com/spreadsheets/d/17GUAk5CP9rW7Dj918biIY5bKyCHY_45iq9iKKqX8yP0/export?format=xlsx"' +
			' > {}/../resources/tests.xlsx'.format(path))

	sh = logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level='WARNING')
	sh.push_application()
	log_file = path + '/../test-results.txt'
	fh = logbook.FileHandler(log_file, format_string=LOG_FORMAT, mode='w', bubble=True, level='INFO')
	fh.push_application()

	SAM.CACHE = True
	solver.SHOW_CN_INFO = False

	def get_samples(sheet, gene, samples_path, profile, threshold):
		if 'CDC' in sheet:
			data = pd.read_excel(path + '/../resources/tests.xlsx', header=[0, 2], sheetname=sheet, parse_dates=False)
			data = data.set_index(['PGXT'])
		else:
			data = pd.read_excel(path + '/../resources/tests.xlsx', header=[0, 1], sheetname=sheet, parse_dates=False)
		samples = []
		for r in data[gene].iterrows():
			sample_id = r[0]
			if isinstance(sample_id, tuple):
				sample_id = sample_id[0].upper()
			p = [sample_id, r[1].Validated]
			if not pd.isnull(r[1].Additional):
				p += r[1].Additional.split(';')
			samples.append(p)

		log.warn('{} [{}]', sheet, gene)
		test_samples(samples, path + samples_path, profile=profile, threshold=threshold, gene=gene.lower())

	# get_samples('PGRNseq-v2 (Old)', gene='CYP2D6', samples_path='/../data/pgrnseq-old/{}.dz', profile='pgrnseq-v2', threshold=.5)
	# get_samples('PGRNseq-v2',       gene='CYP2D6', samples_path='/../data/pgrnseq-96/{}.dz',  profile='pgrnseq-v2', threshold=.5)
	
	# get_samples('Illumina',      gene='CYP2D6', samples_path='/../data/illumina/{}.bam',        profile='illumina', threshold=.5)
	# get_samples('Illumina (IU)', gene='CYP2D6', samples_path='/../data/illumina-milan/{}.bam',  profile='illumina', threshold=.5)
	
	# get_samples('PGRNseq-v1 (CDC)', gene='CYP2D6',  samples_path='/../data/pgrnseq-cdc/{}.dz',  profile='pgrnseq-v1', threshold=.5)
	# get_samples('PGRNseq-v1 (CDC)', gene='CYP2A6',  samples_path='/../data/pgrnseq-cdc/{}.dz',  profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='CYP2C19', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='CYP2C8', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='CYP2C9', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='CYP3A4', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='CYP3A5', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='CYP4F2', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='TPMT', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	get_samples('PGRNseq-v1 (CDC)', gene='DPYD', samples_path='/../data/pgrnseq-cdc/{}.dz', profile='pgrnseq-v1', threshold=.5)
	
	# TODO:
	# get_samples('PGRNseq-v1 (CDC)', gene='CYP2B6',  samples_path='/../data/pgrnseq-cdc/{}.dz',  profile='pgrnseq-v1', threshold=.5)
	# get_samples('PGRNseq-v1 (CDC)', gene='SLCO1B1', samples_path='/../data/pgrnseq-cdc/{}.dz',  profile='pgrnseq-v1', threshold=.5)
