#!/usr/bin/env python
# 786

# Aldy source: __main__.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from __future__ import division

import logbook
import argparse
import os

from . import genotype
from . import protein
from . import sam
from . import lpinterface

from .common import *


def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--threshold', '-T', default=50, help="Cut-off rate for variations (percent per copy). Default is 50%%.")
	parser.add_argument('--profile', '-p', default=None, help='Sequencing profile. Currently, only "illumina", "pgrnseq-v1" and "pgrnseq-v2" are supported.')
	parser.add_argument('--gene', '-g', default='cyp2d6', help='Gene profile. Default is "CYP2D6".')
	parser.add_argument('--verbosity', '-v', default='INFO', help='Logging verbosity. Acceptable values are T (trace), D (debug), I (info) and W (warn). Default is I (info).')
	parser.add_argument('--log', '-l', default=None, help='Location of the output log file (default: [input].[gene].aldylog)')
	parser.add_argument('--output', '-o', default=None, help='Location of the output file (default: [input].[gene].aldy)')
	parser.add_argument('--solver', '-s', default='any', help='IP Solver (default: any available)')
	parser.add_argument('--phase', '-P', default=0, action='store_true', help='Phase reads (default: no; slows down pipeline)')
	parser.add_argument('--license', default=0, action='store_true', help='Show the Aldy license')
	parser.add_argument('--test', default=0, action='store_true', help='Sanity-check Aldy on NA10860 sample')
	parser.add_argument('file', nargs='?', help='SAM or BAM input file')

	# internal parameters for development purposes
	# NOT TO BE USED
	parser.add_argument('--threads', '-t', default=0, help='DO NOT USE / ONLY FOR DEVELOPMENT USE')
	parser.add_argument('--coverage', '-c', default=0, help='Manually set the average coverage. Use only if you are really sure what you\'re doing!')
	parser.add_argument('--cache', '-C', dest='cache', action='store_true', help='DO NOT USE / ONLY FOR DEVELOPMENT USE')
	parser.add_argument('--no-display-cnv', '-N', dest='no_display_cnv', action='store_true', help='DO NOT USE / ONLY FOR DEVELOPMENT USE')
	parser.add_argument('--generate-profile', dest='generate_profile', action='store_true', help='DO NOT USE / ONLY FOR DEVELOPMENT USE')
	parser.add_argument('--profile-factor', dest='profile_factor', default=2.0, help='DO NOT USE / ONLY FOR DEVELOPMENT USE')

	return parser.parse_args()


def print_licence():
	with open(script_path('aldy.resources', 'LICENSE.md')) as f:
		for l in f: print(l.strip())


def run_test():
	log.warn('Aldy Sanity-Check Test')
	log.warn('Expected result is: *1/*4+*4\n')
	result = genotype.genotype(
		script_path('aldy.resources', 'NA10860.bam'), 
		'', os.devnull, 'CYP2D6', 'illumina',
		0.5, 'any', 0
	)
	log.warn('Result{}: ', '' if len(result) == 1 else 's')
	for rd, r in result:
		log.warn('  {:30} ({})', rd, ', '.join(r))


def main(args=None):
	args = get_args()

	level = args.verbosity.lower()
	for k, v in logbook.base._reverse_level_names.items():
		if k.lower().startswith(level):
			level = v
			break

	sh = logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level=level)
	sh.push_application()

	if args.license:
		print_licence()
		exit(0)
	elif args.test:
		run_test()
		exit(0)
	elif args.file is None:
		log.critical('Input file must be specified!')
		exit(1)

	if args.generate_profile:
		sam.SAM.PROFILE = True
		sam.SAM(args.file, None, float(args.profile_factor))

	if args.cache:
		sam.SAM.CACHE = True
	if args.no_display_cnv:
		protein.SHOW_CN_INFO = False
	if args.phase:
		sam.SAM.PHASE = True

	try:
		m = lpinterface.model('init', args.solver)
	except Exception as e:
		log.critical('{}', e.message)


	if args.profile is None:
		log.critical('No profile provided! Please run aldy --help for available profile choices.')
		exit(1)

	result = genotype.genotype(
		args.file, args.output,
		args.log,
		args.gene, args.profile,
		float(args.threshold) / 100.0,
		args.solver,
		float(args.coverage)
	)
	log.warn('Result{}: ', '' if len(result) == 1 else 's')
	for rd, r in result:
		log.warn('  {:30} ({})', rd, ', '.join(r))


if __name__ == "__main__":
	main()
