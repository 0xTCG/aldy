# 786

# Aldy source: genotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from builtins import str

import os
import sys
import traceback
import logbook
import logbook.more

from . import cn
from . import filtering
from . import protein
from . import refiner
from . import sam
from . import diplotype

from .common import *
from .gene import Gene


@timing
def genotype(sample, output, log_output, gene, profile, threshold, solver, coverage=0):
	try:
		with open(sample): # Check does file exist
			pass

		profile = profile.lower()
		if output is None:
			output = '{}.{}.aldy'.format(os.path.splitext(sample)[0], gene.upper())
		if log_output is None:
			log_output = '{}.{}.aldylog'.format(os.path.splitext(sample)[0], gene.upper())
		fh = logbook.FileHandler(log_output, format_string=LOG_FORMAT, mode='w', bubble=True, level='TRACE')
		fh.push_application()

		log.info('*** Aldy v1.1 ***')
		log.info('(c) 2017 SFU, MIT & IUB. All rights reserved.')
		log.info('Arguments:')
		log.info('  Gene:      {}', gene.upper())
		log.info('  Profile:   {}', profile.lower())
		log.info('  Threshold: {:.0f}%', threshold * 100)
		log.info('  Input:     {}', sample)
		log.info('  Output:    {}', output)
		log.info('  Log:       {}', log_output)
		log.info('  Phasing:   {}', sam.SAM.PHASE)

		database_file = script_path('aldy.resources.genes', '{}.yml'.format(gene.lower()))
		gene = Gene(database_file)

		sample = sam.SAM(sample, gene, threshold)
		gene.alleles = filtering.initial_filter(gene, sample)
		cn_sol = cn.estimate_cn(gene, sample, profile, coverage, solver)
		score, init_sol = protein.get_initial_solution(gene, sample, cn_sol, solver)
		score, sol = refiner.get_refined_solution(gene, sample, init_sol, solver)

		sol = diplotype.assign_diplotype(gene, sol, output)

		return sol
	except ValueError as ex:
		log.critical('Input BAM has no index. Please create index by running samtools index.')
		log.debug(ex)
		exit(1)
	except IOError as ex:
		if ex.filename is not None:
			log.critical('File cannot be accessed: {}', ex.filename)
		else:
			log.critical('File cannot be accessed: {}', str(ex))
		log.debug(ex)
		exit(1)
	except SystemExit as ex:
		log.debug(ex)
		exit(ex.code)
	except:
		ex = sys.exc_info()[0]
		log.critical('Unrecoverable error: {}', str(ex))
		traceback.print_exc()
		exit(1)
