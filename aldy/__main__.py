#!/usr/bin/env python
# 786

# Aldy source: __main__.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from __future__ import division

import logbook
import argparse
import os
import sys
import platform
import multiprocessing
import functools
import traceback


from . import genotype
from . import protein
from . import gene
from . import sam
from . import lpinterface
from .version import __version__


from .common import *


def _get_args():
   parser = argparse.ArgumentParser()
   parser.add_argument('--threshold', '-T', default=50, 
      help="Cut-off rate for variations (percent per copy). Default is 50%%.")
   parser.add_argument('--profile', '-p', default=None, 
      help='Sequencing profile. Currently, "illumina", "pgrnseq-v1", ' + 
          '"pgrnseq-v2" and "wxs" are supported out of the box. You can also use SAM/BAM as a profile. Please check documentation for more information.')
   parser.add_argument('--gene', '-g', default='all', 
      help='Gene profile. Default is "all" which calls all supported genes.')
   parser.add_argument('--verbosity', '-v', default='INFO', 
      help='Logging verbosity. Acceptable values are ' + 
          'T (trace), D (debug), I (info) and W (warn). Default is I (info).')
   parser.add_argument('--log', '-l', default=None, 
      help='Location of the output log file (default: [input].[gene].aldylog)')
   parser.add_argument('--output', '-o', default=None, 
      help='Location of the output file (default: [input].[gene].aldy)')
   parser.add_argument('--solver', '-s', default='any', 
      help='IP Solver (default: any available)')
   parser.add_argument('--cn-neutral', '-n', default=None, 
      help='Manually specify copy number-neutral region in the sample to be used for calculating region copy numbers. '
      + 'Format is chromosome:start-end (e.g. chr1:10000-20000). ' +
       'Default is as specified in the gene YML file (typically CYP2D8 region).')
   parser.add_argument('--phase', '-P', default=0, action='store_true', 
      help='Phase reads for better variant calling. May provide neglegible benefits at the cost of significant slowdown. Default is no.')
   parser.add_argument('--license', default=0, action='store_true', 
      help='Show the Aldy license')
   parser.add_argument('--test', default=0, action='store_true', 
      help='Sanity-check Aldy on NA10860 sample')
   parser.add_argument('--remap', default=0, #action='store_true', 
      help='Run realignment via bowtie2')
   parser.add_argument('file', nargs='?', 
      help='SAM or BAM input file')

   parser.add_argument('--show-cn', dest='show_cn', action='store_true',
      help='Shows all copy number configurations supported by a gene (requires -g).')
   parser.add_argument('--cn', '-c', default=None,
      help='Manually set copy number (input: a comma-separated list CN1,CN2,...). ' +
          'For a list of supported configurations, please run --show-cn. For standard diploid case specify -c 1,1')
   parser.add_argument('--generate-profile', dest='generate_profile', action='store_true', 
      help='Generate a sequencing profile for a given SAM/BAM. Please check documentation for more information.')

   # internal parameters for development purposes: please do not use unless instructed
   parser.add_argument('--cache', '-C', dest='cache', action='store_true', help=argparse.SUPPRESS)
   parser.add_argument('--profile-factor', dest='profile_factor', default=2.0, help=argparse.SUPPRESS)

   # WARNING: internally all chromosome names do not have chr prefix 
   return parser.parse_args()


def _print_licence():
   with open(script_path('aldy.resources', 'LICENSE.md')) as f:
      for l in f: print(l.strip())


def run(gene, file, output, log_output, profile, threshold, solver, cn_region, cn, remap):
   try:
      result = genotype.genotype(
         file, output,
         log_output,
         gene, profile,
         float(threshold) / 100.0,
         solver,
         cn_region, cn,
         remap
      )
      log.warn('Result{} for {}: ', '' if len(result) == 1 else 's', gene.upper())
      for rd, r in result:
         log.warn('  {:30} ({})', rd, ', '.join(r))
   except AldyException as ex:
      log.error(ex.message)


def _run_test():
   log.warn('Aldy Sanity-Check Test')
   log.warn('Expected result is: *1/*4+*4')
   run('cyp2d6', script_path('aldy.resources', 'NA10860.bam'),
      '', os.devnull, 'illumina', 0.5, 'any', None)


def main(args=None):
   args = _get_args()

   level = args.verbosity.lower()
   for k, v in logbook.base._reverse_level_names.items():
      if k.lower().startswith(level):
         level = v
         break

   sh = logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level=level)
   sh.push_application()

   log.info('*** Aldy v{} (Python {}) ***', __version__, platform.python_version())
   log.info('(c) 2016-2018 SFU, MIT & IUB. All rights reserved.')
   log.info('Arguments: {}', ' '.join(k+'='+str(v) for k, v in vars(args).items() if k is not None))
   

   if args.license:
      _print_licence()
      exit(0)
   elif args.test:
      _run_test()
      exit(0)
   elif args.show_cn:
      database_file = script_path('aldy.resources.genes', '{}.yml'.format(args.gene.lower()))
      g = gene.Gene(database_file)
      g.print_configurations()
      exit(0)
   elif args.file is None:
      log.critical('Input file must be specified!')
      exit(1)

   if args.generate_profile:
      p = sam.SAM.load_profile(args.file, float(args.profile_factor))
      for a, b, c, d in p:
         print(a, b, c, d)
      exit(0)

   if args.cache:
      sam.SAM.CACHE = True
   if args.phase:
      sam.SAM.PHASE = True

   try:
      m = lpinterface.model('init', args.solver)
   except Exception as ex:
      if hasattr(ex, 'message'):
         log.critical(ex.message)
      else:
         log.critical(ex)

   if args.profile is None:
      log.critical('No profile provided! Please run aldy --help for available profile choices.')
      exit(1)

   if args.gene.lower() == 'all':
      avail_genes = pkg_resources.resource_listdir('aldy.resources', 'genes')
      avail_genes = [i[:-4] for i in avail_genes if len(i) > 4 and i[-4:] == '.yml']
      avail_genes = sorted(avail_genes)
   else:
      avail_genes = [args.gene.lower()]

   log.info('  Profile:   {}', args.profile)
   log.info('  Threshold: {:.0f}%', float(args.threshold))
   log.info('  Input:     {}', args.file)
   
   if args.log is None:
      log_output = '{}.aldylog'.format(os.path.splitext(args.file)[0])
   else:
      log_output = args.log
   fh = logbook.FileHandler(log_output, format_string=LOG_FORMAT, mode='w', bubble=True, level='TRACE')
   fh.push_application()
   log.info('  Log:       {}', log_output)
   log.info('  Phasing:   {}', sam.SAM.PHASE)

   # pool = multiprocessing.Pool(4)
   # pool.map(functools.partial(run, args=args, log_output=log_output), avail_genes)
   try:
      output = args.output
      if output == '-':
         output = sys.stdout
      elif output:
         log.info('  Output:    {}', output)
         output = open(output, 'w')
      else:
         output = '{}.aldy'.format(os.path.splitext(args.file)[0])
         log.info('  Output:    {}', output)
         output = open(output, 'w')
      for gene in avail_genes:
         run(gene, args.file, output, log_output, args.profile, args.threshold, args.solver, args.cn_neutral, args.cn, args.remap)
      if output != sys.stdout:
         output.close()
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
   except Exception as ex:
      if hasattr(ex, 'message'):
         log.critical(ex.message)
      else:
         log.critical(ex)
      log.debug(traceback.format_exc())
      exit(1)
   except:
      ex = sys.exc_info()[0]
      log.critical('Unrecoverable error: {}', str(ex))
      log.debug(traceback.format_exc())
      exit(1)

if __name__ == "__main__":
   main()
