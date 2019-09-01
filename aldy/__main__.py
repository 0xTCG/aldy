#!/usr/bin/env python
# 786

# Aldy source: __main__.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Optional

import logbook
import argparse
import os
import sys
import platform
import datetime
import multiprocessing
import functools
import traceback

from .common import *
from .gene import Gene
from .sam import Sample, DEFAULT_CN_NEUTRAL_REGION
from .genotype import genotype
from .version import __version__


def main():
   """
   Entry point.
   """

   parser, args = _get_args()

   # Set the logging verbosity
   level = args.verbosity.lower()
   level = next(v for k, v in logbook.base._reverse_level_names.items() if k.lower().startswith(level))
   
   # Set the command-line logging
   sh = logbook.more.ColorizedStderrHandler(format_string='{record.message}', level=level)
   sh.push_application()

   log.info('*** Aldy v{} (Python {}) ***', __version__, platform.python_version())
   log.info('*** (c) 2016-{} Aldy Authors & Indiana University Bloomington. All rights reserved.', datetime.datetime.now().year)
   log.info('*** Free for non-commercial/academic use only.')
   log.debug('Arguments: {}', ' '.join(k+'='+str(v) for k, v in vars(args).items() if k is not None))

   try:
      if args.subparser == 'help':
         parser.print_help()
      elif args.subparser == 'license':
         _print_licence()
      elif args.subparser == 'test':
         _run_test()
      elif args.subparser == 'show':
         database_file = script_path('aldy.resources.genes/{}.yml'.format(args.gene.lower()))
         Gene(database_file).print_configurations()
      elif args.subparser == 'profile':
         p = Sample.load_sam_profile(args.file)
         for i in p:
            print(*i)
      elif args.subparser == 'genotype':
         # Prepare the list of genes to be genotypes
         if args.gene.lower() == 'all':
            avail_genes = pkg_resources.resource_listdir('aldy.resources', 'genes')
            avail_genes = [i[:-4] for i in avail_genes if len(i) > 4 and i[-4:] == '.yml']
            avail_genes = sorted(avail_genes)
         else:
            avail_genes = [args.gene.lower()]

         # Prepare the log files
         if args.debug is not None:
            log_output = f'{args.debug}/{os.path.splitext(args.file)[0]}.log'
            fh = logbook.FileHandler(log_output, mode='w', bubble=True, level='TRACE')
            fh.formatter = lambda record, _: '[{}:{}/{}] {}'.format(
               record.level_name[0], os.path.splitext(os.path.basename(record.filename))[0], record.func_name, record.message) 
            fh.push_application()

         # Prepare output file
         output = args.output
         if output == '-':
            output = sys.stdout
         elif output:
            output = open(output, 'w')
         else:
            output = '{}.aldy'.format(os.path.splitext(args.file)[0])
            output = open(output, 'w')
         
         for gene in avail_genes:
            _genotype(gene, output, args) 
         if output != sys.stdout:
            output.close()
      else:
         raise AldyException('Invalid sub-command ' + args.subparser)
   except IOError as ex:
      if ex.filename is not None:
         log.critical('File cannot be accessed: {}', ex.filename)
      else:
         log.critical('File cannot be accessed: {}', str(ex))
      log.debug(ex)
      log.debug(traceback.format_exc())
      exit(1)
   except SystemExit as ex:
      log.debug(ex)
      log.debug(traceback.format_exc())
      exit(ex.code)
   except Exception as ex:
      log.critical(ex)
      log.warn(traceback.format_exc())
      exit(1)
   except:
      exc = sys.exc_info()[0]
      log.critical('Unrecoverable error: {}', str(exc))
      log.warn(traceback.format_exc())
      exit(1)


def _get_args():
   """
   Prepares the command-line arguments.
   """

   parser = argparse.ArgumentParser(prog='aldy', 
      description=td("""
         Allelic decomposition and exact genotyping of highly polymorphic 
         and structurally variant genes"""))

   base = argparse.ArgumentParser(add_help=False)
   base.add_argument('--verbosity', '-v', default='INFO', 
      help=td("""
         Logging verbosity. Acceptable values are:
         - T (trace)
         - D (debug)
         - I (info) and
         - W (warn). 
         Default is "I" (info)."""))
   base.add_argument('--log', '-l', default=None, 
      help='Location of the output log file. Default is [input].[gene].aldylog')
   
   subparsers = parser.add_subparsers(dest="subparser")
   
   genotype_parser = subparsers.add_parser('genotype', 
      help='Genotype a SAM/BAM/CRAM/DeeZ file', parents=[base])
   genotype_parser.add_argument('file', nargs='?', 
      help='Input sample in SAM, BAM, CRAM or DeeZ format.') 
   genotype_parser.add_argument('--gene', '-g', default='all', 
      help='Gene to be genotyped. Default is "all" which attempt to genotype all supported genes.')
   genotype_parser.add_argument('--profile', '-p', 
      required=True,
      help=td("""
         Sequencing profile. Currently, the following profiles are supported out of the box:
         - illumina
         - pgrnseq-v1
         - pgrnseq-v2 and
         - [wxs] (coming soon). 
         You can also provide a SAM/BAM file as a profile. 
         Please check documentation for more information."""))
   genotype_parser.add_argument('--threshold', '-T', default=50, 
      help="Cut-off rate for variations (percent per copy). Default is 50.")
   genotype_parser.add_argument('--reference', '-r', default=None, 
      help="CRAM or DeeZ FASTQ reference")
   genotype_parser.add_argument('--cn-neutral-region', '-n', default='22:42547463-42548249', 
      help=td("""
         Copy-number neutral region. Format of the region is chromosome:start-end 
         (e.g. chr1:10000-20000). 
         Default value is CYP2D8 within hg19 (22:42547463-42548249)."""))
   genotype_parser.add_argument('--output', '-o', default=None, 
      help='Location of the output file (default: [input].[gene].aldy)')
   genotype_parser.add_argument('--solver', '-s', default='any', 
      help=td("""
         ILP Solver. Available solvers:
         - gurobi (Gurobi)
         - scip (SCIP)
         - any (attempts to use Gurobi, and if fails, uses SCIP).
         Default is "any" """))
   genotype_parser.add_argument('--phase', '-P', default=0, action='store_true', 
      help=td("""
      Phase aligned reads for better variant calling. 
      May provide neglegible benefits at the cost of significant slowdown. 
      Default is off."""))
   genotype_parser.add_argument('--remap', default=0, #action='store_true', 
      help='Realign reads for better mutation calling. Requires samtools and bowtie2 in $PATH.')
   genotype_parser.add_argument('--debug', default=None,
      help='Specify the debug directory where the solver info is saved for debugging purposes.')
   genotype_parser.add_argument('--cn', '-c', default=None,
      help=td("""
         Manually set the copy number configuration.
         Input format is a comma-separated list of configuration IDs: e.g. "CN1,CN2".
         For a list of supported configuration IDs, please run aldy --show-cn. 
         For a commonly used diploid case (e.g. 2 copies of the main gene) specify -c 1,1"""))

   _ = subparsers.add_parser('test', parents=[base],
      help='Sanity-check Aldy on NA10860 sample. Recommended prior to the first use')
   
   _ = subparsers.add_parser('license', parents=[base], 
      help='Show Aldy license')
   
   show_parser = subparsers.add_parser('show', parents=[base],
      help='Show all available copy number configurations of a given gene.')
   show_parser.add_argument('--gene', '-g', default='all', 
      help='Gene to be shown.')

   profile_parser = subparsers.add_parser('profile', parents=[base],
      help=td("""
         Generate a sequencing profile for a given alignment file in SAM/BAM/CRAM format. 
         Please check documentation for more information."""))
   profile_parser.add_argument('file', nargs='?', 
      help='Input sample in SAM, BAM, CRAM or DeeZ format.') 

   _ = subparsers.add_parser('help', parents=[base], 
      help='Show a help message and exit.')

   return parser, parser.parse_args()


def _print_licence():
   """
   Prints Aldy license.
   """
   with open(script_path('aldy.resources/LICENSE.md')) as f:
      for l in f: print(l.strip())


def _genotype(gene: str, output: Optional, args) -> None:
   """
   Attempts to genotype a file.

   Args:
      gene (str)
      output (file, optional)
      args: remaining arguments 

   Raises:
      :obj:`aldy.common.AldyException` if ``cn_region`` is invalid.
   """
   
   cn_region = args.cn_neutral_region
   if cn_region is not None:
      r = re.match(r'^(.+?):(\d+)-(\d+)$', cn_region)
      if not r:
         raise AldyException(f'Parameter --cn-neutral={cn_region} is not in the format chr:start-end (where start and end are numbers)')
      ch = r.group(1)
      if ch.startswith('chr'):
         ch = ch[3:]
      cn_region = GRange(ch, int(r.group(2)), int(r.group(3)))
      log.info('Using {} as copy-number neutral region', cn_region)
   
   cn_solution = args.cn
   if cn_solution:
      cn_solution = cn_solution.split(',')

   threshold = float(args.threshold) / 100
   
   if debug:
      os.makedirs(debug, exist_ok=True)
      debug = f'{debug}/{os.path.splitext(args.file)[0]}'
   
   try:
      result = genotype(gene_db=     gene, 
                        sam_path=    args.file,
                        profile=     args.profile,
                        output_file= output,
                        cn_region=   cn_region,
                        cn_solution= cn_solution,
                        threshold=   threshold,
                        solver=      args.solver,
                        phase=       args.phase,
                        reference=   args.reference,
                        debug=       args.debug)
      log.info(colorize('Result{} for {}: '.format('' if len(result) == 1 else 's', gene.upper())))
      for r in result:
         log.info(colorize('  {:30} ({})'.format(r.diplotype, ', '.join(f.major_repr() for f in r.solution))))
   except AldyException as ex:
      log.error(ex)


def _run_test() -> None:
   """
   Runs the full pipeline as a sanity check on NA12878 sample (located within the package).
   """

   log.warn('Aldy Sanity-Check Test')
   log.warn('Expected result is: *1/*4+*4')
   
   # Masquerade args via this helper class
   class DictWrapper:
      def __init__(self, d):
         self.d = d
      def __getattr__(self, key):
         if key in self.d:
            return self.d[key]
         else:
            return None
   _genotype(gene='cyp2d6', 
             output=None,
             args=DictWrapper({'file': script_path('aldy.resources/NA10860.bam'),
                               'profile': 'illumina',
                               'threshold': 50,
                               'solver': 'any',
                               'phase': False}))


if __name__ == "__main__":
   main()
   #import cProfile
   #cProfile.run('main()', filename='aldy.prof', sort='cumulative')
