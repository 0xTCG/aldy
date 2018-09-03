#!/usr/bin/env python
# 786

# Aldy source: __main__.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import logbook
import argparse
import os
import sys
import platform
import multiprocessing
import functools
import traceback

from .common import *
from .gene import Gene
from .sam import Sample
from .genotype import genotype_init
from .lpinterface import model as lp_model
from .version import __version__


def main():
   """
   Entry point
   """

   parser, args = _get_args()


   # Set the logging verbosity
   level = args.verbosity.lower()
   level = next(v for k, v in logbook.base._reverse_level_names.items() if k.lower().startswith(level))
   
   # Set the command-line logging
   sh = logbook.more.ColorizedStderrHandler(format_string=LOG_FORMAT, level=level)
   sh.push_application()
   sh.formatter = lambda record, _: '[{}:{}/{}] {}'.format(
      record.level_name[0], os.path.splitext(os.path.basename(record.filename))[0], record.func_name, record.message) 

   log.info('*** Aldy v{} (Python {}) ***', __version__, platform.python_version())
   log.info('(c) 2016-2018 SFU, MIT & IUB. All rights reserved.')
   log.info('Arguments: {}', ' '.join(k+'='+str(v) for k, v in vars(args).items() if k is not None))

   try:
      if args.subparser == 'help':
         parser.print_help()
      elif args.subparser == 'license':
         _print_licence()
      elif args.subparser == 'test':
         _run_test()
      elif args.subparser == 'show':
         database_file = script_path('aldy.resources.genes', '{}.yml'.format(args.gene.lower()))
         Gene(database_file).print_configurations()
      elif args.subparser == 'profile':
         p = Sample.load_sam_profile(args.file, float(args.profile_factor))
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

         log.info('  Profile:   {}', args.profile)
         log.info('  Threshold: {:.0f}%', float(args.threshold))
         log.info('  Input:     {}', args.file)
      
         # Prepare the log files
         if args.log is None:
            log_output = '{}.aldylog'.format(os.path.splitext(args.file)[0])
         else:
            log_output = args.log
         fh = logbook.FileHandler(log_output, format_string=LOG_FORMAT, mode='w', bubble=True, level='TRACE')
         fh.push_application()
         log.info('  Log:       {}', log_output)
         log.info('  Phasing:   {}', args.phase)

         # Prepare output file
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
            _genotype(gene, args.file, output, log_output, 
                      args.profile, args.threshold, args.solver, 
                      args.cn_neutral, args.cn, args.remap)
         if output != sys.stdout:
            output.close()
      else:
         raise AldyException('Invalid sub-command ' + args.subparser)
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
      log.critical(ex)
      log.debug(traceback.format_exc())
      exit(1)
   except:
      ex = sys.exc_info()[0]
      log.critical('Unrecoverable error: {}', str(ex))
      log.debug(traceback.format_exc())
      exit(1)


def _get_args():
   """
   Prepares the command-line arguments
   """

   parser = argparse.ArgumentParser(prog='aldy', 
      description=td("""
         Allelic decomposition and exact genotyping of highly polymorphic 
         and structurally variant genes"""))
   parser.add_argument('--verbosity', '-v', default='INFO', 
      help=td("""
         Logging verbosity. Acceptable values are:
         - T (trace)
         - D (debug)
         - I (info) and
         - W (warn). 
         Default is "I" (info)."""))
   parser.add_argument('--log', '-l', default=None, 
      help='Location of the output log file. Default is [input].[gene].aldylog')

   subparsers = parser.add_subparsers(dest="subparser")
   
   genotype_parser = subparsers.add_parser('genotype', 
      help='Genotype a SAM/BAM/CRAM/DeeZ file')
   genotype_parser.add_argument('file', nargs='?', 
      help='Input sample in SAM, BAM, CRAM or DeeZ format.') 
   genotype_parser.add_argument('--gene', '-g', default='all', 
      help='Gene to be genotyped. Default is "all" which attempt to genotype all supported genes.')
   genotype_parser.add_argument('--profile', '-p', 
      help=td("""
         Sequencing profile. Currently, the following profiles are supported out of the box:
         - illumina
         - pgrnseq-v1
         - pgrnseq-v2 and
         - wxs. 
         You can also provide a SAM/BAM file as a profile. 
         Please check documentation for more information."""))
   genotype_parser.add_argument('--threshold', '-T', default=50, 
      help="Cut-off rate for variations (percent per copy). Default is 50.")
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
         - any (attempts to use gurobi, and if fails, uses scip).
         Default is "any" """))
   genotype_parser.add_argument('--phase', '-P', default=0, action='store_true', 
      help=td("""
      Phase aligned reads for better variant calling. 
      May provide neglegible benefits at the cost of significant slowdown. 
      Default is off."""))
   genotype_parser.add_argument('--remap', default=0, #action='store_true', 
      help='Realign reads for better mutation calling. Requires samtools and bowtie2 in $PATH.')
   genotype_parser.add_argument('--cn', '-c', default=None,
      help=td("""
         Manually set the copy number configuration.
         Input format is a comma-separated list of configuration IDs: e.g. "CN1,CN2".
         For a list of supported configuration IDs, please run aldy --show-cn. 
         For a commonly used diploid case (e.g. 2 copies of the main gene) specify -c 1,1"""))
   # HACK: Internal parameters for development purposes: please do not use unless instructed
   genotype_parser.add_argument('--cache', dest='cache', action='store_true', help=argparse.SUPPRESS)
   genotype_parser.add_argument('--profile-factor', dest='profile_factor', default=2.0, help=argparse.SUPPRESS)

   _ = subparsers.add_parser('test',
      help='Sanity-check Aldy on NA10860 sample. Recommended prior to the first use')
   
   _ = subparsers.add_parser('license', 
      help='Show Aldy license')
   
   show_parser = subparsers.add_parser('show',
      help='Show all available copy number configurations of a given gene.')
   show_parser.add_argument('--gene', '-g', default='all', 
      help='Gene to be shown.')

   profile_parser = subparsers.add_parser('profile',
      help=td("""
         Generate a sequencing profile for a given alignment file in SAM/BAM/CRAM format. 
         Please check documentation for more information."""))
   profile_parser.add_argument('file', nargs='?', 
      help='Input sample in SAM, BAM, CRAM or DeeZ format.') 

   _ = subparsers.add_parser('help', help='Show a help message and exit.')

   return parser, parser.parse_args()


def _print_licence():
   """Prints Aldy license"""

   with open(script_path('aldy.resources', 'LICENSE.md')) as f:
      for l in f: print(l.strip())


def _genotype(gene, file, output, log_output, profile, threshold, solver, cn_region, cn, remap):
   """Attempts to genotype a file"""
   
   def test_lp_solver(solver):
      _ = lp_model('init', solver)

   try:
      test_lp_solver(solver)
      result = genotype_init(
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
      log.error(ex)


def _run_test():
   """
   Runs the full pipeline as a sanity check on NA12878 sample (located within the package)
   """

   log.warn('Aldy Sanity-Check Test')
   log.warn('Expected result is: *1/*4+*4')
   _genotype(
      gene='cyp2d6', 
      file=script_path('aldy.resources', 'NA10860.bam'),
      output='', 
      log_output=os.devnull, 
      profile='illumina', 
      threshold=0.5, 
      solver='any', 
      cn_region='22:42547463-42548249', 
      cn=None, 
      remap=False)


if __name__ == "__main__":
   main()
