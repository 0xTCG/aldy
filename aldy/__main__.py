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
import tempfile
import pkg_resources
import traceback
import re

from .common import log, script_path, AldyException, td, colorize
from .gene import Gene, GRange
from .cn import LEFT_FUSION_PENALTY
from .sam import Sample
from .genotype import genotype
from .version import __version__


def main():
    """
    The main entry point.
    """

    parser, args = _get_args()

    # Set the logging verbosity
    level = args.verbosity.lower()
    level = next(
        v
        for k, v in logbook.base._reverse_level_names.items()
        if k.lower().startswith(level)
    )

    # Set the command-line logging
    sh = logbook.more.ColorizedStderrHandler(
        format_string="{record.message}", level=level
    )
    sh.push_application()

    log.info("*** Aldy v{} (Python {}) ***", __version__, platform.python_version())
    log.info(
        "*** (c) 2016-{} Aldy Authors & Indiana University Bloomington. "
        + "All rights reserved.",
        datetime.datetime.now().year,
    )
    log.info("*** Free for non-commercial/academic use only.")

    try:
        if args.subparser == "help":
            parser.print_help()
        elif args.subparser == "license":
            _print_licence()
        elif args.subparser == "test":
            _run_test()
        elif args.subparser == "show":
            database_file = script_path(
                "aldy.resources.genes/{}.yml".format(args.gene.lower())
            )
            Gene(database_file).print_configurations()
        elif args.subparser == "profile":
            p = Sample.load_sam_profile(args.file)
            for i in p:
                print(*i)
        elif args.subparser == "genotype":
            # Prepare the list of available genes
            if args.gene.lower() == "all":
                avail_genes = pkg_resources.resource_listdir("aldy.resources", "genes")
                avail_genes = [
                    i[:-4] for i in avail_genes if len(i) > 4 and i[-4:] == ".yml"
                ]
                avail_genes = sorted(avail_genes)
            else:
                avail_genes = [args.gene.lower()]

            # Prepare the output file
            output = args.output
            if output == "-":
                output = sys.stdout
            elif output:
                output = open(output, "w")
            else:
                output = "{}.aldy".format(os.path.splitext(args.file)[0])
                output = open(output, "w")

            for gene in avail_genes:
                _genotype(gene, output, args)
            if output != sys.stdout:
                output.close()
        else:
            raise AldyException("Invalid sub-command " + args.subparser)
    except IOError as ex:
        if ex.filename is not None:
            log.critical("File cannot be accessed: {}", ex.filename)
        else:
            log.critical("File cannot be accessed: {}", str(ex))
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
    except:  # noqa
        exc = sys.exc_info()[0]
        log.critical("Unrecoverable error: {}", str(exc))
        log.warn(traceback.format_exc())
        exit(1)


def _get_args():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        prog="aldy",
        description="Allelic decomposition and exact genotyping of highly polymorphic "
        + "and structurally variant genes",
    )

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument(
        "--verbosity",
        "-v",
        default="INFO",
        help=td(
            """Logging verbosity:
               - T (trace)
               - D (debug)
               - I (info) and
               - W (warn).
               Default is "I" (info)."""
        ),
    )

    subparsers = parser.add_subparsers(dest="subparser")

    genotype_parser = subparsers.add_parser(
        "genotype",
        help="Call the most likely genotype and diplotype within a sample.",
        parents=[base],
    )
    genotype_parser.add_argument(
        "file", nargs="?", help="Input file in SAM, BAM, CRAM or DeeZ format."
    )
    genotype_parser.add_argument(
        "--gene",
        "-g",
        default="all",
        help='Gene whose genotype is to be called. Default is "all" which calls '
        + "genotypes for all supported genes.",
    )
    genotype_parser.add_argument(
        "--profile",
        "-p",
        required=True,
        help=td(
            """Sequencing profile. The following profiles are supported:
               - illumina
               - pgrnseq-v1
               - pgrnseq-v2,
               - pgrnseq-v3 and
               - [wxs] (coming soon).
               You can also provide a SAM/BAM file as a profile.
               Please check documentation for more details."""
        ),
    )
    genotype_parser.add_argument(
        "--threshold",
        "-T",
        default=50,
        help="Cut-off rate for variations (percent per copy). Default is 50.",
    )
    genotype_parser.add_argument(
        "--reference",
        "-r",
        default=None,
        help="Genome reference used for reading CRAM or DeeZ files",
    )
    genotype_parser.add_argument(
        "--cn-neutral-region",
        "-n",
        default="22:42547463-42548249",
        help=td(
            """Copy-number neutral region in the format chromosome:start-end
               (e.g. chr1:10000-20000).
               Default is CYP2D8 region within hg19 (22:42547463-42548249)."""
        ),
    )
    genotype_parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output file location. Default is [input].[gene].aldy.",
    )
    genotype_parser.add_argument(
        "--solver",
        "-s",
        default="any",
        help=td(
            """ILP solver:
               - gurobi (Gurobi)
               - scip (SCIP)
               - cbc (Google OR-Tools/CBC)
               - any (attempts to use Gurobi, then SCIP, then CBC).
               Default is "any"."""
        ),
    )
    genotype_parser.add_argument(
        "--gap",
        "-G",
        default=0,
        help=td(
            """Solver optimality gap.
               Any solution whose score is less than (1+gap) times the optimal solution
               score will be reported.
               Default is 0 (report only optimal solutions)."""
        ),
    )
    # genotype_parser.add_argument('--remap', default=0,
    #   help='Realign reads for better mutation calling.
    # Requires samtools and bowtie2 in $PATH.')
    genotype_parser.add_argument(
        "--debug",
        default=None,
        help="Create a directory that will contain the debug information "
        + "and core dumps.",
    )
    genotype_parser.add_argument(
        "--fusion-penalty",
        "-f",
        default=LEFT_FUSION_PENALTY,
        help="Fusion penalty. Use higher values to avoid fusions. "
        + f"Default is {LEFT_FUSION_PENALTY}.",
    )
    genotype_parser.add_argument(
        "--cn",
        "-c",
        default=None,
        help=td(
            """Manually set the copy number configuration as a list of comma-separated
               configuration IDs (e.g. "CN1,CN2").
               For a list of supported configuration IDs, please run `aldy show`.
               For the most common diploid case that contains no structural variations
               (e.g. two copies of the main gene), use 1,1."""
        ),
    )

    _ = subparsers.add_parser(
        "test",
        parents=[base],
        help="Run Aldy test suite. Recommended prior to the first use",
    )

    _ = subparsers.add_parser("license", parents=[base], help="Show Aldy license")

    show_parser = subparsers.add_parser(
        "show",
        parents=[base],
        help="Show all available copy number configurations for a given gene.",
    )
    show_parser.add_argument(
        "--gene", "-g", default="all", help="Gene whose configurations are to be shown."
    )

    profile_parser = subparsers.add_parser(
        "profile",
        parents=[base],
        help=td(
            """Generate a sequencing profile for a SAM/BAM/CRAM file.
               Please check the documentation for more details."""
        ),
    )
    profile_parser.add_argument("file", nargs="?", help="SAM/BAM/CRAM file.")

    _ = subparsers.add_parser(
        "help", parents=[base], help="Show program usage and exit."
    )

    return parser, parser.parse_args()


def _print_licence():
    """
    Print Aldy license.
    """
    with open(script_path("aldy.resources/LICENSE.md")) as f:
        for l in f:
            print(l.strip())


def _genotype(gene: str, output: Optional[str], args) -> None:
    """
    Genotype a file.

    Args:
        gene (str)
        output (str, optional)
        args: remaining command-line arguments

    Raises:
        :obj:`aldy.common.AldyException` if ``cn_region`` is invalid.
    """

    cn_region = args.cn_neutral_region
    if cn_region is not None:
        r = re.match(r"^(.+?):(\d+)-(\d+)$", cn_region)
        if not r:
            raise AldyException(
                f"Parameter --cn-neutral={cn_region} cannot be parsed. "
                + "Must be chr:start-end (where start and end are numbers)"
            )
        ch = r.group(1)
        if ch.startswith("chr"):
            ch = ch[3:]
        cn_region = GRange(ch, int(r.group(2)), int(r.group(3)))

    cn_solution = args.cn
    if cn_solution:
        cn_solution = cn_solution.split(",")

    threshold = float(args.threshold) / 100

    def run(debug):
        log.debug(
            "\nArguments: {}",
            " ".join(k + "=" + str(v) for k, v in vars(args).items() if k is not None),
        )
        log.info("Genotyping sample {}...", os.path.basename(args.file))
        try:
            result = genotype(
                gene_db=gene,
                sam_path=args.file,
                profile=args.profile,
                output_file=output,
                cn_region=cn_region,
                cn_solution=cn_solution,
                threshold=threshold,
                solver=args.solver,
                phase=False,
                fusion_penalty=float(args.fusion_penalty),
                reference=args.reference,
                gap=float(args.gap),
                debug=debug,
            )
            log.info(colorize(f"{gene.upper()} result{' s'[len(result) > 1]}:"))
            for r in result:
                minors = ", ".join(sorted([f.major_repr() for f in r.solution]))
                log.info(colorize(f"  {r.diplotype:30} ({minors})"))
        except AldyException as ex:
            log.error(ex)

    if args.debug:
        with tempfile.TemporaryDirectory() as tmp:
            try:
                prefix = f"{tmp}/{os.path.splitext(os.path.basename(args.file))[0]}"
                log_output = f"{prefix}.log"
                fh = logbook.FileHandler(
                    log_output, mode="w", bubble=True, level="TRACE"
                )
                fh.formatter = lambda record, _: "[{}:{}/{}] {}".format(
                    record.level_name[0],
                    os.path.splitext(os.path.basename(record.filename))[0],
                    record.func_name,
                    record.message,
                )
                fh.push_application()
                log.debug(f"Using {tmp} as temporary debug directory")
                run(prefix)
            finally:
                log.info("Preparing debug archive...")
                os.system(f"tar czvf {args.debug}.tar.gz2 -C {tmp} .")
    else:
        run(None)


def _run_test() -> None:
    """
    Run the Aldy test suite.
    """

    pass


if __name__ == "__main__":
    main()
