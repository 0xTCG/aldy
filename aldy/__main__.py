#!/usr/bin/env python
# 786
# Aldy source: __main__.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Optional, Any
import logbook
import logbook.more
import logbook.base
import argparse
import os
import sys
import platform
import yaml
import datetime
import tempfile
import traceback
import pytest

from . import common
from .common import log, script_path, AldyException, td, parse_cn_region
from .gene import Gene
from .profile import Profile
from .genotype import genotype
from .query import query
from .version import __version__


def get_version():
    return "{} {}".format(
        platform.system() if platform.system() != "Darwin" else "macOS",
        platform.platform()[6:]
        if platform.system() == "Linux"
        else platform.mac_ver()[0]
        if platform.system() == "Darwin"
        else platform.platform(),
    )


def main(argv):
    """
    The main entry point.
    """

    parser, args = _get_args(argv)

    # Set the logging verbosity
    level = (args.verbosity if "verbosity" in args else "I").lower()
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

    log.info(
        "🐿  Aldy v{} (Python {} on {})",
        __version__,
        platform.python_version(),
        get_version(),
    )
    log.info(
        "   (c) 2016-{} Aldy Authors. All rights reserved.\n"
        + "   Free for non-commercial/academic use only.",
        datetime.datetime.now().year,
    )

    try:
        if not args.subparser or args.subparser == "help":
            parser.print_help()
        elif args.subparser == "license":
            _print_licence()
        elif args.subparser == "test":
            _run_test()
        elif args.subparser in ["query", "q"]:
            q = args.query
            if "*" in q:
                gene, q = q.split("*", maxsplit=1)
            else:
                gene, q = q, ""
            if args.gene:
                gene = args.gene
            db_file = script_path("aldy.resources.genes/{}.yml".format(gene.lower()))
            if not os.path.exists(db_file):
                db_file = gene
            with open(db_file):  # Check if file exists
                pass
            query(Gene(db_file, genome=args.genome), q)
        elif args.subparser == "profile":
            params = {}
            if args.param:
                for pl in args.param:
                    for p in pl:
                        if "=" not in p:
                            raise AldyException(f"Invalid parameter {p}")
                        k, v = p.split("=", 1)
                        params[k.replace("-", "_")] = v
            try:
                p = Profile.get_sam_profile_data(
                    args.file,
                    cn_region=parse_cn_region(args.cn_neutral_region),
                    genome=args.genome,
                    params=params,
                )
                print(yaml.dump(p, default_flow_style=None))
            except AldyException as ex:
                log.error("ERROR: {}", str(ex))
        elif args.subparser in ["genotype", "g"]:
            # Prepare the output file
            output = args.output
            if output == "-":
                output = sys.stdout
            elif output:
                output = open(output, "w")
            _genotype(args.gene, output, args)
            if output and output != sys.stdout:
                output.close()
        else:
            raise AldyException("Invalid sub-command " + args.subparser)
    except IOError as ex:
        if ex.filename is not None:
            log.critical("ERROR: {} cannot be accessed", ex.filename)
        else:
            log.critical("ERROR: {} cannot be accessed", str(ex))
        log.debug(ex)
        log.debug(traceback.format_exc())
        exit(1)
    except SystemExit as ex:
        log.debug(repr(ex))
        log.debug(traceback.format_exc())
        exit(ex.code)
    except Exception as ex:
        gene = "" if "gene" not in args else f"gene= {args.gene}, "
        log.critical(f"ERROR: {gene}file= {args.file if 'file' in args else '-'}")
        log.critical(repr(ex))
        log.warn(traceback.format_exc())
        exit(1)
    except:  # noqa
        gene = "" if "gene" not in args else f"gene= {args.gene}, "
        exc = sys.exc_info()[0]
        log.critical(f"ERROR: {gene}file= {args.file if 'file' in args else '-'}")
        log.critical("Unrecoverable error: {}", repr(exc))
        log.warn(traceback.format_exc())
        exit(1)


def _get_args(argv):
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
        "file", nargs="?", help="Input file in SAM, BAM or CRAM format."
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
        help=td(
            """Sequencing profile. The following profiles are supported:
               - illumina (also: wgs),
               - pgrnseq-v1 (also: pgx1),
               - pgrnseq-v2 (also: pgx2),
               - pgrnseq-v3 (also: pgx3),
               - 10x,
               - pacbio-hifi-targeted,
               - pacbio-hifi-targeted-twist,
               - exome (also: wxs, wes).
               You can also provide a custom SAM/BAM file as a profile.
               Please check documentation for more details."""
        ),
    )
    genotype_parser.add_argument(
        "--reference",
        "-r",
        default=None,
        help="Genome reference used for reading CRAM files",
    )
    genotype_parser.add_argument(
        "--genome",
        default=None,
        help="SAM/BAM reference genome (hg19 or hg38; none for auto-detection)",
    )
    genotype_parser.add_argument(
        "--cn-neutral-region",
        "-n",
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
               - cbc (Google OR-Tools/CBC)
               - any (attempts to use Gurobi and then CBC).
               Default is "any"."""
        ),
    )
    genotype_parser.add_argument(
        "--debug",
        default=None,
        help="Create a directory that will contain the debug information "
        + "and core dumps.",
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
    genotype_parser.add_argument("--log", "-l", default=None, help="Log file location")
    genotype_parser.add_argument(
        "--multiple-warn-level",
        "-W",
        default=1,
        type=int,
        help="Show warning if multiple solutions are found. "
        + "Can be 1 (warn after genotyping) or 2 "
        + "(also warn if there are multiple major solutions)."
        + "Default is 1 (warn after the genotyping).",
    )
    genotype_parser.add_argument(
        "--simple",
        action="store_true",
        default=False,
        help=td("""Print one-line result per gene (debug). Default is off."""),
    )
    genotype_parser.add_argument("--param", action="append", nargs="+")

    _ = subparsers.add_parser(
        "test",
        parents=[base],
        help="Run Aldy test suite. Recommended prior to the first use",
    )

    _ = subparsers.add_parser("license", parents=[base], help="Show Aldy license")

    show_parser = subparsers.add_parser(
        "query",
        aliases=["q"],
        parents=[base],
        help="Query database definitions for a given gene.",
    )
    show_parser.add_argument("--gene", "-g", default=None, help="Gene file.")
    show_parser.add_argument("--genome", default="hg19", help="Reference genome.")
    show_parser.add_argument("query", help="Gene or allele to show.")

    profile_parser = subparsers.add_parser(
        "profile",
        parents=[base],
        help=td(
            """Generate a sequencing profile for a SAM/BAM/CRAM file.
               Please check the documentation for more details."""
        ),
    )
    profile_parser.add_argument("file", nargs="?", help="SAM/BAM/CRAM file.")
    profile_parser.add_argument(
        "--cn-neutral-region",
        "-n",
        default=None,
        help=td(
            """Copy-number neutral region in the format chromosome:start-end
               (e.g. chr1:10000-20000).
               Default is CYP2D8 region within hg19 (22:42547463-42548249)."""
        ),
    )
    profile_parser.add_argument(
        "--genome",
        default=None,
        help="SAM/BAM reference genome (hg19 or hg38; hg19 by default)",
    )
    profile_parser.add_argument("--param", action="append", nargs="+")

    _ = subparsers.add_parser(
        "help", parents=[base], help="Show program usage and exit."
    )

    return parser, parser.parse_args(argv)


def _print_licence():
    """
    Print Aldy license.
    """
    with open(script_path("aldy.resources/LICENSE.rst")) as f:
        for ll in f:
            print(ll.strip())


def _genotype(gene: str, output: Optional[Any], args) -> None:
    """
    Genotype a file.

    :raise: :py:class:`aldy.common.AldyException` if `cn_region` is invalid.
    """

    cn_region = parse_cn_region(args.cn_neutral_region)
    cn_solution = args.cn
    if cn_solution:
        cn_solution = cn_solution.split(",")

    def run(debug):
        log.trace(
            "\n[main] arguments= {}",
            " ".join(k + "=" + str(v) for k, v in vars(args).items() if k is not None),
        )
        log.info("Genotyping sample {}...", os.path.basename(args.file))

        if args.profile in ["exome", "wxs"]:
            log.warn("WARNING: Copy-number calling is not available for exome data.")
            log.warn(
                "WARNING: Aldy will NOT be able to detect gene duplications, "
                + "deletions and fusions."
            )
            log.warn(
                "WARNING: Calling of alleles that are defined by non-exonic mutations "
                + "is not available."
            )
            log.warn("         Results might not be biologically relevant!")

        try:
            params = {
                k: v
                for k, v in vars(args).items()
                if k in ["solver", "reference", "multiple_warn_level", "genome"]
            }
            if args.param:
                for pl in args.param:
                    for p in pl:
                        if "=" not in p:
                            raise AldyException(f"Invalid parameter {p}")
                        k, v = p.split("=", 1)
                        params[k.replace("-", "_")] = v
            _ = genotype(
                gene_db=gene,
                sam_path=args.file,
                profile_name=args.profile,
                output_file=output,
                cn_region=cn_region,
                cn_solution=cn_solution,
                report=True,
                is_simple=args.simple,
                debug=debug,
                **{k: v for k, v in params.items() if v is not None},
            )
        except AldyException as ex:
            log.critical(
                f"ERROR: gene= {gene}, profile= {args.profile}, file= {args.file}"
            )
            log.error(str(ex))

    if args.log:
        fh = logbook.FileHandler(
            args.log, mode="w", bubble=True, level="TRACE"  # type: ignore
        )
        fh.formatter = lambda record, _: record.message  # type: ignore
        fh.push_application()
    if args.debug:
        with tempfile.TemporaryDirectory() as tmp:
            prefix = None
            try:
                prefix = f"{tmp}/{os.path.splitext(os.path.basename(args.file))[0]}"
                log_output = f"{prefix}.log"
                fh = logbook.FileHandler(
                    log_output, mode="w", bubble=True, level="TRACE"  # type: ignore
                )
                fh.formatter = lambda record, _: "[{}:{}/{}] {}".format(  # type: ignore
                    record.level_name[0],
                    os.path.splitext(os.path.basename(record.filename))[0],
                    record.func_name,
                    record.message,
                )
                fh.push_application()
                log.debug(f"Using {tmp} as temporary debug directory")
                run(prefix)
            finally:
                if prefix:
                    log.info("Preparing debug archive...")
                    with open(f"{prefix}.yml", "w") as f:
                        yaml.Dumper.ignore_aliases = lambda *_: True  # type: ignore
                        yaml.dump(common.json, f, default_flow_style=None)
                    os.system(f"tar czf {args.debug}.tar.gz -C {tmp} .")
    else:
        run(None)


def _run_test() -> None:
    """
    Run the Aldy test suite.
    """

    pytest.main(["--pyargs", "aldy"])


def console():
    main(sys.argv[1:])


if __name__ == "__main__":
    console()
