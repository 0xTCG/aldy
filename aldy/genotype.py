# 786
# Aldy source: genotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Optional, Any, Set, Dict
import os
import sys
import pkg_resources
import datetime
import time

from . import sam
from . import cn
from . import major
from . import minor
from . import diplotype
from . import solutions
from .common import (
    colorize,
    log,
    script_path,
    json,
    AldyException,
    SOLUTION_PRECISION,
)
from .profile import Profile
from .gene import Gene, GRange
from .diplotype import OUTPUT_COLS
from .lpinterface import model as lp_model


def genotype(
    gene_db: str,
    sam_path: str,
    profile_name: Optional[str],
    output_file: Optional[Any] = sys.stdout,
    cn_region: Optional[GRange] = None,
    cn_solution: Optional[List[str]] = None,
    solver: str = "any",
    reference: Optional[str] = None,
    debug: Optional[str] = None,
    multiple_warn_level: int = 1,
    report: bool = False,
    genome=None,
    is_simple: bool = False,
    **params,
) -> Dict[str, List[solutions.MinorSolution]]:
    """Genotype a sample.

    :param gene_db: Gene name (if it is shipped with Aldy)
        or the location of the gene database in YAML format.
    :param sam_path: Location of SAM/BAM/CRAM file that is to be analyzed.
    :param profile_name: Coverage profile (e.g. WGS).
        `None` if `cn_solution` is provided.
    :param output_file: Location of the output file. Use `None` for no output.
        Default: `sys.stdout`.
    :param cn_region: Copy-number neutral region.
        Default: None (uses the provided CYP2D8 region).
    :param cn_solution: List of the copy number configurations.
        Copy-number detection will not run if this parameter is set.
        Default: `None`.
    :param solver: ILP solver (see :py:mod:`aldy.lpinterface` for supported solvers).
    :param reference: Reference genome (for reading CRAM files).
        Default: `None`.
    :param debug: Prefix for debug and core dump files.
        Default: `None` (no debug information).
    :param multiple_warn_level: Warning level (1 for optimal solutions, 2 for major
        solutions and 3 for CN solutions).
        Default: 1 (warn only on multiple optimal solutions).
    :param report: If set, write the solution summary to the stderr.
        Default: `False`.
    :param genome: Reference genome (e.g., hg19 or hg38).
        Default: `None` (auto-detect).
    :param is_simple: Use simple output format.
        Default: `False`.
    :param params: Model parameters. See :py:mod:`aldy.profile` for details.
    """

    t1 = time.time()
    log.debug("[genotype] gene={}; start={}", gene_db, datetime.datetime.now())

    # Test the existence of LP solver
    _ = lp_model("init", solver)

    with open(sam_path):  # Check if file exists
        pass
    kind, g = sam.detect_genome(sam_path)
    if genome is None:
        genome = g
        if not genome:
            log.warn("WARNING: Cannot detect genome, defaulting to hg19.")
            genome = "hg19"
        log.debug(f"[genotype] reference= {genome}")
    if genome not in ["hg19", "hg38"]:
        raise AldyException(f"Unknown genome {genome}")

    avail_genes = []
    if gene_db == "all":
        avail_genes = pkg_resources.resource_listdir("aldy.resources", "genes")
        avail_genes = [
            i[:-4]
            for i in avail_genes
            if len(i) > 4 and i.endswith(".yml") and not i.startswith("pharma-")
        ]
        avail_genes = sorted(avail_genes)
    elif gene_db == "pharmacoscan":
        avail_genes = pkg_resources.resource_listdir(
            "aldy.resources.genes", "pharmacoscan"
        )
        avail_genes = [
            f"pharmacoscan/{i[:-4]}" for i in avail_genes if i.endswith(".yml")
        ]
        avail_genes = sorted(avail_genes)
    else:
        avail_genes = gene_db.lower().split(",")
    if len(avail_genes) != 1:
        res: Dict = {}
        for a in avail_genes:
            log.warn("=" * 50)
            log.warn("Gene {}", a.upper())
            try:
                res = {
                    **res,
                    **genotype(
                        a,
                        sam_path,
                        profile_name,
                        output_file,
                        cn_region,
                        cn_solution,
                        solver,
                        reference,
                        debug,
                        multiple_warn_level,
                        report,
                        genome,
                        is_simple,
                        **params,
                    ),
                }
            except AldyException as ex:
                log.error(f"Failed gene {a.upper()}")
                log.error(f"Message: {str(ex)}")
            log.warn("")
        return res
    else:
        db_file = script_path(
            "aldy.resources.genes/{}.yml".format(avail_genes[0].lower())
        )

    # Load the gene specification
    if os.path.exists(db_file):
        gene_db = db_file
    with open(gene_db):  # Check if file exists
        pass
    gene = Gene(gene_db, genome=genome)

    if "/wgs/" in sam_path:
        profile_name = "wgs"
    elif "/10x/" in sam_path:
        profile_name = "10x"
    elif "/pgx3/" in sam_path:
        profile_name = "pgx3"

    if profile_name in ["exome", "wxs", "wes"]:
        gene.do_copy_number = False
        profile_name = "illumina"
    elif profile_name == "wgs":
        profile_name = "illumina"
    elif profile_name == "pgrnseq-v1":
        profile_name = "pgx1"
    elif profile_name == "pgrnseq-v2":
        profile_name = "pgx2"
    elif profile_name == "pgrnseq-v3":
        profile_name = "pgx3"

    if kind == "vcf":
        log.warn("WARNING: Using VCF file. Copy-number calling is not available.")
        profile = Profile("user_provided", cn_solution=["1", "1"], **params)
        sample = sam.Sample(gene, profile, sam_path, debug=debug)
    else:
        if cn_solution:
            profile = Profile("user_provided", cn_solution=cn_solution, **params)
        elif kind != "dump":
            if not profile_name:
                raise AldyException("Profile not provided")
            profile = Profile.load(gene, profile_name, cn_region, **params)
        else:
            profile = None
        sample = sam.Sample(gene, profile, sam_path, reference, debug)
    profile = sample.profile  # if loaded for a dump
    assert profile, "Profile not set"
    if kind == "dump":
        profile.update(params)

    json[gene.name].update({"sample": sample.name})
    is_vcf = output_file and output_file.name.endswith(".vcf")
    is_simple = is_simple or (
        output_file is not None and output_file.name.endswith(".simple")
    )
    is_aldy = output_file and not (is_vcf or is_simple)
    if is_simple:
        print(
            sample.name,
            gene.name,
            sep="\t",
            end="\t",
            file=output_file,
        )
    if kind != "vcf":
        avg_cov = sample.coverage.average_coverage()
        if profile.cn_region and avg_cov < 2:
            if is_simple:
                print(file=output_file)
            raise AldyException(
                f"Average coverage of {avg_cov:.2f} for gene {gene.name} is too low; "
                + f"skipping gene {gene.name}. "
                + f"Please ensure that {gene.name} is present in the input SAM/BAM."
            )
        elif profile.cn_region and avg_cov < 20:
            log.warn(
                f"Average sample coverage is {avg_cov}. "
                + "We recommend at least 20x coverage for the best results."
            )

    # Get copy-number solutions
    cn_sols = cn.estimate_cn(
        gene,
        profile,
        sample.coverage,
        solver=solver,
        debug=debug,
    )

    # Add SLACK to each score to avoid division by zero or other numerical issues
    # when the optimum is close to zero.
    # HACK: Value of 1 (1 gene copy) seems to be reasonable, but is not
    # theoretically guaranteed to be better than any other value.
    SLACK = 1

    # Get major solutions and pick the best one
    if multiple_warn_level >= 3 and len(cn_sols) > 1:
        log.warn("WARNING: multiple gene structures found!")
    log.info(f"Potential {gene.name} gene structures for {sample.name}:")
    major_sols: list = []
    cn_sols = sorted(cn_sols, key=lambda m: (int(1000 * m.score), m._solution_nice()))
    if len(cn_sols) == 0:
        if is_simple:
            print(file=output_file)
        raise AldyException("No solutions found!")
    min_cn_score = min(cn_sols, key=lambda m: m.score).score
    for i, cn_sol in enumerate(cn_sols):
        conf = 100 * (min_cn_score + SLACK) / (cn_sol.score + SLACK)
        log.info(f"  {i + 1:2}: {cn_sol._solution_nice()} (confidence: {conf:.0f}%)")
    log.debug("*" * 80)

    # return []

    for i, cn_sol in enumerate(cn_sols):
        sols = major.estimate_major(
            gene,
            sample.coverage,
            cn_sol,
            solver=solver,
            identifier=i,
            debug=debug,
        )
        for s in sols:
            s.score += cn_sol.score - min_cn_score
        major_sols += sols
    if len(major_sols) == 0:
        if is_simple:
            print(file=output_file)
        raise AldyException("No major solutions found!")

    major_sols = [
        solutions.MajorSolution(
            m.score,  # * ((m.cn_solution.score + SLACK) / (min_cn_score + SLACK)),
            m.solution,
            m.cn_solution,
            m.added,
        )
        for m in major_sols
    ]
    min_major_score = min(major_sols, key=lambda m: m.score).score
    major_sols = sorted(
        [
            m
            for m in major_sols
            if m.score - min_major_score - profile.gap < SOLUTION_PRECISION
        ],
        key=lambda m: (int(1000 * m.score), m._solution_nice()),
    )

    if multiple_warn_level >= 2 and len(major_sols) > 1:
        log.warn("WARNING: multiple major solutions found!")
    log.info(f"Potential major {gene.name} star-alleles for {sample.name}:")
    for i, major_sol in enumerate(major_sols):
        conf = 100 * (min_major_score + SLACK) / (major_sol.score + SLACK)
        log.info(f"  {i + 1:2}: {major_sol._solution_nice()} (confidence: {conf:.0f}%)")
    log.debug("*" * 80)

    minor_sols = []
    for m in minor.estimate_minor(
        gene,
        sample.coverage,
        major_sols,
        solver,
        max_solutions=profile.max_minor_solutions,
    ):
        n = solutions.MinorSolution(
            m.score
            * ((m.major_solution.cn_solution.score + SLACK) / (min_cn_score + SLACK)),
            m.solution,
            m.major_solution,
        )
        n.set_diplotype(m.get_diplotype())
        minor_sols.append(n)

    if len(minor_sols) == 0:
        if is_simple:
            print(file=output_file)
        raise AldyException(
            "Aldy could not phase any major solution.\n"
            + "Possible solutions:\n"
            + " - Check the coverage. Extremely low coverage prevents Aldy from "
            + "calling star-alleles.\n"
            + " - Run with --debug parameter and notify the authors of Aldy."
        )
    min_minor_score = min(minor_sols, key=lambda m: m.score).score
    minor_sols = sorted(
        [
            m
            for m in minor_sols
            if m.score - min_minor_score - profile.gap < SOLUTION_PRECISION
        ],
        key=lambda m: (int(1000 * m.score), m._solution_nice()),
    )
    log.debug("*" * 80)

    if multiple_warn_level >= 1 and len(minor_sols) > 1:
        log.warn("WARNING: multiple optimal solutions found!")
    log.info(
        f"{{}} {gene.name} star-alleles for {sample.name}:",
        "Best" if len(minor_sols) == 1 else "Potential",
    )
    if is_aldy:
        print("#" + "\t".join(OUTPUT_COLS), file=output_file)
    simple = [sample.name, gene.name]
    for i, minor_sol in enumerate(minor_sols):
        conf = 100 * (min_minor_score + SLACK) / (minor_sol.score + SLACK)
        log.info(
            f"  {i + 1:2}: {minor_sol.get_major_diplotype()}"
            + f" (confidence={conf:.0f}%)"
        )
        log.info(f"      Minor alleles: {minor_sol._solution_nice()}")
        simple += [
            minor_sol.get_major_diplotype().replace(" ", ""),
            minor_sol.get_minor_diplotype(legacy=True).replace(" ", ""),
        ]
        if is_aldy:
            print(f"#Solution {i + 1}: {minor_sol._solution_nice()}", file=output_file)
            diplotype.write_decomposition(
                sample.name, gene, sample.coverage, i + 1, minor_sol, output_file
            )
        elif is_simple:
            print(simple[-2], simple[-1], sep="\t", end="\t", file=output_file)
    if is_vcf:
        diplotype.write_vcf(sample.name, gene, sample.coverage, minor_sols, output_file)
    elif is_simple:
        print(file=output_file)
    log.debug("[simple]\t{}", "\t".join(simple))

    if report:
        log.info(colorize(f"{gene.name} results:"))
        reported: Set[str] = set()
        for r in minor_sols:
            st = f"  - {r.get_major_diplotype()}"
            if st not in reported:
                log.info(colorize(st))
                log.info(f"    Minor: {r.get_minor_diplotype()}")
                log.info(f"    Legacy notation: {r.get_minor_diplotype(legacy=True)}")
                reported.add(st)
                # Output the activity
                for m in r.solution:
                    mn = [
                        ma.minors[m.minor]
                        for ma in gene.alleles.values()
                        if m.minor in ma.minors
                    ]
                    msg = f"Estimated activity for {m.major_repr()}"
                    if len(mn) == 1 and mn[0].activity:
                        msg = f"{msg}: {mn[0].activity}"
                        if mn[0].evidence:
                            msg = f"{msg} (evidence: {mn[0].evidence})"
                        if mn[0].pharmvar:
                            msg = f"{msg}; see {mn[0].pharmvar} for details"
                    else:
                        msg = f"{msg}: unknown"
                    log.info(f"    {msg}")
                # Calculate coverage of each mutation in a solution
                if multiple_warn_level >= 2:
                    for m, mc, sc in r.get_mutation_coverages(sample.coverage):
                        if abs(mc - sc) > 0.5:
                            log.warn(
                                "    Warning: Coverage of {} is not in line with "
                                "the prediction (predicted: {:.1f}, observed: {:.1f})",
                                gene.get_rsid(m),
                                mc,
                                sc,
                            )

    if (
        any(len(m.major_solution.added) > 0 for m in minor_sols)
        and not sample.is_long_read
    ):
        novels = {
            gene.get_rsid(mm) for m in minor_sols for mm in m.major_solution.added
        }
        log.warn(
            "WARNING: mutations {} suggest presence of a novel major star-allele.\n"
            "However, exact allele contents cannot be determined without the "
            "long-read phasing data. The reported assignments of novel mutations "
            "are for that reason currently randomly assigned.",
            ", ".join(sorted(novels)),
        )

    log.debug("[genotype] gene={}; took={}", gene_db, time.time() - t1)
    return {gene_db: minor_sols}
