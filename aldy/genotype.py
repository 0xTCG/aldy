# 786
# Aldy source: genotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Optional, Any, Set

import os
import sys

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
from .gene import Gene, GRange
from .diplotype import OUTPUT_COLS
from .lpinterface import model as lp_model


def genotype(
    gene_db: str,
    sam_path: str,
    profile: str,
    output_file: Optional[Any] = sys.stdout,
    cn_region: Optional[GRange] = None,
    cn_solution: Optional[List[str]] = None,
    threshold: float = 0.5,
    fusion_penalty: float = cn.LEFT_FUSION_PENALTY,
    solver: str = "any",
    reference: Optional[str] = None,
    gap: int = 0,
    max_minor_solutions: int = 1,
    debug: Optional[str] = None,
    multiple_warn_level: int = 1,
    phase: Optional[str] = None,
    report: bool = False,
    genome=None,
    min_cov: Optional[str] = None,
) -> List[solutions.MinorSolution]:
    """
    Genotype a sample.

    Returns:
        list[:obj:`aldy.solutions.MinorSolution`]: List of genotype solutions.

    Args:
        gene_db (str):
            Gene name (if it is located in the Aldy's gene database)
            or the location of a gene database in YML format.
        sam_path (str):
            Location of SAM/BAM/CRAM file that is to be analyzed.
        profile (str):
            Coverage profile (e.g. 'illumina').
            Can be ``None`` if ``cn_solution`` is provided.
        cn_region (:obj:`aldy.common.GRange`, optional):
            Copy-number neutral region.
            Can be ``None`` (will use the default CYP2D8 region or ``None``
            if ``cn_solution`` is provided).
        output_file (file, optional):
            Location of an output decomposition file.
            Provide ``None`` for no output.
            Default is ``sys.stdout``.
        cn_solution (list[str], optional):
            User-specified list of the copy number configurations.
            Copy-number detection solver will not run is this parameter is set.
            Default is ``None``.
        threshold (float):
            Mutation filtering threshold.
            Default is 0.5 (for 50%).
        fusion_penalty (float):
            Fusion penalty. Use higher values to avoid fusions.
            Default is 0.1.
        solver (str):
            ILP solver. Check :obj:`aldy.lpinterface` for the list of available solvers.
            Default is ``'any'``.
        phase (str, optional):
            Location of HapCUT or HapTree-compatible phased blocks.
        reference (str, optional):
            A reference genome for reading CRAM files.
            Default is ``None``.
        gap (float):
            Relative optimality gap. Use non-zero values to allow non-optimal solutions.
            Default is 0 (reports only optimal solutions).
        max_minor_solutions (int):
            Maximum number of minor solutions to report for each major solution.
            Default is 1.
        debug (str, optional):
            Prefix for debug information and core dump files.
            ``None`` for no debug information.
            Default is ``None``.
        min_cov (float):
            Minimum mutation read coverage.
            Default is 1.
    Raises:
        :obj:`aldy.common.AldyException` if the average coverage is too low (below 2).
    """

    # Test the existence of LP solver
    _ = lp_model("init", solver)

    sample_name = os.path.splitext(os.path.basename(sam_path))[0]
    with open(sam_path):  # Check if file exists
        pass
    kind, g = sam.detect_genome(sam_path)
    # assert kind in ["vcf", "sam", "dump"]
    if genome is None:
        genome = g
        if not genome:
            log.warn("WARNING: Cannot detect genome, defaulting to hg19.")
            genome = "hg19"
        log.debug(f"[genotype] reference= {genome}")
    if genome not in ["hg19", "hg38"]:
        raise AldyException(f"Unknown genome {genome}")

    # Load the gene specification
    db_file = script_path("aldy.resources.genes/{}.yml".format(gene_db.lower()))
    if os.path.exists(db_file):
        gene_db = db_file
    with open(gene_db):  # Check if file exists
        pass
    gene = Gene(gene_db, genome=genome)

    if not cn_region:
        cn_region = sam.DEFAULT_CN_NEUTRAL_REGION[genome]
    if profile in ["exome", "wxs"]:
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
        cn_region = None
        cn_solution = ["1", "1"]
        profile = "illumina"
        if not min_cov:
            min_cov = "5.0"
    elif profile == "wgs":
        profile = "illumina"
    elif profile == "pgrnseq-v1":
        profile = "pgx1"
    elif profile == "pgrnseq-v2":
        profile = "pgx2"
    elif profile == "pgrnseq-v3":
        profile = "pgx3"

    if kind == "vcf":
        log.warn("WARNING: Using VCF file. Copy-number calling is not available.")
        cn_region = None
        cn_solution = ["1", "1"]
        sample = sam.Sample(gene=gene, vcf_path=sam_path, debug=debug)
    else:
        sample = sam.Sample(
            gene=gene,
            sam_path=sam_path,
            threshold=threshold,
            profile=profile,
            reference=reference,
            cn_region=None if cn_solution else cn_region,
            debug=debug,
            min_cov=float(min_cov) if min_cov else 1.0,
        )
        avg_cov = sample.coverage.average_coverage()
        if cn_region and avg_cov < 2:
            raise AldyException(
                f"Average coverage of {avg_cov:.2f} for gene {gene.name} is too low; "
                + f"skipping gene {gene.name}. "
                + f"Please ensure that {gene.name} is present in the input SAM/BAM."
            )
        elif cn_region and avg_cov < 20:
            log.warn(
                f"Average sample coverage is {avg_cov}. "
                + "We recommend at least 20x coverage for the best results."
            )

    json.update({"sample": os.path.basename(sam_path).split(".")[0], "gene": gene.name})

    # Get copy-number solutions
    cn_sols = cn.estimate_cn(
        gene,
        sample.coverage,
        solver=solver,
        gap=gap,
        user_solution=cn_solution,
        fusion_penalty=fusion_penalty,
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
    log.info(f"Potential {gene.name} gene structures for {sample_name}:")
    major_sols: list = []
    cn_sols = sorted(cn_sols, key=lambda m: (int(1000 * m.score), m._solution_nice()))
    if len(cn_sols) == 0:
        raise AldyException("No solutions found!")
    min_cn_score = min(cn_sols, key=lambda m: m.score).score
    for i, cn_sol in enumerate(cn_sols):
        conf = 100 * (min_cn_score + SLACK) / (cn_sol.score + SLACK)
        log.info(f"  {i + 1:2}: {cn_sol._solution_nice()} (confidence: {conf:.0f}%)")
    log.debug("*" * 80)

    phases = sam.load_phase(gene, phase) if phase else None

    for i, cn_sol in enumerate(cn_sols):
        major_sols += major.estimate_major(
            gene,
            sample.coverage,
            cn_sol,
            solver=solver,
            gap=gap,
            identifier=i,
            debug=debug,
        )
    if len(major_sols) == 0:
        raise AldyException("No major solutions found!")

    major_sols = [
        solutions.MajorSolution(
            m.score * ((m.cn_solution.score + SLACK) / (min_cn_score + SLACK)),
            m.solution,
            m.cn_solution,
            m.added,
        )
        for m in major_sols
    ]
    min_major_score = min(major_sols, key=lambda m: m.score).score
    major_sols = sorted(
        [m for m in major_sols if m.score - min_major_score - gap < SOLUTION_PRECISION],
        key=lambda m: (int(1000 * m.score), m._solution_nice()),
    )

    if multiple_warn_level >= 2 and len(major_sols) > 1:
        log.warn("WARNING: multiple major solutions found!")
    log.info(f"Potential major {gene.name} star-alleles for {sample_name}:")
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
        max_solutions=max_minor_solutions,
        debug=debug,
        phases=phases,
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
        raise AldyException(
            "Aldy could not phase any major solution.\n"
            + "Possible solutions:\n"
            + " - Check the coverage. Extremely low coverage prevents Aldy from "
            + "calling star-alleles.\n"
            + " - Run with --debug parameter and notify the authors of Aldy."
        )
    min_minor_score = min(minor_sols, key=lambda m: m.score).score
    minor_sols = sorted(
        [m for m in minor_sols if m.score - min_minor_score - gap < SOLUTION_PRECISION],
        key=lambda m: (int(1000 * m.score), m._solution_nice()),
    )
    log.debug("*" * 80)

    if multiple_warn_level >= 1 and len(minor_sols) > 1:
        log.warn("WARNING: multiple optimal solutions found!")
    log.info(
        f"{{}} {gene.name} star-alleles for {sample_name}:",
        "Best" if len(minor_sols) == 1 else "Potential",
    )
    is_vcf = output_file and output_file.name[-4:] == ".vcf"
    if output_file and not is_vcf:
        print("#" + "\t".join(OUTPUT_COLS), file=output_file)
    for i, minor_sol in enumerate(minor_sols):
        conf = 100 * (min_minor_score + SLACK) / (minor_sol.score + SLACK)
        log.info(
            f"  {i + 1:2}: {minor_sol.get_major_diplotype()}"
            + f" (confidence={conf:.0f}%)"
        )
        log.info(f"      Minor alleles: {minor_sol._solution_nice()}")
        if output_file and not is_vcf:
            print(f"#Solution {i + 1}: {minor_sol._solution_nice()}", file=output_file)
            diplotype.write_decomposition(
                sample_name, gene, i + 1, minor_sol, output_file
            )
    if is_vcf:
        diplotype.write_vcf(sample_name, gene, sample.coverage, minor_sols, output_file)

    if report:
        log.info(colorize(f"{gene.name} results:"))
        reported: Set[str] = set()
        for r in minor_sols:
            s = f"  - {r.get_major_diplotype()}"
            if s not in reported:
                log.info(colorize(s))
                log.info(f"    Minor: {r.get_minor_diplotype()}")
                log.info(f"    Legacy notation: {r.get_minor_diplotype(legacy=True)}")
                reported.add(s)
                # Calculate coverage of each mutation in a solution
                for m, mc, sc in r.get_mutation_coverages(sample.coverage):
                    if abs(mc - sc) > 0.5:
                        log.warn(
                            "    Warning: Coverage of {} is not in line with "
                            "the prediction\n"
                            "             (predicted: {:.1f}, observed: {:.1f})",
                            gene.get_rsid(m),
                            mc,
                            sc,
                        )
    if any(len(m.major_solution.added) > 0 for m in minor_sols) and not phase:
        novels = {
            gene.get_rsid(mm) for m in minor_sols for mm in m.major_solution.added
        }
        log.warn(
            "WARNING: mutations {} suggest presence of a novel major star-allele."
            + "\nHowever, such alleles cannot be determined without phasing data."
            + "\nPlease provide --phase parameter for Aldy to accurately call novel "
            + "major star-alleles.\nThe above-reported assignments of these mutations "
            + "are random.",
            ", ".join(sorted(novels)),
        )

    return minor_sols
