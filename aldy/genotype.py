# 786
# Aldy source: genotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import List, Optional, Any

import os
import sys
import textwrap

from . import sam
from . import cn
from . import major
from . import minor
from . import diplotype
from . import solutions

from .common import log, script_path, json_print, AldyException, SOLUTION_PRECISION
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
    phase: bool = False,
    reference: Optional[str] = None,
    gap: int = 0,
    max_minor_solutions: int = 1,
    debug: Optional[str] = None,
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
        phase (bool):
            Perform basic rudimentary phasing of the reads
            to aid the genotyping process.
            DEPRECATED: currently does nothing.
            Default is ``False``.
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

    Raises:
        :obj:`aldy.common.AldyException` if the average coverage is too low (below 2).
    """

    # Test the existence of LP solver
    _ = lp_model("init", solver)

    # Load the gene specification
    db_file = script_path("aldy.resources.genes/{}.yml".format(gene_db.lower()))
    if os.path.exists(db_file):
        gene_db = db_file
    with open(gene_db):  # Check if file exists
        pass
    gene = Gene(gene_db)

    sample_name = os.path.splitext(os.path.basename(sam_path))[0]
    with open(sam_path):  # Check if file exists
        pass
    if not cn_region:
        cn_region = sam.DEFAULT_CN_NEUTRAL_REGION
    sample = sam.Sample(
        sam_path=sam_path,
        gene=gene,
        threshold=threshold,
        profile=profile,
        phase=phase,
        reference=reference,
        cn_region=cn_region,
        debug=debug,
    )

    avg_cov = sample.coverage.average_coverage()
    if avg_cov < 2:
        raise AldyException(
            f"Average coverage of {avg_cov:.2f} for gene {gene.name} is too low; "
            + f"skipping gene {gene.name}. "
            + f"Please ensure that {gene.name} is present in the input SAM/BAM."
        )
    elif avg_cov < 20:
        log.warn(
            f"Average sample coverage is {avg_cov}. "
            + "We recommend at least 20x coverage for the best results."
        )

    json_print(debug, f'"{os.path.basename(sam_path).split(".")[0]}": {{')

    # Get copy-number solutions
    json_print(debug, '  "cn": {')
    cn_sols = cn.estimate_cn(
        gene,
        sample.coverage,
        solver=solver,
        gap=gap,
        user_solution=cn_solution,
        fusion_penalty=fusion_penalty,
        debug=debug,
    )
    log.debug("*" * 80)
    json_print(debug, "  },")

    # Add SLACK to each score to avoid division by zero or other numerical issues
    # when the optimum is close to zero.
    # HACK: Value of 1 (1 gene copy) seems to be reasonable, but is not
    # theoretically guaranteed to be better than any other value.
    SLACK = 1

    # Get major solutions and pick the best one
    log.info(f"Potential {gene.name} copy number configurations for {sample_name}:")
    major_sols: list = []
    json_print(debug, '  "major": [', end="")
    cn_sols = sorted(cn_sols, key=lambda m: (int(1000 * m.score), m._solution_nice()))
    if len(cn_sols) == 0:
        raise AldyException("No solutions found!")
    min_cn_score = min(cn_sols, key=lambda m: m.score).score
    for i, cn_sol in enumerate(cn_sols):
        log.info(f"  {i + 1:2}: {cn_sol._solution_nice()}")
        conf = (min_cn_score + SLACK) / (cn_sol.score + SLACK)
        log.info(f"      Confidence: {conf:.2f} (score = {cn_sol.score:.2f})")
    log.info("")

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
    log.debug("*" * 80)
    json_print(debug, "],")
    if len(major_sols) == 0:
        raise AldyException("No major solutions found!")

    major_sols = [
        solutions.MajorSolution(
            m.score * ((m.cn_solution.score + SLACK) / (min_cn_score + SLACK)),
            m.solution,
            m.cn_solution,
        )
        for m in major_sols
    ]
    min_major_score = min(major_sols, key=lambda m: m.score).score
    major_sols = sorted(
        [m for m in major_sols if m.score - min_major_score - gap < SOLUTION_PRECISION],
        key=lambda m: (int(1000 * m.score), m._solution_nice()),
    )

    log.info(f"Potential major {gene.name} star-alleles for {sample_name}:")
    for i, major_sol in enumerate(major_sols):
        log.info(f"  {i + 1:2}: {major_sol._solution_nice()}")
        conf = (min_major_score + SLACK) / (major_sol.score + SLACK)
        log.info(f"      Confidence: {conf:.2f} (score = {major_sol.score:.2f})")
    log.info("")

    json_print(debug, '  "minor": [', end="")
    minor_sols = []
    for m in minor.estimate_minor(
        gene,
        sample.coverage,
        major_sols,
        solver,
        max_solutions=max_minor_solutions,
        debug=debug,
    ):
        n = solutions.MinorSolution(
            m.score
            * ((m.major_solution.cn_solution.score + SLACK) / (min_cn_score + SLACK)),
            m.solution,
            m.major_solution,
        )
        n.diplotype = m.diplotype
        minor_sols.append(n)

    if len(minor_sols) == 0:
        raise AldyException(
            "Aldy could not phase any major solution. "
            + "Please run with  --debug parameter and notify the authors of Aldy."
        )
    min_minor_score = min(minor_sols, key=lambda m: m.score).score
    minor_sols = sorted(
        [m for m in minor_sols if m.score - min_minor_score - gap < SOLUTION_PRECISION],
        key=lambda m: (int(1000 * m.score), m._solution_nice()),
    )
    log.debug("*" * 80)
    json_print(debug, "  ],")

    log.info(f"Best {gene.name} star-alleles for {sample_name}:")
    is_vcf = output_file and output_file.name[-4:] == ".vcf"
    if output_file and not is_vcf:
        print("#" + "\t".join(OUTPUT_COLS), file=output_file)
    for i, minor_sol in enumerate(minor_sols):
        log.info(f"  {i + 1:2}: {minor_sol.diplotype}")
        tt = textwrap.wrap(minor_sol._solution_nice(), width=60, break_long_words=False)
        t = "\n             ".join(tt)
        log.info(f"      Minor: {t}")
        conf = (min_minor_score + SLACK) / (minor_sol.score + SLACK)
        log.info(f"      Confidence: {conf:.2f} (score = {minor_sol.score:.2f})")
        if output_file and not is_vcf:
            print(f"#Solution {i + 1}: {minor_sol._solution_nice()}", file=output_file)
            diplotype.write_decomposition(
                sample_name, gene, i + 1, minor_sol, output_file
            )
    if is_vcf:
        diplotype.write_vcf(sample_name, gene, sample.coverage, minor_sols, output_file)
    json_print(debug, "},")

    return minor_sols
