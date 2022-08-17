API Tutorial
************

Aldy can be used as a standalone application or as a Python module.

Here is a quick overview of the functions that Aldy provides.
Detailed API documentation is available in :ref:`aldy_api`.


One-step genotyping API
=======================

You can get the genotypes of a SAM/BAM with a single call::

  import aldy.genotype

  result = aldy.genotype.genotype(
    "cyp2d6", "/path/to/sample.bam", profile_name="wgs", output_file=None
  )
  print(result)

Please check :py:class:`aldy.solution.SolvedAllele` and :py:mod:`aldy.genotype` for more details.


Loading a gene
==============

We always start by loading a gene database for a gene that is to be genotyped::

  import aldy.gene
  gene = aldy.gene.Gene("path/to/gene.yaml", genome="hg38")


Alternatively, we can use genes shipped with Aldy::

  import aldy.gene
  import aldy.common

  gene_path = aldy.common.script_path("aldy.resources.genes/cyp2d6.yml")
  gene = aldy.gene.Gene(gene_path)  # for CYP2D6


Loading a SAM/BAM or VCF file
=============================

Before working with alignment files, it is necessary to load a profile as follows::

  import aldy.profile
  profile = aldy.profile.Profile.load(gene, "illumina")

If using VCF files, construct the profile as follows::

  profile = Profile("user_provided", cn_solution=["1", "1"])

Various model parameters can be passed to the :py:meth:`aldy.profile.Profile.load` method.
One such parameter is a custom copy-number neutral region that can be passed via the ``cn_region`` parameter
(by default, the *CYP2D8* region is used).
Please consult :py:class:`aldy.profile.Profile` for details about the other parameters.

A SAM/BAM/CRAM or a VCF file can be loaded as follows::

  import aldy.sam
  sample = sample = aldy.sam.Sample(gene, profile, path="my/sample.bam")

Pass a path to the reference genome (``reference`` parameter) when using a CRAM file.

Detecting copy number configurations and fusions
================================================

To detect copy number configurations, run::

  import aldy.cn
  cn_sols = aldy.cn.estimate_cn(gene, profile, sample.coverage, solver="cbc")

You can also use ``solver="gurobi"`` if you have Gurobi installed.

The result is a list of :py:class:`aldy.solutions.CNSolution` objects described in API.


Calling the major and minor star-alleles
========================================

Once you get copy number solutions, you can call valid major star-alleles for each copy number solution::

  import aldy.major

  major_sols = [sol
                for cn_sol in cn_sols
                for sol in aldy.major.estimate_major(gene, sample.coverage, cn_sol, "cbc")]
  # Get the best major star-allele score
  min_score = min(major_sols, key=lambda m: m.score).score
  # Take the best major star-allele calls
  major_sols = sorted([m for m in major_sols if abs(m.score - min_score) < 1e-3],
                     key=lambda m: m.score)

You are pretty much set if you need only major star-alleles.
However, if you have multiple equally likely major star-allele calls, or if you need
to get the complete star-allele decomposition, then do the following::

  import aldy.minor

  minor_sols = aldy.minor.estimate_minor(gene, sample.coverage, major_sols, "cbc")
  # Get the best minor star-allele score
  min_score = min(minor_sols, key=lambda m: m.score).score
  # Get the best minor (and major) star-allele
  minor_sols = [m for m in minor_sols if abs(m.score - min_score) < 1e-3]

``minor_sols`` will contain a list of :py:class:`aldy.solutions.MinorSolution` objects
that point to the optimal major star-alleles and the corresponding copy numbers.

Finally, if you want to get a nice diplotype (e.g., \*1/\*2+\*3), just type::

  minor_solution.get_major_diplotype()
  minor_solution.get_minor_diplotype()  # if interested in minor star-alleles

A more detailed explanation of these functions is available in the :ref:`aldy_api`.
