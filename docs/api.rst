API Tutorial
************

Aldy can be used as a standalone application and as a Python module.

Here is a quick overview of the functions that Aldy provides.
Detailed API documentation is available in :ref:`aldy_api`.


One-step genotyping API
=======================

You can get the genotypes of a SAM/BAM with a single call::

  import aldy.genotype

  result = aldy.genotype.genotype('cyp2d6', '/path/to/sample.bam', profile='pgrnseq-v1')
  print(result)

Please check :py:class:`aldy.major.SolvedAllele` and :py:class:`aldy.genotype` for more details.


Loading a gene
==============

We always start by loading a gene database for a gene that is to be genotyped::

  import aldy.gene
  gene = aldy.gene.Gene('path/to/gene.yaml')


Alternatively we can use genes shipped with Aldy::

  import aldy.gene
  import aldy.common

  gene_path = aldy.common.script_path('aldy.resources.genes/cyp2d6.yml')
  gene = aldy.gene.Gene(gene_path) # for CYP2D6


Loading a SAM/BAM file
======================

A SAM/BAM/CRAM file can be loaded as follows::

  import aldy.sam
  sample = aldy.sam.Sample(sam_path='my/sample.bam',
                           gene=gene,
                           threshold=0.5,
                           profile='illumina',
                           reference=None,
                           cn_region=aldy.sam.DEFAULT_CN_NEUTRAL_REGION)

Pass a path to the reference genome (``reference`` parameter) if using a CRAM file.
Custom copy-number region can be passed via ``cn_region`` parameter (by default it points to *CYP2D8* region).
Threshold should be always set to 0.5.


Detecting copy number configurations and fusions
================================================

To detect copy number configurations, run::

  import aldy.cn
  cn_sols = aldy.cn.estimate_cn(gene, sample.coverage, solver='gurobi')

You can also use ``solver='cbc'``, or alternatively ``solver='scip'`` if you have
`PySCIPOpt <https://github.com/SCIP-Interfaces/PySCIPOpt>`_ installed.

Result is a list of :py:class:`aldy.solutions.CNSolution` objects described in API.


Calling major and minor star-alleles
====================================

Once you get copy number solutions, you can call valid major star-alleles for each copy number solution::

  import aldy.major

  major_sols = [sol
                for cn_sol in cn_sols
                for sol in aldy.major.estimate_major(gene, sample.coverage, cn_sol, 'gurobi')]
  # Get the best major star-allele score
  min_score = min(major_sols, key=lambda m: m.score).score
  # Take the best major star-allele calls
  major_sols = sorted([m for m in major_sols if abs(m.score - min_score) < 1e-3],
                     key=lambda m: m.score)

You are pretty much set if you need only major star-alleles.
However, if you have multiple equally likely major star-allele calls, or if you need
to get the whole decomposition, then do the following::

  import aldy.minor

  minor_sols = aldy.minor.estimate_minor(gene, sample.coverage, major_sols, 'gurobi')
  # Get the best minor star-allele score
  min_score = min(minor_sols, key=lambda m: m.score).score
  # Get the best minor (and major) star-allele
  minor_sols = [m for m in minor_sols if abs(m.score - min_score) < 1e-3]


``minor_sols`` will contain a list of :py:class:`aldy.solutions.MinorSolution` objects that point to the optimal major star-alleles and copy numbers.

Finally, if you want to get a nice diplotype (e.g. \*1/\*2+\*3), just type::

  minor_solution.diplotype

More detailed explanation of these functions is available in the :ref:`aldy_api`.
