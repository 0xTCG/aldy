Aldy
****

.. image:: https://badge.fury.io/py/aldy.svg
    :target: https://badge.fury.io/py/aldy
.. image:: https://travis-ci.com/inumanag/aldy.svg?branch=master
    :target: https://travis-ci.com/inumanag/aldy
.. image:: https://readthedocs.org/projects/aldy/badge/?version=latest
    :target: https://aldy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. image:: https://codecov.io/github/inumanag/aldy/coverage.svg?branch=master
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black


Aldy is a tool for allelic decomposition (haplotype reconstruction) and exact genotyping
of highly polymorphic and structurally variant genes.
More simply, it is a tool which can detect the copy number of a target gene,
and the structure and genotype of each gene copy present in the sample.

Aldy has been published in `Nature Communications <https://www.nature.com/articles/s41467-018-03273-1>`_
(`doi:10.1038/s41467-018-03273-1 <http://doi.org/10.1038/s41467-018-03273-1>`_).
Preprint `is available here <https://github.com/inumanag/aldy/blob/master/docs/preprint.pdf>`_.
Full experimental pipeline `is available here <https://github.com/inumanag/aldy-paper-resources>`_.

Documentation is available `at Read the Docs <https://aldy.readthedocs.io/en/latest/>`_.


Installation
============

Aldy is written in Python, and requires Python 3.6+.
It is intended to be run on POSIX-based systems
(so far, only Linux and macOS have been tested).

The easiest way to install Aldy is to use `pip`::

    pip install aldy

Append ``--user`` to the previous command to install Aldy locally
if you cannot write to the system-wide Python directory.


Prerequisite: ILP solver
------------------------

Aldy requires a mixed integer solver to run.

The following solvers are currently supported:

* `CBC / Google OR-Tools <https://developers.google.com/optimization/>`_:
  a free, open-source MIP solver that is shipped by default with Google's OR-Tools.
  `pip` installs it by default when installing Aldy.

* `Gurobi <http://www.gurobi.com>`_:
  a commercial solver which is free for academic purposes.
  Most thoroughly tested solver: if you encounter any issues with CBC, try Gurobi.
  After installing it, don't forget to install ``gurobipy`` package by going to
  Gurobi's installation directory
  (e.g. ``/opt/gurobi/linux64`` on Linux or ``/Library/gurobi751/mac64/`` on macOS)
  and typing::

      python3 setup.py install

* `SCIP <http://scip.zib.de>`_: another solver which is also free for academic purposes.
  SCIP is easier to install than Gurobi (no registration or activation required).
  However, it might be slower than Gurobi.
  Once you you install SCIP, please install
  `PySCIPPpt <https://github.com/SCIP-Interfaces/PySCIPOpt>`_ module for the Python
  SCIP bindings via `pip`: ``pip install pyscipopt``. If it fails, make sure to set
  `SCIPOPTDIR` environmental variable to point to SCIP's install directory.



Sanity check
============

After installing Aldy and a compatible ILP solver, please make sure to test
the installation by issuing the following command (this should take a few minutes)::

    aldy test

In case everything is set up properly, you should see something like this::

    *** Aldy v2.0 (Python 3.6.6, darwin) ***
    *** (c) 2016-2019 Aldy Authors & Indiana University Bloomington. All rights reserved.
    *** Free for non-commercial/academic use only.
    ========================================== test session starts ==========================================
    platform darwin -- Python 3.6.6, pytest-5.1.2, py-1.8.0, pluggy-0.12.0
    plugins: xdist-1.29.0, forked-1.0.2, cov-2.7.1
    collected 62 items

    tests/test_cn_real.py ........                                                                    [ 12%]
    tests/test_cn_synthetic.py .....                                                                  [ 20%]
    tests/test_diplotype_real.py ....                                                                 [ 27%]
    tests/test_diplotype_synthetic.py ......                                                          [ 37%]
    tests/test_full.py .....                                                                          [ 45%]
    tests/test_gene.py ....                                                                           [ 51%]
    tests/test_major_real.py ...........                                                              [ 69%]
    tests/test_major_synthetic.py .......                                                             [ 80%]
    tests/test_minor_real.py ......                                                                   [ 90%]
    tests/test_minor_synthetic.py .....                                                               [ 98%]
    tests/test_paper.py s                                                                             [100%]

    =============================== 61 passed, 1 skipped in 106.83s (0:01:46) ===============================


Running
=======

Aldy needs a SAM, BAM, CRAM or a DeeZ file for genotyping.
We will be using BAM as an example.

.. attention::
  It is assumed that reads are mapped to hg19 or GRCh37. hg38 is not yet supported.

An index is needed for BAM files. Get one by running::

    samtools index file.bam

Aldy is invoked as::

    aldy genotype -p [profile] -g [gene] file.bam

The ``[gene]`` parameter indicates the name of the gene to be genotyped.
Currently, Aldy supports:

- *CYP2D6*
- *CYP2A6*
- *CYP2C19*
- *CYP2C8*
- *CYP2C9*
- *CYP3A4*
- *CYP3A5*
- *CYP4F2*
- *TPMT* and
- *DPYD*.


Sequencing profile selection
----------------------------

The ``[profile]`` argument refers to the sequencing profile.
The following profiles are available:

- ``illumina`` for Illumina WGS (or any uniform-coverage technology).

.. attention::
  It is highly recommended to use samples with at least 40x coverage.
  Anything lower than 20x will result in tears and agony.

- ``pgrnseq-v1`` for PGRNseq v.1 capture protocol data
- ``pgrnseq-v2`` for PGRNseq v.2 capture protocol data

If you are using different technology (e.g. some home-brewed capture kit),
you can proceed provided that the following requirements are met:

- all samples have the similar coverage distribution
  (i.e. two sequenced samples with the same copy number configuration
  MUST have similar coverage profiles; please consult us if you are not sure about this)
- your panel includes a copy-number neutral region
  (currently, Aldy uses *CYP2D8* as a copy-number neutral region,
  but it can be overridden)

Having said that, you can use a sample BAM that is known to have two copies
of the genes you wish to genotype (without any fusions or copy number alterations)
as a profile as follows::

    aldy genotype -p profile-sample.bam -g [gene] file.bam

Alternatively, you can generate a profile for your panel/technology by running::

    # Get the profile
    aldy profile profile-sample.bam > my-cool-tech.profile
    # Run Aldy
    aldy genotype -p my-cool-tech.profile -g [gene] file.bam


Output
======

Aldy will by default generate the following file: ``file-[gene].aldy``
(default location can be changed via ``-o`` parameter),
The summary of results are shown at the end of the output::

    $ aldy -p pgrnseq-v2 -g cyp2d6 NA19788_x.bam
    *** Aldy v2.0 (Python 3.7.4) ***
    *** (c) 2016-2019 Aldy Authors & Indiana University Bloomington. All rights reserved.
    *** Free for non-commercial/academic use only.
    Genotyping sample NA07048.cram...
    Potential CYP2D6 copy number configurations for NA07048:
      1: 2x*1
          Confidence: 1.00 (score = 3.22)

    Potential major CYP2D6 star-alleles for NA07048:
      1: 1x*1 +42525810:SNP.TC*, 1x*4.b
          Confidence: 1.00 (score = 22.47)
      2: 1x*1, 1x*4.b +42525810:SNP.TC*
          Confidence: 1.00 (score = 22.47)

    Best CYP2D6 star-alleles for NA07048:
      1: *1-like/*4
          Minor: *1 +42525810:SNP.TC*, *4EW
          Confidence: 1.00 (score = 25.73)
      2: *1/*4-like
          Minor: *1, *4EW +42525810:SNP.TC*
          Confidence: 1.00 (score = 25.73)
    CYP2D6 results:
      *1-like/*4                     (*1 +42525810:SNP.TC*, *4.b)
      *1/*4-like                     (*1, *4.b +42525810:SNP.TC*)

In this example, *CYP2D6* genotype is \*1/\*4 as expressed in terms of
major star-alleles.
Minor star-alleles are given after each "best" star-allele (here, \*1 and \*4EW).
Note that there is a novel SNP here (42525810:SNP.TC) that Aldy assigned to \*1
(and \*4 in the second solution). The presence of a novel functional SNP causes Aldy to
report modified allele with the suffix `-like` (e.g. `*1-like`).
Minor alleles might have additional mutations, or might lose some default mutations.
Additions are marked with `+` in front (e.g. `*1 +42525810:SNP.TC*`).
Losses carry `-` in front.

Confidence scores express Aldy's confidence in a solution.
Maximum score is 1.0. By default, Aldy only reports solutions that have the 
confidence score of 1.0. Use `--gap` to report more solutions.

Explicit decomposition is given in the ``file-[gene].aldy``
(in the example above, it is ``NA19788_x.CYP2D6.aldy``).
An example of such file is::

    #Sample	Gene	SolutionID	Major	Minor	Copy	Allele	Location	Type	Coverage	Effect	dbSNP	Code	Status
    #Solution 1: *1 +42528223:SNP.GA, *4AW, *4N -42522391:SNP.GA
    NA10860	CYP2D6	1	*1/*4+*4	1;4AW;4N	0	1	42528223	SNP.GA	-1	NEUTRAL	rs28588594	-1426:C>T
    NA10860	CYP2D6	1	*1/*4+*4	1;4AW;4N	1	4AW	42522391	SNP.GA	-1	NEUTRAL	rs28371738	4401:C>T
    NA10860	CYP2D6	1	*1/*4+*4	1;4AW;4N	1	4AW	42522612	SNP.CG	-1	DISRUPTING	rs1135840	4180:G>C    ...[redacted]...
    ...[redacted]...
    #Solution 2: *1, *4AW +42528223:SNP.GA, *4N -42522391:SNP.GA
    NA10860	CYP2D6	2	*1/*4+*4	1;4AW;4N	0	1
    NA10860	CYP2D6	2	*1/*4+*4	1;4AW;4N	1	4AW	42522391	SNP.GA	-1	NEUTRAL	rs28371738	4401:C>T
    ...[redacted]...

The columns stand for:
- sample name,
- gene name,
- solution count (different solutions have different counts),
- major star-allele call,
- minor star-allele call,
- allele copy identifier (0 for the first allele in the minor column, 1 for the second and so on)
- mutation locus,
- mutation type (SNP or indel),
- mutation coverage,
- mutation functionality:
  - ``DISRUPTING`` for gene-disrupting
  - ``NEUTRAL`` for neutral mutation,
- dbSNP ID (if available),
- traditional Karolinska-style mutation code from CYP allele database, and
- mutation status, which indicates the status of the mutation in the decomposition:

    + ``NORMAL``: mutation is associated with the star-allele in the database, and is found in the sample
    + ``NOVEL``: gene-disrupting mutation is **NOT** associated with the star-allele in the database,
      but is found in the sample (this indicates that Aldy found a novel major star-allele)
    + ``EXTRA``: neutral mutation is **NOT** associated with the star-allele in the database,
      but is found in the sample (this indicates that Aldy found a novel minor star-allele)
    + ``MISSING``: neutral mutation is associated with the star-allele in the database,
      but is **NOT** found in the sample (this also indicates that Aldy found a novel minor star-allele)


Problems & Debugging
--------------------

If you encounter any issues with Aldy, please run Aldy with debug parameter:

   aldy genotype ... --debug debuginfo

This will produce `debuginfo.tar.gz` file that contains sample and LP model dumps.
Please send us this file and we will try to resolve the issue.

This file contains no private information of any kind except for the mutation counts
at the target gene locus and the file name.


Sample datasets
===============

Sample datasets are also available for download. They include:

- `HG00463 <http://cb.csail.mit.edu/cb/aldy/data/HG00463.bam>`_ (PGRNseq v.2), containing *CYP2D6* configuration with multiple copies
- `NA19790 <http://cb.csail.mit.edu/cb/aldy/data/NA19790.bam>`_ (PGRNseq v.2), containing a fusion between *CYP2D6* and *CYP2D7* deletion (\*78 allele)
- `NA24027 <http://cb.csail.mit.edu/cb/aldy/data/NA24027.bam>`_ (PGRNseq v.1), containing novel *DPYD* allele and multiple copies of *CYP2D6*
- `NA10856 <http://cb.csail.mit.edu/cb/aldy/data/NA10856.bam>`_ (PGRNseq v.1), containing *CYP2D6* deletion (\*5 allele)
- `NA10860 <http://cb.csail.mit.edu/cb/aldy/data/NA10860.bam>`_ (Illumina WGS), containing 3 copies of *CYP2D6*. This sample contains only *CYP2D6* region.

Expected results are:

============= ===================== ================ ================= ============ ==============
Gene (`-g`)   HG00463               NA19790          NA24027           NA10856      NA10860
============= ===================== ================ ================= ============ ==============
*CYP2D6*      \*36+\*10/\*36+\*10   \*1/\*78+\*2     \*6/\*2+\*2       \*1/\*5      \*1/\*4+\*4
*CYP2A6*      \*1/\*1               \*1/\*1          \*1/\*35          \*1/\*1
*CYP2C19*     \*1/\*3               \*1/\*1          \*1/\*2           \*1/\*2
*CYP2C8*      \*1/\*1               \*1/\*3          \*1/\*3           \*1/\*1
*CYP2C9*      \*1/\*1               \*1/\*2          \*1/\*2           \*1/\*2
*CYP3A4*      \*1/\*1               \*1/\*1          \*1/\*1           \*1/\*1
*CYP3A5*      \*3/\*3               \*3/\*3          \*1/\*3           \*1/\*3
*CYP4F2*      \*1/\*1               \*3/\*4          \*1/\*1           \*1/\*1
*TPMT*        \*1/\*1               \*1/\*1          \*1/\*1           \*1/\*1
*DPYD*        \*1/\*1               \*1/\*1          \*4/\*5           \*5/\*6
============= ===================== ================ ================= ============ ==============


License
=======

© 2016-2019 Aldy Authors, Indiana University Bloomington. All rights reserved.

**Aldy is NOT free software.**
Complete legal license is available in :ref:`aldy_license`.

For non-legal folks, here is a TL;DR version:

- Aldy can be freely used in academic and non-commercial environments
- Please contact us if you intend to use Aldy for any commercial purpose


Parameters & Usage
==================

**NAME**:
---------

Aldy --- tool for allelic decomposition (haplotype reconstruction) and exact genotyping
         of highly polymorphic and structurally variant genes.

**SYNOPSIS**:
-------------

    aldy [--verbosity VERBOSITY] [--log LOG] command

Commands::

    aldy help
    aldy test
    aldy license
    aldy show [-g/--gene GENE]
    aldy profile [FILE]
    aldy genotype [-h]
                  --profile PROFILE
                  [--verbosity VERBOSITY]
                  [--gene GENE]
                  [--threshold THRESHOLD]
                  [--reference REFERENCE]
                  [--cn-neutral-region CN_NEUTRAL_REGION]
                  [--output OUTPUT]
                  [--solver SOLVER]
                  [--gap GAP]
                  [--debug DEBUG]
                  [--log LOG]
                  [--fusion-penalty FUSION_PENALTY]
                  [--max-minor-solutions MAX_MINOR_SOLUTIONS]
                  [--cn CN]
                  [FILE]

**OPTIONS**:
------------

Global arguments:
^^^^^^^^^^^^^^^^^

* ``-h, --help``

  Show the help message and exit.

* ``-v, --verbosity VERBOSITY``

  Logging verbosity. Acceptable values:

  - ``T`` (trace)
  - ``D`` (debug),
  - ``I`` (info), and
  - ``W`` (warn)

  *Default:* ``I``

* ``-l, --log LOG``

  Location of the output log file.

  *Default:* no log file


Commands:
^^^^^^^^^

* ``help``

  Show the help message and exit.

* ``license``

  Print Aldy license.

* ``test``

  Run Aldy test suite.

* ``show``

  Show a gene description (requires ``--gene``).

  - ``-g, --gene GENE``

    Gene profile.


  Optional parameters:

  - ``-c, --cn-config [CN_CONFIG]``

    Describe the copy number configuration CN_CONFIG.

  - ``-m, --major [MAJOR]``

    Describe the major star-allele MAJOR.

  - ``-M, --minor [MINOR]``

    Describe the minor star-allele MINOR.


* ``profile [FILE]``

  Generate a copy-number profile for a custom sequencing panel and
  print it on the standard output.
  ``FILE`` is a SAM/BAM sample that is known to have two copies of the gene of interest
  (without any fusions or copy number alterations).

* ``genotype``

  Genotype a SAM/BAM sample. Arguments:

  - ``FILE``

    SAM, BAM, CRAM or DeeZ file. CRAM and DeeZ require ``--reference`` as well.

  - ``-T, --threshold THRESHOLD``

    Cut-off rate for variations (percent per copy).
    Any variation with normalized coverage less than the threshold will be ignored.

    *Default:* `50`

  - ``-p, --profile PROFILE``

    Sequencing profile. Supported values are:

    + ``illumina``
    + ``pgrnseq-v1``
    + ``pgrnseq-v2``.

    You can also pass a SAM/BAM file
    (please check the documentation quick-start for more details).
    Also consult ``profile`` command.

  - ``-g, --gene GENE``

    Gene profile.

    *Default:* ``CYP2D6``

  - ``-o, --output OUTPUT``

    Location of the output file.

    *Default:* ``[input].[gene].aldy``

  - ``-s, --solver SOLVER``

    ILP Solver. Currently supported solvers are Gurobi, SCIP and CBC.
    You can also pass ``any`` to let Aldy choose the best (available) solver.

    *Default:* ``any``

  - ``-c, --cn CN``

    Manually specify a copy number configuration.
    Input: a comma-separated list of configurations ``CN1,CN2,...``.
    For a list of supported configurations, please run::

        aldy show --gene [GENE]

  - ``-r, --reference REF``

    FASTA reference for reference-encoded CRAM/DeeZ files.

  - ``-n, --cn-neutral-region CN_NEUTRAL``

    Provide a custom copy-number neutral region.
    Format is ``chr:start-end``.

    *Default:* *CYP2D8* (22:42547463-42548249 in hg19)

  - ``-G, --gap GAP``

    Solution gap.
    By setting this to any positive value, Aldy will also report solutions whose score 
    is less than (1+GAP) times the optimal solution score.
    Useful for exploring the solution space.
    
    *Default:* 0 (only optimal solutions allowed)

  - ``-d, --debug DEBUG``

    Create a DEBUG.tar.gz file that can be shared with the authors for easier debugging.
    Contains no private information except the file name and sample mutation counts in 
    the gene of interest.
    
  - ``-f, --fusion-penalty FUSION_PENALTY``

    Penalize each fusion additional FUSION_PENALTY times.
    Larger values mean lower likelihood of seeing fusions.

    *Default:* 0.1


Acknowledgements
================

The following people made Aldy much better software:

- Michael Ford `@michael-ford <https://github.com/michael-ford>`_
- Farid Rashidi `@faridrashidi <https://github.com/faridrashidi>`_
- David Twesigomwe `@twesigomwedavid <https://github.com/twesigomwedavid>`_
- Lawrence Hon `@lhon <https://github.com/lhon>`_
- Zach Langley `@zlangley <https://github.com/zlangley>`_


Contact & Bug Reports
=====================

`Ibrahim Numanagić <mailto:inumanag.at.uvic.ca>`_

or open a `GitHub issue <https://github.com/inumanag/aldy/issues>`_.

If you have an urgent problem, I suggest using e-mail.
