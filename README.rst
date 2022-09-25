.. raw:: html

   <h1 align="center">
   <img src="https://user-images.githubusercontent.com/10132487/100571499-1ee1fd00-3288-11eb-9760-75c4b0b98d2a.png" alt="Aldy" width=100px/>
   </h1>
   <p align="center">
   <a href="https://badge.fury.io/py/aldy"><img src="https://badge.fury.io/py/aldy.svg" alt="Version"/></a>
   <img src="https://github.com/0xTCG/aldy/workflows/aldy-test/badge.svg" alt="CI Status"/>
   <a href="https://aldy.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/aldy/badge/?version=latest" alt="ReadTheDocs"/></a>
   <a href="https://codecov.io/github/0xTCG/aldy"><img src="https://codecov.io/github/0xTCG/aldy/coverage.svg?branch=master" alt="Code Coverage"/></a>
   <a href="https://github.com/psf/black"><img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black"/></a>
   <a href="https://lgtm.com/projects/g/0xTCG/aldy/context:python"><img src="https://img.shields.io/lgtm/grade/python/g/0xTCG/aldy.svg?logo=lgtm&logoWidth=18" alt="Language grade: Python"/></a>
   <a href="https://www.nature.com/articles/s41467-018-03273-1"><img src="https://img.shields.io/badge/Published%20in-Nature%20Communications-red.svg" alt="Published in Nature Communications" /></a>
  <br/>
  <b><i>A quick and nifty tool for genotyping and phasing popular pharmacogenes.</i></b>
  </p>


Aldy 4 calls genotypes of many highly polymorphic pharmacogenes and reports them in a phased star-allele nomenclature.
It can also call copy number of a given pharmacogene and genotype each copy present in the sample—something that standard
genotype callers like GATK cannot do.

Algorithm details
=================

TL;DR: Aldy 4 uses star-allele databases to guide the process of detecting the most likely genotype.
The optimization is done in three stages via integer linear programming.
See `Gene Support`_ for more details about the supported pharmacogene databases.

More details, together with the API documentation, are available
`at Read the Docs <https://aldy.readthedocs.io/en/latest/>`_.

Experimental data is available `here <paper>`_.

If you are using Aldy, please cite our papers in the
`Nature Communications <https://www.nature.com/articles/s41467-018-03273-1>`_
and `bioRxiv <https://www.biorxiv.org/content/10.1101/2022.08.11.503701v1>`_.

⚠️ Warning
==========

**Please read this carefully** if you are using Aldy in a clinical or commercial environment.

Aldy is a computational tool whose purpose is to *aid the genotype detection process*. It can be of tremendous help in that process. However, it is not perfect, and it can easily make a wrong call if the data is noisy, ambiguous or if the target sample contains a previously unknown allele.

☣️🚨 **Do not use the raw output of Aldy (or any other computational tool for that matter) to diagnose a disease or prescribe a drug!**
**You are responsibe for inspecting and validating the results (ideally) in a wet lab before doing something that can have major consequences.** 🚨☣️

We really mean it.

Finally, note that the allele databases are still a work in progress and that we still do not know the downstream impact of the vast majority of genotypes.

Installation
============

Aldy is written in Python and requires Python 3.7+ to run.
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

       If you have trouble installing `ortools` on a Nix-based Linux distro, try this::

           pip install --platform=manylinux1_x86_64 --only-binary=:all: --target ~/.local/lib/python3.8/site-packages ortools

* `Gurobi <http://www.gurobi.com>`_:
  a commercial solver which is free for academic purposes.
  Most thoroughly tested solver: if you encounter any issues with CBC, try Gurobi.
  After installing it, don't forget to install ``gurobipy`` package by going to
  Gurobi's installation directory
  (e.g., ``/opt/gurobi/linux64`` on Linux or ``/Library/gurobi751/mac64/`` on macOS)
  and typing::

      python3 setup.py install


Sanity check
============

After installing Aldy and a compatible ILP solver, please make sure to test
the installation by issuing the following command (this should take a few minutes)::

    aldy test

In case everything is set up properly, you should see something like this::

    🐿  Aldy v4.0 (Python 3.7.5 on macOS 12.4)
        (c) 2016-2022 Aldy Authors. All rights reserved.
        Free for non-commercial/academic use only.
    ================================ test session starts ================================
    platform darwin -- Python 3.7.5, pytest-5.3.1, py-1.8.0, pluggy-0.13.1
    rootdir: aldy, inifile: setup.cfg
    plugins: anyio-3.6.1, xdist-1.31.0, cov-2.10.1, forked-1.1.3
    collected 76 items
    aldy/tests/test_cn_real.py ........                                            [ 10%]
    aldy/tests/test_cn_synthetic.py .....                                          [ 17%]
    aldy/tests/test_diplotype_real.py ....                                         [ 22%]
    aldy/tests/test_diplotype_synthetic.py ......                                  [ 30%]
    aldy/tests/test_full.py ...........                                            [ 44%]
    aldy/tests/test_gene.py .......                                                [ 53%]
    aldy/tests/test_major_real.py ...........                                      [ 68%]
    aldy/tests/test_major_synthetic.py .......                                     [ 77%]
    aldy/tests/test_minor_real.py .......                                          [ 86%]
    aldy/tests/test_minor_synthetic.py ......                                      [ 94%]
    aldy/tests/test_query.py ....                                                  [100%]
    =========================== 76 passed in 131.10s (0:02:11) ==========================

Running
=======

Aldy needs a SAM, BAM, CRAM or VCF file for genotyping.
We will be using BAM as an example.

.. attention::
  It is assumed that reads are mapped to hg19 (GRCh37) or hg38 (GRCh38). Other reference genomes are not yet supported.

An index is needed for BAM files. Get one by running::

    samtools index file.bam

Aldy is invoked as::

    aldy genotype -p [profile] -g [gene] file.bam

Sequencing profile selection
----------------------------

The ``[profile]`` argument refers to the sequencing profile.
The following profiles are available:

- ``illumina`` or ``wgs`` for the Illumina WGS or exome (WXS) data (or any uniform-coverage technology).

   .. attention::

    It is highly recommended to use samples with at least 40x coverage.
    Anything below 20x might result in noisy copy number calls and missed variants.

- ``pgx1`` for the PGRNseq v.1 capture protocol data
- ``pgx2`` for the PGRNseq v.2 capture protocol data
- ``pgx3`` for the PGRNseq v.3 capture protocol data

- ``10x`` for 10X Genomics data

   .. attention::

    For the best results on the 10X Genomics datasets, use the `EMA aligner <https://github.com/arshajii/ema/>`_,
    especially if doing *CYP2D6* analysis. Aldy will also use the EMA read cloud information for
    improved variant phasing.

- ``exome``, ``wxs``, ``wes`` for the whole-exome sequencing data

   .. attention::

    ⚠️ **Be warned!:** whole-exome data is incomplete *by definition*, and Aldy will not be able to call major star-alleles
    defined by their intronic or upstream variants.
    Aldy also assumes that there are only two (2) gene copies if the `wxs` profile is used, as it cannot call copy number changes nor fusions from exome data.

- ``pacbio-hifi-targeted``, ``pacbio-hifi-targeted-twist`` for PacBio HiFi target capture data

   .. attention::

    The provided PacBio capture profiles are custom and are not standard.
    Please ensure to generate a custom profile if using different PacBio HiFi capture protocols.


If you are using a different technology (e.g., some home-brewed capture kit),
you can proceed provided that the following requirements are met:

- all samples have a similar coverage distribution
  (i.e., two sequenced samples with the same copy number configuration
  **must** have similar coverage profiles; please consult us if you are not sure about this)
- your panel includes a copy-number neutral region
  (currently, Aldy uses *CYP2D8* as a copy-number neutral region, but it can be overridden).

Having said that, you can use a sample BAM that is known to have two copies
of the genes you wish to genotype (without any fusions or copy number alterations)
as a profile as follows::

    aldy genotype -p profile-sample.bam -g [gene] file.bam -n [cn-neutral-region]

Alternatively, you can generate a profile for your panel/technology by running::

    # Get the profile
    aldy profile profile-sample.bam > my-cool-tech.profile
    # Run Aldy
    aldy genotype -p my-cool-tech.profile -g [gene] file.bam


Output
======

By default, Aldy will generate ``file-[gene].aldy``
(the default location can be changed via ``-o`` parameter).
Aldy also supports VCF file output: to enable it, just append `.vcf` to the output file name.
The summary of the calls is shown at the end of the output::

    $ aldy -p pgx2 -g cyp2d6 NA19788.bam
    🐿  Aldy v4.0 (Python 3.8.2 on Linux 3.10.0-1160.71.1.el7.x86_64-x86_64-with-glibc2.2.5)
        (c) 2016-2022 Aldy Authors. All rights reserved.
        Free for non-commercial/academic use only.
    Genotyping sample NA07048.cram...
    Potential CYP2D6 gene structures for NA07048:
      1: 2x*1 (confidence: 100%)
    Potential major CYP2D6 star-alleles for NA07048:
      1: 1x*1, 1x*4.021 (confidence: 100%)
      2: 1x*4, 1x*139 (confidence: 100%)
      3: 1x*4.021.ALDY_2, 1x*74 (confidence: 100%)
    Best CYP2D6 star-alleles for NA07048:
      1: *1 / *4.021 (confidence=100%)
          Minor alleles: *(1.016 +rs112568578 +rs113889384 +rs28371713 +rs28633410), *(4.021 +rs28371729 -rs28371702 -rs28588594)
    CYP2D6 results:
      - *1 / *4.021
        Minor: [*1.016 +rs112568578 +rs113889384 +rs28371713 +rs28633410] / [*4.021 +rs28371729 -rs28371702 -rs28588594]
        Legacy notation: [*1.016 +rs112568578 +rs113889384 +rs28371713 +rs28633410] / [*4.021 +rs28371729 -rs28371702 -rs28588594]

In this example, the *CYP2D6* genotype is \*1/\*4 in terms of major star-alleles.
The minor star-alleles are given after each major star-allele call (here, \*1.016 and \*4.021).
The minor alleles might also have additional or removed mutations.
The additions are marked with `+` in front (e.g., `+rs112568578`), while the losses carry `-` in front (e.g., `-rs28588594`).
In some instances, even the major alleles might contain additions (e.g., `(*1 +rs1234)`).
This indicates the presence of a novel star-allele that has not been cataloged yet.

By default, Aldy only reports solutions with the maximum confidence.
Use `--param gap=XY` (where `XY` is greater than 0) to report less likely solutions.

Explicit decomposition is given in the ``file-[gene].aldy``
(in the example above, it is ``NA19788_x.CYP2D6.aldy``).
An example of such a file is::

    #Sample Gene    SolutionID      Major   Minor   Copy    Allele  Location        Type    Coverage        Effect  dbSNP   Code    Status
    #Solution 1: *1.001, *4, *4.021
    NA10860 CYP2D6  1       *1/*4+*4.021    1.001;4;4.021   0       1.001
    NA10860 CYP2D6  1       *1/*4+*4.021    1.001;4;4.021   1       4       42522612        C>G     15      S486T   rs1135840
    ...[redacted]...
    #Solution 2: *4, *4, *139.001
    NA10860 CYP2D6  2       *4+*4/*139      4;139.001;4     0       4       42522612        C>G     15      S486T   rs1135840
    NA10860 CYP2D6  2       *4+*4/*139      4;139.001;4     0       4       42524946        C>T     32      splicing defect/169frameshift    rs3892097
    ...[redacted]...

The columns are:

- the sample name,
- the gene name,
- the solution count (different solutions have different counts),
- the major star-allele call,
- the minor star-allele call,
- the allele copy identifier (0 for the first allele in the minor column, 1 for the second and so on)
- the mutation location,
- the mutation type (SNP or indel),
- the mutation coverage,
- the mutation functionality:

  - ``DISRUPTING`` for gene-disrupting (functional) mutations, and
  - ``NEUTRAL`` for neutral (silent) mutations

- the dbSNP ID (if available),
- traditional Karolinska-style mutation code from the CYP allele database (if available); and
- the mutation status, which indicates the status of the mutation in the decomposition:

    + ``NORMAL``: mutation is associated with the star-allele in the database and is found in the sample
    + ``NOVEL``: gene-disrupting mutation is **NOT** associated with the star-allele in the database,
      but is found in the sample (this indicates that Aldy found a novel major star-allele)
    + ``EXTRA``: neutral mutation is **NOT** associated with the star-allele in the database,
      but is found in the sample (this indicates that Aldy found a novel minor star-allele)
    + ``MISSING``: neutral mutation is associated with the star-allele in the database,
      but is **NOT** found in the sample (this also indicates that Aldy found a novel minor star-allele)

VCF support
-----------

The output will be a VCF file if the output file extension is `.vcf`.
Aldy will report a VCF sample for each potential solution and the appropriate genotypes.
Aldy will also output tags `MA` and `MI` for major and minor solutions.

  **Note:** VCF is not an optimal format for star-allele reporting. Unless you really need it,
  we recommend using Aldy's default format.


Problems & Debugging
--------------------

If you encounter any issues with Aldy, please run Aldy with debug parameter:

   aldy genotype ... --debug debuginfo

This will produce `debuginfo.tar.gz` file that contains the sample and LP model dumps.
Please send us this file, and we will try to resolve the issue.

This file contains no private information of any kind except for the phasing information
and mutation counts at the target gene locus as well as the file name.


Sample datasets
===============

Sample datasets are also available for download. They include:

- `HG00463 <http://cb.csail.mit.edu/cb/aldy/data/HG00463.bam>`_ (PGRNseq v.2), containing *CYP2D6* configuration with multiple copies
- `NA19790 <http://cb.csail.mit.edu/cb/aldy/data/NA19790.bam>`_ (PGRNseq v.2), containing a fusion between *CYP2D6* and *CYP2D7* deletion (\*78 allele)
- `NA24027 <http://cb.csail.mit.edu/cb/aldy/data/NA24027.bam>`_ (PGRNseq v.1), containing novel *DPYD* allele and multiple copies of *CYP2D6*
- `NA10856 <http://cb.csail.mit.edu/cb/aldy/data/NA10856.bam>`_ (PGRNseq v.1), containing *CYP2D6* deletion (\*5 allele)
- `NA10860 <http://cb.csail.mit.edu/cb/aldy/data/NA10860.bam>`_ (Illumina WGS), containing three copies of *CYP2D6*. This sample contains only the *CYP2D6* region.

The expected results are:

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

© 2016-2022 Aldy Authors, Indiana University Bloomington. All rights reserved.

**Aldy is NOT free software.**
A complete legal license is available in :ref:`aldy_license`.

For non-legal folks, here is a TL;DR version:

- Aldy can be freely used in academic and non-commercial environments
- Please contact us if you intend to use Aldy for any commercial purpose


Parameters & Usage
==================

**NAME**:
---------

Aldy --- a tool for allelic decomposition (haplotype reconstruction) and exact genotyping
         of highly polymorphic and structurally variant genes.

**SYNOPSIS**:
-------------

    aldy [--verbosity VERBOSITY] [--log LOG] command

Commands::

    aldy help
    aldy test
    aldy license
    aldy query (q)
    aldy profile [FILE]
    aldy genotype [-h] [--verbosity VERBOSITY] [--gene GENE] [--profile PROFILE]
                  [--reference REFERENCE] [--genome GENOME] [--cn-neutral-region CN_NEUTRAL_REGION]
                  [--output OUTPUT] [--solver SOLVER] [--debug DEBUG] [--cn CN] [--log LOG]
                  [--multiple-warn-level MULTIPLE_WARN_LEVEL] [--simple]
                  [--param PARAM=VALUE [PARAM2=VALUE2 ...]]
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

* ``query``, ``q``

  Query a gene or an allele.

  You can specify a gene name (e.g. ``aldy query CYP2D6``) or an allele (e.g. ``aldy query 'CYP2D6*121'`` or ``aldy q 'CYP2D6*4C'``).

* ``profile [FILE]``

  Generate a copy-number profile for a custom sequencing panel and
  print it on the standard output.
  ``FILE`` is a SAM/BAM sample that is known to have two copies of the gene of interest
  (without any fusions or copy number alterations).

* ``genotype``

  Genotype a sequencing sample. Arguments:

  - ``FILE``

    A SAM, BAM, CRAM or VCF file. A CRAM file requires ``--reference`` as well.

  - ``-p, --profile PROFILE``

    Sequencing profile. Supported values are:

    + ``illumina`` (or ``wgs``)
    + ``exome`` (or ``wxs`` or ``wes``)
    + ``pgx1`` (or ``pgrnseq-v1``)
    + ``pgx2`` (or ``pgrnseq-v2``)
    + ``pgx3`` (or ``pgrnseq-v3``)
    + ``10x``
    + ``pacbio-hifi-targeted``
    + ``pacbio-hifi-targeted-twist``

    You can also pass a SAM/BAM file as a profile(please check the documentation quick-start for more details).
    Also consult ``profile`` command.

  - ``-g, --gene GENE``

    Gene profile.

    *Default:* ``CYP2D6``

  - ``-o, --output OUTPUT``

    Location of the output file.

    *Default:* ``[input].[gene].aldy``

  - ``-s, --solver SOLVER``

    ILP Solver. Currently supported solvers are Gurobi and CBC.
    You can also pass ``any`` to let Aldy choose the best (available) solver.

    *Default:* ``any`` (uses CBC if available, then Gurobi).

  - ``-c, --cn CN``

    Manually specify a copy number configuration.
    Input: a comma-separated list of configurations ``CN1,CN2,...``.
    For a list of supported configurations, please run::

        aldy query [GENE]

  - ``-r, --reference REF``

    FASTA reference for the reference-encoded CRAM files.

  - ``-n, --cn-neutral-region CN_NEUTRAL``

    Provide a custom copy-number neutral region.
    Format is ``chr:start-end``.

    *Default:* *CYP2D8* (22:42547463-42548249 for hg19)

  - ``-d, --debug DEBUG``

    Create a `DEBUG.tar.gz`` file that can be shared with the authors for easier debugging.
    Contains no private information except the file name and sample mutation counts in
    the gene of interest.

  - ``--multiple-warn-level MULTIPLE_WARN_LEVEL``

    Warning level when multiple optimal solutions are found.

    If set to 1, Aldy will warn if multiple final optimal solutions are found.
    If set to 2, Aldy will also warn if multiple optimal major star-allele solutions are found.
    If set to 3, Aldy will even warn if multiple copy-number configurations are found.

    *Default:* 1

  - ``--param PARAM1=VAL1 [PARAM2=VAL2 ...]``

    Additional model parameters. Please check
    `the parameter documentation <https://aldy.readthedocs.io/en/latest/source/aldy.html#aldy.profile.Profile>`_
    for the list of the available parameters.

Gene Support
============

.. list-table::
   :header-rows: 1

   * - Gene
     - Version
     - Status
     - Notes
   * - *CYP2D6*
     - PharmVar 5.2.3
     - ✅
     - - Copy number and structural variation supported
       - Alleles with the *CYP2D7* exon 9 retention such as \*36, \*57, \*83 and \*141
         can be accurately called only when the copy number detection is enabled
         (i.e., they cannot be called in WES mode)
       - Detection of the non-functional *CYP2D7* intron 1 retention is spotty
   * - *CYP2A6*
     - PharmVar 5.2.3
     - ✅
     - - Copy number and structural variation supported
       - Detection of the *CYP2A7* 3' UTR retention not yet supported
   * - *CYP2B6*
     - PharmVar 5.2.3
     - ✅
     - Some allele calls should be further validated (e.g., \*6/\*9)
   * - *CYP1A1*
     - PharmGKB (Dec 2014) and Pharmacoscan R9
     - ✅
     -
   * - *CYP1A2*
     - PharmGKB (Mar 2014) and Pharmacoscan R9
     - ✅
     -
   * - *CYP2A13*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP2C19*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP2C8*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP2C9*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP2E1*
     - PharmGKB (Nov 2013)
     - ⚠️
     - Thorough testing on the real datasets pending
   * - *CYP2F1*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP2J2*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP2R1*
     - PharmVar 5.2.3
     - ⚠️
     - Thorough testing on the real datasets pending
   * - *CYP2S1*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP2W1*
     - PharmVar 5.2.3
     - ⚠️
     - Thorough testing on the real datasets pending
   * - *CYP3A43*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP3A4*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP3A5*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP3A7*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CYP4F2*
     - PharmVar 5.2.3
     - ✅
     -
   * - *CFTR*
     - PharmGKB (Jun 2020) and Pharmacoscan R9
     - ✅
     -
   * - *COMT*
     - Pharmacoscan R9
     - ✅
     -
   * - *DPYD*
     - PharmVar 5.2.3
     - ✅
     -
   * - *G6PD*
     - PharmGKB and Pharmacoscan R9 (Sep 2018)
     - ⚠️
     - - Thorough testing on the real datasets pending
       - Null allele calling is unstable
   * - *GSTM1*
     - Pharmacoscan R9
     - ✅
     -
   * - *GSTP1*
     - Pharmacoscan R9
     - ✅
     -
   * - *IFNL3*
     - PharmGKB and Pharmacoscan R9
     - ✅
     -
   * - *NAT1*
     - PharmGKB (Mar 2014) and Pharmacoscan R9
     - ✅
     -
   * - *NAT2*
     - PharmGKB (Mar 2014) and Pharmacoscan R9
     - ✅
     -
   * - *NUDT15*
     - PharmVar 5.2.3
     - ✅
     -
   * - *SLCO1B1*
     - PharmVar 5.2.3
     - ✅
     -
   * - *TPMT*
     - PharmGKB (Jun 2020) and Pharmacoscan R9
     - ✅
     -
   * - *UGT1A1*
     - PharmGKB (Feb 2020) and Pharmacoscan R9
     - ⚠️
     - Thorough testing on the real datasets pending
   * - *UGT2B7*
     - pharmacogenomics.pha.ulaval.ca (Apr 2015) / Pharmacoscan R9
     - ⚠️
     - Thorough testing on the real datasets pending
   * - *VKORC1*
     - PharmGKB (Jan 2021) and Pharmacoscan R9
     - ⚠️
     - Thorough testing on the real datasets pending

Change log
==========

- Aldy v4.2 (Sep 25th, 2022)
   - Fix indelpost setup errors
   - Various small fixes

- Aldy v4.1 (Aug 28th, 2022)
   - Output allele's activity and/or impact when available
   - Updated and tested gene definitions
     - Major changes to *NAT1*, *NAT2*, *UGT1A1*, *CYP2E1* and *CYP2A6*

   - Indel realignment support via `indelpost <https://github.com/stjude/indelPost>`_
   - New debug format
   - Various small fixes

- Aldy v4.0 (Aug 17th, 2022)
   - Major model changes
   - Phasing support
   - Long-read sequencing support (PacBio HiFi, 10X Genomics)
   - Support for new pharmacogenes
   - New allele databases
   - New profile format (**⚠️ WARNING:** Please make sure to re-generate custom profiles if using older Aldy profiles.)
   - Major API changes
   - New debug format
   - Various small fixes

- Aldy v3.0 (Nov 30th, 2020)
   - Support for hg38
   - Support for 15+ new pharmacogenes
   - New profile format (**⚠️ WARNING:** Please make sure to re-generate custom profiles if using Aldy v2 profiles.)
   - Better genotype calling models
   - Major API changes

Acknowledgments
===============

The following people made Aldy much better software:

- Ananth Hari
- Qinghui Zhou
- Michael Ford `@michael-ford <https://github.com/michael-ford>`_
- Farid Rashidi `@faridrashidi <https://github.com/faridrashidi>`_
- David Twesigomwe `@twesigomwedavid <https://github.com/twesigomwedavid>`_
- Tyler Shrug `@tshugg <https://github.com/tshugg>`_
- Reynold C. Ly
- Pieter W. Smit
- Lawrence Hon `@lhon <https://github.com/lhon>`_
- Zach Langley `@zlangley <https://github.com/zlangley>`_


Contact & Bug Reports
=====================

`Ibrahim Numanagić <mailto:inumanag.at.uvic.ca>`_

or open a `GitHub issue <https://github.com/inumanag/aldy/issues>`_.

If you have an urgent problem, I suggest using e-mail.
