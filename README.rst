Aldy
****

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

* `Gurobi <http://www.gurobi.com>`_ (**recommended**):
  a commercial solver which is free for academic purposes. 
  Most thoroughly tested solver.
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
  `SCIPOPTDIR` enviromental variable to point to SCIP's install directory.


* `CBC / Google OR-Tools <https://developers.google.com/optimization/>`_: 
  a free, open-source MIP solver that is shipped by default with Google's OR-Tools.
  Install OR-Tools via ``pip install ortools`` to use this solver.
  

Sanity check
============

After installing Aldy and a compatible ILP solver, please make sure to test 
the installation by issuing the following command (this should take around a minute)::

    aldy test

In case everything is set up properly, you should see something like this::

    *** Aldy v1.9 (Python 3.6.6) ***
    (c) 2016-2018 Aldy Authors & Indiana University Bloomington. All rights reserved.
    Aldy Sanity-Check Test
    Expected result is: *1/*4+*4
    Result:
      *1/*4+*4                       (1, 4N, 4AW)

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
  but it can be overriden)

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
    *** Aldy v2.0 ***
    [...]
    Result:
      *2/*78+*2                      (2MW, 2MW, 78/2|2M)

In this example, *CYP2D6* genotype is \*2/\*78+\*2 as expressed in terms of major star-alleles. 
Minor star-alleles are given in the parenthesis 
(in this case, two copies of \*2MW, and one copy of \*78 fusion on the \*2M background).

Explicit decomposition is given in the ``file-[gene].aldy`` (in the example above, it is ``NA19788_x.CYP2D6.aldy``).  
An example of such file is::

    # Aldy v1.0
    # Gene: CYP2D6
    # Number of solutions: 1

    # Solution 0
    # Predicted diplotype: *2/*78+*2
    # Composition: 2MW,2MW,78/2|2M
    Copy   Allele   Location   Type     Coverage  Effect      dbSNP       Code        Status
    0      78/2     42522311   SNP.CT   1760      NEUTRAL     rs12169962  4481:G>A    NORMAL
    0      78/2     42522612   SNP.CG   1287      DISRUPTING  rs1135840   4180:G>C    NORMAL
    ...[redacted]...
    1      2MW      42522311   SNP.CT   1760      NEUTRAL     rs12169962  4481:G>A    NORMAL
    1      2MW      42527541   DEL.TC   0         NEUTRAL     rs536645539 -750:delGA  MISSING
    ...[redacted]...


Each solution is indicated with the **"Solution"** line. 
The first column (copy) shows the ordinary number of the allelic copy (e.g. 0, 1 and 2 for 2MW, 2MW and 78/2M, respectively). 
The following columns indicate:

- star-allele, 
- mutation loci,
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


Logging
-------

Detailed execution log will be located in ``file-[gene].aldylog``. It is used mainly for debugging purposes.
In case you have issues with Aldy, please provide this file as it will greatly help us during the debugging process.


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

© 2016-2018 Aldy Authors, Indiana University Bloomington. All rights reserved.

**Aldy is NOT a free software.** Complete legal license is available in :ref:`aldy_license`. 

For non-legal folks, here is a TL;DR version:

- Aldy can be freely used in academic and non-commercial environments
- Please contact us if you intend to use Aldy for any commercial purpose


Parameters & Usage
==================

**NAME**:
---------

Aldy --- Tool for allelic decomposition and exact genotyping of highly polymorphic and structurally variant genes.

**SYNOPSIS**:
-------------

    aldy [--verbosity VERBOSITY] [--log LOG] command

Commands::

    aldy help
    aldy test
    aldy license
    aldy show [-g/--gene GENE]
    aldy profile [FILE]
    aldy genotype [-T/--threshold THRESHOLD] 
                  [-p/--profile PROFILE] 
                  [-g/--gene GENE] 
                  [-o/--output OUTPUT] 
                  [-n/--cn-neutral-region CN_NEUTRAL]
                  [--solver SOLVER]
                  [-r/--reference REF]
                  [-c/--cn CN] 
                  [FILE]

**OPTIONS**:
------------

Global arguments:
^^^^^^^^^^^^^^^^^

* ``-h, --help`` 

  Show the help message and exit.  

* ``-v, --verbosity VERBOSITY``  

  Logging verbosity. Acceptable values are:

  - ``T`` (trace)
  - ``D`` (debug), 
  - ``I`` (info) and 
  - ``W`` (warn)
    
  *Default:* ``I``

* ``-l, --log LOG``  

  Location of the output log file .  
  
  *Default:* ``[FILE].[GENE].aldylog``


Commands:
^^^^^^^^^

* ``help``
  
  Show the help message and exit.

* ``license`` 

  Print Aldy license.  

* ``test``  

  Sanity-check on NA10860 sample.

* ``show``  

  Show all copy number configurations supported by a gene (requires ``--gene``).

* ``profile [FILE]``

  Generate a copy-number profile for a custom sequencing panel and 
  print it on the standard output.
  ``FILE`` is a SAM/BAM of a sample that is known to have two copies of a target genes 
  (without any fusions or copy number alterations).

* ``genotype``  

  Genotype SAM/BAM sample. Arguments:

  - ``FILE``

    SAM, BAM, CRAM or DeeZ input file. CRAM and DeeZ require ``--reference`` as well.

  - ``-T, --threshold THRESHOLD``
  
    Cut-off rate for variations (percent per copy)  
    
    *Default:* `50`

  - ``-p, --profile PROFILE``
  
    Sequencing profile. Supported values are:

    + ``illumina``
    + ``pgrnseq-v1``
    + ``pgrnseq-v2``. 

    You can also pass a SAM/BAM file 
    (please check documentation quick-start for more information).
    Also check ``profile`` command.

  - ``-g, --gene GENE``
  
    Gene profile.  

    *Default:* ``CYP2D6``

  - ``-o, --output OUTPUT``
   
    Location of the output file.   

    *Default:* ``[input].[gene].aldy``

  - ``-s, --solver SOLVER``
  
    ILP Solver. Currently supported solvers are Gurobi and SCIP.    
    
    *Default:* ``any``

  - ``-c, --cn CN``
   
    Manually set copy number configuration.
    Input: a comma-separated list ``CN1,CN2,...``. 
    For a list of supported configurations, please run::

        aldy show --gene [GENE]

  - ``-r, --reference REF``
   
    Specify FASTA reference for reference-encoded CRAM/DeeZ files.

  - ``-n, --cn-neutral-region CN_NEUTRAL``
   
    Provide a custom copy-number neutral region.
    Format is ``chr:start-end``.

    *Default:* *CYP2D8* (22:42547463-42548249 in hg19)


Contact & Bug Reports
=====================

`Ibrahim Numanagić <mailto:inumanag.at.mit.dot.edu>`_

If you have an urgent question, I suggest using e-mail. 
GitHub issues are not handled as fast as email requests are.
