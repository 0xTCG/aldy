Database Specification
**********************

Aldy gene descriptions are stored in the database filesd in `YAML <https://en.wikipedia.org/wiki/YAML>`_ (or alternatively JSON) format.


Caveats
=======

- All coordinates are zero-indexed.
- All intervals are closed at the beginning and open at the end (e.g. 1-10 includes position 1 but excludes position 10).
- All positions and intervals must be contained within the ``region`` interval.


YML specification
=================

Each gene YML must contain the following fields:

- ``name`` (string): 
  The name of a gene.
  Example: ``CYP2D6``.
- ``region`` (list): 
  The region where a gene is located. 
  List members are ``[chromosome`` (str), ``start,`` (int), ``end]`` (int). 
  Example: ``["1", 1, 10]``
- ``seq`` (string): 
  The canonical sequence in the region ``region``. 
  This sequence contains \*1 allele and *must be* locus-identical to 
  the reference genome (meaning that ``seq`` must not contain any indels with 
  respect to the same region in the reference; substitutions are allowed).
- ``rev_comp`` (int): 
  1 if the gene is reverse-complemented in the reference genome 
  (e.g. the first exon is first on 3' strand). 
  Zero otherwise. 
- ``exons`` (dictionary of int: list): 
  A description of exons coordinates within a gene.
  Each exon is specified as ``exon number``: ``[start, end]`` (all integers).
  Example: ``"exons": {1: [10, 20], 2: [30, 40]}`` defines 2 exons where exon 1 
  spans loci 10-20 while exon 2 spans loci 30-40.  
  This implies that intron 1 spans loci [20, 30].
- ``special_regions`` (dictionary of int: list): 
  A dictionary that defines other regions of importance within a gene.
  Key is the ordinal of the region, where each region has an ordinal number that 
  specifies its order within a structure of a gene (for example, ordinal 9 means
  that a region comes after exon and intron 8; ordinal 0 is given to a region
  before the exon 1).
  Each special region is described as a list in format 
  ``[name`` (string), ``start, end]`` (integers).
  Example: UTR region between loci 50-60 can be specified as ``["utr", 50, 60]``.
- ``pseudogenes`` (dictionary of string: dictionary): 
  A description of pseudogenes. 
  Each pseudogene is associated a name (dictionary key) and 
  a description of its exons and special regions in the same format 
  as described above.
  Can be empty if no pseudogenes are present.
  Example is given below.
- ``unique_regions`` (list of string):
  List of the regions whose coverage will be used during the copy number calling.
  This list typically includes regions that are not sequence-similar between gene
  and a pseudogene (i.e. region where read alignments are not completely ambiguous).
  Each region is defined via its name, where names for introns and exons are in
  format ``(ordinal)(i or e)`` (for example, exon 9 is ``9e``).
  Example: ``["1e", "1i", "utr"]``.
- ``common_tandems`` (list of list):
  List of common allelic tandems. 
  For example, in _CYP2D6_ \*13 is often followed by \*1. 
  This information can be passed to Aldy to aid its diplotype forming heuristics.
  Example: ``[["13", "1"], ["68", "4"]]`` indicates that \*13 and \*1, and \*68 and\*4 are commonly paired together.
- ``alleles``: dictionary of string: dictionary
  Detailed description of gene alleles.
  Key of the dictionary is the name of allele in the format
  ``(alphanum?)*(num)(alphanum?)`` (allowing names such as ``CYP2D6*2``, ``CYP2D6*2X``, 
  ``*3XD`` **but not** ``gene`` or ``*X``).
  Each key points to a dictionary with the following members:
  - ``mutations``: list of mutations describing an allele. 
    Each mutation is a dictionary with the following fields:
      + ``pos`` (int or string): Position within a reference genome. 
        If ``"pseudogene"``, describes a fusion.
      + ``op`` (string): Mutation operation. Can be:
          * ``SNP.AB`` for a SNP that changes ``A`` to ``B`` 
            (e.g. ``15:C>T`` is specified as ``SNP.CT`` with ``pos=15``).
          * ``INS.ACGT`` for an insertion that inserts ``ACGT`` at the given position.
          * ``DEL.ACGT`` for an deletion that removes ``ACGT`` at the given position.
          * ``(region)+`` for a right fusion that stitches the left portion of gene 
            to the right portion of pseudogene.
            For example, ``e9+`` indicates a fusion that takes 
            regions before exon 9 in gene (i.e. exons 1, 2, ..., 8)
            and attaches it to the pseudogene region starting at exon 9 (e.g. exons 9, 10, 11 etc).
          * ``(region)-`` for a left fusion that stitches the left portion of pseudogene 
            to the right portion of gene.
            For example, ``e5-`` indicates a fusion that takes 
            regions before exon 5 in pseudogene (i.e. exons 1, 2, ..., 4)
            and attaches it to the genic region after the exon 5 (e.g. exons 5, 6, 7 etc).
          * ``deletion`` for a whole gene deletion
      + ``functional`` (int): Describes the functionality of mutation as follows:
          * 0 for silent mutations
          * 1 for gene-disrupting SNPs
          * 2 for mutations that affect splicing
          * 3 for indels
        If a mutation has non-zero functionality, it will be considered as a mutation
        that defines a major star-allele.
      + ``dbsnp`` (list of string, optional): 
        List of dnSNP IDs associated with this mutation.
      + ``old`` (string, optional): Karolinska-style mutation description (e.g. ``-1426:C>T``).
  - ``phenotype`` (dictionary of string: string, optional):
    Describes phenotype of an allele. Currently not used by Aldy.


Here is a short example of a valid YAML file that describes a simple gene with 3 exons
and 5 alleles, together with a pseudogene::

  {
    "name": "GENE", 
    "pseudogenes": {
       "PSEUDO": {
          "special_regions": { 0: ["tmp", 0, 10] }, 
          "exons": { 1: [10, 20], 2: [30, 40], 3: [50, 60] }
       } 
    }, 
    "region": ["1", 0, 200], 
    "rev_comp": 0,
    "seq": "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",

    "special_regions": { 0: ["tmp", 100, 110] }, 
    "exons": { 1: [110, 120], 
               2: [130, 140], 
               3: [150, 160] },
    "common_tandems": [ ['1', '4'] ],
    "unique_regions": ["1e", "1i", "2e", "2i", "3e"], 

    "alleles": {
       "GENE*1": { # Normal *1 allele
          "phenotype": {"invivo": "Normal", "invitro": "Normal"}, 
          "mutations": []
       }, 
       "GENE*1B": { 
          "mutations": [ {"pos": 115, "op": "SNP.TA", "functional": 0,
                          "old": "3828:T>A",  
                          "dbsnp": ["rs28371732", "rs28371741"] # Has 2 dbSNP IDs 
                          } ]
       }, 
       "GENE*1C": { 
          "mutations": [ {"pos": 105, "op": "SNP.TA", "functional": 2} ] # Affects splicing
       }, 
       "GENE*2": {
          "mutations": [ {"pos": 111, "op": "DEL.AC", "functional": 3},  
                         {"pos": 119, "op": "INS.TT", "functional": 3} ]
       }, 
       "GENE*3": {
          "mutations": [ {"pos": 151, "op": "SNP.CT", "functional": 1}, 
                         {"pos": 148, "op": "INS.A", "functional": 0} ]
       }, 
       "GENE*4": { # Left fusion at intron 2
          "mutations": [ {"pos": "pseudogene", "op": "i2-" } ]
       }, 
       "GENE*5": { # Right fusion at exon 2
          "mutations": [ {"pos": "pseudogene", "op": "e2+" },
                         {"pos": 111, "op": "DEL.AC", "functional": 3} ]
       }, 
       "GENE*6DEL": { # Whole gene deletion 
          "mutations": [ {"op": "deletion" } ]
       }, 
    }
  }


