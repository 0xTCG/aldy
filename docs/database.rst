Database Specification
**********************

Aldy gene descriptions are stored in the `YAML <https://en.wikipedia.org/wiki/YAML>`_ format.


Caveats
=======

- All coordinates start from 1.
- All intervals are closed at the beginning and open at the end (e.g., 1-10 includes 1 but excludes 10).


YAML specification
==================

Each valid YAML file must contain the following fields:

- ``name`` (string): the name of a gene.
- ``version`` (string): the gene database version.
- ``generated`` (string): date of the last update.
- ``pharmvar`` (int): PharmVar ID (when applicable)
- ``ensembl`` (string): ENSEMBL gene ID (when applicable)
- ``structure`` (dict):

  - ``genes`` (list): a list of genes and pseudogenes
  - ``regions`` (dict): locations of the genic regions **within the reference genomes**.
    Each region should be named and specified as ``name: [start1, end1, start2, end2, ...]``.
    The number of intervals should match the number of the genes.
    Regions between exons (named as ``e1``, ``e2``, etc.) are automatically named as ``i1``, ``i2``
    and so on. For other regions, ensure that there are no gaps between them.
  - ``cn_regions`` (list): a list of the names of the regions whose coverage will be used for copy number calling.
    These typically includes regions whose gene and pseudogene sequence differs significantly
    (i.e., the regions where the read alignments are not ambiguous).
  - ``tandems`` (list): a list of common allelic tandems.
    For example, _CYP2D6\*13_ is often followed by _CYP2D6\*1_.
    This information can be passed to Aldy to aid its diplotype detection heuristics.
    Example: ``["13", "1"]`` indicates that \*13 and \*1 are commonly paired together.

- ``reference`` (dict): the RefSeq description.

  - ``name`` (str): the RefSeq sequence name.
  - ``mappings`` (dict): the alignment of the RefSeq sequence to each of the reference genomes.
    Each alignment is expressed as `[chromosome, start, end, strand, cigar]`.
  - ``exons`` (list): the locations of the main gene's exonic regions **within the RefSeq sequence**.
  - ``seq`` (string): the RefSeq sequence.

- ``alleles`` (dict): a detailed description of gene alleles.
  The name of each allele follows the following format:
  ``GENE_NAME*ALLELE_NAME`` (allowing names such as ``CYP2D6*2``, ``CYP2D6*2X``,
  ``DPYD*hello`` **but not** ``gene`` or ``rs123``).
  Each dictionary value is another dictionary that describes the allele with the following keys:

  - ``pharmvar``: PharmVar allele URL
  - ``activity``: allele activity
  - ``evidence``: the evidence for allele activity
  - ``label``: an alternative allele label
  - ``mutations``: the list of variants that define the allele.
    Each variants is expressed as ``[pos, op, rsID, functionality]``:::

      + ``pos`` (int): the position within the RefSeq sequence.
        Alternatively, a string ``"pseudogene"`` describes a pseudogene fusion.
      + ``op`` (str): variant description::
          * ``A>B`` for an SNP that changes ``A`` to ``B``
            (e.g., ``15.C>T`` is specified as ``C>T`` with ``pos=15``).
          * ``insACGT`` for an insertion that inserts ``ACGT`` at the given position.
          * ``delACGT`` for a deletion that removes ``ACGT`` at the given position.
          * ``(region)+`` for a right fusion that stitches the left portion of the gene
            to the right portion of the pseudogene.
            For example, ``e9+`` indicates a fusion that takes
            regions before exon 9 in gene (i.e., exons 1, 2, ..., 8)
            and attaches it to the pseudogene region starting at exon 9 (e.g., exons 9, 10, 11 and so on).
          * ``(region)-`` for a left fusion that stitches the left portion of the pseudogene
            to the right portion of the gene.
            For example, ``e5-`` indicates a fusion that takes
            regions before exon 5 in pseudogene (i.e., exons 1, 2, ..., 4)
            and attaches it to the genic region after exon 5 (e.g., exons 5, 6, 7 and so on).
          * ``deletion`` for the whole gene deletion
      + ``rsID`` (str): variant rsID. ``-`` if not available.
      + ``functionality`` (str): functionality of the variant.
        Silent variants do not have this field. Any variant that has this field specified
        is considered as a functional variant that defines a major star-allele.


See `<https://github.com/0xTCG/aldy/blob/master/aldy/tests/resources/toy.yml>`_ and `<https://github.com/0xTCG/aldy/blob/master/aldy/resources/genes/cyp2d6.yml>`_ for the complete examples of database files.
