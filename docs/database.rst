Database Specification
**********************

Aldy gene descriptions are stored in the database filesd in `YAML <https://en.wikipedia.org/wiki/YAML>`_ (or alternatively JSON) format.


Caveats
=======

- All coordinates start from 1.
- All intervals are closed at the beginning and open at the end (e.g., 1-10 includes 1 but excludes 10).


YML specification
=================

Each gene YML must contain the following fields:

- ``name`` (string): The name of a gene.
- ``version`` (string): The gene database version.
- ``generated`` (string): Date of the last update.
- ``pharmvar`` (int): PharmVar ID (when applicable)
- ``ensembl`` (string): ENSEMBL gene ID (when applicable)
- ``structure`` (dict):
  - ``genes`` (list): List of genes and pseudogenes
  - ``regions`` (dict of dict):
    Location of the genic regions **within the reference genomes**.
    Each region should be names and specified as ``name: [start1, end1, start2, end2, ...]``.
    The number of intervals should match the number of the genes.
    Regions between exons (named as ``e1``, ``e2``, etc.) are automatically named as ``i1``, ``i2``
    and so on. For other regions you need to ensure that there are no gaps between them.
  - ``cn_regions`` (list):
    List of the regions names whose coverage will be used for copy number calling.
    This typically includes regions whose gene and pseudogene sequence differs significantly
    (i.e., the regions where the read alignments are not ambiguous).
  - ``tandems`` (list):
    List of common allelic tandems.
    For example, _CYP2D6\*13_ is often followed by _CYP2D6\*1_.
    This information can be passed to Aldy to aid its diplotype detection heuristics.
    Example: ``["13", "1"]`` indicates that \*13 and \*1 are commonly paired together.
- ``reference`` (dict): RefSeq description.
  - ``name`` (str): RefSeq sequence name.
  - ``mappings`` (dict):
    Alignment of the RefSeq sequence to each of the reference genomes.
    Each alignment is expressed as `[chromosome, start, end, strand, cigar]`.
  - ``exons`` (list):
    Location of the main gene's exonic regions **within the RefSeq sequence**.
  - ``seq`` (string):
    RefSeq sequence.
- ``alleles``: dict
  Detailed description of gene alleles.
  The name of each allele follows the following format:
  ``GENE_NAME*ALLELE_NAME`` (allowing names such as ``CYP2D6*2``, ``CYP2D6*2X``,
  ``DPYD*hello`` **but not** ``gene`` or ``rs123``).
  Each dictionary value is another dictionary that describes the allele with the following keys:
  - ``pharmvar``: PharmVar allele URL
  - ``activity``: Allele activity
  - ``evidence``: Evidence for allele activity
  - ``label``: Alternative allele label
  - ``mutations``: list of variants that define the allele.
    Each variants is expressed as ``[pos, op, rsID, functionality]``.
      + ``pos`` (int): Position within the RefSeq sequence.
        Alternatively, a string ``"pseudogene"`` describes a pseudogene fusion.
      + ``op`` (string): Variant operation. Can be:
          * ``A>B`` for a SNP that changes ``A`` to ``B``
            (e.g., ``15.C>T`` is specified as ``C>T`` with ``pos=15``).
          * ``insACGT`` for an insertion that inserts ``ACGT`` at the given position.
          * ``delACGT`` for an deletion that removes ``ACGT`` at the given position.
          * ``(region)+`` for a right fusion that stitches the left portion of gene
            to the right portion of pseudogene.
            For example, ``e9+`` indicates a fusion that takes
            regions before exon 9 in gene (i.e., exons 1, 2, ..., 8)
            and attaches it to the pseudogene region starting at exon 9 (e.g., exons 9, 10, 11 etc).
          * ``(region)-`` for a left fusion that stitches the left portion of pseudogene
            to the right portion of gene.
            For example, ``e5-`` indicates a fusion that takes
            regions before exon 5 in pseudogene (i.e., exons 1, 2, ..., 4)
            and attaches it to the genic region after the exon 5 (e.g., exons 5, 6, 7 etc).
          * ``deletion`` for the whole gene deletion
      + ``rsID`` (str): Variant rsID. ``-`` if not available.
      + ``functionality`` (str): Describes the functionality of the variant.
        Silent variants do not have this field. Any variant that has this field specified
        is considered as a functional variant that defines a major star-allele.


See `aldy/tests/resources/toy.yml` and `aldy/resources/genes/cyp2d6.yml`
for complete examples of database files.
