# 786
# Aldy source: common.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Iterable, Any, List
import pkg_resources
import re
import time
import pprint
import logbook
import textwrap
import collections


PROTEINS = {
    "TTT": "F",
    "CTT": "L",
    "ATT": "I",
    "GTT": "V",
    "TTC": "F",
    "CTC": "L",
    "ATC": "I",
    "GTC": "V",
    "TTA": "L",
    "CTA": "L",
    "ATA": "I",
    "GTA": "V",
    "TTG": "L",
    "CTG": "L",
    "ATG": "M",
    "GTG": "V",
    "TCT": "S",
    "CCT": "P",
    "ACT": "T",
    "GCT": "A",
    "TCC": "S",
    "CCC": "P",
    "ACC": "T",
    "GCC": "A",
    "TCA": "S",
    "CCA": "P",
    "ACA": "T",
    "GCA": "A",
    "TCG": "S",
    "CCG": "P",
    "ACG": "T",
    "GCG": "A",
    "TAT": "Y",
    "CAT": "H",
    "AAT": "N",
    "GAT": "D",
    "TAC": "Y",
    "CAC": "H",
    "AAC": "N",
    "GAC": "D",
    "TAA": "X",
    "CAA": "Q",
    "AAA": "K",
    "GAA": "E",
    "TAG": "X",
    "CAG": "Q",
    "AAG": "K",
    "GAG": "E",
    "TGT": "C",
    "CGT": "R",
    "AGT": "S",
    "GGT": "G",
    "TGC": "C",
    "CGC": "R",
    "AGC": "S",
    "GGC": "G",
    "TGA": "X",
    "CGA": "R",
    "AGA": "R",
    "GGA": "G",
    "TGG": "W",
    "CGG": "R",
    "AGG": "R",
    "GGG": "G",
}
"""Codon table (stop codon is X)."""


REV_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
"""Reverse-complement DNA table."""


log = logbook.Logger("Aldy")
"""Default console logger."""


SOLUTION_PRECISION = 1e-2
"""
Solution precision (all values whose absolute difference falls below the specified
precision are considered equal).
"""


class AldyException(Exception):
    """Aldy exception class."""

    pass


class GRange(collections.namedtuple("GRange", ["chr", "start", "end"])):
    """Reference genome range (e.g. `chr22:10-20`). Immutable."""

    def samtools(self, pad_left=500, pad_right=1, prefix="") -> str:
        """Samtools-compatible region representation (e.g. chr1:10-20).

        :param pad_left: Left padding.
        :param pad_right: Right padding.
        :param prefix: Chromosome prefix."""

        return "{}:{}-{}".format(
            prefix + self.chr, self.start - pad_left, self.end + pad_right
        )

    def __str__(self):
        return self.samtools(0, 0, "")


def allele_name(x: str) -> str:
    """:returns: Major allele number of the star-allele name (e.g. `'12A'` -> `12`)."""
    if "*" in x:
        x = x.split("*", maxsplit=1)[1]
    return x.replace("/", "_")


def rev_comp(seq: str) -> str:
    """:returns: Reverse-complemented DNA sequence."""

    return "".join([REV_COMPLEMENT.get(x, x) for x in seq[::-1]])


def seq_to_amino(seq: str) -> str:
    """:returns: Protein sequence formed from the provided DNA sequence."""

    return "".join(
        PROTEINS[seq[i : i + 3]] for i in range(0, len(seq) - len(seq) % 3, 3)
    )


def freezekey(x):
    """Hashing support for dictionaries."""
    return tuple(i[1] for i in sorted(x[0].items())) + tuple(
        i[1] for i in sorted(x[1].items())
    )


def sorted_tuple(x: Iterable) -> tuple:
    """:returns: Sorted tuple."""
    return tuple(sorted(x))


def td(s: str) -> str:
    """
    Abbreviation for textwrap.dedent. Used for stripping indentation in multi-line
    docstrings.
    """
    return textwrap.dedent(s)


class Timing:
    """
    Context manager for timing code blocks. Prints the time spent in the function after
    it is completed.
    """

    def __init__(self, name="Block"):
        self.name = name

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *_):
        self.end = time.time()
        log.debug("{} took {}", self.name, self.end - self.start)


def pp(x) -> str:
    """:returns: Pretty-printed variable string."""

    return pprint.pformat(x)


def script_path(key: str) -> str:
    """
    Obtain the full path of a resource.

    :param key: resource to be extracted.
    :param key: resource to be extracted in `path/file` format
        (e.g., `aldy.resources/test.txt`).
    :returns: Full path of the resource.
    :raises: :py:class:`aldy.common.AldyException` if the resource does not exist.
    """
    components = key.split("/")
    if len(components) < 2:
        raise AldyException(f'"{key}"" is not valid resource name')
    return pkg_resources.resource_filename(components[0], "/".join(components[1:]))


def colorize(text: str, color: str = "green") -> str:
    """:returns: xterm-compatible colorized string with a given color."""

    import logbook._termcolors

    return logbook._termcolors.colorize(color, text)


def parse_cn_region(cn_region):
    """
    :returns: :py:class:`GRange` object that represents the user-provided CN region in
        Samtools format (i.e., `chr1:100-200`).
    :raises: :py:class:`aldy.common.AldyException` if the region is invalid.
    """
    if cn_region is not None:
        r = re.match(r"^(.+?):(\d+)-(\d+)$", cn_region)
        if not r:
            raise AldyException(
                f"Parameter --cn-neutral={cn_region} cannot be parsed. "
                + "Must be chr:start-end (where start and end are numbers)"
            )
        ch = r.group(1)
        if ch.startswith("chr"):
            ch = ch[3:]
        return GRange(ch, int(r.group(2)), int(r.group(3)))
    return None


def chr_prefix(ch: str, chrs: List[str]) -> str:
    """
    Check if a chromosome needs "chr" prefix given the available chromosomes.
    :returns: Chromosome prefix if the chromosome does not have it.
    """
    if ch not in chrs and "chr" + ch in chrs:
        return "chr"
    return ""


class JsonDict(dict):
    """
    Dictionary that adds a dictionary for each missing key. Used to ease handling and
    populating JSON objects.
    """

    def __getitem__(self, key):
        if key not in self:
            self[key] = JsonDict()
        return self.get(key)


json: Any = JsonDict()
