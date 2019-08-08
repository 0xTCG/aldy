#!/usr/bin/env python
# 786

# Aldy source: common.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple, Iterable

import pkg_resources 
import os
import re
import time
import pprint
import logbook
import textwrap
import collections


PROTEINS = {
   'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V', 'TTC': 'F',
   'CTC': 'L', 'ATC': 'I', 'GTC': 'V', 'TTA': 'L', 'CTA': 'L',
   'ATA': 'I', 'GTA': 'V', 'TTG': 'L', 'CTG': 'L', 'ATG': 'M',
   'GTG': 'V', 'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
   'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'TCA': 'S',
   'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'TCG': 'S', 'CCG': 'P',
   'ACG': 'T', 'GCG': 'A', 'TAT': 'Y', 'CAT': 'H', 'AAT': 'N',
   'GAT': 'D', 'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
   'TAA': 'X', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'TAG': 'X',
   'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'TGT': 'C', 'CGT': 'R',
   'AGT': 'S', 'GGT': 'G', 'TGC': 'C', 'CGC': 'R', 'AGC': 'S',
   'GGC': 'G', 'TGA': 'X', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
   'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}
"""dict[str, str]: Codon table (stop codon is X)."""


REV_COMPLEMENT = {
   'A': 'T', 'T': 'A',
   'C': 'G', 'G': 'C',
}
"""dict[str, str]: Reverse-complement nucleotide table."""


log = logbook.Logger('Aldy')
"""Default console logger."""


SOLUTION_PRECISION = 1e-2
"""float: Solution precision (all objectives that differ less than this are considered equal)"""


class AldyException(Exception):
   """
   Base Aldy exception class.
   """
   pass


class GRange(collections.namedtuple('GRange', ['chr', 'start', 'end'])):
   """
   Describes the range in the reference genome (e.g. chr22:10-20).
   Immutable class.

   Attributes:
      chr (str): Chromosome identifier.
      start (int): Starting position of the interval.
      end (int): Ending position of the interval.

   Notes:
      Has custom printer (``__str__``).
   """
   
   def samtools(self, pad_left=500, pad_right=1, prefix='') -> str:
      """
      Samtools-compatible region string (e.g. chr1:10-20).

      Returns:
         str
      """
      return '{}:{}-{}'.format(prefix + self.chr, self.start - 500, self.end + 1)
   
   def __str__(self):
      return self.samtools(0, 0, '')


class GeneRegion(collections.namedtuple('GeneRegion', ['number', 'kind'])):
   """
   Describes the region of a gene.
   
   Attributes:
      number (int): 
         The region number (e.g. for exon 9, number is 9).
      kind (str): 
         Type of the region. Typically 'e' (for **e**\ xon), 'i' (for **i**\ ntron) 
         but can be anything else.

   Notes:
      Has custom printer (``__str__``).
   """

   def __str__(self):
      return 'GR({}.{})'.format(self.number, self.kind)


### Aldy auxiliaries 


def allele_number(x: str) -> str:
   """
   Returns:
      str: Major allele number of the allele name string (e.g. ``'12A'`` -> ``12``).
   """
   p = re.split(r'(\d+)', x)
   return p[1]


def allele_sort_key(x: str) -> Tuple[int, str]:
   """
   Returns: 
      tuple[int, str]: Key for sorting allele names (e.g. ``'13a'`` -> ``(13, 'a')``).
   """
   p = re.split(r'(\d+)', x)
   return (int(p[1]), ''.join(p[2:]))


def rev_comp(seq: str) -> str:
   """
   Returns:
      str: Reverse-complemented DNA sequence.
   """

   return ''.join([REV_COMPLEMENT[x] for x in seq[::-1]])


def seq_to_amino(seq: str) -> str:
   """
   Returns:
      str: Protein sequence formed from a DNA sequence.
   """
   
   return ''.join(PROTEINS[seq[i:i + 3]] for i in range(0, len(seq) - len(seq) % 3, 3))


### Language auxiliaries


def sorted_tuple(x: Iterable) -> tuple:
   """
   Sorts a tuple.
   """
   return tuple(sorted(x))


def td(s: str) -> str:
   """
   Abbreviation for textwrap.dedent. Useful for stripping indentation 
   in multi-line docstrings.
   """
   return textwrap.dedent(s)


def static_vars(**kwargs):
   """
   Decorator that adds static variables to a function.
   
   Usage::
      
      @static_vars(var=init_val)
   """
   def decorate(func):
      for k in kwargs:
         setattr(func, k, kwargs[k])
      return func
   return decorate


def timing(f):
   """
   Decorator for timing a function.
   Prints the time spent in function after the function is completed.
   
   Usage:: 

      @timing 

   (without any parameters).
   """
   def wrap(*args, **kwargs):
      time1 = time.time()
      ret = f(*args, **kwargs)
      time2 = time.time()
      log.warn('Time needed: ({:.1f})', time2 - time1)
      return ret
   return wrap


def pp(x) -> str:
   """
   Returns:
      str: Pretty-printed variable string.
   """
   return pprint.pformat(x)


def pr(x):
   """
   Pretty-prints a variable to stdout.
   """
   pprint.pprint(x)


def script_path(key: str) -> str:
   """
   Args: 
      key (str): a resource to be extracted. 
      Specify with ``path/file`` (e.g. ``aldy.resources/test.txt``).

   Returns:
      str: Full path of a package resource.

   Raises:
      :obj:`aldy.common.AldyException`.
   """
   components = key.split('/')
   if len(components) < 2:
      raise AldyException(f'"{key}"" is not valid resource name')
   return pkg_resources.resource_filename(components[0], '/'.join(components[1:]))


def colorize(text: str, color:str = 'green') -> str:
   """
   Returns:
      str: Colorized string (on xterm-compatible terminals) with a given color.
   """
   return logbook._termcolors.colorize(color, text)


def check_path(cmd: str) -> bool:
   """
   Returns:
      bool: ``True`` if a command ``cmd`` is an executable in ``PATH`` or local directory.
   """

   def is_exe(path): 
      """
      Based on http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028
      """
      return os.path.isfile(path) and os.access(path, os.X_OK)

   if not is_exe(cmd):
      for path in os.environ["PATH"].split(os.pathsep):
         path = path.strip('"')
         if is_exe(os.path.join(path, cmd)):
            return True
      return False
   return True
