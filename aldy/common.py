#!/usr/bin/env python
# 786

# Aldy source: common.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple

import pkg_resources 
import os
import re
import time
import pprint
import logbook
import textwrap
import collections


PROTEINS = { # X is stop
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

REV_COMPLEMENT = {
   'A': 'T', 'T': 'A',
   'C': 'G', 'G': 'C',
}

LOG_FORMAT = '{record.message}'


log = logbook.Logger('Aldy')


class AldyException(Exception):
   pass


class GRange(collections.namedtuple('GRange', ['chr', 'start', 'end'])):
   """
   Describes the range in the reference genome (e.g. chr22:10-20)
   """
   
   def samtools(self, pad_left=500, pad_right=1, prefix='') -> str:
      """Samtools-compatible region string"""
      return '{}:{}-{}'.format(prefix + self.chr, self.start - 500, self.end + 1)
   
   def __str__(self):
      return self.samtools(0, 0, '')


class GeneRegion(collections.namedtuple('GeneRegion', ['number', 'kind'])):
   """
   Describes the region in the gene.
   Members:
      number (int): 
         region number (e.g. for exon 9 it is 9)
      kind (str): 
         type of the region. Can be 'e' (EXON), 'i' (INTRON) or anything else.
   """
   def __repr__(self):
      return 'GR({}.{})'.format(self.number, self.kind)


### Aldy auxiliaries 


def allele_number(x: str) -> str:
   """Returns a major allele number for an allele string (e.g. '12A' -> 12)"""
   p = re.split(r'(\d+)', x)
   return p[1]


def allele_sort_key(x: str) -> Tuple[int, str]:
   """Sort key for allele names (e.g. '13a' -> (13, 'a')). Useful for numeric sorting."""
   p = re.split(r'(\d+)', x)
   return (int(p[1]), ''.join(p[2:]))


def rev_comp(seq: str) -> str:
   """Reverse-complement a DNA sequence"""
   return ''.join([REV_COMPLEMENT[x] for x in seq[::-1]])


def seq_to_amino(seq: str) -> str:
   """Converts DNA sequence to protein sequence"""
   return ''.join(PROTEINS[seq[i:i + 3]] for i in range(0, len(seq) - len(seq) % 3, 3))


### Language auxiliaries


def sorted_tuple(x: tuple) -> tuple:
   """Sorts a tuple"""
   return tuple(sorted(x))


def td(s: str) -> str:
   """Abbreviation for textwrap.dedent (useful for stripping indentation in multi-line strings)"""
   return textwrap.dedent(s)


def nt(*args):
   """Abbreviation for a collections.namedtuple"""
   return collections.namedtuple(*args)


def static_vars(**kwargs):
   """Decorator that adds static variables to a function"""
   def decorate(func):
      for k in kwargs:
         setattr(func, k, kwargs[k])
      return func
   return decorate


def timing(f):
   """Decorator for timing a function"""
   def wrap(*args, **kwargs):
      time1 = time.time()
      ret = f(*args, **kwargs)
      time2 = time.time()
      log.warn('Time needed: ({:.1f})', time2 - time1)
      return ret
   return wrap


@static_vars(pp=pprint.PrettyPrinter(indent=4))
def pp(x):
   """Returns a pretty-printed variable string"""
   return pp.pformat(x)
def pr(x):
   """Pretty-prints a variable to stdout"""
   return pprint.pprint(x)


def script_path(path: str, file: str) -> str:
   """Gives a full path of a package resource"""
   return pkg_resources.resource_filename(path, file)


def colorize(text: str, color:str = 'green') -> str:
   """Colorizes a string (for xterm) with a given color"""
   return logbook._termcolors.colorize(color, text)


def check_path(cmd: str) -> bool:
   """Checks whether `cmd` is in PATH or local directory"""

   def is_exe(path): 
      """Based on http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028"""
      return os.path.isfile(path) and os.access(path, os.X_OK)

   if not is_exe(cmd):
      for path in os.environ["PATH"].split(os.pathsep):
         path = path.strip('"')
         if is_exe(os.path.join(path, cmd)):
            return True
      return False
   return True
