# 786
# Aldy source: test_cn_real.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest
import os
import re

from aldy.gene import GeneRegion
from aldy.__main__ import _genotype
from aldy.common import *


def escape_ansi(line):
   """
   Inspired by https://www.tutorialspoint.com/How-can-I-remove-the-ANSI-escape-sequences-from-a-string-in-python
   """
   return re.compile(r'(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]').sub('', line)


def assert_file(monkeypatch, file, expected, params=None):
   lines = []
   def log_info(*args):
      s = str.format(*args)
      lines.append(s)
   monkeypatch.setattr(log, "info", log_info)
   solver = os.getenv('ALDY_SOLVER', default='gurobi')
   args = DictWrapper( {**{'file': file,
                               'profile': 'illumina',
                               'threshold': 50,
                               'gap': 0,
                               'solver': solver},
                               **(params or {})})
   _genotype('cyp2d6', None, args)
   expected = expected.strip()
   lines = escape_ansi('\n'.join(lines)).strip()
   assert lines == expected


def test_NA10860(monkeypatch):
   file = script_path('aldy.resources/NA10860.bam')
   expected = td("""
      Genotyping sample NA10860.bam...
      Potential CYP2D6 copy number configurations for NA10860:
         1: 2x*1,1x*36
            Confidence: 1.00 (score = 6.73)
         2: 2x*1,1x*61
            Confidence: 1.00 (score = 6.73)
      Potential major CYP2D6 star-alleles for NA10860:
         1: 1x*1, 1x*4, 1x*4.i
            Confidence: 1.00 (score = 1.33)
         2: 1x*4, 1x*4.f, 1x*61
            Confidence: 1.00 (score = 1.33)
         3: 1x*4, 1x*4.f, 1x*83
            Confidence: 1.00 (score = 1.33)
         4: 1x*4, 1x*4.h, 1x*36
            Confidence: 1.00 (score = 1.33)
         5: 1x*4.f, 1x*4.i, 1x*39
            Confidence: 1.00 (score = 1.33)
         6: 1x*4.h, 1x*4.i, 1x*10
            Confidence: 1.00 (score = 1.33)
      Best CYP2D6 star-alleles for NA10860:
         1: *1/*4+*4
            Minor: *1 +42528223:SNP.GA, *4AW, *4N -42522391:SNP.GA
            Confidence: 1.00 (score = 7.77)
      CYP2D6 result:
        *1/*4+*4                       (*1, *4, *4.i)
   """)
   assert_file(monkeypatch, file, expected)


def test_NA10860_gap(monkeypatch):
   file = script_path('aldy.resources/NA10860.bam')
   expected = td("""
      Genotyping sample NA10860.bam...
      Potential CYP2D6 copy number configurations for NA10860:
         1: 2x*1,1x*36
            Confidence: 1.00 (score = 6.73)
         2: 2x*1,1x*61
            Confidence: 1.00 (score = 6.73)
         3: 2x*1,1x*5,1x*36
            Confidence: 0.94 (score = 7.23)
         4: 2x*1,1x*5,1x*61
            Confidence: 0.94 (score = 7.23)
      Potential major CYP2D6 star-alleles for NA10860:
         1: 1x*1, 1x*4, 1x*4.i
            Confidence: 1.00 (score = 1.33)
         2: 1x*4, 1x*4.f, 1x*61
            Confidence: 1.00 (score = 1.33)
         3: 1x*4, 1x*4.f, 1x*83
            Confidence: 1.00 (score = 1.33)
         4: 1x*4, 1x*4.h, 1x*36
            Confidence: 1.00 (score = 1.33)
         5: 1x*4.f, 1x*4.i, 1x*39
            Confidence: 1.00 (score = 1.33)
         6: 1x*4.h, 1x*4.i, 1x*10
            Confidence: 1.00 (score = 1.33)
         7: 1x*1, 1x*4, 1x*4.i, 1x*5
            Confidence: 0.96 (score = 1.42)
         8: 1x*4, 1x*4.f, 1x*5, 1x*61
            Confidence: 0.96 (score = 1.42)
         9: 1x*4, 1x*4.f, 1x*5, 1x*83
            Confidence: 0.96 (score = 1.42)
        10: 1x*4, 1x*4.h, 1x*5, 1x*36
            Confidence: 0.96 (score = 1.42)
        11: 1x*4.f, 1x*4.i, 1x*5, 1x*39
            Confidence: 0.96 (score = 1.42)
        12: 1x*4.h, 1x*4.i, 1x*5, 1x*10
            Confidence: 0.96 (score = 1.42)
      Best CYP2D6 star-alleles for NA10860:
         1: *1/*4+*4
            Minor: *1 +42528223:SNP.GA, *4AW, *4N -42522391:SNP.GA
            Confidence: 1.00 (score = 7.77)
      CYP2D6 result:
        *1/*4+*4                       (*1, *4, *4.i)
   """)
   assert_file(monkeypatch, file, expected, {'gap': 0.10})


def test_hard(monkeypatch):
   file = script_path('aldy.resources/HARD.dump')
   expected = td("""
      Genotyping sample HARD.dump...
      Potential CYP2D6 copy number configurations for HARD:
         1: 4x*1,2x*80
            Confidence: 1.00 (score = 25.81)
      Potential major CYP2D6 star-alleles for HARD:
         1: 2x*2, 1x*4, 1x*10, 2x*80/4
            Confidence: 1.00 (score = 3.21)
         2: 2x*2, 1x*4, 1x*4.b, 1x*80/10, 1x*80/4
            Confidence: 1.00 (score = 3.21)
      Best CYP2D6 star-alleles for HARD:
         1: *2+*2+*80+*80/*4+*4
            Minor: *2MW -42527541:DEL.TC, *2MW -42527541:DEL.TC, *4AW, *4AW,
                   *4EW +42522311:SNP.CT +42523002:SNP.GA, *88 +42522311:SNP.CT
            Confidence: 1.00 (score = 20.72)
      CYP2D6 result:
        *2+*2+*80+*80/*4+*4            (*2, *2, *4, *4.b, *80/10, *80/4)
   """)
   assert_file(monkeypatch, file, expected)
