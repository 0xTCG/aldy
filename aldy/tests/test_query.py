# 786
# Aldy source: test_cn_real.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.

# flake8: noqa

import pytest  # noqa


from aldy.__main__ import main
from aldy.common import script_path, log
from aldy.version import __version__
from .test_full import HEADER, escape_ansi


def assert_show(monkeypatch, expected, query=None, gene="aldy.tests.resources/toy.yml"):
    lines = []

    def log_info(*args):
        s = str.format(*args)
        lines.append(s)  # type: ignore

    monkeypatch.setattr(log, "info", log_info)

    args = ["q", script_path(gene)]
    if query:
        args[1] += "*" + query
    main(args)
    expected = expected.strip()
    lines = escape_ansi("\n".join(lines)).strip()
    assert lines == expected


EXPECTED_SHOW = f"""
{HEADER}
--------------------------------------------------------------------------------
Gene TOY
--------------------------------------------------------------------------------
hg19 genome locus: 20:100000000-100000199 (NG_TEST)
5'
Pseudogenes: TOYP (ID 1)
Aminoacid:
  VRTCTYVYVR
Regions of interest:
     Region:                        TOY                      TOYP
         tmp:    20:100000100-100000110    20:100000000-100000010
          e1:    20:100000110-100000120    20:100000010-100000020
          i1:    20:100000120-100000130    20:100000020-100000030
          e2:    20:100000130-100000140    20:100000030-100000040
          i2:    20:100000140-100000150    20:100000040-100000050
          e3:    20:100000150-100000160    20:100000050-100000060
        down:    20:100000160-100000200    20:100000060-100000100
Structural alleles (deletions, conservations and fusions):
       *1: Standard copy-number configuration
       *4: TOYP fusion until i2
       *5: TOYP conservation after e2
       *6: TOY deletion
Major star-alleles:
  *1:
    Key mutations: none
    Minor alleles: *1.001 (*1),
                   *1.002 (*1B)
  *1C:
    Key mutations: 20:100000105.T>A (-5T>A, NG_TEST:105T>A, functional)
    Minor alleles: *1.003 (*1C)
  *2:
    Key mutations: 20:100000111.delAC (2delAC, NG_TEST:111delAC, frameshift),
                   20:100000119.insTT (10insTT, NG_TEST:119insTT, frameshift)
    Minor alleles: *2.001 (*2)
  *3:
    Key mutations: 20:100000151.C>T (42C>T, NG_TEST:151C>T, functional)
    Minor alleles: *3.001 (*3)
  *5:
    Key mutations: 20:100000111.delAC (2delAC, NG_TEST:111delAC, frameshift)
    Minor alleles: *5.001 (*5)
  *6:
    Key mutations: none
    Minor alleles: *6.001 (*6DEL)
"""


def test_show(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW)


EXPECTED_SHOW_MAJOR = f"""
{HEADER}
Gene TOY, major star-allele TOY*2:
  Structure: 1
  Key mutations: 20:100000111.delAC (2delAC, NG_TEST:111delAC, frameshift),
                 20:100000119.insTT (10insTT, NG_TEST:119insTT, frameshift)
  Minor star-alleles:
    *2.001:
      Legacy name: *2
      Silent mutations: none
"""

EXPECTED_SHOW_MAJOR_2 = f"""
{HEADER}
Gene TOY, major star-allele TOY*1C:
  Structure: 1
  Key mutations: 20:100000105.T>A (-5T>A, NG_TEST:105T>A, functional)
  Minor star-alleles:
    *1.003:
      Legacy name: *1C
      Silent mutations: none
"""


def test_show_major(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW_MAJOR, "2")
    assert_show(monkeypatch, EXPECTED_SHOW_MAJOR_2, "1C")


EXPECTED_SHOW_MINOR = f"""
{HEADER}

Gene TOY, minor star-allele TOY*1.002:
  Structure: 1
  Major star-allele: TOY*1
  Legacy name: *1B
  Key mutations: none
  Silent mutations: 20:100000115.T>A (rs28371732, 6T>A, NG_TEST:115T>A)

"""  # noqa


def test_show_minor(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW_MINOR, "1B")


EXPECTED_SHOW_CN = f"""
{HEADER}
Gene CYP2D6, structural allele CYP2D6*13:
  Structure: CYP2D6   CYP2D7
         up: 0        1
       utr5: 0        1
         e1: 0        1
         i1: 1        0
         e2: 1        0
         i2: 1        0
         e3: 1        0
         i3: 1        0
         e4: 1        0
         i4: 1        0
         e5: 1        0
         i5: 1        0
         e6: 1        0
         i6: 1        0
         e7: 1        0
         i7: 1        0
         e8: 1        0
         i8: 1        0
         e9: 1        0
       utr3: 1        0
        ins: 1        0
        pce: 0        0
        rep: 1        0
"""


def test_show_cn(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW_CN, "13", "aldy.resources.genes/cyp2d6.yml")
