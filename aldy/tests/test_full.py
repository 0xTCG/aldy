# 786
# Aldy source: test_cn_real.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa
import os
import re
import subprocess

from tempfile import NamedTemporaryFile as tmpfile

from aldy.__main__ import _genotype
from aldy.common import DictWrapper, script_path, td, log
from aldy.cn import LEFT_FUSION_PENALTY


def escape_ansi(line):
    """
    Inspired by 
    https://www.tutorialspoint.com/How-can-I-remove-the-ANSI-escape-sequences-from-a-string-in-python
    """  # noqa
    return re.compile(r"(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]").sub("", line)


def assert_file(monkeypatch, file, expected, params=None, output=None):
    lines = []

    def log_info(*args):
        s = str.format(*args)
        lines.append(s)

    monkeypatch.setattr(log, "info", log_info)
    solver = os.getenv("ALDY_SOLVER", default="gurobi")
    args = DictWrapper(
        {
            **{
                "file": file,
                "profile": "illumina",
                "threshold": 50,
                "gap": 0,
                "solver": solver,
                "fusion_penalty": LEFT_FUSION_PENALTY,
                "max_minor_solutions": 10,
            },
            **(params or {}),
        }
    )
    _genotype("cyp2d6", output, args)
    expected = expected.strip()
    lines = escape_ansi("\n".join(lines)).strip()
    assert lines == expected


def test_NA10860(monkeypatch):
    file = script_path("aldy.tests.resources/NA10860.bam")
    expected = td(
        """
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
           2: *1/*4+*4
              Minor: *1, *4AW +42528223:SNP.GA, *4N -42522391:SNP.GA
              Confidence: 1.00 (score = 7.77)
        CYP2D6 results:
          *1/*4+*4                       (*1, *4, *4.i)"""
    )
    assert_file(monkeypatch, file, expected)


def test_NA10860_gap(monkeypatch):
    file = script_path("aldy.tests.resources/NA10860.bam")
    expected = td(
        """
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
           2: *1/*4+*4
              Minor: *1, *4AW +42528223:SNP.GA, *4N -42522391:SNP.GA
              Confidence: 1.00 (score = 7.77)
        CYP2D6 results:
          *1/*4+*4                       (*1, *4, *4.i)"""
    )
    assert_file(monkeypatch, file, expected, {"gap": 0.10})


def test_hard(monkeypatch):
    file = script_path("aldy.tests.resources/HARD.dump")
    expected = td(
        """
        Genotyping sample HARD.dump...
        Potential CYP2D6 copy number configurations for HARD:
           1: 3x*1,1x*68,1x*80
              Confidence: 1.00 (score = 10.07)

        Potential major CYP2D6 star-alleles for HARD:
           1: 2x*2, 1x*4, 1x*68, 1x*80/4
              Confidence: 1.00 (score = 2.46)

        Best CYP2D6 star-alleles for HARD:
           1: *2+*2/*68+*4+*80
              Minor: *2MW -42527541:DEL.TC, *2MW -42527541:DEL.TC, *4AW
                     +42523002:SNP.GA, *4AW, *68 -42526483:SNP.CA
                     -42528223:SNP.GA
              Confidence: 1.00 (score = 15.49)
           2: *2+*2/*68+*4+*80
              Minor: *2MW -42527541:DEL.TC, *2MW -42527541:DEL.TC, *4AW, *4AW
                     +42523002:SNP.GA, *68 -42526483:SNP.CA -42528223:SNP.GA
              Confidence: 1.00 (score = 15.49)
        CYP2D6 results:
          *2+*2/*68+*4+*80               (*2, *2, *4, *68, *80/4)"""
    )
    assert_file(monkeypatch, file, expected, {"profile": "pgrnseq-v1"})


def test_hard_fusion(monkeypatch):
    file = script_path("aldy.tests.resources/HARD.dump")
    expected = td(
        """
        Genotyping sample HARD.dump...
        Potential CYP2D6 copy number configurations for HARD:
           1: 4x*1
              Confidence: 1.00 (score = 11.03)

        Potential major CYP2D6 star-alleles for HARD:
           1: 2x*2, 1x*4, 1x*4.b
              Confidence: 1.00 (score = 2.84)

        Best CYP2D6 star-alleles for HARD:
           1: *2+*2/*4+*4
              Minor: *2MW -42527541:DEL.TC, *2MW -42527541:DEL.TC, *4AW
                     +42523002:SNP.GA, *4EW
              Confidence: 1.00 (score = 18.22)
           2: *2+*2/*4+*4
              Minor: *2MW -42527541:DEL.TC, *2MW -42527541:DEL.TC, *4AW, *4EW
                     +42523002:SNP.GA
              Confidence: 1.00 (score = 18.22)
        CYP2D6 results:
          *2+*2/*4+*4                    (*2, *2, *4, *4.b)"""
    )
    assert_file(
        monkeypatch, file, expected, {"profile": "pgrnseq-v1", "fusion_penalty": 5}
    )


def test_NA10860_debug(monkeypatch):
    file = script_path("aldy.tests.resources/NA10860.bam")
    expected = td(
        """
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
           2: *1/*4+*4
              Minor: *1, *4AW +42528223:SNP.GA, *4N -42522391:SNP.GA
              Confidence: 1.00 (score = 7.77)
        CYP2D6 results:
          *1/*4+*4                       (*1, *4, *4.i)
        Preparing debug archive..."""
    )
    with tmpfile(suffix=".tar.gz") as tmp:
        with tmpfile(mode="w") as out, tmpfile(mode="w") as out_log:
            assert_file(
                monkeypatch,
                file,
                expected,
                {"debug": tmp.name[:-7], "log": out_log.name},
                out,
            )
            out.flush()
            out_log.flush()

            with open(script_path("aldy.tests.resources/NA10860.out.expected")) as f:
                expected = f.read()
            with open(out.name) as f:
                produced = f.read()
            assert produced == expected

            with open(out_log.name) as f:
                log = f.read()
            s = "Major solver: MajorSol[1.33; sol=(1x*4, 1x*4.f, 1x*61); "
            s += "cn=CNSol[6.73; sol=(2x*1,1x*61); "
            s += "cn=3333333333333322222_2|333333333333334444444]"
            assert s in log

        out = subprocess.check_output(f"tar tzf {tmp.name}", shell=True).decode("utf-8")
        out = "\n".join(sorted(out.strip().split("\n")))
        expected = td(
            """
            ./
            ./NA10860.cn.lp
            ./NA10860.dump
            ./NA10860.json
            ./NA10860.log
            ./NA10860.major0.lp
            ./NA10860.major1.lp
            ./NA10860.minor0.lp
            ./NA10860.minor1.lp
            ./NA10860.minor2.lp
            ./NA10860.minor3.lp
            ./NA10860.minor4.lp
            ./NA10860.minor5.lp"""
        ).strip()
        assert out == expected
