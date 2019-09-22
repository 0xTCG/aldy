# 786
# Aldy source: test_cn_real.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa
import re
import subprocess
import platform
import sys
import datetime

from tempfile import NamedTemporaryFile as tmpfile

from aldy.__main__ import main
from aldy.common import script_path, log
from aldy.cn import LEFT_FUSION_PENALTY
from aldy.version import __version__


YEAR = datetime.datetime.now().year
HEADER = f"""
*** Aldy v{__version__} (Python {platform.python_version()}, {sys.platform}) ***
*** (c) 2016-{YEAR} Aldy Authors & Indiana University Bloomington. All rights reserved.
*** Free for non-commercial/academic use only.
""".strip()


def escape_ansi(line):
    """
    Inspired by 
    https://www.tutorialspoint.com/How-can-I-remove-the-ANSI-escape-sequences-from-a-string-in-python
    """  # noqa
    return re.compile(r"(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]").sub("", line)


def assert_file(monkeypatch, file, solver, expected, params=None, output=None):
    lines = []

    def log_info(*args):
        s = str.format(*args)
        lines.append(s)

    monkeypatch.setattr(log, "info", log_info)

    args = {
        **{
            "--gene": "CYP2D6",
            "--profile": "illumina",
            "--threshold": "50",
            "--gap": "0",
            "--solver": solver,
            "--fusion-penalty": f"{LEFT_FUSION_PENALTY}",
            "--max-minor-solutions": "10",
        },
        **(params or {}),
    }
    main(["genotype", file] + [i for k, v in args.items() for i in [k, v]])
    expected = expected.strip()
    lines = escape_ansi("\n".join(lines)).strip()
    assert lines == expected


def assert_show(monkeypatch, expected, params=None):
    lines = []

    def log_info(*args):
        s = str.format(*args)
        lines.append(s)

    monkeypatch.setattr(log, "info", log_info)

    args = ["show", "--gene", script_path("aldy.tests.resources/toy.yml")]
    main(args + (params or []))
    expected = expected.strip()
    lines = escape_ansi("\n".join(lines)).strip()
    assert lines == expected


EXPECTED_NA10860 = f"""
{HEADER}
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
"""


def test_NA10860(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, EXPECTED_NA10860)


EXPECTED_NA10860_GAP = f"""
{HEADER}
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
  *1/*4+*4                       (*1, *4, *4.i)
"""


def test_NA10860_gap(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, EXPECTED_NA10860_GAP, {"--gap": "0.1"})


EXPECTED_HARD = f"""
{HEADER}
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
  *2+*2/*68+*4+*80               (*2, *2, *4, *68, *80/4)
"""


def test_hard(monkeypatch, solver):
    file = script_path("aldy.tests.resources/HARD.dump")
    assert_file(monkeypatch, file, solver, EXPECTED_HARD, {"--profile": "pgrnseq-v1"})


EXPECTED_HARD_FUSION = f"""
{HEADER}
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
  *2+*2/*4+*4                    (*2, *2, *4, *4.b)
"""


def test_hard_fusion(monkeypatch, solver):
    file = script_path("aldy.tests.resources/HARD.dump")
    assert_file(
        monkeypatch,
        file,
        solver,
        EXPECTED_HARD_FUSION,
        {"--profile": "pgrnseq-v1", "--fusion-penalty": "5"},
    )


EXPECTED_NA10860_DEBUG_TAR = """
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
./NA10860.minor5.lp
"""


def test_NA10860_debug(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    with tmpfile(suffix=".tar.gz") as tmp:
        with tmpfile(mode="w") as out, tmpfile(mode="w") as out_log:
            out.close()
            out_log.close()
            assert_file(
                monkeypatch,
                file,
                solver,
                EXPECTED_NA10860 + "Preparing debug archive...",
                {"--debug": tmp.name[:-7], "--log": out_log.name, "--output": out.name},
            )
            with open(script_path("aldy.tests.resources/NA10860.out.expected")) as f:
                expected = f.read()
            with open(out.name) as f:
                produced = f.read()
            assert produced == expected
            # Check logs
            with open(out_log.name) as f:
                log = f.read()
            s = "Major solver: MajorSol[1.33; sol=(1x*4, 1x*4.f, 1x*61); "
            s += "cn=CNSol[6.73; sol=(2x*1,1x*61); "
            s += "cn=3333333333333322222_2|333333333333334444444]"
            assert s in log

        out = subprocess.check_output(f"tar tzf {tmp.name}", shell=True).decode("utf-8")
        out = "\n".join(sorted(out.strip().split("\n")))
        assert out == EXPECTED_NA10860_DEBUG_TAR.strip()


EXPECTED_NA10860_CN = f"""
{HEADER}
Genotyping sample NA10860.bam...
Potential CYP2D6 copy number configurations for NA10860:
   1: 2x*1
      Confidence: 1.00 (score = 0.00)

Potential major CYP2D6 star-alleles for NA10860:
   1: 1x*1, 1x*4
      Confidence: 1.00 (score = 2.69)
   2: 1x*4.f, 1x*39
      Confidence: 1.00 (score = 2.69)
   3: 1x*4.h, 1x*10
      Confidence: 1.00 (score = 2.69)

Best CYP2D6 star-alleles for NA10860:
   1: *1/*4
      Minor: *1 +42528223:SNP.GA, *4AW
      Confidence: 1.00 (score = 11.88)
   2: *1/*4
      Minor: *1, *4AW +42528223:SNP.GA
      Confidence: 1.00 (score = 11.88)
CYP2D6 results:
  *1/*4                          (*1, *4)"""


def test_NA10860_cn(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, EXPECTED_NA10860_CN, {"--cn": "1,1"})


def test_NA10860_vcf(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    with tmpfile(suffix=".vcf", mode="w") as out:
        out.close()
        assert_file(monkeypatch, file, solver, EXPECTED_NA10860, {"--output": out.name})
        with open(script_path("aldy.tests.resources/NA10860.vcf.expected")) as f:
            expected = f.read()
        with open(out.name) as f:
            produced = f.read()
        assert produced == expected.replace("aldy-v2.2", f"aldy-v{__version__}")


EXPECTED_INS = f"""
{HEADER}
Genotyping sample INS.dump...
Potential CYP2D6 copy number configurations for INS:
   1: 2x*1
      Confidence: 1.00 (score = 4.33)

Potential major CYP2D6 star-alleles for INS:
   1: 1x*1, 1x*40
      Confidence: 1.00 (score = 0.51)

Best CYP2D6 star-alleles for INS:
   1: *1/*40
      Minor: *1, *40
      Confidence: 1.00 (score = 0.76)
CYP2D6 results:
  *1/*40                         (*1, *40)"""


def test_fix_insertions(monkeypatch, solver):
    file = script_path("aldy.tests.resources/INS.dump")
    assert_file(
        monkeypatch,
        file,
        solver,
        EXPECTED_INS,
        {"--profile": "pgrnseq-v1", "--max-minor-solutions": "1"},
    )


EXPECTED_PROFILE = f"""
{HEADER}
Generating profile for DPYD (1:97541297-98388616)
Generating profile for CYP2C19 (10:96444999-96615001)
Generating profile for CYP2C9 (10:96690999-96754001)
Generating profile for CYP2C8 (10:96795999-96830001)
Generating profile for CYP4F2 (19:15618999-16009501)
Generating profile for CYP2A6 (19:41347499-41400001)
Generating profile for CYP2D6 (22:42518899-42553001)
Generating profile for TPMT (6:18126540-18157375)
Generating profile for CYP3A5 (7:99244999-99278001)
Generating profile for CYP3A4 (7:99353999-99465001)
"""


def test_profile(monkeypatch, capsys):
    lines = []

    def log_info(*args):
        s = str.format(*args)
        lines.append(s)

    monkeypatch.setattr(log, "info", log_info)

    main(["profile", script_path("aldy.tests.resources/NA10860.bam")])
    lines = escape_ansi("\n".join(lines)).strip()
    assert lines == EXPECTED_PROFILE.strip()

    captured = capsys.readouterr()
    with open(script_path("aldy.tests.resources/NA10860.profile")) as f:
        expected = f.read()
    assert captured.out == expected


EXPECTED_SHOW = f"""
{HEADER}
Gene TOY
  Reference genome locus: 1:0-200
  Pseudogenes: PSEUDO (ID 1)
  Aminoacid:
    VRTCTYVYVR
  Genic regions:
    Region                        Gene                Pseudogene
     0 tmp :                 1:100-110                    1:0-10
     1 e   :                 1:110-120                   1:10-20
     1 i   :                 1:120-130                   1:20-30
     2 e   :                 1:130-140                   1:30-40
     2 i   :                 1:140-150                   1:40-50
     3 e   :                 1:150-160                   1:50-60
  Copy number configurations:
    1, 4, 5, 6
  Major star-alleles:
    *1, *1.a, *2, *3, *4/1, *4/3, *5, *6
"""


def test_show(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW)


EXPECTED_SHOW_CN = f"""
{HEADER}
Gene TOY
  Copy number configuration 1:
  Normal allelic configuration: all regions present in both gene and pseudogene
  Major star-alleles:
    *1, *1.a, *2, *3
  Copy number profile:
               Gene Pseudo
     0 tmp :      1      1
     1 e   :      1      1
     1 i   :      1      1
     2 e   :      1      1
     2 i   :      1      1
     3 e   :      1      1
"""


def test_show_cn(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW_CN, ["--cn", "1"])


EXPECTED_SHOW_MAJOR = f"""
{HEADER}
Gene TOY
Major star-allele 2:
  Copy number configuration 1
  Functional mutations:
    1:111       DEL.AC              
    1:119       INS.TT              
  Minor star-alleles:
    *2
"""  # noqa


def test_show_major(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW_MAJOR, ["--major", "2"])


EXPECTED_SHOW_MINOR = f"""
{HEADER}
Gene TOY

Minor star-allele 1C:
  Major star-allele: 1.a:
  Copy number configuration: 1
  Functional mutations:
    1:105       SNP.TA              
  Silent mutations:
"""  # noqa


def test_show_minor(monkeypatch):
    assert_show(monkeypatch, EXPECTED_SHOW_MINOR, ["--minor", "1C"])
    with pytest.raises(SystemExit):
        assert_show(monkeypatch, EXPECTED_SHOW_MINOR, ["--minor", "1c"])
