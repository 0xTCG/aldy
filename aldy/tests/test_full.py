# 786
# Aldy source: test_cn_real.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.

# flake8: noqa

import pytest  # noqa
import re
import subprocess
import platform
import datetime

from tempfile import NamedTemporaryFile as tmpfile

from aldy.__main__ import get_version, main
from aldy.common import script_path, log
from aldy.cn import LEFT_FUSION_PENALTY
from aldy.version import __version__


YEAR = datetime.datetime.now().year
HEADER = f"""
🐿  Aldy v{__version__} (Python {platform.python_version()} on {get_version()})
   (c) 2016-{datetime.datetime.now().year} Aldy Authors. All rights reserved.
   Free for non-commercial/academic use only.
""".strip()


def escape_ansi(line):
    """
    Inspired by
    https://www.tutorialspoint.com/How-can-I-remove-the-ANSI-escape-sequences-from-a-string-in-python
    """  # noqa
    return re.compile(r"(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]").sub("", line)


def assert_file(monkeypatch, file, solver, expected, params=None, warn=False):
    lines = []

    def log_info(*args):
        s = str.format(*args)
        lines.append(s)

    monkeypatch.setattr(log, "info", log_info)
    if warn:
        monkeypatch.setattr(log, "warn", log_info)

    args = {
        **{
            "--gene": "CYP2D6",
            "--profile": "illumina",
            "--threshold": "50",
            "--gap": "0",
            "--solver": solver,
            "--fusion-penalty": f"{LEFT_FUSION_PENALTY}",
            "--max-minor-solutions": "3",
        },
        **(params or {}),
    }
    main(["genotype", file] + [i for k, v in args.items() for i in [k, v]])
    expected = expected.strip()
    lines = escape_ansi("\n".join(lines)).strip()
    assert lines == expected


EXPECTED_NA10860 = f"""
{HEADER}
Genotyping sample NA10860.bam...
Potential CYP2D6 gene structures for NA10860:
   1: 2x*1,1x*36.ALDY (confidence: 100%)
   2: 2x*1,1x*61 (confidence: 100%)
Potential major CYP2D6 star-alleles for NA10860:
   1: 1x*1, 1x*4.021, 1x*4N.ALDY (confidence: 100%)
   2: 1x*4, 1x*4N.ALDY, 1x*139 (confidence: 100%)
   3: 1x*4.021, 1x*4J, 1x*61 (confidence: 100%)
   4: 1x*4.021, 1x*4J, 1x*83.ALDY (confidence: 100%)
   5: 1x*4.021, 1x*4M, 1x*36.ALDY (confidence: 100%)
Best CYP2D6 star-alleles for NA10860:
   1: *1 / *4.021 + *4N.ALDY (confidence=100%)
      Minor alleles: *1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713, *4.021, *4.1013 -rs28371738
CYP2D6 results:
  - *1 / *4.021 + *4N.ALDY
    Minor: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713] / [*4.021] + [*4.1013 -rs28371738]
    Legacy notation: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713] / [*4.021] + [*4N.ALDY -rs28371738]
"""


def test_NA10860(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, EXPECTED_NA10860)


def test_NA10860_hg38(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860_hg38.bam")
    assert_file(
        monkeypatch, file, solver, EXPECTED_NA10860.replace("NA10860", "NA10860_hg38")
    )


EXPECTED_NA10860_GAP = f"""
{HEADER}
Genotyping sample NA10860.bam...
Potential CYP2D6 gene structures for NA10860:
   1: 2x*1,1x*36.ALDY (confidence: 100%)
   2: 2x*1,1x*61 (confidence: 100%)
   3: 2x*1,1x*5,1x*36.ALDY (confidence: 94%)
   4: 2x*1,1x*5,1x*61 (confidence: 94%)
Potential major CYP2D6 star-alleles for NA10860:
   1: 1x*1, 1x*4.021, 1x*4N.ALDY (confidence: 100%)
   2: 1x*4, 1x*4N.ALDY, 1x*139 (confidence: 100%)
   3: 1x*4.021, 1x*4J, 1x*61 (confidence: 100%)
   4: 1x*4.021, 1x*4J, 1x*83.ALDY (confidence: 100%)
   5: 1x*4.021, 1x*4M, 1x*36.ALDY (confidence: 100%)
   6: 1x*1, 1x*4.021, 1x*4N.ALDY, 1x*5 (confidence: 96%)
   7: 1x*4, 1x*4N.ALDY, 1x*5, 1x*139 (confidence: 96%)
   8: 1x*4.021, 1x*4J, 1x*5, 1x*61 (confidence: 96%)
   9: 1x*4.021, 1x*4J, 1x*5, 1x*83.ALDY (confidence: 96%)
  10: 1x*4.021, 1x*4M, 1x*5, 1x*36.ALDY (confidence: 96%)
Best CYP2D6 star-alleles for NA10860:
   1: *1 / *4.021 + *4N.ALDY (confidence=100%)
      Minor alleles: *1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713, *4.021, *4.1013 -rs28371738
CYP2D6 results:
  - *1 / *4.021 + *4N.ALDY
    Minor: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713] / [*4.021] + [*4.1013 -rs28371738]
    Legacy notation: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713] / [*4.021] + [*4N.ALDY -rs28371738]
"""


def test_NA10860_gap(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, EXPECTED_NA10860_GAP, {"--gap": "0.1"})


EXPECTED_HARD = f"""
{HEADER}
Genotyping sample HARD.dump...
Potential CYP2D6 gene structures for HARD:
   1: 3x*1,1x*68,1x*80 (confidence: 100%)
Potential major CYP2D6 star-alleles for HARD:
   1: 1x*2, 1x*4, 1x*39, 1x*68, 1x*80#4K (confidence: 100%)
   2: 2x*2, 1x*4, 1x*68, 1x*80#4 (confidence: 100%)
Best CYP2D6 star-alleles for HARD:
   1: *2 + *2 / *68 + *4 + *80 (confidence=100%)
      Minor alleles: *2.001, *2.001, *4.001 +rs4987144 -rs28588594, *68.001 -rs28371699 -rs28588594, *80#4.001
CYP2D6 results:
  - *2 + *2 / *68 + *4 + *80
    Minor: [*2.001] + [*2.001] / [*68.001 -rs28371699 -rs28588594] + [*4.001 +rs4987144 -rs28588594] + [*80#4.001]
    Legacy notation: [*2A] + [*2A] / [*68 -rs28371699 -rs28588594] + [*4A +rs4987144 -rs28588594] + [*80#4.001]
"""


def test_hard(monkeypatch, solver):
    file = script_path("aldy.tests.resources/HARD.dump")
    assert_file(monkeypatch, file, solver, EXPECTED_HARD, {"--profile": "pgrnseq-v1"})


EXPECTED_HARD_FUSION = f"""
{HEADER}
Genotyping sample HARD.dump...
Potential CYP2D6 gene structures for HARD:
   1: 4x*1 (confidence: 100%)
Potential major CYP2D6 star-alleles for HARD:
   1: 1x*2, 1x*4, 1x*4K, 1x*39 (confidence: 100%)
   2: 2x*2, 1x*4, 1x*4C (confidence: 100%)
Best CYP2D6 star-alleles for HARD:
   1: *2 + *2 / *4 + *4C (confidence=100%)
      Minor alleles: *2.001 +rs28371738 +rs2004511 +rs2267447 +rs1080989, *2.001, *4.001 +rs4987144 -rs28588594, *4.011 +rs1985842 +rs28371702 +rs28371701 +rs28371699 +rs28735595
CYP2D6 results:
  - *2 + *2 / *4 + *4C
    Minor: [*2.001 +rs28371738 +rs2004511 +rs2267447 +rs1080989] + [*2.001] / [*4.001 +rs4987144 -rs28588594] + [*4.011 +rs1985842 +rs28371702 +rs28371701 +rs28371699 +rs28735595]
    Legacy notation: [*2A +rs28371738 +rs2004511 +rs2267447 +rs1080989] + [*2A] / [*4A +rs4987144 -rs28588594] + [*4L +rs1985842 +rs28371702 +rs28371701 +rs28371699 +rs28735595]
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
./NA10860.log
./NA10860.major0.lp
./NA10860.major1.lp
./NA10860.minor0.lp
./NA10860.minor1.lp
./NA10860.minor2.lp
./NA10860.minor3.lp
./NA10860.yml
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
            s = "    rs1058172    42523528.C>T    3267G>A    "
            s += "(cov=  21, cn= 0.9; impact=R365H)\n"
            s = "[major] status= optimal; opt= 1.48; solution= 1x*4.021, 1x*4J, 1x*61\n"
            assert s in log

        out = subprocess.check_output(f"tar tzf {tmp.name}", shell=True).decode("utf-8")
        out = "\n".join(sorted(out.strip().split("\n")))
        assert out == EXPECTED_NA10860_DEBUG_TAR.strip()


EXPECTED_NA10860_CN = f"""
{HEADER}
Genotyping sample NA10860.bam...
Potential CYP2D6 gene structures for NA10860:
   1: 2x*1 (confidence: 100%)
Potential major CYP2D6 star-alleles for NA10860:
   1: 1x*1, 1x*4.021 (confidence: 100%)
   2: 1x*4, 1x*139 (confidence: 100%)
Best CYP2D6 star-alleles for NA10860:
   1: *1 / *4.021 (confidence=100%)
      Minor alleles: *1.018 +rs28371713, *4.021
CYP2D6 results:
  - *1 / *4.021
    Minor: [*1.018 +rs28371713] / [*4.021]
    Legacy notation: [*1.018 +rs28371713] / [*4.021]
"""


def test_NA10860_cn(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, EXPECTED_NA10860_CN, {"--cn": "1,1"})


def test_NA10860_vcf_out(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    with tmpfile(suffix=".vcf", mode="w") as out:
        out.close()
        assert_file(monkeypatch, file, solver, EXPECTED_NA10860, {"--output": out.name})
        with open(script_path("aldy.tests.resources/NA10860.vcf.expected")) as f:
            expected = f.read()
        with open(out.name) as f:
            produced = f.read()
        assert produced == expected.replace("aldy-v2.2", f"aldy-v{__version__}")


EXPECTED_NA07000 = f"""
{HEADER}
Genotyping sample NA07000_SLCO1B1.vcf.gz...
WARNING: Cannot detect genome, defaulting to hg19.
WARNING: Using VCF file. Copy-number calling is not available.
Using VCF sample NA07000
Potential SLCO1B1 gene structures for NA07000_SLCO1B1.vcf:
   1: 2x*1 (confidence: 100%)
Potential major SLCO1B1 star-alleles for NA07000_SLCO1B1.vcf:
   1: 1x*1, 1x*15 & rs4149057, rs2291075 (confidence: 100%)
   2: 1x*1B, 1x*5 & rs4149057, rs2291075 (confidence: 100%)
WARNING: multiple optimal solutions found!
Potential SLCO1B1 star-alleles for NA07000_SLCO1B1.vcf:
   1: *1+rs4149057+rs2291075 / *15 (confidence=100%)
      Minor alleles: *1.001 +rs4149057 +rs2291075, *15.001
   2: *1B+rs4149057+rs2291075 / *5 (confidence=100%)
      Minor alleles: *1.002 +rs4149057 +rs2291075, *5.001
SLCO1B1 results:
  - *1+rs4149057+rs2291075 / *15
    Minor: [*1.001 +rs4149057 +rs2291075] / [*15.001]
    Legacy notation: [*1A +rs4149057 +rs2291075] / [*15]
  - *1B+rs4149057+rs2291075 / *5
    Minor: [*1.002 +rs4149057 +rs2291075] / [*5.001]
    Legacy notation: [*1B +rs4149057 +rs2291075] / [*5]
WARNING: mutations rs2291075, rs4149057 suggest presence of a novel major star-allele.
However, such alleles cannot be determined without phasing data.
Please provide --phase parameter for Aldy to accurately call novel major star-alleles.
The above-reported assignments of these mutations are random.
"""


def test_NA07000_vcf_in(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA07000_SLCO1B1.vcf.gz")
    assert_file(
        monkeypatch,
        file,
        solver,
        EXPECTED_NA07000,
        {"--max-minor-solutions": "1", "--gene": "SLCO1B1"},
        warn=True,
    )


EXPECTED_NA07000_PHASE = f"""
{HEADER}
Genotyping sample NA07000_SLCO1B1.vcf.gz...
WARNING: Cannot detect genome, defaulting to hg19.
WARNING: Using VCF file. Copy-number calling is not available.
Using VCF sample NA07000
Potential SLCO1B1 gene structures for NA07000_SLCO1B1.vcf:
   1: 2x*1 (confidence: 100%)
Potential major SLCO1B1 star-alleles for NA07000_SLCO1B1.vcf:
   1: 1x*1, 1x*15 & rs4149057, rs2291075 (confidence: 100%)
   2: 1x*1B, 1x*5 & rs4149057, rs2291075 (confidence: 100%)
Using phasing information
Using phasing information
Best SLCO1B1 star-alleles for NA07000_SLCO1B1.vcf:
   1: *1+rs4149057 / *15+rs2291075 (confidence=100%)
      Minor alleles: *1.001 +rs4149057, *15.001 +rs2291075
SLCO1B1 results:
  - *1+rs4149057 / *15+rs2291075
    Minor: [*1.001 +rs4149057] / [*15.001 +rs2291075]
    Legacy notation: [*1A +rs4149057] / [*15 +rs2291075]
"""


def test_NA07000_phase(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA07000_SLCO1B1.vcf.gz")
    assert_file(
        monkeypatch,
        file,
        solver,
        EXPECTED_NA07000_PHASE,
        {
            "--max-minor-solutions": "1",
            "--gene": "SLCO1B1",
            "--phase": script_path("aldy.tests.resources/NA07000_SLCO1B1.phase"),
        },
        warn=True,
    )


EXPECTED_INS = f"""
{HEADER}
Genotyping sample INS.dump...
Potential CYP2D6 gene structures for INS:
   1: 2x*1 (confidence: 100%)
Potential major CYP2D6 star-alleles for INS:
   1: 1x*1, 1x*40 (confidence: 100%)
Best CYP2D6 star-alleles for INS:
   1: *1 / *40 (confidence=100%)
      Minor alleles: *1.006, *40.001
CYP2D6 results:
  - *1 / *40
    Minor: [*1.006] / [*40.001]
    Legacy notation: [*1.006] / [*40]
"""


def test_fix_insertions(monkeypatch, solver):
    file = script_path("aldy.tests.resources/INS.dump")
    assert_file(
        monkeypatch,
        file,
        solver,
        EXPECTED_INS,
        {"--profile": "pgrnseq-v3", "--max-minor-solutions": "1"},
    )


EXPECTED_PROFILE = f"""
{HEADER}
Scanning 1:60355979-98392615...
Scanning 6:18125541-18161374...
Scanning 7:1016834-99467173...
Scanning 10:96516437-135355620...
Scanning 11:14896554-14919751...
Scanning 12:21278127-21395730...
Scanning 13:48605702-48629878...
Scanning 15:75008882-75051948...
Scanning 19:15985833-41716444...
Scanning 22:42518175-42549249...
Scanning X:153756605-153781787...
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
