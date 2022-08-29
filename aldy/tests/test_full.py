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
        nonlocal lines

        s = str.format(*args)
        lines += s.split("\n")  # type: ignore

    monkeypatch.setattr(log, "info", log_info)
    if warn:
        monkeypatch.setattr(log, "warn", log_info)

    args = {
        **{
            "--gene": "cyp2d6",
            "--profile": "illumina",
            "--param": "threshold=0.5",
            "--param": "gap=0",
            "--param": "max-minor-solutions=3",
            "--param": "minor_phase_vars=10",
            "--solver": solver,
        },
        **(params or {}),
    }
    main(["genotype", file] + [i for k, v in args.items() for i in [k, v]])
    expected = "\n".join(e.strip() for e in expected.strip().split("\n"))
    lines = "\n".join(escape_ansi(l).strip() for l in lines).strip()
    assert lines == expected


EXPECTED_NA10860 = f"""
{HEADER}
Genotyping sample NA10860.bam...
Potential CYP2D6 gene structures for NA10860:
   1: 2x*1,1x*141.1001 (confidence: 100%)
   2: 2x*1,1x*61 (confidence: 100%)
Potential major CYP2D6 star-alleles for NA10860:
   1: 1x*1, 1x*4, 1x*4.021.ALDY (confidence: 100%)
   2: 1x*1, 1x*4.021, 1x*4N.ALDY (confidence: 100%)
   3: 1x*4, 1x*4N.ALDY, 1x*139 (confidence: 100%)
   4: 1x*4.021, 1x*4J, 1x*61 (confidence: 100%)
   5: 1x*4.021, 1x*4J, 1x*83.ALDY (confidence: 100%)
   6: 1x*4.021, 1x*4M, 1x*36.ALDY (confidence: 100%)
   7: 1x*4.021.ALDY, 1x*4J, 1x*39 (confidence: 100%)
   8: 1x*4.021.ALDY, 1x*4M, 1x*10 (confidence: 100%)
   9: 1x*4.021.ALDY_2, 1x*4N.ALDY, 1x*74 (confidence: 100%)
Potential CYP2D6 star-alleles for NA10860:
   1: *1 / *4 + *4.021.ALDY (confidence=100%)
      Minor alleles: *(1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701), *4.001, *(4.1021 -rs28371738)
   2: *1 / *4.021 + *4N.ALDY (confidence=100%)
      Minor alleles: *(1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701), *4.021, *(4.1013 -rs28371738)
CYP2D6 results:
  - *1 / *4 + *4.021.ALDY
    Minor: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4.001] + [*4.1021 -rs28371738]
    Legacy notation: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4A] + [*4.021.ALDY -rs28371738]
    Estimated activity for *1: normal function (evidence: D); see https://www.pharmvar.org/haplotype/662 for details
    Estimated activity for *4: no function (evidence: D); see https://www.pharmvar.org/haplotype/235 for details
    Estimated activity for *4.021.ALDY: no function (evidence: D)
  - *1 / *4.021 + *4N.ALDY
    Minor: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4.021] + [*4.1013 -rs28371738]
    Legacy notation: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4.021] + [*4N.ALDY -rs28371738]
    Estimated activity for *1: normal function (evidence: D); see https://www.pharmvar.org/haplotype/662 for details
    Estimated activity for *4.021: no function (evidence: D); see https://www.pharmvar.org/haplotype/652 for details
    Estimated activity for *4N.ALDY: no function (evidence: D)
"""


def test_NA10860(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, EXPECTED_NA10860)


def test_NA10860_debug(monkeypatch, solver):
    expected_tar = """
    ./
    ./NA10860.CYP2D6.cn.lp
    ./NA10860.CYP2D6.dump
    ./NA10860.CYP2D6.genome
    ./NA10860.CYP2D6.major0.lp
    ./NA10860.CYP2D6.major1.lp
    ./NA10860.log
    ./NA10860.yml
    """
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
            s += "(cov=  17, cn= 0.8; region=e7; impact=R365H; )\n"
            s = "[major] status= optimal; opt= 1.47; solution= 1x*1, 1x*4, 1x*4.021"
            assert s in log
        out = subprocess.check_output(f"tar tzf {tmp.name}", shell=True).decode("utf-8")
        out = "\n".join(i.strip() for i in sorted(out.strip().split("\n")))
        assert out == "\n".join(i.strip() for i in expected_tar.strip().split("\n"))


def test_NA10860_hg38(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860_hg38.bam")
    assert_file(
        monkeypatch, file, solver, EXPECTED_NA10860.replace("NA10860", "NA10860_hg38")
    )


def test_NA10860_gap(monkeypatch, solver):
    expected = f"""
    {HEADER}
    Genotyping sample NA10860.bam...
    Potential CYP2D6 gene structures for NA10860:
    1: 2x*1,1x*141.1001 (confidence: 100%)
    2: 2x*1,1x*61 (confidence: 100%)
    3: 1x*1,1x*61,1x*141.1001 (confidence: 90%)
    4: 1x*1,2x*141.1001 (confidence: 90%)
    5: 1x*1,2x*61 (confidence: 90%)
    6: 3x*1 (confidence: 90%)
    Potential major CYP2D6 star-alleles for NA10860:
    1: 1x*1, 1x*4, 1x*4.021.ALDY (confidence: 100%)
    2: 1x*1, 1x*4.021, 1x*4N.ALDY (confidence: 100%)
    3: 1x*4, 1x*4N.ALDY, 1x*139 (confidence: 100%)
    4: 1x*4.021, 1x*4J, 1x*61 (confidence: 100%)
    5: 1x*4.021, 1x*4J, 1x*83.ALDY (confidence: 100%)
    6: 1x*4.021, 1x*4M, 1x*36.ALDY (confidence: 100%)
    7: 1x*4.021.ALDY, 1x*4J, 1x*39 (confidence: 100%)
    8: 1x*4.021.ALDY, 1x*4M, 1x*10 (confidence: 100%)
    9: 1x*4.021.ALDY_2, 1x*4N.ALDY, 1x*74 (confidence: 100%)
    Potential CYP2D6 star-alleles for NA10860:
    1: *1 / *4 + *4.021.ALDY (confidence=100%)
        Minor alleles: *(1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701), *4.001, *(4.1021 -rs28371738)
    2: *1 / *4.021 + *4N.ALDY (confidence=100%)
        Minor alleles: *(1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701), *4.021, *(4.1013 -rs28371738)
    CYP2D6 results:
    - *1 / *4 + *4.021.ALDY
        Minor: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4.001] + [*4.1021 -rs28371738]
        Legacy notation: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4A] + [*4.021.ALDY -rs28371738]
        Estimated activity for *1: normal function (evidence: D); see https://www.pharmvar.org/haplotype/662 for details
        Estimated activity for *4: no function (evidence: D); see https://www.pharmvar.org/haplotype/235 for details
        Estimated activity for *4.021.ALDY: no function (evidence: D)
    - *1 / *4.021 + *4N.ALDY
        Minor: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4.021] + [*4.1013 -rs28371738]
        Legacy notation: [*1.018 +rs111564371 +rs112568578 +rs113889384 +rs28371713 +rs28371701] / [*4.021] + [*4N.ALDY -rs28371738]
        Estimated activity for *1: normal function (evidence: D); see https://www.pharmvar.org/haplotype/662 for details
        Estimated activity for *4.021: no function (evidence: D); see https://www.pharmvar.org/haplotype/652 for details
        Estimated activity for *4N.ALDY: no function (evidence: D)
    """
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, expected, {"--param": "gap=0.3"})


def test_NA10860_cn(monkeypatch, solver):
    expected = f"""
    {HEADER}
    Genotyping sample NA10860.bam...
    Potential CYP2D6 gene structures for NA10860:
    1: 2x*1 (confidence: 100%)
    Potential major CYP2D6 star-alleles for NA10860:
    1: 1x*1, 1x*4.021 (confidence: 100%)
    2: 1x*4, 1x*139 (confidence: 100%)
    3: 1x*4.021.ALDY_2, 1x*74 (confidence: 100%)
    Best CYP2D6 star-alleles for NA10860:
    1: *1 / *4.021 (confidence=100%)
        Minor alleles: *(1.018 +rs113889384 +rs28371713), *4.021
    CYP2D6 results:
    - *1 / *4.021
        Minor: [*1.018 +rs113889384 +rs28371713] / [*4.021]
        Legacy notation: [*1.018 +rs113889384 +rs28371713] / [*4.021]
        Estimated activity for *1: normal function (evidence: D); see https://www.pharmvar.org/haplotype/662 for details
        Estimated activity for *4.021: no function (evidence: D); see https://www.pharmvar.org/haplotype/652 for details
    """
    file = script_path("aldy.tests.resources/NA10860.bam")
    assert_file(monkeypatch, file, solver, expected, {"--cn": "1,1"})


def test_NA10860_vcf_out(monkeypatch, solver):
    file = script_path("aldy.tests.resources/NA10860.bam")
    with tmpfile(suffix=".vcf", mode="w") as out:
        out.close()
        assert_file(monkeypatch, file, solver, EXPECTED_NA10860, {"--output": out.name})
        with open(script_path("aldy.tests.resources/NA10860.vcf.expected")) as f:
            expected = f.read()
        with open(out.name) as f:
            produced = f.read()
        assert produced == expected.replace("aldy-v4.0", f"aldy-v{__version__}")


def test_fusion(monkeypatch, solver):
    # NA19785 pgrnseq-v1
    expected = f"""
    {HEADER}
    Genotyping sample HARD.dump.tar.gz...
    Potential CYP2D6 gene structures for NA19785:
    1: 2x*1,1x*79 (confidence: 100%)
    Potential major CYP2D6 star-alleles for NA19785:
    1: 1x*1, 1x*2, 1x*79#2 (confidence: 100%)
    2: 1x*2, 1x*34, 1x*79#10 (confidence: 100%)
    3: 2x*2, 1x*79#1 (confidence: 100%)
    Potential CYP2D6 star-alleles for NA19785:
    1: *1 / *79 + *2 (confidence=100%)
        Minor alleles: *1, *2, *79#11
    2: *2 / *79 + *2 (confidence=100%)
        Minor alleles: *2, *2, *79#1
    3: *34 / *79 + *2 (confidence=100%)
        Minor alleles: *2, *34, *79#10
    CYP2D6 results:
    - *1 / *79 + *2
        Minor: [*1] / [*79#11] + [*2]
        Legacy notation: [*1] / [*79#11] + [*2]
        Estimated activity for *1: unknown
        Estimated activity for *2: unknown
        Estimated activity for *79#2: unknown
    - *2 / *79 + *2
        Minor: [*2] / [*79#1] + [*2]
        Legacy notation: [*2] / [*79#1] + [*2]
        Estimated activity for *2: unknown
        Estimated activity for *79#1: unknown
        Estimated activity for *2: unknown
    - *34 / *79 + *2
        Minor: [*34] / [*79#10] + [*2]
        Legacy notation: [*34] / [*79#10] + [*2]
        Estimated activity for *2: unknown
        Estimated activity for *34: unknown
        Estimated activity for *79#10: unknown
    """
    file = script_path("aldy.tests.resources/HARD.dump.tar.gz")
    assert_file(monkeypatch, file, solver, expected, {"--gene": "pharmacoscan/cyp2d6"})


def test_fusion_off(monkeypatch, solver):
    expected = f"""
    {HEADER}
    Genotyping sample HARD.dump.tar.gz...
    Potential CYP2D6 gene structures for NA19785:
    1: 3x*1 (confidence: 100%)
    Potential major CYP2D6 star-alleles for NA19785:
    1: 1x*1, 2x*2 (confidence: 100%)
    Best CYP2D6 star-alleles for NA19785:
    1: *1 / *2 + *2 (confidence=100%)
        Minor alleles: *1, *2, *2
    CYP2D6 results:
    - *1 / *2 + *2
        Minor: [*1] / [*2] + [*2]
        Legacy notation: [*1] / [*2] + [*2]
        Estimated activity for *1: unknown
        Estimated activity for *2: unknown
        Estimated activity for *2: unknown
    """
    file = script_path("aldy.tests.resources/HARD.dump.tar.gz")
    assert_file(
        monkeypatch,
        file,
        solver,
        expected,
        {"--param": "cn-fusion-left=10", "--gene": "pharmacoscan/cyp2d6"},
    )


def test_NA07000_vcf_in(monkeypatch, solver):
    expected = f"""
    {HEADER}
    Genotyping sample NA07000_SLCO1B1.vcf.gz...
    WARNING: Cannot detect genome, defaulting to hg19.
    WARNING: Using VCF file. Copy-number calling is not available.
    Using VCF sample NA07000
    Potential SLCO1B1 gene structures for NA07000:
    1: 2x*1 (confidence: 100%)
    Potential major SLCO1B1 star-alleles for NA07000:
    1: 1x*1, 1x*15 (confidence: 100%)
    2: 1x*5, 1x*37 (confidence: 100%)
    Best SLCO1B1 star-alleles for NA07000:
    1: *1 / *15 (confidence=100%)
        Minor alleles: *1.003, *15.001
    SLCO1B1 results:
    - *1 / *15
        Minor: [*1.003] / [*15.001]
        Legacy notation: [*1.003] / [*15A, SLCO1B1*15B, SLCO1B1*17]
        Estimated activity for *1: normal function (evidence: D); see https://www.pharmvar.org/haplotype/1713 for details
        Estimated activity for *15: function not assigned (evidence: D); see https://www.pharmvar.org/haplotype/1704 for details
    """
    file = script_path("aldy.tests.resources/NA07000_SLCO1B1.vcf.gz")
    assert_file(
        monkeypatch,
        file,
        solver,
        expected,
        {"--param": "max-minor-solutions=1", "--gene": "SLCO1B1"},
        warn=True,
    )


def test_fix_insertions(monkeypatch, solver):
    expected = f"""
    {HEADER}
    Genotyping sample INS.dump.tar.gz...
    Potential CYP2D6 gene structures for NA17102:
    1: 2x*1 (confidence: 100%)
    Potential major CYP2D6 star-alleles for NA17102:
    1: 1x*1, 1x*40 (confidence: 100%)
    Best CYP2D6 star-alleles for NA17102:
    1: *1 / *40 (confidence=100%)
        Minor alleles: *1, *40
    CYP2D6 results:
    - *1 / *40
        Minor: [*1] / [*40]
        Legacy notation: [*1] / [*40]
        Estimated activity for *1: unknown
        Estimated activity for *40: unknown
    """
    file = script_path("aldy.tests.resources/INS.dump.tar.gz")
    assert_file(
        monkeypatch,
        file,
        solver,
        expected,
        {"--gene": "pharmacoscan/cyp2d6", "--param": "max-minor-solutions=1"},
    )


def test_profile(monkeypatch, capsys):
    lines = []

    def log_info(*args):
        nonlocal lines

        s = str.format(*args)
        lines += s.split("\n")  # type: ignore

    monkeypatch.setattr(log, "info", log_info)

    expected = f"""
    {HEADER}
    Scanning 1:60355979-110239367...
    Scanning 2:234492389-234684945...
    Scanning 4:69916092-69980005...
    Scanning 6:18125541-18161374...
    Scanning 7:1016834-117359025...
    Scanning 8:18021970-18261723...
    Scanning 10:96516437-135355620...
    Scanning 11:14896554-67357124...
    Scanning 12:21278127-21395730...
    Scanning 13:48605702-48629878...
    Scanning 15:75008882-75051948...
    Scanning 16:31099174-31112276...
    Scanning 19:15985833-41716444...
    Scanning 22:19923262-42549249...
    Scanning X:153756605-153781787...
    """
    main(["profile", script_path("aldy.tests.resources/NA10860.bam")])
    lines = "\n".join(escape_ansi(l).strip() for l in lines).strip()
    expected = "\n".join(e.strip() for e in expected.strip().split("\n"))
    assert lines == expected

    captured = capsys.readouterr()
    with open(script_path("aldy.tests.resources/NA10860.profile")) as f:
        expected = f.read()
    assert captured.out == expected

    lines = []
    expected = f"""
    {HEADER}
    Scanning 1:59890307-109696745...
    Scanning 2:233583743-233776299...
    Scanning 4:69050374-69114287...
    Scanning 6:18125310-18161143...
    Scanning 7:977198-117718971...
    Scanning 8:18164461-18404213...
    Scanning 10:94756680-133542116...
    Scanning 11:14875008-67589653...
    Scanning 12:21125193-21242796...
    Scanning 13:48031566-48055742...
    Scanning 15:74716541-74759607...
    Scanning 16:31087853-31100955...
    Scanning 19:15875023-41210539...
    Scanning 22:19935739-42153258...
    Scanning X:154528390-154553572...
    """
    main(
        ["profile", script_path("aldy.tests.resources/NA10860.bam"), "--genome", "hg38"]
    )
    lines = "\n".join(escape_ansi(l).strip() for l in lines).strip()
    expected = "\n".join(e.strip() for e in expected.strip().split("\n"))
    assert lines == expected

    captured = capsys.readouterr()
    with open(script_path("aldy.tests.resources/NA10860.profile.hg38")) as f:
        expected = f.read()
    assert captured.out == expected


def test_pacbio(monkeypatch, solver):
    expected = f"""
    {HEADER}
    Genotyping sample HG03166.pb.bam...
    Potential CYP2D6 gene structures for HG03166:
    1: 2x*1 (confidence: 100%)
    Potential major CYP2D6 star-alleles for HG03166:
    1: 1x*2, 1x*40 (confidence: 100%)
    Best CYP2D6 star-alleles for HG03166:
    1: *2 / *40 (confidence=100%)
        Minor alleles: *2.023, *40.001
    CYP2D6 results:
    - *2 / *40
        Minor: [*2.023] / [*40.001]
        Legacy notation: [*2.023] / [*40]
        Estimated activity for *2: normal function (evidence: D); see https://www.pharmvar.org/haplotype/1370 for details
        Estimated activity for *40: no function (evidence: D); see https://www.pharmvar.org/haplotype/231 for details
    """
    file = script_path("aldy.tests.resources/HG03166.pb.bam")
    assert_file(
        monkeypatch, file, solver, expected, {"--profile": "pacbio-hifi-targeted"}
    )
