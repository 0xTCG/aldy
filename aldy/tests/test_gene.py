# 786
# Aldy source: test_gene.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa

from aldy.gene import Mutation, CNConfig, CNConfigType, MajorAllele, MinorAllele
from aldy.common import GRange


def test_gene_basic(toy_gene):
    assert toy_gene.name == "TOY"
    assert toy_gene.version == "test-2.0 (2020-11-17)"
    assert toy_gene.pharmvar is None
    assert toy_gene.ensembl is None

    assert toy_gene.refseq == "NG_TEST"
    assert toy_gene.seq == "ACGT" * 50, "Sequence"


def test_gene_regions(toy_gene):
    assert toy_gene.genome == "hg19"
    assert toy_gene.chr == "20"
    assert toy_gene.strand == 1
    assert toy_gene.chr_to_ref == {100_000_000 + i: i for i in range(200)}
    assert toy_gene.ref_to_chr == {i: 100_000_000 + i for i in range(200)}

    # Test gene regions
    assert toy_gene.regions == [
        {
            "tmp": GRange("20", 100_000_100, 100_000_110),
            "e1": GRange("20", 100_000_110, 100_000_120),
            "i1": GRange("20", 100_000_120, 100_000_130),
            "e2": GRange("20", 100_000_130, 100_000_140),
            "i2": GRange("20", 100_000_140, 100_000_150),
            "e3": GRange("20", 100_000_150, 100_000_160),
            "down": GRange("20", 100_000_160, 100_000_200),
        },
        {
            "tmp": GRange("20", 100_000_000, 100_000_010),
            "e1": GRange("20", 100_000_010, 100_000_020),
            "i1": GRange("20", 100_000_020, 100_000_030),
            "e2": GRange("20", 100_000_030, 100_000_040),
            "i2": GRange("20", 100_000_040, 100_000_050),
            "e3": GRange("20", 100_000_050, 100_000_060),
            "down": GRange("20", 100_000_060, 100_000_100),
        },
    ], "Regions"
    assert toy_gene.get_wide_region() == ("20", 100_000_000, 100_000_200)
    assert toy_gene.pseudogenes == ["TOYP"], "Pseudogenes"

    coding_seq = "".join(toy_gene.seq[s:e] for [s, e] in toy_gene.exons)
    assert coding_seq == "GTACGTACGTGTACGTACGTGTACGTACGT", "Coding sequence"
    assert toy_gene.aminoacid == "VRTCTYVYVR", "Protein"

    assert toy_gene[100_000_040:100_000_050] == "ACGTACGTAC"
    assert toy_gene[100_000_050] == "G"
    assert toy_gene[100_000_196:100_000_210] == "ACGT" + "N" * 10


def test_gene_alleles(toy_gene):
    def cnify(s):
        s = s.replace(" ", "").replace("|", ",").split(",")
        assert len(s) == len(toy_gene.regions[0]) + len(toy_gene.regions[1])
        i = 0
        cn = [{}, {}]
        for g in range(2):
            for r in toy_gene.regions[g]:
                cn[g][r] = float(s[i])
                i += 1
        return cn

    # Test copy number configurations
    assert toy_gene.do_copy_number
    assert toy_gene.cn_configs == {
        "1": CNConfig(
            cnify("1,  1,1,  1,1,  1,1 | 1,  1,1,  1,1,  1,1"),
            kind=CNConfigType.DEFAULT,
            alleles=set(["1", "1C", "2", "3"]),
            description="Standard copy-number configuration",
        ),
        "4": CNConfig(
            cnify("0,  0,0,  0,1,  1,1 | 1,  1,1,  1,0,  0,0"),
            kind=CNConfigType.LEFT_FUSION,
            alleles=set(["4#1", "4#3"]),
            description="TOYP fusion until i2",
        ),
        "5": CNConfig(
            cnify("1,  1,1,  0,0,  0,0 | 1,  1,1,  2,2,  2,2"),
            kind=CNConfigType.RIGHT_FUSION,
            alleles=set(["5"]),
            description="TOYP conservation after e2",
        ),
        "6": CNConfig(
            cnify("0,  0,0,  0,0,  0,0 | 1,  1,1,  1,1,  1,1"),
            kind=CNConfigType.DELETION,
            alleles=set(["6"]),
            description="TOY deletion",
        ),
    }, "CN configs"
    assert toy_gene.unique_regions == ["e1", "i1", "e2", "i2", "e3"], "Unique regions"

    # Test allele configurations
    assert sorted(toy_gene.alleles.keys()) == [
        "1",
        "1C",
        "2",
        "3",
        "4#1",
        "4#3",
        "5",
        "6",
    ], "Major alleles"
    assert toy_gene.common_tandems == [("1", "4")], "Common tandems"
    assert toy_gene.alleles["1"] == MajorAllele(
        "1",
        cn_config="1",
        minors={
            "1.001": MinorAllele(
                "1.001", "1", activity="normal function", evidence="D"
            ),
            "1.002": MinorAllele(
                "1.002", "1B", neutral_muts=set([Mutation(100_000_114, "T>A")])
            ),
        },
    )
    assert toy_gene.alleles["1C"] == MajorAllele(
        "1C",
        cn_config="1",
        func_muts=set([Mutation(100_000_104, "T>A")]),
        minors={
            "1.003": MinorAllele(
                "1.003", "1C", pharmvar="https://www.pharmvar.org/haplotype/NA"
            )
        },
    )
    assert toy_gene.alleles["2"] == MajorAllele(
        "2",
        cn_config="1",
        func_muts=set([Mutation(100_000_110, "delAC"), Mutation(100_000_118, "insTT")]),
        minors={"2.001": MinorAllele("2.001", "2")},
    )
    assert toy_gene.alleles["3"] == MajorAllele(
        "3",
        cn_config="1",
        func_muts=set([Mutation(100_000_150, "C>T")]),
        minors={
            "3.001": MinorAllele(
                "3.001", "3", neutral_muts=set([Mutation(100_000_147, "insA")])
            )
        },
    )
    assert toy_gene.alleles["4#1"] == MajorAllele(
        "4#1",
        cn_config="4",
        minors={"4#1.001": MinorAllele("4#1.001")},
    )
    assert toy_gene.alleles["4#3"] == MajorAllele(
        "4#3",
        cn_config="4",
        func_muts=set([Mutation(100_000_150, "C>T")]),
        minors={
            "4#3.001": MinorAllele(
                "4#3.001", neutral_muts=set([Mutation(100_000_147, "insA")])
            )
        },
    )
    assert toy_gene.alleles["5"] == MajorAllele(
        "5",
        cn_config="5",
        func_muts=set([Mutation(100_000_110, "delAC")]),
        minors={"5.001": MinorAllele("5.001", "5")},
    )
    assert toy_gene.alleles["6"] == MajorAllele(
        "6",
        cn_config="6",
        minors={"6.001": MinorAllele("6.001", "6DEL")},
    )
    assert toy_gene.deletion_allele() == "6"

    # Test mutations
    assert toy_gene.mutations == {
        (100_000_000 + 105 - 1, "T>A"): ("functional", "-", 105 - 1, 105 - 1, "T>A"),
        (100_000_000 + 111 - 1, "delAC"): (
            "frameshift",
            "-",
            111 - 1,
            111 - 1,
            "delAC",
        ),
        (100_000_000 + 115 - 1, "T>A"): (None, "rs28371732", 115 - 1, 115 - 1, "T>A"),
        (100_000_000 + 119 - 1, "insTT"): (
            "frameshift",
            "-",
            119 - 1,
            119 - 1,
            "insTT",
        ),
        (100_000_000 + 148 - 1, "insA"): (None, "-", 148 - 1, 148 - 1, "insA"),
        (100_000_000 + 151 - 1, "C>T"): ("functional", "-", 151 - 1, 151 - 1, "C>T"),
    }, "Mutations"
    assert toy_gene.has_coverage("1", 100_000_147)
    assert not toy_gene.has_coverage("5", 100_000_147)

    assert toy_gene.get_rsid((100_000_104, "T>A"))


def test_region_at(toy_gene):
    # Test gene auxiliary functions
    assert toy_gene.region_at(100_000_105) == (0, "tmp"), "region_at test"
    assert toy_gene.region_at(100_000_127) == (0, "i1"), "region_at test"
    assert toy_gene.region_at(100_000_155) == (0, "e3"), "region_at test"
    assert toy_gene.region_at(100_000_165) == (0, "down"), "region_at test"
    assert toy_gene.region_at(100_000_265) is None, "region_at test"
    assert toy_gene.region_at(100_000_005) == (1, "tmp"), "region_at test"
    assert toy_gene.region_at(100_000_027) == (1, "i1"), "region_at test"
    assert toy_gene.region_at(100_000_055) == (1, "e3"), "region_at test"
    assert toy_gene.region_at(100_000_080) == (1, "down"), "region_at test"
    assert toy_gene.region_at(90_000_000) is None, "region_at test"


def test_get_functional(toy_gene):
    assert toy_gene.get_functional((100_000_118, "insTT")) == "frameshift"
    assert toy_gene.get_functional((100_000_111, "A>C"), infer=False) is None
    assert toy_gene.get_functional((100_000_111, "A>C"), infer=True) == "V1A"
    assert not toy_gene.is_functional((100_000_111, "A>C"), infer=False)
    assert toy_gene.is_functional((100_000_104, "T>A"))


def test_get_rsid(toy_gene):
    assert toy_gene.get_rsid(100_000_114, "T>A") == "rs28371732"
    assert toy_gene.get_rsid(Mutation(100_000_118, "insTT"), default=False) == "-"
    assert toy_gene.get_rsid(Mutation(100_000_118, "insTT")) == "100000119.insTT"


def test_get_refseq(toy_gene):
    assert toy_gene.get_refseq(100_000_114, "T>A") == "115T>A"
    assert toy_gene.get_refseq(100_000_114, "T>A", from_atg=True) == "6T>A"
