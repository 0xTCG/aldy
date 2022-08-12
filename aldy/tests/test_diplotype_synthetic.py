# 786
# Aldy source: test_diplotype_synthetic.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from collections import Counter
import pytest  # noqa

from aldy.solutions import SolvedAllele, MinorSolution, MajorSolution, CNSolution
from aldy.diplotype import estimate_diplotype
from aldy.gene import Mutation


def assert_diplotype(gene, test, majors):
    sols = []
    cns = []
    for ma in majors:
        if isinstance(ma, tuple):
            sols.append(SolvedAllele(gene, ma[0], added=[Mutation(*m) for m in ma[1]]))
            cns.append(ma[0])
        else:
            sols.append(SolvedAllele(gene, ma))
    cns = [c for m in cns for c, cc in gene.cn_configs.items() if m in cc.alleles]
    minor = MinorSolution(
        0, sols, MajorSolution(0, Counter(sols), CNSolution(gene, 0, cns), [])
    )
    estimate_diplotype(gene, minor)
    assert test == minor.get_major_diplotype(), "Diplotype not equal"


def test_basic(toy_gene):
    assert_diplotype(toy_gene, "*1 / *1C", ["1", "1C"])


def test_tandem(toy_gene):
    assert_diplotype(toy_gene, "*1C + *4 / *3", ["4", "1C", "3"])


def test_multi(toy_gene):
    assert_diplotype(toy_gene, "*1 + *4 / *3 + *3", ["4", "1", "3", "3"])
    assert_diplotype(
        toy_gene, "*1C + *4 / *3 + *3 + *3 + *3", ["4", "1C", "3", "3", "3", "3"]
    )


def test_combination(toy_gene):
    assert_diplotype(toy_gene, "*1 + *4 / *1C + *4 + *3", ["4", "1C", "1", "3", "4"])


def test_novel(toy_gene):
    assert_diplotype(
        toy_gene,
        "*1 + *4 / *1C+100000151.C>T + *4 + *3",
        ["4", ("1C", [(100_000_150, "C>T")]), "1", "3", "4"],
    )


def test_deletion(toy_gene):
    assert_diplotype(toy_gene, "*2 / *6", ["6", "2"])
