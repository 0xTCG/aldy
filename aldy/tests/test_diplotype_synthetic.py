# 786
# Aldy source: test_diplotype_synthetic.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa

from aldy.solutions import SolvedAllele, MinorSolution, MajorSolution
from aldy.diplotype import estimate_diplotype
from aldy.gene import Mutation


def assert_diplotype(gene, test, majors):
    sols = []
    for ma in majors:
        if isinstance(ma, tuple):
            sols.append(
                SolvedAllele(
                    ma[0], None, [Mutation(m[0], m[1], True) for m in ma[1]], []
                )
            )
        else:
            sols.append(SolvedAllele(ma, None, [], []))
    minor = MinorSolution(0, sols, MajorSolution(0, sols, None))
    res = estimate_diplotype(gene, minor)
    assert test == res, "Diplotype not equal"
    assert test == minor.diplotype, "Diplotype not set"


def test_basic(toy_gene):
    assert_diplotype(toy_gene, "*1/*1", ["1", "1.a"])


def test_tandem(toy_gene):
    assert_diplotype(toy_gene, "*1+*4/*3", ["4", "1.a", "3"])


def test_multi(toy_gene):
    assert_diplotype(toy_gene, "*1+*4/*3+*3", ["4", "1", "3", "3"])
    assert_diplotype(toy_gene, "*1+*4/*3+*3+*3+*3", ["4", "1", "3", "3", "3", "3"])


def test_combination(toy_gene):
    assert_diplotype(toy_gene, "*1+*4/*1+*4+*3", ["4", "1.a", "1", "3", "4"])


def test_novel(toy_gene):
    assert_diplotype(
        toy_gene, "*1+*4+*3/*1-like+*4", ["4", ("1", [(151, "SNP.CT")]), "1", "3", "4"]
    )


def test_deletion(toy_gene):
    assert_diplotype(toy_gene, "*2/*6", ["6", "2"])
