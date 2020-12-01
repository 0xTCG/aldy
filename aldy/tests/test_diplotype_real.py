# 786
# Aldy source: test_diplotype_real.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa

from .test_diplotype_synthetic import assert_diplotype


def test_basic(real_gene):
    assert_diplotype(real_gene, "*1 / *1", ["1", "1"])


def test_tandem(real_gene):
    assert_diplotype(real_gene, "*2 + *2 / *71", ["2", "2", "71"])
    assert_diplotype(real_gene, "*4.021 + *4C / *41", ["4C", "41", "4.021"])
    assert_diplotype(real_gene, "*4 + *4.019 / *4J", ["4.019", "4", "4J"])
    assert_diplotype(real_gene, "*3 / *68 + *4", ["3", "4", "68"])


def test_36(real_gene):
    assert_diplotype(real_gene, "*36 + *10 / *36 + *41", ["10", "36", "41", "36"])
    assert_diplotype(real_gene, "*1 + *36 / *36 + *10", ["1", "10", "36", "36"])
    assert_diplotype(real_gene, "*10 / *36 + *10", ["10", "36", "10"])
    assert_diplotype(real_gene, "*36 + *10 / *36 + *10", ["10", "36", "36", "10"])
    assert_diplotype(real_gene, "*10 + *10 / *36 + *10", ["10", "36", "10", "10"])


def test_fusion(real_gene):
    assert_diplotype(real_gene, "*2 + *2 / *68 + *4", ["2", "4", "68", "2"])
    assert_diplotype(real_gene, "*1 / *79 + *2", ["1", "2", "79#2"])
    assert_diplotype(real_gene, "*2 / *78 + *2", ["2", "78#2", "2"])
    assert_diplotype(real_gene, "*1 / *78 + *2", ["2", "1", "78#2"])
