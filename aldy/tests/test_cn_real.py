# 786
# Aldy source: test_cn_real.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa
import re

from aldy.gene import GeneRegion
from .test_cn_synthetic import assert_cn


def make_coverage(d):
    cov = {}
    for k, v in d.items():
        k = re.split(r"(\d+)", k)
        cov[GeneRegion(int(k[1]), k[2])] = v
    return cov


def test_many_copies_multiple_solutions(real_gene, solver):
    # HG00465
    assert_cn(
        real_gene,
        solver,
        [
            {"1": 2, "36": 2},
            {"1": 2, "61": 2},
            {"1": 2, "63": 2},
            {"1": 2, "36": 1, "61": 1},
            {"1": 2, "36": 1, "63": 1},
            {"1": 2, "61": 1, "63": 1},
        ],
        make_coverage(
            {
                "1e": (3.8, 2.1),
                "1i": (3.3, 2.5),
                "2e": (3.7, 2.1),
                "2i": (4.0, 1.9),
                "5e": (3.9, 2.0),
                "5i": (4.1, 2.0),
                "6e": (4.0, 2.0),
                "6i": (3.9, 1.8),
                "3e": (4.5, 1.8),
                "9e": (2.5, 3.6),
                "11pce": (0, 4.1),
            }
        ),
    )


def test_right_fusion(real_gene, solver):
    # HG01190
    assert_cn(
        real_gene,
        solver,
        [{"1": 1, "68": 1}],
        make_coverage(
            {
                "1e": (1.8, 2.0),
                "1i": (1.9, 2.0),
                "2e": (0.9, 2.8),
                "2i": (1.0, 3.0),
                "5e": (1.0, 3.3),
                "5i": (1.1, 3.2),
                "6e": (1.0, 3.1),
                "6i": (1.1, 2.9),
                "3e": (1.1, 2.7),
                "9e": (1.3, 2.7),
                "11pce": (0, 3.1),
            }
        ),
    )


def test_normal(real_gene, solver):
    # HG02260
    assert_cn(
        real_gene,
        solver,
        [{"1": 2}],
        make_coverage(
            {
                "1e": (1.8, 2.0),
                "1i": (1.9, 2.0),
                "2e": (2.0, 1.9),
                "2i": (1.9, 1.8),
                "5e": (1.8, 2.1),
                "5i": (1.9, 2.0),
                "6e": (2.0, 2.1),
                "6i": (2.0, 1.9),
                "3e": (2.0, 1.8),
                "9e": (2.1, 1.8),
                "11pce": (0, 2.0),
            }
        ),
    )


def test_deletion(real_gene, solver):
    # NA12336
    assert_cn(
        real_gene,
        solver,
        [{"1": 1, "5": 1}],
        make_coverage(
            {
                "1e": (1.0, 1.9),
                "1i": (1.3, 1.6),
                "2e": (0.9, 1.6),
                "2i": (1.0, 1.9),
                "5e": (0.9, 2.2),
                "5i": (0.9, 2.0),
                "6e": (0.9, 1.9),
                "6i": (0.9, 1.8),
                "3e": (1.1, 1.7),
                "9e": (1.0, 1.8),
                "11pce": (0, 1.9),
            }
        ),
    )


def test_right_fusion_with_copy(real_gene, solver):
    # NA12878
    assert_cn(
        real_gene,
        solver,
        [{"1": 2, "68": 1}],
        make_coverage(
            {
                "1e": (2.6, 2.0),
                "1i": (2.4, 2.3),
                "2e": (1.8, 2.9),
                "2i": (2.0, 2.9),
                "5e": (1.9, 3.0),
                "5i": (1.9, 3.0),
                "6e": (1.8, 2.9),
                "6i": (1.9, 2.9),
                "3e": (2.2, 2.5),
                "9e": (2.1, 2.6),
                "11pce": (0, 3.0),
            }
        ),
    )


def test_normal2(real_gene, solver):
    # NA19239
    assert_cn(
        real_gene,
        solver,
        [{"1": 2}],
        make_coverage(
            {
                "1e": (1.6, 2.2),
                "1i": (2.0, 2.0),
                "2e": (2.0, 1.9),
                "2i": (2.0, 2.1),
                "5e": (1.9, 2.2),
                "5i": (2.0, 2.0),
                "6e": (1.9, 2.0),
                "6i": (1.9, 2.0),
                "3e": (2.1, 1.9),
                "9e": (2.1, 2.0),
                "11pce": (0, 2.1),
            }
        ),
    )


def test_left_fusion(real_gene, solver):
    # NA19790
    assert_cn(
        real_gene,
        solver,
        [{"1": 2, "78": 1}, {"1": 2, "67": 1}],
        make_coverage(
            {
                "1e": (1.9, 2.0),
                "1i": (2.0, 1.9),
                "2e": (2.0, 2.0),
                "2i": (2.0, 2.0),
                "5e": (2.8, 1.2),
                "5i": (2.9, 1.0),
                "6e": (2.7, 0.9),
                "6i": (2.8, 0.9),
                "3e": (2.1, 1.9),
                "9e": (3.1, 1.0),
                "11pce": (0, 1.0),
            }
        ),
    )


def test_gap(real_gene, solver):
    data = make_coverage(
        {
            "1e": (1.0, 1.9),
            "1i": (1.3, 1.6),
            "2e": (0.9, 1.6),
            "2i": (1.0, 1.9),
            "5e": (0.9, 2.2),
            "5i": (0.9, 2.0),
            "6e": (0.9, 1.9),
            "6i": (0.9, 1.8),
            "3e": (1.1, 1.7),
            "9e": (1.0, 1.8),
            "11pce": (0, 1.9),
        }
    )
    assert_cn(real_gene, solver, [{"1": 1, "5": 1}], data, gap=0.0)
    assert_cn(
        real_gene,
        solver,
        [
            {"1": 1, "5": 1},
            {"13": 1, "68": 1},
            {"16": 1, "36": 1},
            {"16": 1, "61": 1},
            {"16": 1, "63": 1},
            {"36": 1, "76": 1},
            {"36": 1, "77": 1},
            {"61": 1, "76": 1},
            {"61": 1, "77": 1},
            {"63": 1, "76": 1},
            {"63": 1, "77": 1},
        ],
        data,
        gap=0.1,
    )
    assert_cn(
        real_gene,
        solver,
        [
            {"1": 1, "5": 1},
            {"13": 1, "68": 1},
            {"16": 1, "36": 1},
            {"16": 1, "61": 1},
            {"16": 1, "63": 1},
            {"36": 1, "76": 1},
            {"36": 1, "77": 1},
            {"61": 1, "76": 1},
            {"61": 1, "77": 1},
            {"63": 1, "76": 1},
            {"63": 1, "77": 1},
            {"68": 1, "79": 1},
        ],
        data,
        gap=0.65,
    )
