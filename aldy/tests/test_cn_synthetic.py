# 786
# Aldy source: test_cn_synthetic.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa

from aldy.cn import PARSIMONY_PENALTY, LEFT_FUSION_PENALTY, solve_cn_model
from aldy.common import SOLUTION_PRECISION


def assert_cn(gene, solver, expected, cov, expected_obj=None, gap=0):
    sols = solve_cn_model(
        gene,
        cn_configs=gene.cn_configs,
        max_cn=20,
        region_coverage=cov,
        solver=solver,
        gap=gap,
    )
    if expected_obj:
        for s in sols:
            assert abs(s.score - expected_obj) < SOLUTION_PRECISION, "Score mismatch"
    best_sol_score = min(sols, key=lambda s: s.score).score
    for s in sols:
        assert (
            s.score - best_sol_score * (1 + gap)
        ) < SOLUTION_PRECISION, "Gap violation"
    sols = sorted([sorted(s.solution.items()) for s in sols])
    expected = sorted([sorted(s.items()) for s in expected])
    assert sols == expected, "Solution mismatch"


def make_coverage(gene, lst):
    cov = {}
    for r in sorted(gene.unique_regions):
        cov[r] = next(lst)
    return cov


def test_basic(toy_gene, solver):
    # Test two copies of *1
    assert_cn(
        toy_gene,
        solver,
        [{"1": 2}],
        make_coverage(toy_gene, zip([2, 2, 2, 2, 2], [2, 2, 2, 2, 2])),
        2 * PARSIMONY_PENALTY,
    )
    # Test two copies of *1 with slightly perturbed coverage
    assert_cn(
        toy_gene,
        solver,
        [{"1": 2}],
        make_coverage(
            toy_gene, zip([1.8, 2.2, 2.1, 1.9, 1.7], [2.05, 1.95, 2, 2, 2.7])
        ),
        2 * PARSIMONY_PENALTY + (0.2 * 2 + 0.1 * 2 + 0.3 + 0.05 * 2 + 0.7),
    )


def test_deletion(toy_gene, solver):
    # Test a single copy of *1 (*6 is deletion allele)
    assert_cn(
        toy_gene,
        solver,
        [{"1": 1, "6": 1}],
        make_coverage(toy_gene, zip([1, 1, 1, 1, 1], [2, 2, 2, 2, 2])),
        2 * PARSIMONY_PENALTY,
    )
    # Test whole gene deletion
    assert_cn(
        toy_gene,
        solver,
        [{"6": 2}],
        make_coverage(toy_gene, zip([0, 0, 0, 0, 0], [2, 2, 2, 2, 2])),
        2 * PARSIMONY_PENALTY,
    )
    # TODO: test 2 deletions with no coverage
    # (needs special handling as deletions imply 2D7 coverage)
    # assert_cn(toy_gene,
    #           [{'6': 2}],
    #           ke_coverage(zip([0,0, 0,0, 0], [0,0, 0,0, 0])),
    #           * PARSIMONY_PENALTY)


def test_left_fusion(toy_gene, solver):
    # Test two fused copies (*4 is defined as 00011|11100)
    assert_cn(
        toy_gene,
        solver,
        [{"4": 2}],
        make_coverage(toy_gene, zip([0, 0, 0, 2, 2], [2, 2, 2, 0, 0])),
        2 * PARSIMONY_PENALTY + 2 * LEFT_FUSION_PENALTY,
    )
    # Test one fused and one normal (*1) allele
    # Note: each left fusion operates on the whole genic region;
    #       thus, the maximum number of left fusions is 2
    assert_cn(
        toy_gene,
        solver,
        [{"4": 2, "1": 1}],
        make_coverage(toy_gene, zip([1, 1, 1, 3, 3], [2, 2, 2, 0, 0])),
        3 * PARSIMONY_PENALTY + 2 * LEFT_FUSION_PENALTY,
    )


def test_right_fusion(toy_gene, solver):
    # Test one fused and one normal (*1) allele (*5 is defined as 11000|11222)
    assert_cn(
        toy_gene,
        solver,
        [{"1": 1, "5": 1}],
        make_coverage(toy_gene, zip([2, 2, 1, 1, 1], [2, 2, 3, 3, 3])),
        2 * PARSIMONY_PENALTY,
    )


def test_multiplication(toy_gene, solver):
    # Test twelve copies of *1
    assert_cn(
        toy_gene,
        solver,
        [{"1": 12}],
        make_coverage(toy_gene, zip([12, 12, 12, 12, 12], [2, 2, 2, 2, 2])),
        12 * PARSIMONY_PENALTY,
    )
    # Test seven copies of *1 and one fused *5 copy
    assert_cn(
        toy_gene,
        solver,
        [{"1": 7, "5": 1}],
        make_coverage(toy_gene, zip([8, 8, 7, 7, 7], [2, 2, 3, 3, 3])),
        8 * PARSIMONY_PENALTY,
    )
