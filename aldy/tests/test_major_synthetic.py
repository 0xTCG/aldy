# 786
# Aldy source: test_major_synthetic.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Any, List
import pytest  # noqa
import collections


from aldy.profile import Profile
from aldy.major import estimate_major
from aldy.solutions import CNSolution
from aldy.coverage import Coverage
from aldy.common import SOLUTION_PRECISION


profile = Profile("test")


def assert_major(gene, solver, major):
    cn_sol = CNSolution(gene, 0, list(collections.Counter(major["cn"]).elements()))

    cov = collections.defaultdict(dict)
    for (pos, op), c in major["data"].items():
        cov[pos][op] = [(60, 60)] * c
    cov = Coverage(gene, profile, None, cov, None, {})
    sols = estimate_major(gene, cov, cn_sol, solver)

    if "score" in major:
        for s in sols:
            assert abs(major["score"] - s.score) < SOLUTION_PRECISION, "Score"

    sols_parsed: List[Any] = [
        dict(
            collections.Counter(
                [i.major for i, v in s.solution.items() for _ in range(v)]
            )
        )
        for s in sols
    ]
    for si, s in enumerate(sols):
        if s.added:
            sols_parsed[si] = (sols_parsed[si], *[tuple(m) for m in s.added])

    def k(x):
        if isinstance(x, tuple):
            x = x[0]
        return tuple(x.keys())

    assert sorted(major["sol"], key=k) == sorted(sols_parsed, key=k)


def test_basic(toy_gene, solver):
    # Test two copies of *1
    assert_major(
        toy_gene, solver, {"cn": {"1": 2}, "data": {}, "sol": [{"1": 2}], "score": 0}
    )


def test_deletion(toy_gene, solver):
    # Test a single copy of *1C
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 1, "6": 1},
            "data": {(100_000_114, "_"): 0, (100_000_114, "T>A"): 10},
            "sol": [{"1": 1, "6": 1}],
            "score": 0,
        },
    )


def test_two_copies(toy_gene, solver):
    # Test two copies (*1/*2)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_110, "_"): 10,
                (100_000_110, "delAC"): 10,
                (100_000_118, "_"): 20,
                (100_000_118, "insTT"): 10,
            },
            "sol": [{"1": 1, "2": 1}],
            "score": 0,
        },
    )
    # Test two copies (*2/*3)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_110, "_"): 10,
                (100_000_110, "delAC"): 10,
                (100_000_118, "_"): 20,
                (100_000_118, "insTT"): 10,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "sol": [{"2": 1, "3": 1}],
            "score": 0,
        },
    )
    # Test slightly perturbed two copies (*2/*3)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_110, "_"): 11,
                (100_000_110, "delAC"): 9,
                (100_000_118, "_"): 22,
                (100_000_118, "insTT"): 8,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 9,
            },
            "sol": [{"2": 1, "3": 1}],
            "score": 2 / 10 + 3 / 11 + 1 / 10,
        },
    )


def test_multiple_copies(toy_gene, solver):
    # Test four copies (*1+*1/*2+*2)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 4},
            "data": {
                (100_000_110, "_"): 20,
                (100_000_110, "delAC"): 20,
                (100_000_118, "_"): 40,
                (100_000_118, "insTT"): 20,
            },
            "sol": [{"1": 2, "2": 2}],
            "score": 0,
        },
    )


def test_left_fusion(toy_gene, solver):
    # Test left fusion that has no SNPs (*4[*1]/*1+*2)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2, "4": 1},
            "data": {
                (100_000_110, "_"): 10,
                (100_000_110, "delAC"): 10,
                (100_000_118, "_"): 20,
                (100_000_118, "insTT"): 10,
            },
            "sol": [{"1": 1, "2": 1, "4#1": 1}],
            "score": 0,
        },
    )
    # Test left fusion that has SNPs (*4[*3]/*1+*3)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2, "4": 1},
            "data": {(100_000_150, "_"): 10, (100_000_150, "C>T"): 20},
            "sol": [{"1": 1, "3": 1, "4#3": 1}, {"3": 2, "4#1": 1}],
            "score": 0,
        },
    )
    # Test left fusion combination (*4[*3]/*1C+*3)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2, "4": 1},
            "data": {
                (100_000_104, "_"): 10,
                (100_000_104, "T>A"): 10,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 20,
            },
            "sol": [{"1C": 1, "3": 1, "4#3": 1}],
            "score": 0,
        },
    )
    # Test left fusion combination (*4[*3]/*4[*1]+*3)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 1, "4": 2},
            "data": {
                (100_000_104, "_"): 0,
                (100_000_104, "T>A"): 10,
                (100_000_150, "_"): 20,
                (100_000_150, "C>T"): 10,
            },
            "sol": [{"1C": 1, "4#3": 1, "4#1": 1}],
            "score": 0,
        },
    )


def test_right_fusion(toy_gene, solver):
    # Test right fusion (*5/*1)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 1, "5": 1},
            "data": {(100_000_110, "_"): 10, (100_000_110, "delAC"): 10},
            "sol": [{"1": 1, "5": 1}],
            "score": 0,
        },
    )
    # Test right fusion (*5/*2)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 1, "5": 1},
            "data": {
                (100_000_110, "_"): 0,
                (100_000_110, "delAC"): 20,
                (100_000_118, "_"): 20,
                (100_000_118, "insTT"): 10,
            },
            "sol": [{"2": 1, "5": 1}],
            "score": 0,
        },
    )
    # Test right fusion (*5/*3)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 1, "5": 1},
            "data": {
                (100_000_110, "_"): 10,
                (100_000_110, "delAC"): 10,
                (100_000_150, "_"): 0,
                (100_000_150, "C>T"): 10,
            },
            "sol": [{"3": 1, "5": 1}],
            "score": 0,
        },
    )


def test_novel_mutations(toy_gene, solver):
    # Test novel mutations within a single gene (other is deleted)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 1, "6": 1},
            "data": {(100_000_110, "_"): 0, (100_000_110, "delAC"): 20},
            "sol": [({"1": 1, "6": 1}, (100_000_110, "delAC"))],
            "score": profile.major_novel + 0.1 + 1,
            # Last 1 is needed as the model right now considers novel mutations
            # as "independent" alleles and still (incorrectly) counts the '_'
            # for both *1 here.
        },
    )
    # Test novel mutations
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {(100_000_110, "_"): 10, (100_000_110, "delAC"): 10},
            "sol": [({"1": 2}, (100_000_110, "delAC"))],
            "score": profile.major_novel + 0.1 + 1,
        },
    )
