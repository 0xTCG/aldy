# 786
# Aldy source: test_major_synthetic.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa
import collections

from aldy.major import estimate_major
from aldy.solutions import CNSolution
from aldy.coverage import Coverage
from aldy.major import NOVEL_MUTATION_PENAL
from aldy.common import SOLUTION_PRECISION


def assert_major(gene, solver, major):
    cn_sol = CNSolution(0, list(collections.Counter(major["cn"]).elements()), gene)

    cov = collections.defaultdict(dict)
    for (pos, op), c in major["data"].items():
        cov[pos][op] = c
    cov = Coverage(cov, 0.5, {})
    sols = estimate_major(gene, cov, cn_sol, solver)

    if "score" in major:
        for s in sols:
            assert abs(major["score"] - s.score) < SOLUTION_PRECISION, "Score"
    sols_expected = [
        sorted(collections.Counter(c).elements(), key=str) for c in major["sol"]
    ]
    sols_parsed = [
        sorted(
            collections.Counter(
                {
                    tuple([k.major] + [(m[0], m[1]) for m in k.added])
                    if len(k.added) > 0
                    else k.major: v
                    for k, v in s.solution.items()
                }
            ).elements(),
            key=str,
        )
        for s in sols
    ]
    assert sorted(sols_expected) == sorted(sols_parsed)


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
            "data": {(105, "_"): 0, (105, "SNP.TA"): 10},
            "sol": [{"1.a": 1, "6": 1}],
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
                (111, "_"): 10,
                (111, "DEL.AC"): 10,
                (119, "_"): 20,
                (119, "INS.TT"): 10,
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
                (111, "_"): 10,
                (111, "DEL.AC"): 10,
                (119, "_"): 20,
                (119, "INS.TT"): 10,
                (151, "_"): 10,
                (151, "SNP.CT"): 10,
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
                (111, "_"): 11,
                (111, "DEL.AC"): 9,
                (119, "_"): 22,
                (119, "INS.TT"): 8,
                (151, "_"): 10.5,
                (151, "SNP.CT"): 9.5,
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
                (111, "_"): 20,
                (111, "DEL.AC"): 20,
                (119, "_"): 40,
                (119, "INS.TT"): 20,
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
                (111, "_"): 10,
                (111, "DEL.AC"): 10,
                (119, "_"): 20,
                (119, "INS.TT"): 10,
            },
            "sol": [{"1": 1, "2": 1, "4/1": 1}],
            "score": 0,
        },
    )
    # Test left fusion that has SNPs (*4[*3]/*1+*3)
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2, "4": 1},
            "data": {(151, "_"): 10, (151, "SNP.CT"): 20},
            "sol": [{"1": 1, "3": 1, "4/3": 1}, {"3": 2, "4/1": 1}],
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
                (105, "_"): 10,
                (105, "SNP.TA"): 10,
                (151, "_"): 10,
                (151, "SNP.CT"): 20,
            },
            "sol": [{"1.a": 1, "3": 1, "4/3": 1}],
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
                (105, "_"): 0,
                (105, "SNP.TA"): 10,
                (151, "_"): 20,
                (151, "SNP.CT"): 10,
            },
            "sol": [{"1.a": 1, "4/3": 1, "4/1": 1}],
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
            "data": {(111, "_"): 10, (111, "DEL.AC"): 10},
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
                (111, "_"): 0,
                (111, "DEL.AC"): 20,
                (119, "_"): 20,
                (119, "INS.TT"): 10,
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
                (111, "_"): 10,
                (111, "DEL.AC"): 10,
                (151, "_"): 0,
                (151, "SNP.CT"): 10,
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
            "data": {(111, "_"): 0, (111, "DEL.AC"): 20},
            "sol": [{("1", (111, "DEL.AC")): 1, "6": 1}],
            "score": NOVEL_MUTATION_PENAL,
        },
    )
    # Test novel mutations
    assert_major(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {(111, "_"): 10, (111, "DEL.AC"): 10},
            "sol": [{"1": 1, ("1", (111, "DEL.AC")): 1}],
            "score": NOVEL_MUTATION_PENAL,
        },
    )
