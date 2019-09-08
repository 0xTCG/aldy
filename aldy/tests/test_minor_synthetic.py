# 786
# Aldy source: test_minor_synthetic.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa
import os
import collections

from aldy.solutions import CNSolution, SolvedAllele, MajorSolution
from aldy.coverage import Coverage
from aldy.gene import Mutation
from aldy.major import NOVEL_MUTATION_PENAL
from aldy.minor import estimate_minor, ADD_PENALTY_FACTOR, MISS_PENALTY_FACTOR
from aldy.common import SOLUTION_PRECISION, sorted_tuple


def assert_minor(gene, data, shallow=False):
    solver = os.getenv("ALDY_SOLVER", default="gurobi")
    cn_sol = CNSolution(0, list(collections.Counter(data["cn"]).elements()), gene)

    cov = collections.defaultdict(dict)
    for (pos, op), c in data["data"].items():
        cov[pos][op] = c
    cov = Coverage(cov, 0.5, {})

    major_solved = {
        SolvedAllele(maj, None, tuple(), tuple())
        if not isinstance(maj, tuple)
        else SolvedAllele(
            maj[0], None, tuple(Mutation(mp, mo, True) for mp, mo in maj[1:]), tuple()
        ): cnt
        for maj, cnt in data["major"].items()
    }
    major = MajorSolution(0, major_solved, cn_sol)
    sols = estimate_minor(gene, cov, [major], solver)

    if "score" in data:
        for s in sols:
            assert abs(data["score"] - s.score) < SOLUTION_PRECISION, "Score"

    if not shallow:
        sols_expected = [
            sorted_tuple(
                (i[0], sorted_tuple(i[1]), sorted_tuple(i[2])) for i in data["sol"]
            )
        ]
        sols_parsed = [
            sorted_tuple(
                (
                    i.minor,
                    sorted_tuple((m.pos, m.op) for m in i.missing),
                    sorted_tuple((m.pos, m.op) for m in i.added),
                )
                for i in s.solution
            )
            for s in sols
        ]
        assert sorted(sols_expected) == sorted(sols_parsed)
    else:
        # As assignments can vary between multiple optimal solutions,
        # just test minor allele assignments
        eall, emiss, enew = set(), set(), set()
        for i in data["sol"]:
            eall.add(i[0])
            emiss |= set(i[1])
            enew |= set(i[2])
        assert len(sols) == 1, "Single minor solution"
        pall, pmiss, pnew = set(), set(), set()
        for i in sols[0].solution:
            pall.add(i.minor)
            pmiss |= set((m.pos, m.op) for m in i.missing)
            pnew |= set((m.pos, m.op) for m in i.added)

        assert eall == pall, "Alleles"
        assert emiss == pmiss, "Missing mutations"
        assert enew == pnew, "Novel mutations"
    return sols[0].score if len(sols) > 0 else 0


def test_basic(toy_gene):
    # Test two copies of *1
    assert_minor(
        toy_gene,
        {
            "cn": {"1": 2},
            "data": {(115, "_"): 20},
            "major": {"1": 2},
            "sol": [("1", [], []), ("1", [], [])],
            "score": 0,
        },
    )


def test_minor(toy_gene):
    assert_minor(
        toy_gene,
        {
            "cn": {"1": 2},
            "data": {(115, "_"): 9, (115, "SNP.TA"): 11},
            "major": {"1": 2},
            "sol": [("1", [], []), ("1B", [], [])],
            "score": 0.2,
        },
    )


def test_miss(toy_gene):
    assert_minor(
        toy_gene,
        {
            "cn": {"1": 2},
            "data": {
                (115, "_"): 10,
                (115, "SNP.TA"): 10,
                (148, "_"): 20,
                (151, "_"): 10,
                (151, "SNP.CT"): 10,
            },
            "major": {"1": 1, "3": 1},
            "sol": [("1B", [], []), ("3", [(148, "INS.A")], [])],
            "score": MISS_PENALTY_FACTOR,
        },
    )


def test_add(toy_gene):
    assert_minor(
        toy_gene,
        {
            "cn": {"1": 2},
            "data": {
                (115, "_"): 0,
                (115, "SNP.TA"): 20,
                (148, "_"): 20,
                (148, "INS.A"): 10,
                (151, "_"): 10,
                (151, "SNP.CT"): 10,
            },
            "major": {"1": 1, "3": 1},
            "sol": [("1B", [], []), ("3", [], [(115, "SNP.TA")])],
            "score": ADD_PENALTY_FACTOR,
        },
    )


def test_major_novel(toy_gene):
    assert_minor(
        toy_gene,
        {
            "cn": {"1": 2, "6": 1},
            "data": {
                (111, "_"): 10,
                (111, "DEL.AC"): 10,
                (115, "_"): 10,
                (115, "SNP.TA"): 10,
                (148, "_"): 20,
                (151, "_"): 10,
                (151, "SNP.CT"): 10,
            },
            "major": {("1", (111, "DEL.AC")): 1, "3": 1, "6": 1},
            "sol": [
                ("1B", [], [(111, "DEL.AC")]),
                ("3", [(148, "INS.A")], []),
                ("6DEL", [], []),
            ],
            "score": MISS_PENALTY_FACTOR + NOVEL_MUTATION_PENAL,
        },
    )
