# 786
# Aldy source: test_minor_synthetic.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa
import collections

from aldy.profile import Profile
from aldy.sam import Sample
from aldy.solutions import CNSolution, SolvedAllele, MajorSolution
from aldy.coverage import Coverage
from aldy.gene import Mutation
from aldy.minor import estimate_minor
from aldy.common import SOLUTION_PRECISION, sorted_tuple


profile = Profile("test")


def assert_minor(gene, solver, data, shallow=False, skip_check=False):
    cn_sol = CNSolution(gene, 0, list(collections.Counter(data["cn"]).elements()))

    cov = collections.defaultdict(dict)
    for (pos, op), c in data["data"].items():
        cov[pos][op] = [(60, 60)] * c

    profile = Profile("test")
    cov = Coverage(gene, profile, None, cov, None, {})
    if "phase" in data:
        cov.sam = Sample.__new__(Sample)
        cov.sam.phases = {
            r: {k: v for k, v in ph.items()} for r, ph in data["phase"].items()
        }

    if isinstance(data["major"], tuple):
        major = MajorSolution(
            0,
            {SolvedAllele(gene, m): c for m, c in data["major"][0].items()},
            cn_sol,
            [Mutation(*m) for m in data["major"][1:]],
        )
    else:
        major = MajorSolution(
            0,
            {SolvedAllele(gene, m): c for m, c in data["major"].items()},
            cn_sol,
            [],
        )

    sols = estimate_minor(gene, cov, [major], solver)
    if skip_check:
        return sols

    if "score" in data:
        for s in sols:
            print(s, data)
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


def test_basic(toy_gene, solver):
    # Test two copies of *1
    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {(100_000_114, "_"): 20},
            "major": {"1": 2},
            "sol": [("1.001", [], []), ("1.001", [], [])],
            "score": 0,
        },
    )


def test_minor(toy_gene, solver):
    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {(100_000_114, "_"): 9, (100_000_114, "T>A"): 11},
            "major": {"1": 2},
            "sol": [("1.001", [], []), ("1.002", [], [])],
            "score": 0.2,
        },
    )


def test_miss(toy_gene, solver):
    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_114, "_"): 10,
                (100_000_114, "T>A"): 10,
                (100_000_147, "_"): 20,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "major": {"1": 1, "3": 1},
            "sol": [("1.002", [], []), ("3.001", [(100_000_147, "insA")], [])],
            "score": profile.minor_miss,
        },
    )


def test_add(toy_gene, solver):
    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_114, "_"): 0,
                (100_000_114, "T>A"): 20,
                (100_000_147, "_"): 20,
                (100_000_147, "insA"): 10,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "major": {"1": 1, "3": 1},
            "sol": [("1.002", [], []), ("3.001", [], [(100_000_114, "T>A")])],
            "score": profile.minor_add,
        },
    )


def test_phase(toy_gene, solver):
    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_114, "_"): 10,
                (100_000_114, "T>A"): 10,
                (100_000_147, "_"): 20,
                (100_000_147, "insA"): 10,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "major": {"1": 1, "3": 1},
            "sol": [("1.002", [], []), ("3.001", [], [])],
            "score": 0,
        },
    )

    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_114, "_"): 10,
                (100_000_114, "T>A"): 10,
                (100_000_147, "_"): 20,
                (100_000_147, "insA"): 10,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "major": {"1": 1, "3": 1},
            "sol": [("1.001", [], []), ("3.001", [], [(100_000_114, "T>A")])],
            "score": profile.minor_add,
            "phase": {
                "r1": {100_000_114: "T>A", 100_000_150: "C>T"},
                "r2": {100_000_114: "T>A", 100_000_150: "C>T"},
                "r3": {100_000_114: "T>A", 100_000_150: "C>T"},
            },
        },
    )

    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_104, "_"): 10,
                (100_000_104, "T>A"): 10,
                (100_000_147, "_"): 20,
                (100_000_147, "insA"): 10,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "major": {"1": 1, "3": 1},
            "sol": [("1.001", [], []), ("3.001", [], [])],
            "score": 2,
            # ignored because 104T>A is not part of major solution
            "phase": {
                "r1": {100_000_104: "T>A", 100_000_150: "C>T"},
                "r2": {100_000_104: "T>A", 100_000_150: "C>T"},
                "r3": {100_000_104: "T>A", 100_000_150: "C>T"},
            },
        },
    )

    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2},
            "data": {
                (100_000_110, "_"): 10,
                (100_000_110, "delAC"): 10,
                (100_000_147, "_"): 20,
                (100_000_147, "insA"): 10,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "major": ({"1": 1, "3": 1}, (100_000_110, "delAC")),
            "sol": [("1.001", [], []), ("3.001", [], [(100_000_110, "delAC")])],
            "score": profile.minor_add * 1.5 + 2,
            "phase": {"r1": {100_000_110: "delAC", 100_000_150: "C>T"}},
        },
    )


def test_major_novel(toy_gene, solver):
    assert_minor(
        toy_gene,
        solver,
        {
            "cn": {"1": 2, "6": 1},
            "data": {
                (100_000_110, "_"): 10,
                (100_000_110, "delAC"): 10,
                (100_000_114, "_"): 10,
                (100_000_114, "T>A"): 10,
                (100_000_147, "_"): 20,
                (100_000_150, "_"): 10,
                (100_000_150, "C>T"): 10,
            },
            "major": ({"1": 1, "3": 1, "6": 1}, (100_000_110, "delAC")),
            "sol": [
                ("1.002", [], [(100_000_110, "delAC")]),
                ("3.001", [(100_000_147, "insA")], []),
                ("6.001", [], []),
            ],
            "score": profile.minor_miss + (profile.minor_add * 1.5),
        },
    )
