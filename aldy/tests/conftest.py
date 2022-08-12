# 786
# Aldy source: conftest.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest

from aldy.gene import Gene
from aldy.common import script_path


@pytest.fixture
def toy_gene():
    return Gene(script_path("aldy.tests.resources/toy.yml"))


@pytest.fixture
def real_gene():
    return Gene(script_path("aldy.resources.genes/cyp2d6.yml"))


def pytest_addoption(parser):
    parser.addoption("--solvers", action="store", default=None)


def pytest_generate_tests(metafunc):
    solvers = metafunc.config.getoption("solvers")
    if solvers is None:
        solvers = "any"
    solvers = solvers.split(",")

    if "solver" in metafunc.fixturenames:
        metafunc.parametrize("solver", solvers)
