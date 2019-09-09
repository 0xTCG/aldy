# 786
# Aldy source: conftest.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest
import ast

from aldy.gene import Gene
from aldy.common import script_path


@pytest.fixture
def toy_gene():
    return Gene(script_path("aldy.tests/toy.yml"))


@pytest.fixture
def real_gene():
    return Gene(script_path("aldy.resources.genes/cyp2d6.yml"))


def pytest_addoption(parser):
    parser.addoption("--samples", action="store", default=None)
    parser.addoption("--solvers", action="store", default=None)


def pytest_generate_tests(metafunc):
    path = metafunc.config.getoption("samples")
    solvers = metafunc.config.getoption("solvers")
    if solvers is None:
        solvers = "gurobi"
    solvers = solvers.split(",")

    if "solver" in metafunc.fixturenames:
        metafunc.parametrize("solver", solvers)
    if "real_sample" in metafunc.fixturenames and not path:
        metafunc.parametrize("real_sample", [])
    elif "real_sample" in metafunc.fixturenames and path:

        def read_data(path):
            with open(script_path(path)) as f:
                data = f.read()
                data = ast.literal_eval(data)
                return data

        samples = []
        for file, data in read_data("aldy.tests.resources/data-pgx1.json").items():
            samples.append(
                (
                    data,
                    f"{path}/cdc/pgrnseq-v1/bams/{file}.cram",
                    "pgrnseq-v1",
                    f"{path}/cram-genome.fa",
                )
            )
        for file, data in read_data("aldy.tests.resources/data-pgx2.json").items():
            samples.append(
                (
                    data,
                    f"{path}/baylor/pgrnseq-v2/bams/{file}.cram",
                    "pgrnseq-v2",
                    f"{path}/cram-genome.fa",
                )
            )
        for file, data in read_data("aldy.tests.resources/data-illumina.json").items():
            samples.append(
                (data, f"{path}/1000genomes-illumina/bams/{file}.bam", "illumina", None)
            )
        metafunc.parametrize("real_sample", samples)
