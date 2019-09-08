#!/usr/bin/env python
# 786

# Aldy source: test_cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest
import os

import aldy.genotype


def test_samples(real_sample):
    data, path, profile, reference = real_sample
    solver = os.getenv("ALDY_SOLVER", default="gurobi")
    sols = aldy.genotype.genotype(
        "cyp2d6",
        sam_path=path,
        profile=profile,
        output_file=None,
        reference=reference,
        solver=solver,
    )
    assert sorted(data) == sorted([s.diplotype for s in sols])
