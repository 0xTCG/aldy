#!/usr/bin/env python
# 786

# Aldy source: test_cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import pytest  # noqa

import aldy.genotype


def test_samples(real_sample, solver):
    data, path, profile, reference = real_sample
    sols = aldy.genotype.genotype(
        "cyp2d6",
        sam_path=path,
        profile=profile,
        output_file=None,
        reference=reference,
        solver=solver,
    )
    assert sorted(data) == sorted([s.diplotype for s in sols])
