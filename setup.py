#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

exec(open("aldy/version.py").read())

setup(
    name="aldy",
    version=__version__,
    description="A tool for allelic decomposition (haplotype reconstruction) "
    + "and exact genotyping of highly polymorphic and structurally variant genes",
    url="http://aldy.csail.mit.edu/",
    author="Ibrahim Numanagić",
    author_email="inumanag@mit.edu",
    download_url="https://github.com/inumanag/aldy/tarball/master",
    license="Aldy/IURTC License. Free for academic/non-commercial use.",
    keywords=["cyp2d6", "adme", "genotyping", "illumina", "pgrnseq", "getrm", "allele"],
    install_requires=[
        "pyyaml",
        "logbook",
        "pysam",
        "pytest",
        "ortools",
        "natsort",
        "mappy",
        "indelpost @ git+https://github.com/0xTCG/indelpost",
    ],
    entry_points={"console_scripts": ["aldy = aldy.__main__:console"]},
    packages=find_packages(),
    package_data={
        "aldy.resources": ["*.rst", "aldy/resources/*.rst"],
        "aldy.resources.genes": ["*.yml", "aldy/resources/genes/*.yml"],
        "aldy.resources.profiles": ["*.yml", "aldy/resources/profiles/*.yml"],
        "aldy.tests.resources": [
            "*.json",
            "aldy/tests/resources/*.json",
            "*.bai",
            "aldy/tests/resources/*.bai",
            "*.bam",
            "aldy/tests/resources/*.bam",
            "*.expected",
            "aldy/tests/resources/*.expected",
            "*.yml",
            "aldy/tests/resources/*.yml",
            "*.profile",
            "aldy/tests/resources/*.profile",
            "*.hg38",
            "aldy/tests/resources/*.hg38",
            "*.vcf",
            "aldy/tests/resources/*.vcf",
            "*.gz",
            "aldy/tests/resources/*.gz",
            "*.tbi",
            "aldy/tests/resources/*.tbi",
        ],
    },
    test_suite="pytest-runner",
    tests_require=["pytest"],
)
