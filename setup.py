#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

exec(open('aldy/version.py').read())

setup(
    name='aldy',
    entry_points={
        'console_scripts': [
            'aldy = aldy.__main__:main'
        ]
    },
    packages=find_packages(),
    package_data={
        'aldy.resources': [
            '*.md', 'aldy/resources/*.md', 
            '*.bam', 'aldy/resources/*.bam',
            '*.bai', 'aldy/resources/*.bai'],
    	'aldy.resources.genes': ['*.yml', 'aldy/resources/genes/*.yml'],
    	'aldy.resources.profiles': ['*.profile', 'aldy/resources/profiles/*.profile']
    },
    install_requires=['pyyaml', 'logbook', 'six', 'pysam', 'future'],
    version=__version__,
    description='A tool for allelic decomposition and exact genotyping of highly polymorphic and structurally variant genes',
    author='Ibrahim Numanagić',
    author_email='inumanag@mit.edu',
    url='http://aldy.csail.mit.edu/',
    download_url='https://github.com/inumanag/aldy/tarball/master',
    keywords=['cyp2d6', 'adme', 'genotyping', 'illumina', 'pgrnseq', 'getrm', 'allele']
)
