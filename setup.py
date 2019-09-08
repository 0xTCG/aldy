#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

exec(open('aldy/version.py').read())

setup(name='aldy',
      version=__version__,
      description='A tool for allelic decomposition and exact genotyping of highly polymorphic and structurally variant genes',
      url='http://aldy.csail.mit.edu/',
      author='Ibrahim Numanagić',
      author_email='inumanag@mit.edu',
      download_url='https://github.com/inumanag/aldy/tarball/master',
      license='Aldy/IURTC License. Free for academic/non-commercial use.',
      keywords=['cyp2d6', 'adme', 'genotyping', 'illumina', 'pgrnseq', 'getrm', 'allele'],
      
      install_requires=['pyyaml', 'logbook', 'pysam', 'pytest', 'ortools'],
      entry_points={
         'console_scripts': ['aldy = aldy.__main__:main' ]
      },
      
      packages=find_packages(),
      package_data={
        'aldy.resources': [
            '*.md',  'aldy/resources/*.md', 
            '*.bam', 'aldy/resources/*.bam',
            '*.bai', 'aldy/resources/*.bai'],
         'aldy.resources.genes': ['*.yml', 'aldy/resources/genes/*.yml'],
         'aldy.resources.profiles': ['*.profile', 'aldy/resources/profiles/*.profile']
      },

      test_suite='nose.collector',
      tests_require=['nose'],
)
