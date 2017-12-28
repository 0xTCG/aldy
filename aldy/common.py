#!/usr/bin/env python
# 786

# Aldy source: common.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from builtins import range

import pkg_resources 
import os
import re
import time
import pprint
import logbook


PROTEINS = { # X is stop
	'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V', 'TTC': 'F',
	'CTC': 'L', 'ATC': 'I', 'GTC': 'V', 'TTA': 'L', 'CTA': 'L',
	'ATA': 'I', 'GTA': 'V', 'TTG': 'L', 'CTG': 'L', 'ATG': 'M',
	'GTG': 'V', 'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
	'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'TCA': 'S',
	'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'TCG': 'S', 'CCG': 'P',
	'ACG': 'T', 'GCG': 'A', 'TAT': 'Y', 'CAT': 'H', 'AAT': 'N',
	'GAT': 'D', 'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
	'TAA': 'X', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'TAG': 'X',
	'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'TGT': 'C', 'CGT': 'R',
	'AGT': 'S', 'GGT': 'G', 'TGC': 'C', 'CGC': 'R', 'AGC': 'S',
	'GGC': 'G', 'TGA': 'X', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
	'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
}

REV_COMPLEMENT = {
	'A': 'T', 'T': 'A',
	'C': 'G', 'G': 'C',
}


def inspect():
	import code
	code.interact(local=locals())


def sorted_tuple(x):
	return tuple(sorted(x))


def allele_key(x):
	p = re.split(r'(\d+)', x)
	return p[1]


def region_sort_key(x):
	return tuple(int(z) if z.isdigit() else z for z in re.split('(\d+)', x))


def sort_key(x):
	p = re.split(r'(\d+)', x)
	return (int(p[1]),) + tuple(p[2:])


def var_to_allele(x):
	return re.split(r'[\^_]', x)[0]


def rev_comp(seq):
	return ''.join([REV_COMPLEMENT[x] for x in seq[::-1]])


def seq_to_amino(seq):
	return ''.join(PROTEINS[seq[i:i + 3]] for i in range(0, len(seq) - len(seq) % 3, 3))


def check_functional(gene, m):
	if m.pos not in gene.coding_region:
		return False
	if m.op[:3] != 'SNP':
		return True
	seq = ''.join(gene.coding_region[i] if i != m.pos else m.op[5] for i in sorted(gene.coding_region.keys()))
	amino = seq_to_amino(rev_comp(seq) if gene.rev_comp else seq)
	return amino != gene.aminoacid


def static_vars(**kwargs):
	def decorate(func):
		for k in kwargs:
			setattr(func, k, kwargs[k])
		return func
	return decorate


def timing(f):
	def wrap(*args):
		time1 = time.time()
		ret = f(*args)
		time2 = time.time()
		log.debug('Time needed: ({:.1f})', time2 - time1)
		return ret
	return wrap


@static_vars(pp=pprint.PrettyPrinter(indent=4))
def pr(x):
	return pr.pp.pformat(x)


def script_path(path, file):
	return pkg_resources.resource_filename(path, file)


def colorize(text, color='green'):
	return logbook._termcolors.colorize(color, text)


log = logbook.Logger('Aldy')
LOG_FORMAT = '{record.message}'
