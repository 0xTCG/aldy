from __future__ import print_function
# 786

# Aldy source: diplotype.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from builtins import map
from collections import defaultdict, Counter

import sys

from .common import *


def write_decomposition(sol_id, decomposition, genotype, f):
	print('# Solution {}'.format(sol_id), file=f)
	print('# Predicted diplotype: {}'.format(genotype[0]), file=f)
	print('# Composition: {}'.format(','.join(genotype[1])), file=f)
	print('\t'.join('Copy Allele Location Type Coverage Effect dbSNP Code Status'.split()), file=f)
	for i, (allele, mutations) in enumerate(decomposition):
		for t, m, c in sorted(mutations, key=lambda x: x[1].pos):
			print('\t'.join(map(str, [
				i, allele[0] if '/' in allele[0] else allele[1],
				m.pos, m.op, c,
				['NEUTRAL', 'DISRUPTING'][0 if m.functional == 0 else 1],
				m.aux['dbsnp'], m.aux['old'],
				t
			])), file=f)
	print(file=f)


def assign_diplotype(gene, solutions, f):
	print('# Aldy v1.0', file=f)
	print('# Gene: {}'.format(gene.name), file=f)
	print('# Number of solutions: {}'.format(len(solutions)), file=f)
	print(file=f)

	results = []
	for sol_id, (orig_sol, decomposition) in enumerate(solutions):
		# solution is the array of (major, minor) tuples
		sol = [allele_key(x[0]) for x in orig_sol]
		if len(sol) == 1:
			res = '*{}/*{}'.format(gene.deletion_allele, *sol)
		elif len(sol) == 2:
			res = '*{}/*{}'.format(*sol)
		else:
			sol = defaultdict(int, Counter(sol))
			diplotype, dc = [[], []], 0
			# Handle tandems
			for ta, tb in gene.common_tandems:
				while sol[ta] > 0 and sol[tb] > 0:
					diplotype[dc % 2] += [ta, tb]
					dc += 1
					sol[ta] -= 1
					sol[tb] -= 1
			# Handle duplicates
			for allele, count in sol.items():
				if count > 0:
					diplotype[dc % 2] += count * [allele]
					dc += 1
			if len(diplotype[1]) == 0:
				diplotype[1] = [diplotype[0][-1]]
				diplotype[0] = diplotype[0][:-1]
			if len(diplotype[0]) > len(diplotype[1]):
				diplotype[0], diplotype[1] = diplotype[1], diplotype[0]
			res = '/'.join(['+'.join(['*' + y for y in x]) for x in diplotype])

		minor_sol = []
		for al in orig_sol:
			if '/' in al[0]:
				al = '{}|{}'.format(*al)
			elif '+' in al[0]:
				al = al[1] + '+'
			else:
				al = al[1]
			minor_sol.append(al)
		if len(minor_sol) == 1:
			minor_sol.append(gene.deletion_allele)

		results.append((res, minor_sol))
		write_decomposition(sol_id, decomposition, (res, minor_sol), f)

	return results
