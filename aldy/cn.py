#!/usr/bin/env python
# 786

# Aldy source: cn.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from __future__ import division
from builtins import zip
from builtins import range

import collections
import copy
import os

from . import lpinterface

from .common import *


def get_difference(coverage, unique_regions):
	if '7.pce' in coverage and '6.pce' not in coverage:
		coverage['6.pce'] = 0
	if '7.1e' in coverage:
		wo = [(coverage['6.' + r], coverage['7.' + r], coverage['6.' + r] - coverage['7.' + r]) for r in unique_regions]
	else:
		wo = [(coverage['6.' + r], 0, coverage['6.' + r]) for r in unique_regions]
	return dict(list(zip(unique_regions, wo)))


def load_profile(profile_path):
	result = collections.defaultdict(lambda: collections.defaultdict(int))
	with open(profile_path) as f:
		for line in f:
			if line[0] == '#':
				continue
			line = line.strip().split()[1:]
			if not line[0].startswith('chr'):
				line[0] = 'chr' + line[0]
			result[line[0]][int(line[1])] = float(line[2])
	return result


def normalize_coverage(gene, sam, profile):
	sam_ref = sum(sam.cnv_coverage[i] for i in range(gene.cnv_region[1], gene.cnv_region[2]))
	cnv_ref = sum(profile[gene.cnv_region[0]][i] for i in range(gene.cnv_region[1], gene.cnv_region[2]))

	if round(sam_ref, 2) < 1e-2:
		raise Exception('CN-neutral region ({}) coverage is too low ({:.1f})'.format(gene.cnv_region, sam_ref))
	if round(cnv_ref, 2) < 1e-2:
		raise Exception('CN-neutral region ({}) profile coverage is too low ({:.1f})'.format(gene.cnv_region, cnv_ref))
	cnv_ratio = float(cnv_ref) / sam_ref
	log.debug('CNV factor: {} ({})', cnv_ratio, 1.0 / cnv_ratio)

	sam.location_cnv = {}
	sam.region_coverage = {}
	for k, r in sorted(gene.regions.items()):
		s = sum(sam.unfiltered_coverage[i] for i in range(r[0], r[1] + 1))
		p = sum(profile[gene.region[0]][i] for i in range(r[0], r[1] + 1))
		sam.location_cnv.update({
			i: cnv_ratio * float(sam.unfiltered_coverage[i]) / profile[gene.region[0]][i]
			if profile[gene.region[0]][i] > 0 else 0.0
			for i in range(r[0], r[1] + 1)
		})
		sam.region_coverage[k] = (cnv_ratio * float(s) / p) if p != 0 else 0.0

def estimate_cn(gene, sam, profile, solution, solver):
	log.debug('== Copy Number Estimator ==')

	if solution is not None:
		solution = solution.split(',')
		solution = collections.Counter(solution)
		result = []
		res_cn = collections.defaultdict(int)
		log.debug('  Solutions (user-provided):')
		for sol in solution:
			if sol not in gene.cnv_configurations:
				raise Exception(
					('Given copy number solution contains unknown copy number configuration {}. ' + 
					 'Please run aldy --show-cn to list the valid configurations').format(sol)
				)
			for i in range(solution[sol]):
				result.append((sol, i))
				log.debug('    {}', (sol, i))
				for r in gene.cnv_configurations[sol]:
					res_cn[r] += gene.cnv_configurations[sol][r]
		return [(tuple(result), res_cn)]

	profile_path = script_path('aldy.resources.profiles', '{}.profile'.format(profile.lower()))
	if os.path.exists(profile_path):
		profile = load_profile(profile_path)
	else:
		profile = load_profile(profile)
	normalize_coverage(gene, sam, profile)

	coverage_diff = get_difference(sam.region_coverage, gene.unique_regions)

	log.debug('Coverage per region')
	log.debug('  {:5}: {} {}', 'region', ' 2D6', ' 2D7')

	for r in sorted(sam.region_coverage, key=region_sort_key):
		if r[0] != '6':
			continue
		if '7' + r[1:] in sam.region_coverage:
			log.debug(
				'  {:5}: {:5.2f} {:5.2f} {}', r[2:],
				sam.region_coverage['6' + r[1:]],
				sam.region_coverage['7' + r[1:]],
				'= {:5.2f}'.format(coverage_diff[r[2:]][2]) if r[2:] in coverage_diff else ''
			)
		else:
			log.debug('  {:5}: {:5.2f}', r[2:], sam.region_coverage['6' + r[1:]])

	structures = {(name, 0): structure for name, structure in gene.cnv_configurations.items()}

	# create special allele for any possible copy
	max_observed_cn = int(max(round(y) for r, y in sam.region_coverage.items())) + 1
	log.debug('Max. CN = {}', max_observed_cn)

	for a, _ in list(structures.keys()): # (a, -1) = a^ in old notation
		structures[(a, -1)] = structures[(a, 0)]
	for a, ai in list(structures.keys()):
		if ai == -1:
			continue
		if a == gene.deletion_allele:
			continue
		for i in range(1, max_observed_cn):
			structures[(a, i)] = copy.deepcopy(structures[(a, 0)])
			for r in gene.regions:
				if r.startswith('7.'):
					structures[(a, i)][r] = structures[(a, i)][r] - 1

	c = lpinterface.model('aldy_cnv', solver)

	# add one binary variable to model for any allele copy
	A = {(a, ai): c.addVar(vtype='B', name='A_{}_{}'.format(a, ai)) for a, ai in structures}

	pseudo_inducing = c.quicksum(A[a] for a in A if a[1] <= 0)
	c.addConstr(pseudo_inducing == 2)

	for a, ai in structures:
		if ai == 0 and (a, -1) in structures:
			log.trace('M  Sentinel contraint for {}: A_{}_-1 <= A_{}', a, a, a)
			c.addConstr(A[(a, -1)] <= A[(a, 0)])

	for a, ai in structures:
		if ai != 0:
			continue
		if a == gene.deletion_allele:
			continue
		for i in range(1, max_observed_cn):
			log.trace('M  Sentinel contraint for {}: A_{} <= A_{}', a, i, i - 1)
			c.addConstr(A[(a, i)] <= A[(a, i - 1)])

	E = {}
	for r, exp_cov in coverage_diff.items():
		exp_cov6, exp_cov7, exp_cov_diff = exp_cov
		expr = 0
		for a in structures:
			expr += A[a] * structures[a]['6.' + r]
			if '7.' + r in structures[a]:
				expr -= A[a] * structures[a]['7.' + r]
		E[r] = c.addVar(name=r, lb=-10, ub=10)

		log.trace('M  Contraint for {}: {} == E_{} + {}', r, exp_cov_diff, r, expr)
		c.addConstr(expr + E[r] == exp_cov_diff)

	objective = c.abssum(E.values(), coeffs={'pce': 1.5}) + .5 * sum(A.values())

	fusion_cost = c.quicksum(
		A[a] for a in A
		if a[0] in gene.fusions_left and a[0] != gene.deletion_allele
	)
	objective += .1 * fusion_cost
	log.trace('M  Objective: {}', objective)

	status, opt, solutions = c.solveAll(objective, A)
	log.debug('CN Solver status: {}, opt: {}', status, opt)

	log.debug('  Solutions:')
	for s in solutions:
		log.debug('    {}', s)
	for si, s in enumerate(solutions):
		res = collections.defaultdict(int)
		for sx in s:
			for r in structures[sx]:
				res[r] += structures[sx][r]
		solutions[si] = (s, res)

	return solutions
