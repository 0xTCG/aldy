#%%
import collections
import os, sys, os.path, glob
import statistics
from collections import *
from pprint import pprint
from unittest import result
from natsort import natsorted
from imp import reload
from multiprocessing import Pool
import tqdm
from tqdm.contrib.concurrent import process_map
from aldy import coverage
import aldy.gene, aldy.coverage, aldy.sam, aldy.common
# reload(sys.modules['aldy.gene'])
# reload(sys.modules['aldy.coverage'])
# reload(sys.modules['aldy.sam'])
# reload(sys.modules['aldy.common'])
from aldy.gene import CNConfig, CNConfigType, Gene
from aldy.common import GRange, log
from aldy.coverage import Coverage
from aldy.sam import Sample, DEFAULT_CN_NEUTRAL_REGION
import logbook, logbook.more
import pysam, bisect, tempfile, subprocess as sp
import multiprocess
from contextlib import contextmanager
from time import time
import aldy.cn, aldy.lpinterface, copy
from tabulate import tabulate
from statistics import mean, median
import functools
import mappy

@contextmanager
def timing(description: str) -> None:
    start = time()
    yield
    ellapsed_time = time() - start

    print(f"{description}: {ellapsed_time}")


gene = Gene('aldy/resources/genes/cyp2d6.yml', genome='hg38')
aldy.sam.DEFAULT_CN_NEUTRAL_REGION['hg38']
cn_region = GRange('1', 169511951, 169586630)
gene.get_wide_region()

# HG00332', 'ABCG2
# chr4	88090264	88231626

def idx(l, i):
    try:
        return l.index(i)
    except ValueError:
        return -1
def greedy_set_cover(subsets): # subsets: {id: {elems}, ...}
  # https://stackoverflow.com/questions/21973126/set-cover-or-hitting-set-numpy-least-element-combinations-to-make-up-full-set

  import heapq

  parent_set = set(j for i in subsets.values() for j in i)
  max = len(parent_set)
  # create the initial heap. Note 'subsets' can be unsorted,
  # so this is independent of whether remove_redunant_subsets is used.
  heap = []
  for se, s in subsets.items():
    # Python's heapq lets you pop the *smallest* value, so we
    # want to use max-len(s) as a score, not len(s).
    # len(heap) is just proving a unique number to each subset,
    # used to tiebreak equal scores.
    heapq.heappush(heap, [max-len(s), se, set(s)])
  results = []
  result_set = set()
  while result_set < parent_set:
    best = []
    unused = []
    while heap:
      score, count, s = heapq.heappop(heap)
      if not best:
        best = [max-len(s - result_set), count, s]
        continue
      if score >= best[0]:
        # because subset scores only get worse as the resultset
        # gets bigger, we know that the rest of the heap cannot beat
        # the best score. So push the subset back on the heap, and
        # stop this iteration.
        heapq.heappush(heap, [score, count, s])
        break
      score = max-len(s - result_set)
      if score >= best[0]:
        unused.append([score, count, s])
      else:
        unused.append(best)
        best = [score, count, s]
    add_set = best[2]
    results.append((best[1], add_set))
    result_set.update(add_set)
    while unused:
      heapq.heappush(heap, unused.pop())
  return results


#%% (1) Construct new reference
seq = next(mappy.fastx_read("temp/chr22_small.fa"))[1]
off = 42_000_000
genes = {
    '2D8': ( # 42149886-42155001
        seq[42149886 - off:42155001 - off],
        {'2D8': (-1, 0, 42155001-42149886)},
        None
    )
}
for gi, g in enumerate(gene.regions):
    st = min(r.start for r in g.values())
    ed = max(r.end for r in g.values())
    genes[str(gi)] = (
        seq[st - off:ed - off],
        {rn: (gi, r.start - st, r.end - r.start) for rn, r in sorted(g.items(), key=lambda x: x[1].start)},
        None
    )
fusion_signatures = {}
for n, cn in gene.cn_configs.items():
    if cn.kind == aldy.gene.CNConfigType.LEFT_FUSION:
        sig, pos = 0, next(i for i, j in cn.cn[0].items() if j == 1)
    elif cn.kind == aldy.gene.CNConfigType.RIGHT_FUSION:
        sig, pos = 1, next(i for i, j in cn.cn[0].items() if j == 0)
    else:
        continue
    sig = ((1 - sig, list(cn.cn[0])[list(cn.cn[0]).index(pos) - 1]), (sig, pos))
    fusion_signatures[sig] = n

    gi = sig[1][0]
    n_seq = []
    n_rn = {}
    n_sz = 0
    for rn in genes['0'][1]:
        if rn == sig[0][1]: # switch!
            gi = sig[0][0]
        s, rr, _ = genes[str(gi)]
        n_seq.append(s[rr[rn][1]:rr[rn][1] + rr[rn][2]])
        n_rn[rn] = (gi, n_sz, rr[rn][2])
        n_sz += rr[rn][2]
    genes[n] = (''.join(n_seq), n_rn, sig[0][1])
with tempfile.TemporaryDirectory() as tmp:
    with open(f"{tmp}/test.fa", "w") as fo:
        for n, (s, _, _) in genes.items():
            print(f">{n}\n{s}", file=fo)
    index = mappy.Aligner(f"{tmp}/test.fa", preset="map-hifi", best_n=len(genes) + 1)


#%% (2) Remap a file
def get_regions(r_start, _, cigar):
    regs = []
    start, s_start = r_start, 0
    for size, op in cigar:
        if op == 1 or op == 4: s_start += size
        elif op == 2 or op == 3: start += size
        elif op == 0 or op == 7 or op == 8:
            for i in range(size):
                rg = gene.region_at(start + i)
                if not rg and 42_149_886 <= start + i <= 42_155_001:
                    rg = (2, '2D8')
                if rg:
                    if not regs or regs[-1][0] != rg:
                        regs.append([rg, s_start + i, 0])
                    regs[-1][2] += 1
            start += size
            s_start += size
    return regs
def foo(x): # Print read: regions
    dc = []
    for (g, i), _, _ in x:
        if not dc or dc[-1][0] != g:
            dc.append((g, []))
        dc[-1][1].append(i)
    return '; '.join(f"{k}@[" + ','.join(v) + "]" for k, v in dc)
reload(sys.modules['aldy.sam'])
def remap_mappy(*args):
    def best_map(idx, seq):
        for hit in idx.map(seq):
            if hit.is_primary:
                return hit
        return None
    def adjust(seq, p): # convert fusion gene coordinate to 6 or 7 coordinate
        if seq == '2D8': return p+42149886
        for n, (g, st, sz) in genes[seq][1].items():
            if st <= p < st + sz:
                off = gene.regions[g][n].start
                return off + (p - st)
        assert False
    gene, file = args[0]
    sample = os.path.basename(file).split('.')[0]
    slices = []
    for i in range(2):
        slices.append(min(gene.regions[i].values(), key=lambda x: x.start).start)
        slices.append(max(gene.regions[i].values(), key=lambda x: x.end).end)
    reads = {}
    with pysam.AlignmentFile(file) as sam, timing(sample):
        iter = sam.fetch(region=GRange(gene.get_wide_region()[0], slices[0], 42_155_000).samtools(prefix='chr'))
        for read in iter:
            if not read.cigartuples or "H" in read.cigarstring:  # only valid alignments
                continue
            if read.reference_start >= slices[-1] or read.reference_end < slices[0]:
                continue
            seq = read.query_sequence
            mapped = 0
            pieces = []
            max_len = 2500
            while mapped < len(seq):
                h = best_map(index, seq[mapped:mapped+max_len])
                if not h:
                    break
                if not genes[h.ctg][2]:
                    pieces.append(( adjust(h.ctg, h.r_st), seq[mapped + h.q_st:mapped + h.q_en], h.cigar ))
                else:
                    brk = genes[h.ctg][1][genes[h.ctg][2]][1] # breakpoint
                    splits = [[adjust(h.ctg, h.r_st), [], mapped + h.q_st, 0]]  # ref start, cigar, que start, que len
                    ref_start, seq_start = h.r_st, h.q_st
                    for sz, op in h.cigar:
                        if op == 1 or op == 4:
                            seq_start += sz
                            splits[-1][3] += sz
                        elif op == 2 or op == 3:
                            if ref_start <= brk < ref_start + sz:
                                splits[-1][1].append((brk - ref_start, op))
                                splits.append([adjust(h.ctg, brk),
                                    [(ref_start + sz - brk, op)], mapped + seq_start, 0])
                            else:
                                splits[-1][1].append((sz, op))
                            ref_start += sz
                        elif op == 0 or op == 7 or op == 8:
                            if ref_start <= brk < ref_start + sz:
                                splits[-1][1].append((brk - ref_start, op))
                                splits[-1][3] += brk - ref_start
                                splits.append([adjust(h.ctg, brk),
                                    [(ref_start + sz - brk, op)], mapped + seq_start + brk - ref_start, ref_start + sz - brk])
                            else:
                                splits[-1][1].append((sz, op))
                                splits[-1][3] += sz
                            ref_start += sz
                            seq_start += sz
                    pieces += [
                        (st, seq[ss:ss+se], cigar)
                        for [st, cigar, ss, se] in splits
                    ]
                mapped += h.q_en
            if pieces:
                reads[read.query_name] = pieces
        s = aldy.sam.Sample(
            gene,
            file,
            minimap=reads,
            sample_name=sample,
            cn_region=cn_region
        )
        snorm = aldy.sam.Sample(
            gene,
            file,
            profile='illumina',
            cn_region=cn_region
        )
    cov = {}
    for p in range(cn_region.start, cn_region.end):
        cov[(-1, 'cn'), p] = s.coverage._cnv_coverage[p]
    for gi, gr in enumerate(gene.regions):
        for region, rng in gr.items():
            for p in range(rng.start, rng.end):
                cov[(gi, region), p] = int(s.coverage.total(p))
    return sample, (cov, reads, s, snorm)

pool = multiprocess.Pool(8)
data = [(gene, i) for i in sorted(glob.glob("temp/pacbio/???????.bam"))]
with timing("Map"):
    results = dict(pool.map(remap_mappy, data))
    # r, s, dx = remap_mappy((gene, "temp/pacbio/HG00437.bam"))
print('done')

#%%
from matplotlib import pyplot as plt
from functools import partial
def bin_sample(ref, sample):
    # profile_cn = sum(results[sample][0].get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    # sample_cn = sum(results[sample][0].get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    # cn_ratio = 2 * float(profile_cn) / sample_cn
    obins = [0 for i in range(11)]
    nbins = [0 for i in range(11)]
    def filter_fns(mut, cov, total, thres, mc):
        z = Coverage.basic_filter( mut, cov, total, thres / 10, mc )
        return z
    ocov = results[sample][3].coverage.filtered(partial(filter_fns, mc=results[sample][3].coverage.min_cov))
    ncov = results[sample][2].coverage.filtered(partial(filter_fns, mc=results[sample][2].coverage.min_cov))
    for (pos, _) in gene.mutations:
        if pos in ocov._coverage:
            otot = ocov.total(pos)
            # factor = cn_ratio * tot / results[sample][2].coverage.total(pos)
            for i in ocov._coverage[pos].values():
                obins[int(10*i/otot)]+=1
            ntot = ncov.total(pos)
            for i in ncov._coverage[pos].values():
                nbins[int(10*i/ntot)]+=1
                # bins[ int(10*factor*i/tot) ] += 1
            # pd[pos] = sorted([ ], reverse=True)
    return obins[:10], nbins[:10]
for ri, r in enumerate(results):
    o,n = bin_sample('HG00332', r)
    ax = plt.subplot(len(results), 2, ri*2+1)
    plt.plot([i/10 for i in range(len(o))], o)
    ax.title.set_text(r)
    plt.subplot(len(results), 2, ri*2+2)
    plt.plot([i/10 for i in range(len(n))], n)
plt.gcf().set_figwidth(8)
plt.gcf().set_figheight(40)
# plt.ylim(0, 50)
# HG00332


#%% Fusion calculation & Copy number calculation
# HG01089 *1/*4 / 1/5
# HG02649 *1/*2
# NA12762 *1/*3
# NA19466 *1/*17
# HG02087 *36+*10/*36+*10  36/36+10
import aldy.major
reload(sys.modules['aldy.major'])
def get_fusions(reads):
    fusions = {fn: [0,0] for fn in fusion_signatures.values()}
    for rn, rp in reads.items():
        regs = [rg[0] for rr in rp for rg in get_regions(*rr)]
        for (s1, s2), fn in fusion_signatures.items():
            i = idx(regs, s2)
            if i != -1 and i + 1 != len(regs) and regs[i + 1] == s1:
                fusions[fn][0] += 1
                continue
            if i != -1 and i + 1 != len(regs) and regs[i + 1] == (1 - s1[0], s1[1]):
                fusions[fn][1] += 1
                continue
            i = idx(regs, (1 - s2[0], s2[1]))
            if i != -1 and i + 1 != len(regs) and regs[i + 1] == s1:
                fusions[fn][1] += 1
    return fusions
truth = {}
with open("temp/genotypes.dat") as f:
    for l in f:
        l = l.split()
        truth[l[0]] = l[2]
def plot(gene, ref, sample, plot=True):
    profile_cn = sum(ref.get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    sample_cn = sum(sample.get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    cn_ratio = 2 * float(profile_cn) / sample_cn
    data = {}
    for gi, gr in enumerate(gene.regions):
        for region, rng in gr.items():
            s = sum(sample.get(((gi, region), i), 0) for i in range(rng.start, rng.end))
            p = sum(ref.get(((gi, region), i), 0) for i in range(rng.start, rng.end))
            data[gi, region] = (cn_ratio * float(s) / p) if p != 0 else 0.0
    h = sorted(set(k for _, k in data), key=lambda x: -gene.regions[0][x].start)
    r = [[f'*{k}*' if k in gene.unique_regions else k for k in h]] + [[round(data.get((j, i)), 1) for i in h] for j in range(2)]
    r += [[round(r[1][i] - r[2][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
    r += [[round(r[1][i] / r[2][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
    tab = tabulate(r, tablefmt="plain", stralign="right", numalign="right")
    if plot: print(tab)
    return data

def model(rs):
    (ref, _, _), (sample, sample_rd, sample_sam) = results[rs[0]], results[rs[1]]
    data = plot(gene, ref, sample, False)
    max_cn = 1 + max(int(round(data.get((gi, r), 0))) for gi, g in enumerate(gene.regions) for r in g)
    f = get_fusions(sample_rd)
    # print({n: round(f[0] / f[1],2) for n, f in f.items() if f[0] / f[1] > .1})
    structures = {(name, 0): structure
        for name, structure in gene.cn_configs.items()
        if name not in f
        or f[name][0] / f[name][1] > .1
    }
    for a, ai in list(structures.keys()):
        structures[a, -1] = copy.deepcopy(structures[a, 0])
        if a not in gene.cn_configs or str(gene.cn_configs[a].kind) != str(aldy.gene.CNConfigType.DEFAULT):
            continue
        for i in range(1, max_cn):
            structures[a, i] = copy.deepcopy(structures[a, 0])
            for g in range(1, len(structures[a, i].cn)):
                structures[a, i].cn[g] = { r: v - 1 for r, v in structures[a, i].cn[g].items() }
    for i in range(max_cn):
        structures['2D7', i+1] = aldy.gene.CNConfig(
            copy.deepcopy(structures['5', 0].cn), aldy.gene.CNConfigType.DELETION,
            {}, "extra 2D7 copy"
        )
    region_cov = {
        r: (data.get((0, r), 0), data.get((1, r), 0) if len(gene.regions) > 1 else 0)
        for r in gene.regions[1]
    }
    ## BEGIN cn code
    VERR, VCN = {}, {}
    def setup_model():
        model = aldy.lpinterface.model("AldyCN", 'gurobi')
        VCN.clear()
        VCN.update({ (a, ai): model.addVar(vtype="B", name=f"CN_{a}_{ai}") for a, ai in structures })
        diplo_inducing = model.quicksum(VCN[a] for a in VCN if a[1] <= 0)
        model.addConstr(diplo_inducing <= 2, name="CDIPLO")
        model.addConstr(diplo_inducing >= 2, name="CDIPLO")
        del_allele = gene.deletion_allele()
        if del_allele:
            for (a, ai), v in VCN.items():
                if a != del_allele:
                    model.addConstr(v + VCN[del_allele, -1] <= 1, name=f"CDEL_{a}_{ai}")
        for a, ai in structures:
            if ai == -1:
                model.addConstr(VCN[a, ai] <= VCN[a, 0], name=f"CORD_{a}_{ai}")
            elif ai > 1:
                model.addConstr(VCN[a, ai] <= VCN[a, ai - 1], name=f"CORD_{a}_{ai}")
        return model
    def diff_model(model):
        VERR.clear()
        VEPR = {}

        CYP6_COEFF = 1 # model.addVar(name=f"C6C", lb=0.75, ub=1.25)
        for r, (exp_cov0, exp_cov1) in region_cov.items():
            expr = 0
            expr6 = 0
            for s, structure in structures.items():
                if r in structure.cn[0]:
                    expr += structure.cn[0][r] * VCN[s]
                    expr6 += structure.cn[0][r] * VCN[s]
                if len(structure.cn) > 1 and r in structure.cn[1]:
                    expr -= structure.cn[1][r] * VCN[s]

            if r not in gene.unique_regions: continue

            VEPR[r] = model.addVar(name=f"EX_{r}", lb=-aldy.cn.MAX_CN_ERROR, ub=aldy.cn.MAX_CN_ERROR)
            model.addConstr(expr6 + VEPR[r] <= CYP6_COEFF * exp_cov0, name=f"C6COV_{r}")
            model.addConstr(expr6 + VEPR[r] >= CYP6_COEFF * exp_cov0, name=f"C6COV_{r}")

            VERR[r] = model.addVar(name=f"E_{r}", lb=-aldy.cn.MAX_CN_ERROR, ub=aldy.cn.MAX_CN_ERROR)
            model.addConstr(expr + VERR[r] <= (exp_cov0 - exp_cov1), name=f"CCOV_{r}")
            model.addConstr(expr + VERR[r] >= (exp_cov0 - exp_cov1), name=f"CCOV_{r}")

        PCE_COEFF = 1.5
        DIFF_COEFF = 1
        objective = DIFF_COEFF * model.abssum(VERR.values(), coeffs={"E_pce": PCE_COEFF})
        FIT_COEFF = 0.5
        objective += FIT_COEFF * model.abssum(VEPR.values())

        # CYP6_ABS = model.addVar(lb=0, name=f"ABS_C6C")
        # model.addConstr(CYP6_ABS + (1 - CYP6_COEFF) >= 0, name=f"CABSL_C6C")
        # model.addConstr(CYP6_ABS - (1 - CYP6_COEFF) >= 0, name=f"CABSR_C6C")
        # objective += 5 * CYP6_ABS

        penalty = {s: aldy.cn.PARSIMONY_PENALTY for s, _ in VCN}
        for n, s in gene.cn_configs.items():
            if n in penalty and str(s.kind) == str(aldy.gene.CNConfigType.LEFT_FUSION):
                penalty[n] += aldy.cn.LEFT_FUSION_PENALTY
        penalty['2D7'] += aldy.cn.LEFT_FUSION_PENALTY
        # penalty
        objective += model.quicksum(penalty[s] * v for (s, _), v in VCN.items())
        model.setObjective(objective)
    # Solve the model
    model = setup_model()
    diff_model(model)
    # model.addConstr(VCN['1', 0] >= 1)
    # model.addConstr(VCN['1', -1] >= 1)
    # model.addConstr(VCN['2D7', 0] >= 1)
    # model.addConstr(VCN['36.ALDY', 0] >= 1)
    # model.addConstr(VCN['36.ALDY', -1] >= 1)
    lookup = {model.varName(v): a for (a, ai), v in VCN.items()}
    result: dict = {}
    found = {}
    best_opt = None
    for status, opt, sol in model.solutions(0.1): # gap = 0
        if best_opt is not None and abs(opt-best_opt) > 0.5: continue
        if best_opt is None: best_opt = opt
        sol_tuple = tuple(sorted(lookup[v] for v in sol if lookup[v] not in ['5', '2D7']))
        if sol_tuple not in result:
            result[sol_tuple] = aldy.cn.CNSolution(gene, opt, list(sol_tuple))
            s = tuple(natsorted( '36' if s == '36.ALDY' else s for s in sol_tuple ))
            found[result[sol_tuple]._solution_nice()] = opt - best_opt
    if False:
        q = [[c.region_cn[g][r] for r in c.gene.regions[0]] for g, _ in enumerate(c.region_cn)]
        h = sorted(set(k for _, k in data), key=lambda x: -gene.regions[0][x].start)
        r = [[f'*{k}*' if k in gene.unique_regions else k for k in h]] + \
            [[round(data.get((j, i)), 1) for i in h] for j in range(2)]
        r += [[round(r[1][i] - r[2][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
        r += [[round(r[1][i] / r[2][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
        r += [["---" for i in range(len(r[0]))]]
        r += [[c.region_cn[g].get(r, 0) for r in h] for g, _ in enumerate(c.region_cn)]
        r += [[round(r[-2][i] / r[-1][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
        r += [["---" for i in range(len(r[0]))],
                [round(abs(r[-1][i] - r[-5][i]),1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
        tab = tabulate(r, tablefmt="plain", stralign="right", numalign="right")
        print(s,opt); print(tab)
    print(rs[1], truth[rs[1]], 'vs', found, opt)
    # plot(gene, ref, sample, True)

    major_sols = []
    for i, cn_sol in enumerate(result.values()):
        major_sols += aldy.major.estimate_major(
            gene,
            sample_sam.coverage,
            cn_sol,
            solver='gurobi',
            identifier=i,
        )
    for m in major_sols:
        print(m.cn_solution._solution_nice(), found[m.cn_solution._solution_nice()], m.score, m._solution_nice())

    return found
# print(model(results['HG00332'], results['HG01089']))
# for i in sorted(results):
for i in ['HG00130']:
# for i in "HG01089 HG02087 HG02649 HG03624 NA12762 NA19452 NA19466".split():
    model(('HG00332', i))
    print()
# HG01089 *1/*4 vs {('1', '2D7'), ('1', '5')} 15.422513442593278



#%%
import matplotlib.pyplot as plt
def plot(gene, ref, sample, plot=True):
    profile_cn = sum(ref.get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    sample_cn = sum(sample.get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    cn_ratio = 2 * float(profile_cn) / sample_cn
    data = {}
    for gi, gr in enumerate(gene.regions):
        for region, rng in gr.items():
            s = sum(sample.get(((gi, region), i), 0) for i in range(rng.start, rng.end))
            p = sum(ref.get(((gi, region), i), 0) for i in range(rng.start, rng.end))
            data[gi, region] = (cn_ratio * float(s) / p) if p != 0 else 0.0
    h = sorted(set(k for _, k in data), key=lambda x: -gene.regions[0][x].start)
    r = [[f'*{k}*' if k in gene.unique_regions else k for k in h]] + [[round(data.get((j, i)), 1) for i in h] for j in range(2)]
    r += [[round(r[1][i] - r[2][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
    r += [[round(r[1][i] / r[2][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
    tab = tabulate(r, tablefmt="plain", stralign="right", numalign="right")
    if plot: print(tab)
    return r
r = plot(gene, results['HG00332'][0], results['NA19466'][0])

f = plt.figure()
f.set_figwidth(12)
f.set_figheight(5)
plt.plot(r[0], r[1], color='blue')
plt.plot(r[0], r[2], color='red')
plt.plot(r[0], [2 for _ in r[1]], linestyle='dashed')

#%%
# Remap the file

#%%
fusions = {}
for rn, rp in r.items():
    regs = [rg[0] for rr in rp for rg in get_regions(*rr)]
    for (s1, s2), fn in fusion_signatures.items():
        try:
            z = regs.index(s2)
            if z+1 < len(regs) and regs[z+1] == s1:
                fusions.setdefault(fn, set()).add(rn)
        except:
            continue
for f in fusions:
    print(f, len(fusions[f]))
print()

#%%



#%%

for fl in sorted(glob.glob("temp/pacbio/???????.bam")):
    with timing(f"Remap {fl}"):
        # rds = remap_mappy((gene, "temp/pacbio/HG00437.bam"))
        # rds, _, dx = remap_mappy((gene, "temp/pacbio/HG00185.bam"))
        rds, _, dx = remap_mappy((gene, fl))
        fus, ms = cluster(rds)
    fs = {y:x for x,y in fusion_signatures.items()}
    zz = dict(greedy_set_cover(fus))
    for n, s in zz.items():
        print(n, fs[n], round(100.0*len(s)/ms[n],0))
    print()




# %%
