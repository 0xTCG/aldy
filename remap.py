#%%
import os, sys, os.path, glob
from pickle import FALSE
from collections import *
from pprint import pprint
from natsort import natsorted
from imp import reload
import tqdm
from tqdm.contrib.concurrent import process_map
import aldy.gene, aldy.coverage, aldy.sam, aldy.common, aldy.cn, aldy.lpinterface
import logbook, logbook.more
import pysam, bisect, tempfile, subprocess as sp
import multiprocess
from contextlib import contextmanager
from time import time
import copy
from tabulate import tabulate
from statistics import mean, median
import functools
import mappy
from matplotlib import pyplot as plt

@contextmanager
def timing(description: str) -> None:
    start = time()
    yield
    ellapsed_time = time() - start

    print(f"{description}: {ellapsed_time}")

gene = aldy.gene.Gene('aldy/resources/genes/cyp2d6.yml', genome='hg38')
cn_region = aldy.common.GRange('1', 169511951, 169586630)

def best_map(idx, seq):
    for hit in idx.map(seq):
        if hit.is_primary:
            return hit
    return None

#%% (1) Construct new reference
full_index = mappy.Aligner("temp/chr22_small.fa", preset="map-hifi")
off = 42_000_000
index_gene = []
for gi, g in enumerate(gene.regions):
    st = min(r.start for r in g.values())
    ed = max(r.end for r in g.values())
    index_gene.append((
        mappy.Aligner(seq=full_index.seq(full_index.seq_names[0])[st - off:ed - off], preset="map-hifi"), st
    ))
fusion_signatures = {}
for n, cn in gene.cn_configs.items():
    if cn.kind == aldy.gene.CNConfigType.LEFT_FUSION:
        sig, pos = 0, next(i for i, j in cn.cn[0].items() if j == 1)
    elif cn.kind == aldy.gene.CNConfigType.RIGHT_FUSION:
        sig, pos = 1, next(i for i, j in cn.cn[0].items() if j == 0)
    else:
        continue
    sig = ((1 - sig, list(cn.cn[0])[list(cn.cn[0]).index(pos) - 1]), (sig, pos))
    if sig[0][1] in gene.unique_regions or sig[1][1] in gene.unique_regions:
        fusion_signatures[sig] = n
    else:
        print('skip', n)
fusion_signatures

#%% (2) Remap a file
reload(sys.modules['aldy.sam'])
def remap_mappy(*args):
    gene, file = args[0]
    sample = os.path.basename(file).split('.')[0]
    def get_regions(r_start, cigar):
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
    def fix_piece(seq, h, cnt, cnt2):
        regs = get_regions(h.r_st + off - 1, h.cigar)
        if regs:
            fx= []
            for fs, fn in fusion_signatures.items():
                # DO BOTH sides and see does it match at all...
                br = [r for r in regs if r[0][1] == fs[1][1]]
                if not br: continue
                cnt2[fn]+=1
                br = br[0]
                lg, br, rg = fs[1][0], br[1]+br[2], fs[0][0]
                # print(regs)
                lh = best_map(index_gene[lg][0], seq[:br])
                rh = best_map(index_gene[rg][0], seq[br:])
                if lh and rh:
                    unmap = len(seq) - (lh.q_en-lh.q_st) - (rh.q_en-rh.q_st)
                    unmap += lh.NM + rh.NM
                    fx.append([fn, unmap, br, lh, rh, lg, rg])
            fx.sort(key=lambda x: x[1])
            if fx and fx[0][1] <= h.NM:
                fn, _, br, lh, rh, lg, rg = fx[0]
                cnt[fn] += 1
                return [
                    (lh.r_st + index_gene[lg][1]-1, seq[lh.q_st:lh.q_en], lh.cigar),
                    (rh.r_st + index_gene[rg][1]-1, seq[br + rh.q_st:br + rh.q_en], rh.cigar)
                ], max(lh.q_en, br + rh.q_en)
        cnt['NORM'] += 1
        return [(h.r_st + off - 1, seq[h.q_st:h.q_en], h.cigar)], h.q_en
    reads = {}
    cnt = defaultdict(int)
    cnt2 = defaultdict(int)
    with pysam.AlignmentFile(file) as sam, timing(sample):
        r = [*gene.get_wide_region()]
        r = aldy.common.GRange(r[0], r[1], 42_155_000)
        iter = sam.fetch(region=r.samtools(prefix='chr'))
        for read in iter:
            if not read.cigartuples or "H" in read.cigarstring:  # only valid alignments
                continue
            if read.reference_start >= r.end or read.reference_end < r.start:
                continue
            seq = read.query_sequence
            pieces = []
            mapped = 0
            while mapped < len(seq):
                s = seq[mapped:]
                h = best_map(full_index, s)
                if h and h.r_st <= 134_991 < h.r_en and 134_991 - h.r_st > 100:
                    s = seq[mapped:mapped + (134_991 - h.r_st)]
                    h = best_map(full_index, s)
                elif h and h.r_st <= 148_730 < h.r_en and 148_730 - h.r_st > 100:
                    s = seq[mapped:mapped + (148_730 - h.r_st)]
                    h = best_map(full_index, s)
                if not h:
                    mapped += 100
                else:
                    pcs, ed = fix_piece(s, h, cnt,cnt2)
                    pieces += pcs; mapped += ed
            if pieces:
                reads[read.query_name] = (pieces, len(seq))
        s = aldy.sam.Sample(
            gene,
            file,
            minimap=reads,
            sample_name=sample,
            cn_region=cn_region
        )
    cov = {}
    for p in range(cn_region.start, cn_region.end):
        cov[(-1, 'cn'), p] = s.coverage._cnv_coverage[p]
    for gi, gr in enumerate(gene.regions):
        for region, rng in gr.items():
            for p in range(rng.start, rng.end):
                cov[(gi, region), p] = int(s.coverage.total(p))
    return sample, (cov, reads, s, (cnt,cnt2))
# zp = remap_mappy((gene, "temp/pacbio/HG00437.bam"))[1]
# results = {'HG00437': zp}
pool = multiprocess.Pool(8)
data = [(gene, i) for i in sorted(glob.glob("temp/pacbio/???????.bam"))]
with timing("Map"):
    results = dict(pool.map(remap_mappy, data))
print('done')
# zp=results['HG02087']
# d = []
# for gi, g in enumerate(gene.regions):
#     st = min(r.start for r in g.values())
#     ed = max(r.end for r in g.values())
#     d.append([zp[2].coverage.total(i) for i in range(st,ed)])
# for gi, g in enumerate(gene.regions):
#     st = min(r.start for r in g.values())
#     ed = max(r.end for r in g.values())
#     d.append([zp[3].coverage.total(i) for i in range(st,ed)])
# d[0] = [0]*1600+d[0]
# d[2] = [2]*1600+d[2]
# for dx in d:
#     plt.plot(dx)

#%%
truth = {}
with open("temp/genotypes.dat") as f:
    for l in f:
        l = l.split()
        truth[l[0]] = l[2]
def bin_sample(ref, sample):
    # profile_cn = sum(results[sample][0].get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    # sample_cn = sum(results[sample][0].get(((-1, 'cn'), i), 0) for i in range(cn_region.start, cn_region.end))
    # cn_ratio = 2 * float(profile_cn) / sample_cn
    obins = [0 for i in range(11)]
    nbins = [0 for i in range(11)]
    def filter_fns(mut, cov, total, thres, mc):
        z = aldy.coverage.Coverage.basic_filter( mut, cov, total, thres / 10, mc )
        return z
    ocov = results[sample][2].coverage.filtered(functools.partial(filter_fns, mc=results[sample][2].coverage.min_cov))
    ncov = results[sample][2].coverage.filtered(functools.partial(filter_fns, mc=results[sample][2].coverage.min_cov))
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
rs = results # {'HG00437': results['HG00437']}
for ri, r in enumerate(rs):
    o,n = bin_sample('HG00332', r)
    ax = plt.subplot(len(rs), 2, ri*2+1)
    ax.title.set_text(r)
    plt.plot([i/10 for i in range(len(o))], o)
    ax = plt.subplot(len(rs), 2, ri*2+2)
    ax.title.set_text(truth[r])
    plt.plot([i/10 for i in range(len(n))], n)
plt.gcf().set_figwidth(8)
plt.gcf().set_figheight(3 * len(rs))

#%%
def hoho(sx):
    structures = {(name, 0): structure for name, structure in gene.cn_configs.items()}
    for a, ai in list(structures.keys()):
        structures[a, -1] = copy.deepcopy(structures[a, 0])
        if a not in gene.cn_configs or str(gene.cn_configs[a].kind) != str(aldy.gene.CNConfigType.DEFAULT):
            continue
        for i in range(1, 6):
            structures[a, i] = copy.deepcopy(structures[a, 0])
            for g in range(1, len(structures[a, i].cn)):
                structures[a, i].cn[g] = { r: v - 1 for r, v in structures[a, i].cn[g].items() }
    VCN = {}
    model = aldy.lpinterface.model("AldyCN", 'gurobi')
    VCN.update({ (a, ai): model.addVar(vtype="B", name=f"CN_{a}_{ai}") for a, ai in structures if a == '1' })
    for a, ai in VCN:
        if ai == -1:
            model.addConstr(VCN[a, ai] <= VCN[a, 0], name=f"CORD_{a}_{ai}")
        elif ai > 1:
            model.addConstr(VCN[a, ai] <= VCN[a, ai - 1], name=f"CORD_{a}_{ai}")
    VVAF = {}
    ocov = results[sx][2].coverage.filtered(
        lambda mut, cov, total, thres: aldy.coverage.Coverage.basic_filter(
            mut, cov, total, thres / 10, results[sx][2].coverage.min_cov
        )
    )
    for (pos, _) in gene.mutations:
        reg = gene.region_at(pos)
        if reg and pos in ocov._coverage and len(ocov._coverage[pos]) > 1:
            tot = ocov.total(pos)
            # ratios = sorted([round(i/tot,2) for i in ocov._coverage[pos].values()])
            # print(ratios, gene.region_at(pos))
            for ii, i in enumerate(ocov._coverage[pos].values()):
                al =  model.quicksum(
                    v for s, v in VCN.items()
                    if s in structures and structures[s].cn[reg[0]][reg[1]] == 1
                )
                model.addConstr(al>=1)
                model.addConstr(al>=1)
                als = (i / tot) * al
                v = model.addVar(name=f"VAF_{pos}_{ii}", vtype='I')
                model.addConstr(als - 0.5 <= v)
                model.addConstr(als + 0.5 >= v)

                VVAF[pos,ii] = model.addVar(lb=0, name=f"ABSVAF_{pos}_{ii}")
                model.addConstr(VVAF[pos,ii] + (v - als) >= 0)
                model.addConstr(VVAF[pos,ii] - (v - als) >= 0)
    objective = model.quicksum(VVAF.values())
    model.setObjective(objective)
    for status, opt, sol in model.solutions(0):
        print(sx, status, opt, sol)
        break
for i in results:
    hoho(i)

#%% Fusion calculation & Copy number calculation
# HG01089 *1/*4 / 1/5
# HG02649 *1/*2
# NA12762 *1/*3
# NA19466 *1/*17
# HG02087 *36+*10/*36+*10  36/36+10
import aldy.major
reload(sys.modules['aldy.major'])
truth = {}
cn_truth = {}
with open("temp/genotypes.dat") as f:
    for l in f:
        l = l.split()
        truth[l[0]] = l[2]
        cn_truth[l[0]] = l[3].strip().split(',')
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
    r += [[round((r[1][i] - r[2][i])/(max(r[1][i], r[2][i]) + 1),1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
    # r += [[round(r[1][i] / r[2][i],1) if r[0][i].startswith('*') else '' for i in range(len(r[0]))]]
    tab = tabulate(r, tablefmt="plain", stralign="right", numalign="right")
    if plot: print(tab)
    return data
def model(rs):
    import math
    (ref, _, _, _), (sample, _, sample_sam, cx) = results[rs[0]], results[rs[1]]
    data = plot(gene, ref, sample, False)
    max_cn = 1 + max(int(math.ceil(data.get((gi, r)))) for gi, g in enumerate(gene.regions) for r in g)

    fus_sup = {n: cx[0][n] / cx[1][n] for n in cx[0] if cx[1][n]>0}
    print(max_cn, fus_sup)
    structures = {(name, 0): structure
        for name, structure in gene.cn_configs.items()
        if name in ['1', '5']
        or (name in fus_sup and fus_sup[name] >= 1/(2*max_cn))
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

    # def diff_model(model):
    VERR = {}
    VEPR = {}

    CYP6_COEFF = 1 # model.addVar(name=f"C6C", lb=0.75, ub=1.25)
    regs = 0
    for r, (exp_cov0, exp_cov1) in region_cov.items():
        expr = 0
        expr6 = 0
        scale = max(exp_cov0, exp_cov1) + 1
        for s, structure in structures.items():
            if r in structure.cn[0]:
                expr += structure.cn[0][r] * VCN[s]
                expr6 += structure.cn[0][r] * VCN[s]
            if len(structure.cn) > 1 and r in structure.cn[1]:
                expr -= structure.cn[1][r] * VCN[s]

        if r not in gene.unique_regions: continue
        regs += 1

        VEPR[r] = model.addVar(name=f"EX_{r}", lb=-aldy.cn.MAX_CN_ERROR, ub=aldy.cn.MAX_CN_ERROR)
        model.addConstr(expr6 + VEPR[r] <= CYP6_COEFF * exp_cov0, name=f"C6COV_{r}")
        model.addConstr(expr6 + VEPR[r] >= CYP6_COEFF * exp_cov0, name=f"C6COV_{r}")

        VERR[r] = model.addVar(name=f"E_{r}", lb=-aldy.cn.MAX_CN_ERROR, ub=aldy.cn.MAX_CN_ERROR)
        model.addConstr(expr / scale + VERR[r] <= (exp_cov0 - exp_cov1) / scale, name=f"CCOV_{r}")
        model.addConstr(expr / scale + VERR[r] >= (exp_cov0 - exp_cov1) / scale, name=f"CCOV_{r}")

    PCE_COEFF = 2
    DIFF_COEFF = 10 / regs
    o1 = DIFF_COEFF * model.abssum(VERR.values(), coeffs={"E_pce": PCE_COEFF})  # 0-10
    objective = 0
    objective += o1

    FIT_COEFF = DIFF_COEFF/3
    o2 = FIT_COEFF * model.abssum(VEPR.values())
    objective += o2

    # CYP6_ABS = model.addVar(lb=0, name=f"ABS_C6C")
    # model.addConstr(CYP6_ABS + (1 - CYP6_COEFF) >= 0, name=f"CABSL_C6C")
    # model.addConstr(CYP6_ABS - (1 - CYP6_COEFF) >= 0, name=f"CABSR_C6C")
    # objective += 5 * CYP6_ABS

    PARS = 10 / regs
    PARS /= 2
    penalty = {s: PARS for s, _ in VCN}
    for n, s in gene.cn_configs.items():
        if n in penalty and str(s.kind) == str(aldy.gene.CNConfigType.LEFT_FUSION):
            penalty[n] += PARS/2
    penalty['2D7'] += PARS/2
    # penalty
    o3 = model.quicksum(penalty[s] * v for (s, _), v in VCN.items())
    objective += o3
        # VVAF = {}
        # ocov = sample_sam.coverage.filtered(
        #     lambda mut, cov, total, thres: aldy.coverage.Coverage.basic_filter(
        #         mut, cov, total, thres / 10, results['HG03624'][2].coverage.min_cov
        #     )
        # )
        # for (pos, _) in gene.mutations:
        #     if pos in ocov._coverage:
        #         reg = gene.region_at(pos)
        #         if reg and len(ocov._coverage[pos]) > 1:
        #             tot = ocov.total(pos)
        #             for ii, i in enumerate(ocov._coverage[pos].values()):
        #                 als = (i / tot) * model.quicksum(
        #                     v for s, v in VCN.items()
        #                     if s in structures and structures[s].cn[reg[0]][reg[1]] == 1
        #                 )
        #                 v = model.addVar(name=f"VAF_{pos}_{ii}", vtype='I')
        #                 model.addConstr(als - 0.5 <= v)
        #                 model.addConstr(als + 0.5 >= v)

        #                 VVAF[pos,ii] = model.addVar(lb=0, name=f"ABSVAF_{pos}_{ii}")
        #                 model.addConstr(VVAF[pos,ii] + (v - als) >= 0)
        #                 model.addConstr(VVAF[pos,ii] - (v - als) >= 0)
        #                 print('adjusting', round(i/tot, 2))
        # objective += 0.1 * model.quicksum(VVAF.values())
    model.setObjective(objective)
    # Solve the model
    # model.addConstr(VCN['1', 0] >= 1)
    # model.addConstr(VCN['1', -1] >= 1)
    # model.addConstr(VCN['2D7', 0] >= 1)
    # model.addConstr(VCN['36.ALDY', 0] >= 1)
    # model.addConstr(VCN['36.ALDY', -1] >= 1)
    lookup = {model.varName(v): a for (a, ai), v in VCN.items()}
    result: dict = {}
    found = {}
    best_opt = None
    for status, opt, sol in model.solutions(0.2): # gap = 0
        if best_opt is not None and abs(opt-best_opt) > 0.5: continue
        if best_opt is None: best_opt = opt
        sol_tuple = tuple(sorted(lookup[v] for v in sol if lookup[v] not in ['5', '2D7']))
        sol_tuple_x = tuple(sorted(lookup[v] for v in sol))
        if sol_tuple not in result:
            result[sol_tuple] = aldy.cn.CNSolution(gene, opt, list(sol_tuple))
            s = tuple(natsorted( '36' if s == '36.ALDY' else s for s in sol_tuple_x ))
            found[s] = (round(opt - best_opt,2), round(o1.getValue(),2), round(o2.getValue(),2), round(o3.getValue(),2))
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
    print(rs[1], truth[rs[1]], cn_truth[rs[1]], 'vs', found, round(best_opt,2))
    print(tuple(cn_truth[rs[1]]) in found, {i: (cx[0][i],cx[1][i]) for i in cx[0]})
    plot(gene, ref, sample, True)
    # return tuple(cn_truth[rs[1]]) in found

    major_sols = []
    for i, cn_sol in enumerate(result.values()):
        major_sols += aldy.major.estimate_major(
            gene,
            sample_sam.coverage,
            cn_sol,
            solver='gurobi',
            identifier=i,
        )
    major_sols.sort(key=lambda x:x.score)
    for m in major_sols:
        print(m.cn_solution._solution_nice(), m.score, m._solution_nice())

    return found

# HG00185 HG02087? HG03624?

# print(model(results['HG00332'], results['HG01089']))
# for r in 'HG02649 NA19315 NA19466'.split():
r = 'NA19315'
for i in sorted(results):
# for i in ['HG00437', 'HG03624']:
# for i in "HG01089 HG02087 HG02649 HG03624 NA12762 NA19452 NA19466".split():
# for i in "HG00118 HG00437 HG02087 HG03624".split():
        # model(('HG00332', i))
    model((r, i))
    print()
# HG01089 *1/*4 vs {('1', '2D7'), ('1', '5')} 15.422513442593278


#%%
# majors
func_muts = { m[0] for m in gene.mutations if gene.is_functional(m) }
len(func_muts)

# for a in alleles:
#         for m in alleles[a].func_muts:
#             constraints[m] += VA[a]

#%%
alleles = {}
for an, a in gene.alleles.items():
    # print(a)
    va = {p: 'X' for p in func_muts}
    for p, op in a.func_muts:
        if op[1] == '>':
            va[p] = op[2]
        elif op.startswith('ins'):
            pass
            # va[p] = '+' + op[3:]
        elif op.startswith('del'):
            for l in range(len(op) - 3):
                va[p] = '-'
    alleles[an] = va


ch = results['HG01845'][2].choices
ch_i = list(ch.keys())
for gi, g in results['HG01845'][2].grid.items():
    l = {ch_i[i]: sig for i, sig in enumerate(g) if ch_i[i] in func_muts if sig != ''}
    sup = set()
    if l:
        for an, aa in alleles.items():
            ok = True
            for p, o in l.items():
                if p in aa and o != aa[p]:
                    ok = False
                    break
            if ok: sup.add(an)
        print(sup)
        break


#%%
import matplotlib.pyplot as plt

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
