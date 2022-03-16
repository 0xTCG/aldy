#%%
import os, sys, os.path, glob
import statistics
from collections import *
from pprint import pprint
from natsort import natsorted
from imp import reload
from multiprocessing import Pool
import tqdm
from tqdm.contrib.concurrent import process_map
import aldy.gene, aldy.coverage, aldy.sam, aldy.common
# reload(sys.modules['aldy.gene'])
# reload(sys.modules['aldy.coverage'])
# reload(sys.modules['aldy.sam'])
# reload(sys.modules['aldy.common'])
from aldy.gene import Gene
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

# sh = logbook.more.ColorizedStderrHandler(
#     format_string="{record.message}", level='DEBUG'
# )
# sh.push_application()

gene = Gene('aldy/resources/genes/cyp2d6.yml', genome='hg38')
cn_region = aldy.sam.DEFAULT_CN_NEUTRAL_REGION['hg38']
gene.get_wide_region()

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


#%%

def remap(gf):
    gene, file = gf
    slices = [gene.get_wide_region().start]
    slices += sorted({
        x
        for g in gene.regions
        for r, gr in g.items()
        # if r in gene.unique_regions
        if gr.end - gr.start
        for x in (gr.start, gr.end)
    })
    slices += [gene.get_wide_region().end]
    slices += [42_155_000]
    with timing(file), pysam.AlignmentFile(file) as sam, tempfile.TemporaryDirectory() as tmp:
        skip = {'hard': 0, 'sup': 0}
        iter = sam.fetch(region=GRange(gene.get_wide_region()[0], slices[0], slices[-1]).samtools(prefix='chr'))
        with open(f"{tmp}/reads.fq", "w") as fq:
            for read in iter:
                if not read.cigartuples:  # only valid alignments
                    continue
                if read.reference_start >= slices[-1] or read.reference_end < slices[0]:
                    continue
                if "H" in read.cigarstring:  # avoid hard-clipped reads
                    skip['hard'] += 1
                seq, qual = read.get_forward_sequence(), read.get_forward_qualities()

                # seq, qual =  read.query_alignment_sequence, read.query_alignment_qualities
                # pairs = [r for r in read.get_aligned_pairs() if r[1]]
                # read_pos = [r for r, _ in pairs]
                # ref_pos = [r for _, r in pairs]
                # l = [0]
                # for st in slices:
                #     if ref_pos[0] <= st < ref_pos[-1]:
                #         i = bisect.bisect_left(ref_pos, st)
                #         if i < len(ref_pos) and (not l or l[-1] != read_pos[i]):
                #             l.append(read_pos[i])
                # l.append(read.query_alignment_length)
                if False: # simple 500bp reads
                    l = list(range(0, len(seq), 500))
                    if len(seq) - l[-1] < 300:
                        l[-1] = len(seq)
                    else:
                        l.append(len(seq))
                    for i in range(1, len(l)):
                        print(f'@{read.query_name}${i}', file=fq)
                        print(seq[l[i - 1]:l[i]], file=fq)
                        print('+', file=fq)
                        print(''.join(chr(q + 33) for q in qual[l[i - 1]:l[i]]), file=fq)
                if True: # paired-end 100bp reads w/ insert size 300
                    for i in range(0, len(seq), 100):
                        if i + 400 >= len(seq): continue
                        j = i + 400
                        if j + 100 >= len(seq):
                            j = len(seq) - 100
                        for ki, k in enumerate([i, j]):
                            print(f'@{read.query_name}${i}', file=fq)
                            print(seq[k:k+100], file=fq)
                            print('+', file=fq)
                            print(''.join(chr(q + 33) for q in qual[k:k+100]), file=fq)
        # print(file, skip)
        base = os.path.splitext(os.path.basename(file))[0]

        sp.run([
            'bwa', 'mem', '-p', "-x", "pacbio", "-o", f"{tmp}/remap.sam", "-v", "1",
            '/Users/inumanag/Projekti/aldy/pacbio/temp/chr22.fa', f"{tmp}/reads.fq"
        ])
        sp.run([
            'samtools', "sort", f"{tmp}/remap.sam", "-o", f"{tmp}/remap.bam"
        ])

        # sp.run([
        #     'temp/pbmm2', 'align', '-j', '1', '--preset', 'CCS', '--sort', '-c', '0',
        #     '--log-level', 'FATAL',
        #     '-y', '70', 'temp/chr22.fa', f"{tmp}/reads.fq",
        #     f"{tmp}/remap.bam"
        # ])
        sp.run([
            "samtools", "view", "-b", file, "chr4", "-o", f"{tmp}/cn.bam",
        ])
        out = f"{os.path.dirname(file)}/{base}.remap.bam"
        sp.run(["samtools", "merge", "-f", "-o", out, f"{tmp}/cn.bam", f"{tmp}/remap.bam"])
        sp.run(["samtools", "index", out])

        # D KOPI NAMBR RIDZHN
        s = aldy.sam.Sample(gene, out, profile="illumina", cn_region=cn_region)
        # profile_cn = sum(profile[i] for i in range(cn_region.start, cn_region.end))
        # cn_ratio = float(profile_cn) / sam_ref

    dx = {}
    for p in range(cn_region.start, cn_region.end):
        dx[(-1, 'cn'), p] = s.coverage._cnv_coverage[p]
    for gi, gr in enumerate(gene.regions):
        for region, rng in gr.items():
            for p in range(rng.start, rng.end):
                dx[(gi, region), p] = int(s.coverage.total(p))
    return (os.path.basename(file).split('.')[0], dx)

coverages = {}
pool = multiprocess.Pool(8)
data = [(gene, i) for i in sorted(glob.glob("temp/pacbio/???????.bam"))]
with timing("Map"):
    cov = dict(pool.map(remap, data))

#%%
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

# _ = plot(gene, cov['HG00332'], cov['HG00437'])

truth = {}
with open("temp/genotypes.dat") as f:
    for l in f:
        l = l.split()
        truth[l[0]] = l[2]

def model(ref, sample):
    data = plot(gene, ref, sample, True)
    return
    max_cn = 1 + max(int(round(data.get((gi, r), 0))) for gi, g in enumerate(gene.regions) for r in g)
    structures = {(name, 0): structure for name, structure in gene.cn_configs.items()}
    for a, ai in list(structures.keys()):
        structures[a, -1] = copy.deepcopy(structures[a, 0])
        if str(gene.cn_configs[a].kind) != str(aldy.gene.CNConfigType.DEFAULT):
            continue
        for i in range(1, max_cn):
            structures[a, i] = copy.deepcopy(structures[a, 0])
            for g in range(1, len(structures[a, i].cn)):
                structures[a, i].cn[g] = { r: v - 1 for r, v in structures[a, i].cn[g].items() }
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
        for r, (exp_cov0, exp_cov1) in region_cov.items():
            if r not in gene.unique_regions: continue
            expr = 0
            for s, structure in structures.items():
                if r in structure.cn[0]:
                    expr += structure.cn[0][r] * VCN[s]
                if len(structure.cn) > 1 and r in structure.cn[1]:
                    expr -= structure.cn[1][r] * VCN[s]
            VERR[r] = model.addVar(name=f"E_{r}", lb=-aldy.cn.MAX_CN_ERROR, ub=aldy.cn.MAX_CN_ERROR)
            model.addConstr(expr + VERR[r] <= (exp_cov0 - exp_cov1), name=f"CCOV_{r}")
            model.addConstr(expr + VERR[r] >= (exp_cov0 - exp_cov1), name=f"CCOV_{r}")
        objective = model.abssum(VERR.values(), coeffs={"E_pce": aldy.cn.PCE_PENALTY_COEFF})
        objective += model.quicksum(
            (aldy.cn.LEFT_FUSION_PENALTY + aldy.cn.PARSIMONY_PENALTY) * v
            if gene.cn_configs[s].kind == aldy.gene.CNConfigType.LEFT_FUSION
            else aldy.cn.PARSIMONY_PENALTY * v
            for (s, _), v in VCN.items()
        )
        model.setObjective(objective)
    def div_model(model):
        VERR_DIV = {}
        VERR_SUM = {}
        for r, (c0, c1) in region_cov.items():
            expr = 0
            # VERR_SUM[r] = model.addVar(name=f"ES_{r}", lb=-aldy.cn.MAX_CN_ERROR, ub=aldy.cn.MAX_CN_ERROR)
            # for s, structure in structures.items():
            #     if r in structure.cn[0]:
            #         expr += structure.cn[0][r] * VCN[s]
            #     if len(structure.cn) > 1 and r in structure.cn[1]:
            #         expr += structure.cn[1][r] * VCN[s]
            # model.addConstr(expr + VERR_SUM[r] <= c0 + c1, name=f"CCOV_{r}")
            # model.addConstr(expr + VERR_SUM[r] >= c0 + c1, name=f"CCOV_{r}")
            if r not in gene.unique_regions: continue
            expr = 0
            VERR_DIV[r] = model.addVar(name=f"E_{r}", lb=-aldy.cn.MAX_CN_ERROR, ub=aldy.cn.MAX_CN_ERROR)
            if r == 'pce':
                for s, structure in structures.items():
                    if r in structure.cn[0]:
                        expr += structure.cn[0][r] * VCN[s]
                    if len(structure.cn) > 1 and r in structure.cn[1]:
                        expr -= structure.cn[1][r] * VCN[s]
                expr += VERR_DIV[r]
                expr -= (c0 - c1)
            else:
                for s, structure in structures.items():
                    if r in structure.cn[0]:
                        expr += c1 * structure.cn[0][r] * VCN[s]
                    if len(structure.cn) > 1 and r in structure.cn[1]:
                        expr += c1 * VERR_DIV[r] * structure.cn[1][r] * VCN[s]
                        expr -= c0 * structure.cn[1][r] * VCN[s]
            model.addConstr(expr <= 0, name=f"CCOV_{r}")
            model.addConstr(expr >= 0, name=f"CCOV_{r}")
        objective = (
            # 1.0 * model.abssum(VERR_SUM.values()) +
            1.0 * model.abssum(VERR_DIV.values(), coeffs={"E_pce": aldy.cn.PCE_PENALTY_COEFF})
        )
        objective += model.quicksum(
            (aldy.cn.LEFT_FUSION_PENALTY + aldy.cn.PARSIMONY_PENALTY) * v
            if gene.cn_configs[s].kind == aldy.gene.CNConfigType.LEFT_FUSION
            else aldy.cn.PARSIMONY_PENALTY * v
            for (s, _), v in VCN.items()
        )
        model.setObjective(objective)
    # Solve the model
    model = setup_model()
    diff_model(model)
    # model.addConstr(VCN['1', 0] >= 1)
    # model.addConstr(VCN['1', -1] >= 1)
    # model.addConstr(VCN['36.ALDY', 0] >= 1)
    # model.addConstr(VCN['36.ALDY', -1] >= 1)
    lookup = {model.varName(v): a for (a, ai), v in VCN.items()}
    result: dict = {}
    found = set()
    for status, opt, sol in model.solutions(0): # gap = 0
        sol_tuple = tuple(sorted(lookup[v] for v in sol))
        if sol_tuple not in result:
            result[sol_tuple] = c = aldy.cn.CNSolution(gene, opt, list(sol_tuple))
            s = tuple(natsorted( '36' if s == '36.ALDY' else s for s in sol_tuple ))
            found.add(s)
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
    return found
for i in sorted(cov):
    print(i, truth[i])
    print(model(cov['HG00332'], cov[i]))
    print()


#%% DO MAP ANALYSIS
reload(sys.modules['aldy.sam'])
sam = aldy.sam.Sample(gene, "temp/pacbio/NA19315.remap.bam", profile="illumina")


#%%

#%%
#%%
regions = list(gene.regions[0].keys())
def after(p, c):
    if p == c: return 1
    try:
        pi, ci = regions.index(p), regions.index(c)
        return 0 if pi + 1 == ci else -1
    except ValueError:
        return -1

modes = []
for read, chunks in sam.reads.items():
    chunks.sort(key=lambda x: x[0])
    # print(read, chunks)
    # arranged meeting
    p = []
    for _, s in chunks:
        s = s[::-1]
        if not s: continue
        if not p: p.append(s)
        else:
            a = after(p[-1][-1], s[0])
            if a == -1: p.append(s)
            else: p[-1] += s[a:]
    p = [(i[0], len(i)) for i in p if len(i) > 1 or i[0][1] != 'rep']
    # s = [i for _, s in chunks for i in s if i]
    modes.append(tuple(p))

#%%
# Counter(modes).most_common()

#%% --------------
#  set best_n?
minimap2_index = mappy.Aligner("temp/chr22_small.fa", preset="map-hifi", best_n=1)
s6 = minimap2_index.seq(
    "chr22:42000000-43000000",
    gene.regions[0]['rep'].start - 42_000_000, gene.regions[0]['up'].end - 42_000_000
)
s7 = minimap2_index.seq(
    "chr22:42000000-43000000",
    gene.regions[1]['rep'].start - 42_000_000, gene.regions[1]['up'].end - 42_000_000,
)
gene_index = [
    mappy.Aligner(seq=s6, preset="map-hifi", best_n=1),
    mappy.Aligner(seq=s7, preset="map-hifi", best_n=1)
]

#%% remap whole, split map later
# reload(sys.modules['aldy.sam'])
def remap_mappy(*args):
    gene, file = args[0]
    slices = []
    for i in range(2):
        slices.append(min(gene.regions[i].values(), key=lambda x: x.start).start)
        slices.append(max(gene.regions[i].values(), key=lambda x: x.end).end)
    reads = {}
    max_read = 0
    with pysam.AlignmentFile(file) as sam, tempfile.TemporaryDirectory() as tmp:
        skip = {'hard': 0, 'sup': 0}
        iter = sam.fetch(region=GRange(gene.get_wide_region()[0], slices[0], 42_155_000).samtools(prefix='chr'))
        with open(f"{tmp}/reads.fq", "w") as fq:
            for read in iter:
                if not read.cigartuples:  # only valid alignments
                    continue
                if read.reference_start >= slices[-1] or read.reference_end < slices[0]:
                    continue
                if "H" in read.cigarstring:  # avoid hard-clipped reads
                    skip['hard'] += 1
                    continue
                seq = read.query_sequence

                span = slices[1] - slices[0] - 1000
                l = [0]
                if read.reference_start < slices[1]:
                    l.append(slices[1] - read.reference_start)
                else:
                    l.append(slices[3] - read.reference_start)
                l += list(range(l[-1] + span, len(seq), span))
                if len(seq) - l[-1] > 500:
                    l.append(len(seq))
                else:
                    l[-1] = len(seq)
                # SPLIT
                pieces = []
                for i in range(1, len(l)):
                    max_read = max(max_read, l[i] - l[i - 1])
                    s = seq[l[i - 1]:l[i]]
                    for h in minimap2_index.map(s):
                        if h.is_primary and h.strand == 1:
                            # assert h.strand == 1  # BAM already fixed those
                            pieces.append((
                                s, h.q_st, h.q_en, h.r_st, h.r_en, list(h.cigar)
                            ))
                            # reads.setdefault(read.query_name, []).append( (i, h, s) )
                            break
                if not pieces:
                    continue
                nseq, ed, cigar = [], pieces[0][3], []
                for s, qst, qed, rst, red, cig in pieces:
                    nseq.append(s)
                    if qst != 0: cig = [(qst, 4)] + cig # S
                    if rst != ed: cig = [(rst - ed, 3)] + cig # N
                    if qed != len(s): cig += [(len(s) - qed, 4)]
                    cigar += cig
                    ed = red
                reads[read.query_name] = (''.join(nseq), pieces[0][3], cigar)
        s = aldy.sam.Sample(gene, minimap=reads, sample_name=os.path.basename(file))
    dx = {}
    for p in range(cn_region.start, cn_region.end):
        dx[(-1, 'cn'), p] = s.coverage._cnv_coverage[p]
    for gi, gr in enumerate(gene.regions):
        for region, rng in gr.items():
            for p in range(rng.start, rng.end):
                dx[(gi, region), p] = int(s.coverage.total(p))
    # print(os.path.basename(file), max_read)
    return reads, os.path.basename(file).split('.')[0], dx

fusion_signatures = {}
for n, cn in gene.cn_configs.items():
    if cn.kind == aldy.gene.CNConfigType.LEFT_FUSION:
        sig, pos = 0, next(i for i, j in cn.cn[0].items() if j == 1)
    elif cn.kind == aldy.gene.CNConfigType.RIGHT_FUSION:
        sig, pos = 1, next(i for i, j in cn.cn[0].items() if j == 0)
    else:
        continue
    sig = ((1 - sig, list(cn.cn[0])[list(cn.cn[0]).index(pos) - 1]), (sig, pos))
    if sig[0][1] not in gene.unique_regions and sig[1][1] not in gene.unique_regions:
        continue
    fusion_signatures[sig] = n
def get_regions(seq, r_start, cigar):
    regs = []
    start, s_start = r_start + 42_000_000, 0
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
# Print read: regions
def foo(x):
    dc = []
    for (g, i), _, _ in x:
        if not dc or dc[-1][0] != g:
            dc.append((g, []))
        dc[-1][1].append(i)
    return '; '.join(f"{k}@[" + ','.join(v) + "]" for k, v in dc)
def best_map(idx, seq):
    maps = [hit for hit in idx.map(seq) if hit.is_primary]
    return None if not maps else maps[0]
def get_left(regs, li):
    while True:
        left = regs[li][1]
        if regs[li][0][1] == 'e9':
            break
        li -= 1
        if li < 0:
            left = 0
            break
    return left
def get_right(regs, li):
    while True:
        right = regs[li][1] + regs[li][2]
        if regs[li][0][1] == 'e1':
            break
        li += 1
        if li == len(regs):
            break
    return right
def cluster(rds):
    potential_fusions = {}
    missed_fusions = {fs: 0 for fs in fusion_signatures.values()}
    for name, (seq, r_start, cigar) in rds.items():
        regs = get_regions(seq, r_start, cigar)
        breaks = []
        for li in range(1, len(regs)):
            for p, n in fusion_signatures:
                if regs[li][0][1] == p[1] and regs[li - 1][0][1] == n[1]:
                    breaks.append((li, p, n))
        for _, p, n in breaks:
            missed_fusions[fusion_signatures[p, n]] += 1
        for li, p, n in breaks:
            (gene, _), breakpoint, _ = regs[li]

            if gene == p[0]:
                left = get_left(regs, li)
                s = seq[left:breakpoint]
                idx, alt_idx = p[0], 1 - p[0]
            else:
                right = get_right(regs, li)
                s = seq[breakpoint:right]
                idx, alt_idx = 1 - p[0], p[0]

            mp, alt_mp = best_map(gene_index[idx], s), best_map(gene_index[alt_idx], s)
            if mp and alt_mp and alt_mp.NM < mp.NM:
                potential_fusions.setdefault(fusion_signatures[p, n], {})[name] = (breakpoint, mp.NM, alt_mp.NM)
    # for f, x in potential_fusions.items():
    #     wt = round(mean([a-b for a,b in x.values()]), 1)
    #     print(f, wt, len(x), missed_fusions[f], sep='  ')
    return potential_fusions, missed_fusions
        # print(f, fusion_signatures[f], len(x), missed_fusions[f])
# with timing("Remap"):
    # rds = remap_mappy((gene, "temp/pacbio/HG00437.bam"))
    # rds, _, dx = remap_mappy((gene, "temp/pacbio/HG00437.bam"))

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
    heapq.heappush(heap, [max-len(s), se, set(s.keys())])
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

# with open("_x.txt", "w") as fo:
#     for name in zz['13']:
#         print(name, file=fo)
#         rd= rds[name]
#         print('\t', len(rd[0]), foo(get_regions(*rd)), fus['13'][name], sep="\t", file=fo)

#%% stats: perc of non
for (nn, ni, _), aln, seq in rds:
    if aln.q_st != 0 or aln.q_en != len(seq):
        q = 100.0 * (aln.q_st + len(seq) - aln.q_en) / len(seq)
        if q >= 10:
            print(len(seq), q)

#%%

# coverages = {}
# pool = multiprocess.Pool(8)
# data = [(gene, i) for i in sorted(glob.glob("temp/pacbio/???????.bam"))]
# with timing("Map"):
#     cov = dict(pool.map(remap_mappy, data))
# cf = remap_mappy(gene, 'temp/pacbio/HG00437.bam')

#%%
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
#%%
# Remap the file
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
# Print read: regions
def foo(x):
    dc = []
    for (g, i), _, _ in x:
        if not dc or dc[-1][0] != g:
            dc.append((g, []))
        dc[-1][1].append(i)
    return '; '.join(f"{k}@[" + ','.join(v) + "]" for k, v in dc)
reload(sys.modules['aldy.sam'])
def remap_mappy(*args):
    gene, file = args[0]
    slices = []
    for i in range(2):
        slices.append(min(gene.regions[i].values(), key=lambda x: x.start).start)
        slices.append(max(gene.regions[i].values(), key=lambda x: x.end).end)
    reads = {}
    with pysam.AlignmentFile(file) as sam, tempfile.TemporaryDirectory() as tmp, open("_x", "w") as fo:
        iter = sam.fetch(region=GRange(gene.get_wide_region()[0], slices[0], 42_155_000).samtools(prefix='chr'))
        for read in iter:
            if not read.cigartuples or "H" in read.cigarstring:  # only valid alignments
                continue
            if read.reference_start >= slices[-1] or read.reference_end < slices[0]:
                continue
            seq = read.query_sequence

            mapped = 0
            pieces = []
            broken=[]
            # print
            while mapped < len(seq):
                h = best_map(index, seq[mapped:])
                if not h:
                    break
                if not genes[h.ctg][2]:
                    pieces.append((
                        adjust(h.ctg, h.r_st), seq[mapped + h.q_st:mapped + h.q_en], h.cigar
                    ))
                else:
                    brk = genes[h.ctg][1][genes[h.ctg][2]][1] # breakpoint
                    if h.r_st <= brk < h.r_en:
                        broken.append(h.ctg)
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
                # if broken:
                #     print(read.query_name, len(seq), broken, file=fo)
                #     for p in pieces:
                #         print('', foo(get_regions(*p)), file=fo)
                #     print(file=fo)
        s = aldy.sam.Sample(gene, minimap=reads, sample_name=os.path.basename(file),
                cn_region=aldy.sam.DEFAULT_CN_NEUTRAL_REGION['hg38'])
    dx = {}
    for p in range(cn_region.start, cn_region.end):
        dx[(-1, 'cn'), p] = s.coverage._cnv_coverage[p]
    for gi, gr in enumerate(gene.regions):
        for region, rng in gr.items():
            for p in range(rng.start, rng.end):
                dx[(gi, region), p] = int(s.coverage.total(p))
    # print(os.path.basename(file), max_read)
    return reads, os.path.basename(file).split('.')[0], dx
with timing('foo'):
    r, s, dx = remap_mappy((gene, "temp/pacbio/HG00437.bam"))
print('done')

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
