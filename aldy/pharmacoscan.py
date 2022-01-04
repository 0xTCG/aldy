# 786
# Aldy source: pharmacoscan.py

# %%
import collections
import pandas as pd
import numpy as np
import sys
from io import StringIO
from aldy.gene import Gene
from aldy.common import script_path

def compare(*args):
    pd.set_option("display.max_columns", None)

    # %%

    dp = pd.read_csv(sys.argv[1], sep="\t", skiprows=5)


    def aldy(x):
        r = x.Ref_Allele
        o = set(x.iloc[1].split("/")) - set([r])
        if not o or len(o) >= 2:
            if len(o) >= 2:
                print(f'Ignoring {x.dbSNP_RS_ID} from Pharmacoscan', file=sys.stderr)
            return ""
        o = o.pop()
        if r == "-":
            return f"del{o}"
        if o == "-":
            return f"ins{r}"
        return f"{r}>{o}"


    dp["Aldy"] = dp.apply(aldy, axis=1)
    fc = dp.columns[1]
    cols = {
        "Chr_id": "chr",
        "Start": "pos",
        "Aldy": "op",
        "dbSNP_RS_ID": "rsid",
        fc: "call",
        "Trivial_Name": "effect",
    }
    dp.head()  # hg38 +1

    with open(sys.argv[2]) as f:
        s = "".join(l for i, l in enumerate(f.readlines()) if not i or (l and l[0] != "#"))
    da = pd.read_csv(StringIO(s), sep="\t")
    da.head()  # hg38 +0

    # %%
    header = (
        "gene chr pos op pharmacoscan_rsid pharmacoscan_call pharmacoscan_comment "
        "aldy_rsid aldy_call aldy_comment reason"
    )
    print(*header.split(), sep="\t")
    for gene in da.Gene.unique():
        gg = Gene(
            script_path("aldy.resources.genes/{}.yml".format(gene.lower())), genome="hg38"
        )
        ch, st, ed = gg.get_wide_region()
        ch, st, ed
        dp_filt = dp[
            (dp.Chr_id == ch) & (dp.Start >= st) & (dp.Stop <= ed) & (dp.Aldy.str.len() > 0)
        ][cols.keys()].rename(columns=cols)
        dg = da[da.Gene == gene]
        cn = len(dg["Copy"].unique())
        aldy_mut = collections.defaultdict(lambda: [False] * cn)
        for gid, g in da[da.Gene == gene].groupby("Copy"):
            for _, r in g.iterrows():
                if not np.isnan(r.Location):
                    aldy_mut[int(r.Location + 1), r.Type, r.dbSNP, r.Effect][gid] = True
        rows = []
        for (pos, code, rsid, effect), call in aldy_mut.items():
            if code.startswith("ins"):
                ref, alt = "-", code[3:]
            elif code.startswith("del"):
                ref, alt = code[3:], "-"
            else:
                ref, alt = code[0], code[2]
            rows.append(
                [ch, pos, code, rsid, "/".join(alt if c else ref for c in call), effect]
            )
        da_filt = pd.DataFrame.from_records(
            rows, columns="chr pos op rsid call effect".split()
        )
        da_filt

        m = dp_filt.merge(
            da_filt,
            on=["chr", "pos", "op"],
            suffixes=["_pharmacoscan", "_aldy"],
            how="outer",
        )

        for _, r in m[m.call_pharmacoscan != m.call_aldy].iterrows():
            r = r.fillna("-")
            if sorted(r.call_aldy.split("/")) == sorted(r.call_pharmacoscan.split("/")):
                continue
            print(
                gene,
                r.chr,
                r.pos,
                r.op,
                r.rsid_pharmacoscan,
                r.call_pharmacoscan,
                r.effect_pharmacoscan,
                r.rsid_aldy,
                r.call_aldy,
                r.effect_aldy,
                end="\t",
                sep="\t",
            )
            comment = ""
            if r.call_aldy == "-":
                if (r.pos - 1, r.op) not in gg.mutations:
                    comment = "not present in Aldy database"
                else:
                    comment = "not called by Aldy"
            elif r.call_pharmacoscan == "-":
                comment = "not called by Pharmacoscan"
            else:
                comment = "mismatch"
            print(comment)
