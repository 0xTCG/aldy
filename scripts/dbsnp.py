# 786

import os, sys
import myvariant, collections
from common import log


def get_dbsnp(data, region, force=False):
    mv = myvariant.MyVariantInfo()
    q = mv.query(
        "_exists_:dbsnp AND _exists_:hg19 AND {}:{}-{}".format(*region),
        fields="dbsnp",
        fetch_all=True,
    )
    snps = list(q)

    # VCF, dbSNP and myVariant use 1-based indexing
    dbsnp = collections.defaultdict(dict)
    for snp in snps:
        pos, ref, alt, rs = (
            snp["dbsnp"]["hg19"]["start"] - 1,
            snp["dbsnp"]["ref"],
            snp["dbsnp"]["alt"],
            snp["dbsnp"]["rsid"],
        )
        if len(ref) > 1 or len(alt) > 1:
            assert ref[0] == alt[0]
        if len(ref) > 1:
            op = "DEL.{}".format(ref[1:])
        elif len(alt) > 1:
            op = "INS.{}".format(alt[1:].lower())
        else:
            op = "SNP.{}{}".format(ref, alt)
        dbsnp[pos][op] = rs

    mutations = {}
    for a in sorted(data):
        for m in data[a]["mutations"]:
            if m["pos"] == "pseudogene":
                continue
            if m["dbsnp"] not in ["", "*"]:
                m["dbsnp"] = [m["dbsnp"]]
            else:
                m["dbsnp"] = []
            pos, op = m["pos"], m["op"]

            # check reversed SNP
            if op in dbsnp[pos]:
                rsid = str(dbsnp[pos][op])
                if rsid not in m["dbsnp"]:
                    if len(m["dbsnp"]) > 0:
                        m["dbsnp"][0] += "(k)"
                    m["dbsnp"].append(rsid)
                    log.debug("dbSNP: Variant {} assigned to {}:{}", rsid, pos, op)
                else:
                    log.debug(
                        "dbSNP: Variant {} matches the Karolinska's prediction", rsid
                    )
            elif len(dbsnp[pos]) > 0 and (
                op[:3] == "SNP" and op[:4] + op[4:6][::-1] in dbsnp[pos]
            ):
                op = op[:4] + op[4:6][::-1]
                rsid = str(dbsnp[pos][op])
                if rsid not in m["dbsnp"]:
                    if len(m["dbsnp"]) > 0:
                        m["dbsnp"][0] += "(k)"
                    m["dbsnp"].append(rsid)
                    log.debug("dbSNP: Variant {} assigned to {}:{}", rsid, pos, op)
                else:
                    log.debug(
                        "dbSNP: Variant {} matches the Karolinska's prediction", rsid
                    )
            elif len(dbsnp[pos]) != 0:
                log.trace("How about {} for {}:{} ({})", dbsnp[pos], pos, op, m["old"])
    return data
