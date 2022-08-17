# 786
# Aldy source: query.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


import sys
import textwrap
from natsort import natsorted

from .common import AldyException, log
from .gene import Gene


def query(gene: Gene, query: str):
    """Query a gene instance."""

    minors = {m: a for a, al in gene.alleles.items() for m in al.minors}
    alts = {
        mal.alt_name: (a, m)
        for a, al in gene.alleles.items()
        for m, mal in al.minors.items()
        if mal.alt_name
    }
    if query:
        if query.lower() in map(str.lower, gene.cn_configs):
            query = next(q for q in gene.cn_configs if query.lower() == q.lower())
            print_cn(gene, query)
        elif query.lower() in map(str.lower, gene.alleles):
            query = next(q for q in gene.alleles if query.lower() == q.lower())
            print_majors(gene, query)
        elif query.lower() in map(str.lower, minors):
            query = next(q for q in minors if query.lower() == q.lower())
            print_minors(gene, minors[query], query)
        elif query.lower() in map(str.lower, alts):
            query = next(q for q in alts if query.lower() == q.lower())
            print_minors(gene, *alts[query])
        else:
            log.error(f"ERROR: Cannot parse query '{query}'")
            sys.exit(1)
        return

    name = [f"Gene {gene.name}"]
    if gene.ensembl:
        name.append(f"ENSEMBL: {gene.ensembl}")
    if gene.pharmvar:
        name.append(f"PharmVar ID: {gene.pharmvar}")
    log.info("-" * 80)
    log.info("{}{}", name[0], "" if len(name) == 1 else f" ({', '.join(name[1:])})")
    log.info("-" * 80)

    st, ed = gene.ref_to_chr[0], gene.ref_to_chr[len(gene.seq) - 1]
    log.info(f"{gene.genome} genome locus: {gene.chr}:{st}-{ed} ({gene.refseq})")
    log.info("Strand: " + "3' (reverse)" if gene.strand < 0 else "5'")
    pseudo = ", ".join(f"{p} (ID {i + 1})" for i, p in enumerate(gene.pseudogenes))
    log.info(f"Pseudogenes: {pseudo if pseudo else 'none'}")
    amino = "\n  ".join(textwrap.wrap(gene.aminoacid, width=78))
    log.info(f"Aminoacid:\n  {amino}")

    log.info("Regions of interest:")
    s = [f"  {'Region:':>10}  {gene.name:>25}"]
    if gene.pseudogenes:
        s[-1] += f" {gene.pseudogenes[0]:>25}"
    for r in gene.regions[0]:
        lg = str(gene.regions[0].get(r, "-"))
        lp = str((gene.regions[1] if len(gene.regions) > 1 else {}).get(r, ""))
        s.append(f"  {r:>10}: {lg:>25} {lp:>25}")
    log.info("\n".join(s))

    log.info("Structural alleles (deletions, conservations and fusions):")
    for c, cn in gene.cn_configs.items():
        log.info(f"  {'*'+c+':':>8} {cn.description}")

    log.info("Major star-alleles:")
    for a, allele in natsorted(gene.alleles.items()):
        if "#" in a:  # avoid fusions here
            continue
        log.info(f"  *{a}:")
        m = ",\n                   ".join(
            _print_mutation(gene, m) for m in sorted(allele.func_muts)
        )
        log.info(f"    Key mutations: {m if m else 'none'}")

        mins = ",\n                   ".join(
            f"*{a}" + (f" (*{mi.alt_name})" if mi.alt_name else "")
            for a, mi in natsorted(allele.minors.items())
        )
        log.info(f"    Minor alleles: {mins if mins else 'none'}")


def print_cn(gene: Gene, major: str):
    """Pretty-print a copy number configuration."""
    if major not in gene.cn_configs:
        raise AldyException(f"Structural allele {major} not found in {gene.name}")

    config = gene.cn_configs[major]
    log.info(f"Gene {gene.name}, structural allele {gene.name}*{major}:")

    lg = " ".join(f"{i:<8}" for i in [gene.name] + gene.pseudogenes)
    log.info(f"  Structure: {lg}".rstrip())
    for r in gene.regions[0]:
        lg = " ".join(f"{config.cn[i][r]:<8}" for i, _ in enumerate(gene.regions))
        log.info(f"      {r:>5}: {lg}".rstrip())

    if major in gene.alleles:
        print_majors(gene, major)


def print_majors(gene: Gene, major: str):
    """Pretty-print a major star-allele."""

    if major not in gene.alleles:
        raise AldyException(f"Major star-allele {major} not found in {gene.name}")

    allele = gene.alleles[major]
    log.info(f"Gene {gene.name}, major star-allele {gene.name}*{allele.name}:")
    log.info(f"  Structure: {allele.cn_config}")

    activities = {m.activity for m in allele.minors.values() if m.activity}
    if len(activities) == 1:
        log.info(f"  Activity: {activities.pop()}")
    elif len(activities) > 1:
        log.info(f"  Activity: {' or '.join(activities)} (please check minor alleles)")

    ms = ",\n                 ".join(
        _print_mutation(gene, m) for m in sorted(allele.func_muts)
    )
    log.info(f"  Key mutations: {ms if ms else 'none'}")

    log.info("  Minor star-alleles:")
    for a, minor in natsorted(allele.minors.items()):
        log.info(f"    *{a}:")
        m = ",\n                        ".join(
            _print_mutation(gene, m) for m in sorted(minor.neutral_muts)
        )
        if minor.alt_name:
            log.info(f"      Legacy name: *{minor.alt_name}")
        log.info(f"      Silent mutations: {m if m else 'none'}")


def print_minors(gene: Gene, major: str, minor: str):
    """Pretty-print a minor star-allele."""

    allele = gene.alleles[major]
    if minor not in allele.minors:
        raise AldyException(f"Minor star-allele {minor} not found in {gene.name}")

    log.info(f"\nGene {gene.name}, minor star-allele {gene.name}*{minor}:")
    log.info(f"  Structure: {allele.cn_config}")
    log.info(f"  Major star-allele: {gene.name}*{allele.name}")
    mn = allele.minors[minor]
    if mn.activity:
        log.info(f"  Activity: {mn.activity}")
    if mn.evidence:
        evidence = {"L": "limited", "D": "definitive"}.get(mn.evidence, mn.evidence)
        log.info(f"  Evidence: {evidence}")
    if mn.pharmvar:
        log.info(f"  PharmVar details: {mn.pharmvar}")
    if mn.alt_name:
        log.info(f"  Legacy name: *{mn.alt_name}")

    ms = ",\n                 ".join(_print_mutation(gene, m) for m in allele.func_muts)
    log.info(f"  Key mutations: {ms if ms else 'none'}")

    m = ",\n                    ".join(
        _print_mutation(gene, m) for m in sorted(mn.neutral_muts)
    )
    log.info(f"  Silent mutations: {m if m else 'none'}")


def _print_mutation(gene: Gene, m):
    """Pretty-print a mutation."""

    fields = [
        gene.get_rsid(m, default=False),
        gene.get_refseq(m, from_atg=True),
        f"{gene.refseq}:" + gene.get_refseq(m),
    ]
    if m in gene.mutations and gene.mutations[m][0]:
        fields.append(gene.mutations[m][0])  # type: ignore
    if fields[0] == "-":
        fields = fields[1:]
    return f"{gene.chr}:{m} ({', '.join(fields)})"
