# 786
# Aldy source: gene.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple, Dict, List, Optional

import os
import re
import enum
import yaml
import functools
import collections
import textwrap
import itertools

from natsort import natsorted
from pprint import pprint

from .common import (
    GRange,
    AldyException,
    allele_number,
    sorted_tuple,
    seq_to_amino,
    rev_comp,
    log,
    allele_sort_key,
)


EXON = "e"
"""str: Abbreviation for exon."""

INTRON = "i"
"""str: Abbreviation for intron."""


class MajorAllele(
    collections.namedtuple("MajorAllele", ["name", "cn_config", "func_muts", "minors"])
):
    """
    Major allele description.
    Immutable.

    Attributes:
        name (str):
            Major allele name.
        cn_config (str):
            Copy-number configuration.
        func_muts (set[:obj:`Mutation`]):
            Functional mutations that describe the major allele.
        minors (dict[str, :obj:`MinorAllele`]):
            Dictionary of minor alleles that are based upon the major allele.
    """

    def minor_mutations(self, minor: str):
        """
        Sequence of mutations that describe a minor allele.
        """
        for m in self.func_muts:
            yield m
        for m in self.minors[minor].neutral_muts:
            yield m


class MinorAllele(
    collections.namedtuple("MinorAllele", ["name", "alt_name", "neutral_muts"])
):
    """
    Minor allele description.
    Immutable.

    Attributes:
        name (str):
            Minor allale name.
        alt_name (Optional[str]):
            List of alternative names.
        neutral_muts (set[:obj:`Mutation`]):
            Netural mutations that describe the minor allele.

    Notes:
        Has custom printer (``__repr__``).
    """

    def __str__(self):
        return "Minor({}; [{}])".format(
            self.name,
            ", ".join(map(str, self.neutral_muts)),
        )


class Mutation(collections.namedtuple("Mutation", ["pos", "op"])):
    """
    Mutation description.
    Immutable.

    Attributes:
        pos (int): Reference genome position.
        op (str): Mutation description that can be:

            - ``SNP.AB``:  a SNP from A to B (A, B are in ``[ACGT]``)
            - ``INS.xx``:  an insertion of xx (xx is ``[acgt]+``)
            - ``DEL.xx``:  a deletion of xx (xx is ``[ACGT]+``)

        is_functional (int): Functionality of the mutation:
            - 0 for non-functional (silent) mutations
            - 1 for functional (gene-disrupting) mutations
            - 2 for splice-disrupting mutations

        aux (dict[str, str]): Auxiliary information.
            Currently used for storing dbSNP rsIDs (key: ``dbsnp``)
            and old Karolinska-style mutation IDs (key: ``old``).

    Notes:
        Has custom printer (``__str__``).
        Comparable and hashable via ``(pos, op)`` tuple.
        Implements ``total_ordering``.
    """

    def __str__(self):
        return f"{self.pos+1}.{self.op}"


class CNConfig(
    collections.namedtuple("CNConfig", ["cn", "kind", "alleles", "description"])
):
    """
    Copy-number (CN) configuration description.
    Immutable.

    Attributes:
        cn (dict[int, dict[str, int]]):
            Value of the expected region copy number in each gene.
            For example, ``cn[0]["e1"] == 1`` means that the exon 1 of
            the main gene (ID 0) has one copy (and thus should be present)
            in the configuration `cn`.
        kind (:obj:`CNConfigType`):
            Type of the copy-number configuration. See :obj:`CNConfigType` for details.
        alleles (set[str]):
            Allele IDs that have this CN configuration.
        description (str):
            Human-readable description of the configuration (e.g. "deletion").

    Notes:
        Has custom printer (``__str__``).
    """

    CNConfigType = enum.Enum(
        "CNConfigType", "DEFAULT_CN LEFT_FUSION RIGHT_FUSION DELETION"
    )
    """
    Enumeration describing the type of a copy-number configuration:

        - ``LEFT_FUSION``
        - ``RIGHT_FUSION``
        - ``DELETION``
        - ``DEFAULT_CN``
    """

    def __str__(self):
        return "CNConfig({}; vector={}; alleles=[{}])".format(
            str(self.kind)[13:],
            "|".join(
                "".join(str(self.cn[g][r]) for r in self.cn[0])
                for g in self.cn
            ),
            " ".join(sorted(self.alleles)),
        )


class Gene:
    """
    Gene (and associated pseudogenes) description.

    Attributes:
        name (str):
            Gene name (e.g. CYP2D6).
        seq (str):
            Wild-type reference sequence that describes \*1 allele.
            Should be location-compatible with hg19 (must contain no indels when
            aligned to the reference genome).
        region (:obj:`aldy.common.GRange`):
            Wild-type reference sequence coordinates within the reference genome.
        regions (dict[int, dict[str, :obj:`aldy.common.GRange`]]):
            Exonic, intronic and special (e.g. UTR) regions in a gene.
            Maps a gene ID to a dictionary that maps gene region IDs
            (e.g. "e9" for exon 9) to a :obj:`aldy.common.GRange` (e.g. chr1:10-20)
            within the reference genome.
            Gene 0 is the main gene.
        pseudogenes (list[str]):
            List of pseudogene names (genes whose ID is greater than 0).
        common_tandems (list[tuple[str, str]]):
            List of common allele tandems. Used in diplotype assignment heuristics.
            For example, the fact that \*13 is always followed by \*1
            (encoded as ``('13', '1')``) will be used to group \*13 and \*1 together
            within the same haplotype (e.g. \*13+\*1).
        cn_configs (dict[str, :obj:`CNConfig`]):
            Copy-number configurations associated with the gene.
            `1` (akin to \*1) is the default CN config
            (no structural variations of any kind).
        alleles (dict[str, :obj:`Allele`]):
            Major star-alleles in the gene.
            Accessible via major star-allele name (e.g. ``'1'``).
        coding_region (:obj:`CodingRegion`):
            Coding region of the gene (includes the protein sequence and
            the reference genome strand).
        mutations (dict[tuple[int, str], :obj:`Mutation`]):
            A hashtable that provides links ``(position, mutation_type)`` to
            a corresponding :obj:`Mutation`.
            Useful for querying functional mutations and auxiliary information.
    """

    CodingRegion = collections.namedtuple("CodingRegion", ["rev_comp", "aminoacid"])
    """
    Immutable class that describes the coding region of a gene.
    Attributes: aminoacid sequence (``aminoacid``) and reference strand (``rev_comp``).
    """

    def __init__(
        self,
        path: Optional[str],
        name: Optional[str] = None,
        yml: Optional[str] = None,
        genome: str = "hg19",
    ) -> None:
        """
        Initialize the Gene class with the database description
        specified in YML file ``path``.

        Args:
            path (str, optional):
                Location of YML file.
            name (str, optional):
                Gene name.

        Raises:
            :obj:`aldy.common.AldyException` if a path is not set.
        """

        if path and os.path.exists(path):
            with open(path) as file:
                yml = file.read()
            if not name:
                name = os.path.split(path)[-1].split(".")[0].upper()

        if not yml or not name:
            raise AldyException("Either a path or a name should be given")

        yml = yaml.safe_load(yml)
        self.genome = genome
        self._parse_yml(name, yml)

        # import sys
        # for i in self.alleles:
        #     print(i, '->', natsorted([x.name for x in self.alleles[i].minors.values()]))
        # sys.exit(0)

    def _parse_yml(self, gene_name: str, yml) -> None:
        """
        Initialize the gene structure from YML data.
        """

        self._init_basic(gene_name, yml)
        self._init_regions(yml)
        self._init_alleles(yml)
        self._init_partials()
        self._init_structure(yml)

    def _init_basic(self, gene: str, yml) -> None:
        """
        Read basic gene properties (``name``, ``seq`` and ``region``).

        All database YAML coordinates are indexed starting from 1.
        Aldy internally uses 0-based indexing.
        """
        self.name = yml["name"]
        self.seq = yml["reference"]["seq"].replace("\n", "")

        chr, start, end, strand, cigar = yml["reference"]["mappings"][self.genome]
        self.chr_to_ref = {}
        self.ref_to_chr = {}
        self.strand = 1 if strand == '+' else -1
        pos_ref = 0 if self.strand > 0 else (len(self.seq) - 1)
        pos_chr = start - 1
        for i in cigar.split():
            op, sz = i[0], int(i[1:])
            if op == 'M':
                for idx in range(sz):
                    self.chr_to_ref[pos_chr + idx] = pos_ref + idx * self.strand
                    self.ref_to_chr[pos_ref + idx * self.strand] = pos_chr + idx
                pos_chr += sz; pos_ref += sz * self.strand
            elif op == 'I':
                pos_ref += sz * self.strand
            elif op == 'D':
                pos_chr += sz
            else:
                raise AldyException("Invalid CIGAR string")
        self.region = GRange(chr, start, end)  # will be updated later
        self.version = yml["version"]
        # self.access_ids = tuple([yml["ensembl"], yml["pharmvar"], yml["reference"]["name"]])

    def _init_regions(self, yml) -> None:
        """
        Calculate the genic regions and pseudogenes
        (``regions``, ``unique_regions`` and ``pseudogenes``).
        Prepare ``region_at()`` call.
        """

        self.regions = {}
        # Each pseudogene is associated with an index > 0
        self.pseudogenes: List[str] = list()
        for i, g in enumerate(yml["structure"]["genes"]):
            if i > 0:
                self.pseudogenes.append(g)
            regions = {}
            num_exons = 0
            for ri, (name, coord) in enumerate(yml["structure"]["regions"][self.genome].items()):
                if name[0] == 'e' and name[1:].isdigit():  # exon
                    num_exons += 1
                regions[name] = GRange(self.region.chr, coord[i * 2] - 1, coord[i * 2 + 1] - 1)
            for e in range(1, num_exons):  # Fill introns
                r1, r2 = regions[f"e{e}"], regions[f"e{e + 1}"]
                if self.strand < 0:
                    regions[f"i{e}"] = GRange(self.region.chr, r2[2], r1[1])
                else:
                    regions[f"i{e}"] = GRange(self.region.chr, r1[2], r2[1])
            regions = sorted(regions.items(), key=lambda x: (x[1].start, x[1].end))
            self.regions[i] = dict(regions if self.strand > 0 else regions[::-1])

        #: dict[int, (int, str)]:
        #: reverse lookup (gene, region) of gene regions given a location
        #: within the reference genome.
        self._region_at = {
            i: (g, r)
            for g, d in self.regions.items()
            for r, rng in d.items()
            for i in range(rng.start, rng.end)
        }

        self.unique_regions = yml["structure"]["cn_regions"]

    def _init_alleles(self, yml) -> None:
        """
        Initialize allele (``alleles``, ``common_tandems``)
        and copy number configurations (``cn_configs``).
        Initialize old notation lookup table (``old_notation``).

        Raises:
            :obj:`aldy.common.AldyException` if allele name cannot be found.
        """

        alleles = {}
        fusions_left = {}
        fusions_right = {}

        # Human-readable CN descriptions
        descriptions = {
            "1": "Normal allelic configuration: "
            + "all regions present in both gene and pseudogene"
        }
        allele_metadata = {}

        # allele ID of the deletion allele (i.e. whole gene is missing).
        deletion_allele = None

        self.mutation_info = {}
        for allele_name, allele in yml["alleles"].items():
            # if not re.match(r"^[A-Z0-9]+\*[A-Z0-9]+$", allele_name):
            #     raise AldyException(
            #         "Allele names must be in format (alphanum+)*(alphanum+)"
            #         + "(e.g. CYP21*2A, DPYD*NEW)"
            #     )
            if '*' in allele_name:
                allele_name = allele_name.split('*')[1]
            mutations: List[Mutation] = []
            if [self.name, "deletion"] in allele["mutations"]:
                deletion_allele = allele_name
                descriptions[allele_name] = "Gene deletion"
                mutations = []
            else:
                for m in allele["mutations"]:
                    if m[0] in self.pseudogenes:  # has only one mutation indicator
                        assert m[0] == self.pseudogenes[0]  # TODO: relax
                        op = m[1]
                        if op[-1] == '-':
                            fusions_left[allele_name] = op[:-1]
                            descriptions[allele_name] = (
                                "Fusion: pseudogene until "
                                + op[:-1][::-1]
                                + " followed by the gene"
                            )
                        else:
                            if op[-1] == "+":
                                op = op[:-1]
                            fusions_right[allele_name] = op
                            descriptions[allele_name] = (
                                f"Conservation: Pseudogene retention after "
                                + op[::-1]
                                + " within the gene"
                            )
                    else:
                        (pos, op, rsid), function = m[:3], m[3] if len(m) > 3 else None
                        if self.strand < 0:
                            if op[1] == '>':
                                op = f'{rev_comp(op[0])}>{rev_comp(op[2])}'
                            elif op[:3] == 'ins':
                                ins = op[3:]
                                while self.seq[pos - len(ins):pos] == ins:
                                    pos -= len(ins)
                                pos += 1
                                op = f'ins{rev_comp(op[3:])}'
                            elif op[:3] == 'del':
                                pos = pos + len(op) - 4
                        if op[:3] == 'del':
                            op = 'del' + 'N' * (len(op) - 3)
                        pos -= 1  # Cast to 0-based index
                        if pos not in self.ref_to_chr:
                            log.warn('___', allele_name, m)
                        else:
                            mutations.append(Mutation(self.ref_to_chr[pos], op))
                            self.mutation_info.setdefault((self.ref_to_chr[pos], op), (function, rsid))
                            r = self.region_at(self.ref_to_chr[pos])
                            if r[1] == '-':
                                print(self.ref_to_chr[pos])
                                print(m)
                                assert False
            alleles[allele_name] = MinorAllele(allele_name, allele.get('label', None), mutations)

        # TODO:
        self.common_tandems: List[tuple] = []
        if "common_tandems" in yml:
            self.common_tandems = [tuple(y) for y in yml["common_tandems"]]

        # Set copy number configurations
        # TODO: currently fusions only cover the space between the main gene and
        # the first pseudogene. Multi-pseudogene fusions are not supported.

        self.cn_configs: Dict[str, CNConfig] = dict()

        def freezekey(x):  # hashing for dictionaries
            return tuple(i[1] for i in sorted(x[0].items())) + tuple(
                i[1] for i in sorted(x[1].items())
            )


        for i in range(len(self.pseudogenes)):
            if list(self.regions[i + 1].keys()) != list(self.regions[0].keys()):
                raise AldyException('Invalid database structure')
        region_rank = {a: ai for ai, a in enumerate(self.regions[0])}
        inverse_cn: Dict[tuple, str] = dict()
        # Left fusions are PSEUDOGENE + GENE fusions
        for a, brk in fusions_left.items():
            cn = dict()
            cn[0] = { r: int(region_rank[r] >= region_rank[brk]) for r in self.regions[0] }
            cn[1] = { r: int(region_rank[r] < region_rank[brk]) for r in self.regions[1] }

            key = freezekey(cn)
            if key not in inverse_cn:
                self.cn_configs[a] = CNConfig(
                    cn, CNConfig.CNConfigType.LEFT_FUSION, {a}, descriptions[a]
                )
                inverse_cn[key] = a
            else:
                self.cn_configs[inverse_cn[key]].alleles.add(a)

        # Deletion is a special kind of left fusion
        if deletion_allele is not None:
            cn = { 0: {r: 0 for r in self.regions[0]} }
            if len(self.pseudogenes) > 0:
                cn[1] = {r: 1 for r in self.regions[1]}
            self.cn_configs[deletion_allele] = CNConfig(
                cn,
                CNConfig.CNConfigType.DELETION,
                {deletion_allele},
                descriptions[deletion_allele],
            )

        # Right fusions GENE + PSEUDOGENE + whole copy of PSEUDOGENE fusions
        for a, brk in fusions_right.items():
            cn = dict()
            cn[0] = {r: int(region_rank[r] < region_rank[brk]) for r in self.regions[0]}
            cn[1] = {r: 1 + int(region_rank[r] >= region_rank[brk]) for r in self.regions[1]}
            key = freezekey(cn)
            if key not in inverse_cn:
                self.cn_configs[a] = CNConfig(
                    cn, CNConfig.CNConfigType.RIGHT_FUSION, {a}, descriptions[a]
                )
                inverse_cn[key] = a
            else:
                self.cn_configs[inverse_cn[key]].alleles.add(a)

        # Normal CN case
        used_alleles = {
            a for _, (_, _, alleles, _) in self.cn_configs.items() for a in alleles
        }
        has_pseudogenes = len(self.pseudogenes) > 0
        default_cn = {
            g: {r: 1 for r in self.regions[g]} for g in range(1 + has_pseudogenes)
        }

        self.cn_configs = {
            allele_number(min(v.alleles)): v for k, v in self.cn_configs.items()
        }
        self.cn_configs["1"] = CNConfig(
            default_cn,
            CNConfig.CNConfigType.DEFAULT_CN,
            set(alleles.keys()) - used_alleles,
            descriptions["1"],
        )

        # Remove "empty" regions
        for c, cn in self.cn_configs.items():
            for g, r in cn.cn.items():
                for rg in r:
                    if self.regions[g][rg].end - self.regions[g][rg].start <= 0:
                        r[rg] = 0

        # Set up major and minor allele structures
        alleles_inverse: Dict[tuple, set] = collections.defaultdict(set)
        for a in alleles:
            cn_config = next(
                cn for cn, conf in self.cn_configs.items() if a in conf.alleles
            )
            fn_muts = sorted(
                m for m in alleles[a].neutral_muts if self.is_functional(m)
            )
            alleles_inverse[(cn_config, tuple(fn_muts))].add(a)

        self.alleles: Dict[str, MajorAllele] = dict()
        # Karolinska DB has a lot of ambiguities, and two alleles with the same 'number'
        # can correspond to two different major alleles.
        # This step fixes this ambiguity by prepending a different letter to each unique
        # major allele that has number collision.
        used_names: Dict[str, int] = {}
        for key, minors in alleles_inverse.items():
            an = min(minors)
            a = MajorAllele(name="", cn_config=key[0], func_muts=key[1], minors=minors)
            name = an.split('.')[0] # allele_number(an)
            if name in used_names:  # Append letter
                used_names[name] += 1
                # log.warn(f"{name} -> {name}:{used_names[name]}")
                name += f":{used_names[name]}"
                used_names[name] = 1
            else:
                used_names[name] = 1
            self.alleles[name] = MajorAllele(
                name=name,
                cn_config=a.cn_config,
                func_muts=set(a.func_muts),
                minors={
                    sa: MinorAllele(
                        name=sa,
                        alt_name=alleles[sa].alt_name,
                        neutral_muts=set(alleles[sa].neutral_muts) - set(a.func_muts),
                    )
                    for sa in a.minors
                },
            )

    def _init_partials(self) -> None:
        """
        Construct "partial" major alleles.
        If a major allele is cut in half by a fusion, we will create a "new" major
        allele that contains the functional mutations that have survived
        the fusion event.

        Note:
            This is currently supported only for left fusions.
        """

        def preserved_mutations(f, m):
            def filter_f(m):
                gene, region = self.region_at(m.pos)
                return self.cn_configs[f].cn[gene][region] > 0

            return set(filter(filter_f, m))

        for f in filter(
            lambda x: self.cn_configs[x].kind == CNConfig.CNConfigType.LEFT_FUSION,
            self.cn_configs,
        ):
            # Do not extend left fusions that are already defined
            # by some functional mutations
            # HACK: This is hacky but works well for now
            if len(self.alleles[f].func_muts) > 0:
                continue
            add: Dict[tuple, MajorAllele] = {}
            for an, a in self.alleles.items():
                # We only make partial major alleles from non-fused major alleles
                if a.cn_config != "1":
                    continue
                new_name = "{}/{}".format(f, an)
                new_muts = preserved_mutations(f, a.func_muts)
                key = sorted_tuple(new_muts)
                if key in add:
                    add[key].minors.update(
                        {
                            san: MinorAllele(
                                san, [], preserved_mutations(f, sa.neutral_muts)
                            )
                            for san, sa in a.minors.items()
                        }
                    )
                else:
                    add[key] = MajorAllele(
                        name=new_name,
                        cn_config=f,
                        func_muts=new_muts,
                        minors={
                            san: MinorAllele(
                                name=san,
                                alt_name=[],
                                neutral_muts=preserved_mutations(f, sa.neutral_muts),
                            )
                            for san, sa in a.minors.items()
                        },
                    )

            # Remove fusion (will be replaced at least by allele "1/{f}")
            del self.alleles[f]
            self.alleles.update({a.name: a for a in add.values()})

        for an, a in self.alleles.items():
            # Clean up minor alleles (as many might be identical after a fusion).
            # Put a reference to the cleaned-up alleles in ``alt_name`` field.
            minors: Dict[tuple, List[str]] = collections.defaultdict(list)
            for s in a.minors:
                key = sorted_tuple(a.minors[s].neutral_muts)
                minors[key].append(s)
            self.alleles[an] = MajorAllele(
                self.alleles[an].name,
                self.alleles[an].cn_config,
                self.alleles[an].func_muts,
                {
                    min(sa): MinorAllele(
                        min(sa),
                        alt_name=a.minors[min(sa)].alt_name, #sorted(list(set(sa) - {min(sa)})),
                        neutral_muts=set(nm),
                    )
                    for nm, sa in minors.items()
                },
            )

        # TODO: prune identical post-fusion alleles
        # (e.g. *2 after *13A and *13B is the same--- no need to have 2 items).

        for _, v in self.cn_configs.items():
            v.alleles.clear()
        for a in self.alleles.values():
            self.cn_configs[a.cn_config].alleles.add(a.name)

    def _init_structure(self, yml) -> None:
        """
        Initialize ``mutations`` lookup table
        and coding region structure ``coding_region``.
        """

        # is_rev_comp = yml["revcomp"]
        # lookup = {
        #     i: self.seq[i - self.region.start]
        #     for (_, kind), (_, start, end) in self.regions[0].items()
        #     if kind == EXON
        #     for i in range(start, end)
        # }
        self.exons = yml["reference"]["exons"]
        seq = "".join(self.seq[s - 1:e - 1] for [s, e] in self.exons)
        aminoacid = seq_to_amino(rev_comp(seq) if self.strand < 0 else seq)

        # Set up a coding region structure for aminoacid calculation
        self.coding_region = Gene.CodingRegion(self.strand < 0, aminoacid)

        self.mutations: Dict[Tuple[int, str], Mutation] = {}
        for _, a in self.alleles.items():
            self.mutations.update({(m.pos, m.op): m for m in a.func_muts})
            for _, sa in a.minors.items():
                self.mutations.update({(m.pos, m.op): m for m in sa.neutral_muts})

    # -----------------------------------------------------------------------------------

    def region_at(self, pos: int) -> Tuple[int, str]:
        """
        Returns:
            (int, str): Tuple consisting of a gene ID and
            a region that harbours the given position in the gene.

        Args:
            pos (int): Position in the reference genome.
        """
        return self._region_at.get(pos, (0, "-"))

    def is_functional(self, mut, infer=True) -> bool:
        """
        Returns:
            bool: ``True`` if a mutation is functional
            (i.e. does it affect the underlying aminoacid or not).
        """
        pos, op = mut
        if (pos, op) in self.mutation_info:
            return (self.mutation_info[pos, op][0] is not None)

        # Calculate based on aminoacid change
        if infer and any(self.ref_to_chr[s] <= pos < self.ref_to_chr[e] for s, e in self.exons):
            if op[:3] != "SNP":
                return True
            # o = self.region.start
            # seq = "".join(
            #     self.seq[s: pos] + op[5] + self.seq[pos - o + 1 : e - o]
            #     if self.ref_to_chr[s] <= pos < self.ref_to_chr[e]
            #     else self.seq[s - o : e - o]
            #     for s, e in self.exons
            # )
            # amino = seq_to_amino(rev_comp(seq) if self.coding_region.rev_comp else seq)
            # return amino != self.coding_region.aminoacid
        return False

    def get_dbsnp(self, *args) -> bool:
        """
        Returns:
            bool: ``True`` if a mutation is functional
            (i.e. does it affect the underlying aminoacid or not).
        """
        assert 1 <= len(args) <= 2
        if len(args) == 2:
            pos, op = args
        else:
            pos, op = args[0]
        if (pos, op) in self.mutation_info:
            return self.mutation_info[pos, op][1]
        return "-"

    def deletion_allele(self) -> Optional[str]:
        """
        Returns:
            str, optional: The deletion allele ID. Can be ``None``
            if gene has no deletion allele.
        """
        try:
            return next(
                a
                for a, cn in self.cn_configs.items()
                if cn.kind == CNConfig.CNConfigType.DELETION
            )
        except StopIteration:
            return None

    def has_coverage(self, a: str, pos: int):
        """
        Returns: ``True`` if a major allele `a` covers the mutation `m`.
        """
        m_gene, m_region = self.region_at(pos)
        return self.cn_configs[self.alleles[a].cn_config].cn[m_gene][m_region] > 0

    def get_wide_region(self):
        mi = min(r.start for g in self.regions.values() for r in g.values())
        ma = max(r.end for g in self.regions.values() for r in g.values())
        return GRange(self.region.chr, mi, ma)

    def __contains__(self, i: int):
        return i in self.chr_to_ref

    def __getitem__(self, i: int):
        s = self.seq[self.chr_to_ref[i]]
        return rev_comp(s) if self.strand < 0 else s

    def print_summary(self, full=True):
        log.info(f"Gene {self.name}")
        log.info(f"  {self.genome} genome locus: {self.region}")
        pseudo = ", ".join(f"{p} (ID {i + 1})" for i, p in enumerate(self.pseudogenes))
        log.info(f"  Pseudogenes: {pseudo if pseudo else 'none'}")
        amino = "\n    ".join(textwrap.wrap(self.coding_region.aminoacid))
        log.info(f"  Aminoacid:\n    {amino}")
        log.info("  Regions:")

        regions = set(j for _, i in self.regions.items() for j in i)
        log.info("    {:>10}  {:>25} {:>25}", "Region", "Gene", "Pseudogene")

        for r in sorted(
            regions, key=lambda r: self.regions[0].get(r, (0, 99999999))[1]
        ):
            lg = str(self.regions[0].get(r, "-"))
            lp = str(self.regions.get(1, {}).get(r, "-"))
            log.info(f"    {r:>10}: {lg:>25} {lp:>25}")

        log.info(f"  Copy number configurations:")
        a = "\n    ".join(
            textwrap.wrap(", ".join(sorted(self.cn_configs, key=allele_sort_key)))
        )
        log.info(f"    {a}")

        log.info("  Major star-alleles:")
        a = "\n    ".join(
            textwrap.wrap(
                ", ".join(f"*{a}" for a in sorted(self.alleles, key=allele_sort_key))
            )
        )
        log.info(f"    {a}")

        if full:
            for c in self.cn_configs:
                self.print_cns(c)
            for major in self.alleles:
                self.print_majors(major, full)
            # ...

    def print_cns(self, cn_config):
        if cn_config not in self.cn_configs:
            raise AldyException(f"Copy number configuration {cn_config} not found")
        config = self.cn_configs[cn_config]
        log.info(f"Gene {self.name}")
        log.info(f"  Copy number configuration {cn_config}:\n  {config.description}")

        alleles = sorted(config.alleles, key=allele_sort_key)
        log.info("  Major star-alleles:")
        a = "\n    ".join(textwrap.wrap(", ".join(f"*{a}" for a in alleles)))
        log.info(f"    {a}")

        log.info("  Copy number profile:")
        regions = set(j for i in config.cn for j in config.cn[i])
        log.info("    {:10}    Gene Pseudo", " ")
        for r in sorted(
            regions, key=lambda r: self.regions[0].get(r, (0, 99999999))[1]
        ):
            if r.kind == EXON:
                name = f"exon {r.number}"
            elif r.kind == INTRON:
                name = f"intron {r.number}"
            else:
                name = r.kind
            lg = config.cn[0].get(r, 0.0)
            lp = config.cn.get(1, {}).get(r, 0.0)
            lg = "-" if lg == 0 else f"{lg:.0f}"
            lp = "-" if lp == 0 else f"{lp:.0f}"
            log.info(f"    {name:>10}: {lg:>6} {lp:>6}")

    def print_majors(self, major, full=False):
        if major not in self.alleles:
            raise AldyException(f"Major star-allele {major} not found")
        allele = self.alleles[major]
        # MajorAllele
        log.info(f"Gene {self.name}, major star-allele {self.name}*{allele.name}:")
        log.info(f"  Copy number configuration {allele.cn_config}")

        text = [
            f"{self.region.chr}:{m.pos}:{m.op} ({self.get_dbsnp(m)})"
            for m in sorted(allele.func_muts)
        ]
        log.info(
            f"  Functional mutations: {', '.join(text).replace(' (-)', '') if text else 'none'}"
        )

        log.info("  Minor star-alleles:")
        a = "\n    ".join(
            textwrap.wrap(", ".join(sorted(f"*{a}" for a in allele.minors)))
        )
        log.info(f"    {a}")

        if full:
            for a in allele.minors:
                self.print_minors(a)
            log.info("-----------------")

    def print_minors(self, minor):
        found = False
        for major, allele in sorted(self.alleles.items()):
            if minor not in allele.minors:
                continue
            found = True
            log.info(f"\nGene {self.name}, minor star-allele {self.name}*{minor}:")
            log.info(f"  Major star-allele: {self.name}*{allele.name}")
            log.info(f"  Copy number configuration: {allele.cn_config}")
            if len(allele.minors[minor].alt_name) > 0:
                alt = ", ".join(sorted(allele.minors[minor].alt_name))
                log.info(f"  Alternative names: {alt}")

            text = [
                f"{self.region.chr}:{m.pos}:{m.op} ({self.get_dbsnp(m)})"
                for m in sorted(allele.func_muts)
            ]
            log.info(
                f"  Functional mutations: {', '.join(text).replace(' (-)', '') if text else 'none'}"
            )
            text = [
                f"{self.region.chr}:{m.pos}:{m.op} ({self.get_dbsnp(m)})"
                for m in sorted(allele.minors[minor].neutral_muts)
            ]
            log.info(
                f"  Silent mutations: {', '.join(text).replace(' (-)', '') if text else 'none'}"
            )

        if not found:
            raise AldyException(f"Minor star-allele {minor} not found")
