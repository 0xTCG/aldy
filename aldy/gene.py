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

from .common import (
    GeneRegion,
    GRange,
    AldyException,
    allele_number,
    sorted_tuple,
    seq_to_amino,
    rev_comp,
)


EXON = "e"
"""str: Abbreviation for exon."""

INTRON = "i"
"""str: Abbreviation for intron."""


class Allele(
    collections.namedtuple("Allele", ["name", "cn_config", "func_muts", "minors"])
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
        minors (dict[str, :obj:`Suballele`]):
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


class Suballele(
    collections.namedtuple("Suballele", ["name", "alt_names", "neutral_muts"])
):
    """
    Minor allele description.
    Immutable.

    Attributes:
        name (str):
            Minor allale name.
        alt_names (list[str]):
            List of alternative names.
        neutral_muts (set[:obj:`Mutation`]):
            Netural mutations that describe the minor allele.

    Notes:
        Has custom printer (``__repr__``).
    """

    def __str__(self):
        return "Sub({}; [{}])".format(
            "|".join([self.name] + self.alt_names),
            ", ".join(map(str, self.neutral_muts)),
        )


@functools.total_ordering
class Mutation(
    collections.namedtuple("Mutation", ["pos", "op", "is_functional", "aux"])
):
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

    def __new__(self, pos: int, op: str, is_functional=0, aux=None):
        return super(Mutation, self).__new__(
            self, pos, op, is_functional, aux if aux else dict()
        )

    def __str__(self):
        return "{}:{}{}".format(self.pos, self.op, "*" if self.is_functional else "")

    def __eq__(self, other):
        return (self.pos, self.op) == (other.pos, other.op)

    def __lt__(self, other):
        return (self.pos, self.op) < (other.pos, other.op)

    def __hash__(self):
        return (self.pos, self.op).__hash__()


class CNConfig(
    collections.namedtuple("CNConfig", ["cn", "kind", "alleles", "description"])
):
    """
    Copy-number (CN) configuration description.
    Immutable.

    Attributes:
        cn (dict[int, dict[:obj:`aldy.common.GeneRegion`, int]]):
            Value of the expected region copy number in each gene.
            For example, ``cn[0][GeneRegion(0, 1, EXON)] == 1`` means that the exon 1 of
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
        regions = sorted(set(r for g in self.cn for r in self.cn[g]))
        return "CNConfig({}; vector={}; alleles=[{}])".format(
            str(self.kind)[13:],
            "|".join(
                "".join(
                    ("{:.0f}".format(self.cn[g][r]) if self.cn[g][r] != 0.5 else "½")
                    if r in self.cn[g]
                    else "_"
                    for r in regions
                )
                for g in sorted(self.cn)
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
        regions (dict[int, dict[:obj:`aldy.common.GeneRegion`,
                :obj:`aldy.common.GRange`]]):
            Exonic, intronic and special (e.g. UTR) regions in a gene.
            Maps a gene ID to a dictionary that maps :obj:`aldy.common.GeneRegion`
            (e.g. exon 9) to a :obj:`aldy.common.GRange` (e.g. chr1:10-20) within
            the reference genome.
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

    CodingRegion = collections.namedtuple(
        "CodingRegion", ["rev_comp", "aminoacid", "lookup"]
    )
    """
    Immutable class that describes the coding region of a gene.
    Attributes: aminoacid sequence (``aminoacid``) and reference strand (``rev_comp``).
    """

    def __init__(
        self, path: Optional[str], name: Optional[str] = None, yml: Optional[str] = None
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
        self.__parse_yml(name, yml)

    def __parse_yml(self, gene_name: str, yml) -> None:
        """
        Initialize the gene structure from YML data.
        """

        self.__init_basic(gene_name, yml)
        self.__init_regions(yml)
        self.__init_alleles(yml)
        self.__init_partials()
        self.__init_structure(yml)

    def __init_basic(self, gene: str, yml) -> None:
        """
        Read basic gene properties (``name``, ``seq`` and ``region``).
        """
        self.name = gene.upper()
        self.seq = yml["seq"]
        self.region = GRange(*yml["region"])

    def __init_regions(self, yml) -> None:
        """
        Calculate the genic regions and pseudogenes
        (``regions``, ``unique_regions`` and ``pseudogenes``).
        Prepare ``region_at()`` call.
        """
        assert hasattr(self, "region")

        def calculate_regions(gene, exons, special):
            """
            Given the list of exons, calculate intronic regions and
            return a dictionary of all regions.
            """
            regions = {
                GeneRegion(e, EXON): GRange(self.region.chr, *ex)
                for e, ex in exons.items()
            }
            # Fill introns
            for i in range(1, len(regions)):
                r1, r2 = regions[GeneRegion(i + 1, EXON)], regions[GeneRegion(i, EXON)]
                if r1 > r2:
                    r1, r2 = r2, r1
                regions[GeneRegion(i, INTRON)] = GRange(self.region.chr, r1[2], r2[1])
            for n, r in special.items():
                regions[GeneRegion(n, r[0])] = GRange(self.region.chr, r[1], r[2])
            return regions

        # Gene 0 is the main gene (key is the gene ID)
        self.regions = {0: calculate_regions(0, yml["exons"], yml["special_regions"])}

        # Each pseudogene is associated with an index > 0
        self.pseudogenes: List[str] = list()
        if "pseudogenes" in yml:
            for gi, g in enumerate(sorted(yml["pseudogenes"])):
                self.pseudogenes.append(g)
                self.regions[gi + 1] = calculate_regions(
                    gi + 1,
                    yml["pseudogenes"][g]["exons"],
                    yml["pseudogenes"][g]["special_regions"],
                )

        #: dict[int, (int, `GeneRegion`]):
        #: reverse lookup (gene, region) of gene regions given a location
        #: within the reference genome.
        self._region_at = {
            i: (g, r)
            for g, d in self.regions.items()
            for r, rng in d.items()
            for i in range(rng.start, rng.end)
        }

        def parse_unique(u):
            u = list(filter(None, re.split(r"^(\d+)", u)))
            if len(u) > 1:
                return GeneRegion(int(u[0]), u[1])
            else:  # find the matching kind number (legacy support)
                for g in self.regions:
                    for r in self.regions[g]:
                        if u[0] == r.kind:
                            return r
            raise KeyError

        self.unique_regions = [parse_unique(y) for y in yml["unique_regions"]]

    def __init_alleles(self, yml) -> None:
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

        #: dict(:obj:`Mutation`, str):
        #: old Karolinska notation (e.g. 32C>T) for a :obj:`Mutation`.
        self._old_notation: Dict[Mutation, str] = collections.defaultdict(str)

        # allele ID of the deletion allele (i.e. whole gene is missing).
        deletion_allele = None

        for allele_name, allele in yml["alleles"].items():
            if not re.match(r"^[A-Z0-9]+\*[0-9]+[A-Z0-9]*$", allele_name):
                raise AldyException(
                    "Allele names must be in format (alphanum?)*(number)(alphanum?) "
                    + "(e.g. CYP21*2A, *3)"
                )
            allele_name = allele_name.split("*")[1]
            mutations: List[Mutation] = []
            if {"op": "deletion"} in allele["mutations"]:
                deletion_allele = allele_name
                descriptions[allele_name] = "Gene deletion"
                mutations = []
            else:
                for m in allele["mutations"]:
                    if m["pos"] == "pseudogene":  # has only one mutation indicator
                        if "-" in m["op"]:
                            fusions_left[allele_name] = (
                                int(m["op"][1:-1]),
                                (EXON if m["op"][0] == "e" else INTRON),
                            )
                            descriptions[allele_name] = (
                                "Fusion: pseudogene until "
                                + m["op"][:-1][::-1]
                                + " followed by the gene"
                            )
                        else:
                            if m["op"][-1] == "+":
                                m["op"] = m["op"][:-1]
                            fusions_right[allele_name] = (
                                int(m["op"][1:]),
                                (EXON if m["op"][0] == "e" else INTRON),
                            )
                            descriptions[allele_name] = (
                                f"Conservation: Pseudogene retention after "
                                + m["op"][::-1]
                                + " within the gene"
                            )
                    else:
                        m["aux"] = {}
                        if "dbsnp" in m:
                            m["aux"]["dbsnp"] = " or ".join(m.get("dbsnp", []))
                        if "old" in m:
                            m["aux"]["old"] = m.get("old", "")
                        mut = Mutation(m["pos"], m["op"], m["functional"], m["aux"])
                        if "old" in m:
                            self._old_notation[mut] = m["old"]
                        mutations.append(mut)
            alleles[allele_name] = Suballele(allele_name, [], mutations)

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

        inverse_cn: Dict[tuple, str] = dict()

        # Left fusions are PSEUDOGENE + GENE fusions
        for a, brk in fusions_left.items():
            cn = dict()
            cn[0] = {
                r: float(0 if (r.number, r.kind) < brk else 1) for r in self.regions[0]
            }
            cn[1] = {
                r: float(1 if (r.number, r.kind) < brk else 0) for r in self.regions[1]
            }

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
            cn = {
                0: {r: 0 for r in self.regions[0]},
                1: {r: 1 for r in self.regions[1]},
            }
            self.cn_configs[deletion_allele] = CNConfig(
                cn,
                CNConfig.CNConfigType.DELETION,
                {deletion_allele},
                descriptions[deletion_allele],
            )

        # Right fusions GENE + PSEUDOGENE + whole copy of PSEUDOGENE fusions
        for a, brk in fusions_right.items():
            cn = dict()
            cn[0] = {r: (1 if (r.number, r.kind) < brk else 0) for r in self.regions[0]}
            cn[1] = {r: (1 if (r.number, r.kind) < brk else 2) for r in self.regions[1]}
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

        # Set up major and minor allele structures
        alleles_inverse: Dict[tuple, set] = collections.defaultdict(set)
        for a in alleles:
            cn_config = next(
                cn for cn, conf in self.cn_configs.items() if a in conf.alleles
            )
            fn_muts = sorted(m for m in alleles[a].neutral_muts if m.is_functional)
            alleles_inverse[(cn_config, tuple(fn_muts))].add(a)

        self.alleles: Dict[str, Allele] = dict()
        # Karolinska DB has a lot of ambiguities, and two alleles with the same 'number'
        # can correspond to two different major alleles.
        # This step fixes this ambiguity by prepending a different letter to each unique
        # major allele that has number collision.
        used_names: Dict[str, int] = {}
        for key, minors in alleles_inverse.items():
            an = min(minors)
            a = Allele(name="", cn_config=key[0], func_muts=key[1], minors=minors)
            name = allele_number(an)
            if name in used_names:  # Append letter
                used_names[name] += 1
                name += "." + chr(used_names[name] + ord("a"))
            else:
                used_names[name] = -1
            self.alleles[name] = Allele(
                name=name,
                cn_config=a.cn_config,
                func_muts=set(a.func_muts),
                minors={
                    sa: Suballele(
                        name=sa,
                        alt_names=[],
                        neutral_muts=set(alleles[sa].neutral_muts) - set(a.func_muts),
                    )
                    for sa in a.minors
                },
            )

    def __init_partials(self) -> None:
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
            add: Dict[tuple, Allele] = {}
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
                            san: Suballele(
                                san, [], preserved_mutations(f, sa.neutral_muts)
                            )
                            for san, sa in a.minors.items()
                        }
                    )
                else:
                    add[key] = Allele(
                        name=new_name,
                        cn_config=f,
                        func_muts=new_muts,
                        minors={
                            san: Suballele(
                                name=san,
                                alt_names=[],
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
            # Put a reference to the cleaned-up alleles in ``alt_names`` field.
            minors: Dict[tuple, List[str]] = collections.defaultdict(list)
            for s in a.minors:
                key = sorted_tuple(a.minors[s].neutral_muts)
                minors[key].append(s)
            self.alleles[an] = Allele(
                self.alleles[an].name,
                self.alleles[an].cn_config,
                self.alleles[an].func_muts,
                {
                    min(sa): Suballele(
                        min(sa),
                        alt_names=sorted(list(set(sa) - {min(sa)})),
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

    def __init_structure(self, yml) -> None:
        """
        Initialize ``mutations`` lookup table
        and coding region structure ``coding_region``.
        """

        is_rev_comp = True if yml["rev_comp"] == 1 else False
        lookup = {
            i: self.seq[i - self.region.start]
            for (_, kind), (_, start, end) in self.regions[0].items()
            if kind == EXON
            for i in range(start, end)
        }
        seq = "".join(lookup[i] for i in sorted(lookup.keys()))
        aminoacid = seq_to_amino(rev_comp(seq) if is_rev_comp else seq)

        # Set up a coding region structure for aminoacid calculation
        self.coding_region = Gene.CodingRegion(is_rev_comp, aminoacid, lookup)

        self.mutations: Dict[Tuple[int, str], Mutation] = {}
        for _, a in self.alleles.items():
            self.mutations.update({(m.pos, m.op): m for m in a.func_muts})
            for _, sa in a.minors.items():
                self.mutations.update({(m.pos, m.op): m for m in sa.neutral_muts})

    # -----------------------------------------------------------------------------------

    def region_at(self, pos: int) -> Tuple[int, GeneRegion]:
        """
        Returns:
            (int, :obj:`aldy.common.GeneRegion`): Tuple consisting of a gene ID and
            a region that harbours the given position in the gene.

        Args:
            pos (int): Position in the reference genome.
        """
        return self._region_at[pos]

    def old_notation(self, mut: Mutation) -> str:
        """
        Returns:
            str: Karolinska-style notation (e.g. ``32C>T``) for a :obj:`Mutation`.

        Args:
            mut (:obj:`Mutation`)
        """
        return self._old_notation[mut]

    def sequence(self, gene: int, pad=0) -> str:
        """
        Returns:
            str: Genomic sequence of a gene.

        Args:
            id (int): Gene ID (0 for the main gene, 1+ for pseudogenes).
            pad (int): Extra padding on both sides (measured in bps). Default is 0.
        """
        regs = [(st, ed) for _, (_, st, ed) in self.regions[gene].items()]
        mi, ma = min(r[0] for r in regs), max(r[1] for r in regs)
        off = self.region.start
        return self.seq[max(0, mi - off - pad) : min(ma - off + pad, len(self.seq))]

    def is_functional(self, m: Mutation) -> bool:
        """
        Returns:
            bool: ``True`` if a mutation is functional
            (i.e. does it affect the underlying aminoacid or not).
        """
        if m.pos not in self.coding_region.lookup:
            return False
        if m.op[:3] != "SNP":
            return True
        seq = "".join(
            self.coding_region.lookup[i] if i != m.pos else m.op[5]
            for i in sorted(self.coding_region.lookup)
        )
        amino = seq_to_amino(rev_comp(seq) if self.coding_region.rev_comp else seq)
        return amino != self.coding_region.aminoacid

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

    def print_configurations(self):
        raise NotImplementedError

        # TODO: do this
        # log.error('Gene {}\n', self.name)
        # for name, config in sorted(self.cnv_configurations.items()):
        #    c6 = [c[2:] for c in config if c.startswith('6.')]
        #    c7 = [c[2:] for c in config if c.startswith('7.')]
        #    labels = sorted(list(set(c6) | set(c7)))
        #    max_width = str(max(3, max(len(x) for x in labels)) + 1)

        #    description = '{}\n'.format(self.descriptions[name]) if name in self.descriptions else ''
        #    alleles = sorted(set('*' + a.split('/')[0] for a in self.alleles if self.alleles[a].cnv_configuration == name))
        #    description += 'Alleles: '
        #    for i in range(0, len(alleles), 12):
        #       if i != 0: description += '\n         ';
        #       description += ', '.join(alleles[i:i + 12])

        #    log.warn('Name: {}\n{}'.format(name, description))
        #    print('            ', end='')
        #    for c in labels:
        #       print(('{:>' + max_width + '}').format(c), end='')
        #    print(' ')
        #    print('Gene:       ', end='')
        #    for c in labels:
        #       if '6.' + c in config:
        #          print(('{:' + max_width + '}').format(config['6.' + c]), end='')
        #    print(' ')
        #    print('Pseudogene: ', end='')
        #    for c in labels:
        #       if '7.' + c in config:
        #          print(('{:' + max_width + '}').format(config['7.' + c]), end='')
        #    print('\n')
