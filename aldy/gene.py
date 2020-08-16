# 786
# Aldy source: gene.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple, Dict, List, Optional

import os
import sys
import enum
import yaml
import collections
import textwrap

from collections import namedtuple as nt
from natsort import natsorted

from .common import (
    GRange,
    AldyException,
    allele_name,
    sorted_tuple,
    seq_to_amino,
    rev_comp,
    log,
    freezekey,
)


class MajorAllele(nt("MajorAllele", ["name", "cn_config", "func_muts", "minors"])):

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

    def get_minor_mutations(self, minor: str):
        """
        Sequence of mutations that describe a minor allele.
        """
        for m in self.func_muts:
            yield m
        for m in self.minors[minor].neutral_muts:
            yield m


class MinorAllele(
    nt(
        "MinorAllele",
        ["name", "alt_name", "neutral_muts", "activity", "evidence", "pharmvar"],
    )
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
        activity (Optional[str])
        evidence (Optional[str])
        pharmvar (Optional[str])

    Notes:
        Has custom printer (``__repr__``).
    """

    def __str__(self):
        muts = ", ".join(map(str, self.neutral_muts))
        return f"Minor({self.name}; [{muts}])"


class Mutation(nt("Mutation", ["pos", "op"])):
    """
    Mutation description.
    Immutable.

    Attributes:
        pos (int): Reference genome position (0-based).
        op (str): Variant operation in HGVS-like format:

            - ``X>Y``:  a SNP from X to Y
            - ``insX``:  an insertion of X
            - ``delX``:  a deletion of X (X, Y are both of format `[ACGT]+`)

    Notes:
        Has custom printer (``__str__``).
        Comparable and hashable via ``(pos, op)`` tuple.
        Implements ``total_ordering``.
    """

    def __str__(self):
        return f"{self.pos + 1}.{self.op}"


class CNConfig(nt("CNConfig", ["cn", "kind", "alleles", "description"])):
    """
    Copy-number (CN) configuration description.
    Immutable.

    Attributes:
        cn (list[dict[str, int]]):
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
        "CNConfigType", "DEFAULT LEFT_FUSION RIGHT_FUSION DELETION"
    )
    """
    Enumeration describing the type of a copy-number configuration:

        - ``LEFT_FUSION``
        - ``RIGHT_FUSION``
        - ``DELETION``
        - ``DEFAULT``
    """

    def __str__(self):
        cov = "|".join("".join(str(self.cn[g][r]) for r in self.cn[0]) for g in self.cn)
        alleles = " ".join(natsorted(self.alleles))
        return f"CNConfig({str(self.kind)[13:]}; vector={cov}; alleles=[{alleles}])"


class Gene:
    r"""
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
        self.version = f"{yml['version']} ({yml['generated']})"
        self.refseq = yml["reference"]["name"]
        self.pharmvar = yml.get("pharmvar", None)
        self.ensembl = yml.get("ensembl", None)

        self.seq = yml["reference"]["seq"].replace("\n", "")
        if "patches" in yml["reference"]:
            self.seq = list(self.seq)
            for [pos, nuc] in yml["reference"]["patches"]:
                self.seq[pos - 1] = nuc
            self.seq = "".join(self.seq)

        self.chr, start, end, strand, cigar = yml["reference"]["mappings"][self.genome]
        self.chr_to_ref = {}
        self.ref_to_chr = {}
        self.strand = 1 if strand == "+" else -1
        pos_ref = 0 if self.strand > 0 else (len(self.seq) - 1)
        pos_chr = start - 1
        for i in cigar.split():
            op, sz = i[0], int(i[1:])
            if op == "M":
                for idx in range(sz):
                    self.chr_to_ref[pos_chr + idx] = pos_ref + idx * self.strand
                    self.ref_to_chr[pos_ref + idx * self.strand] = pos_chr + idx
                pos_chr += sz
                pos_ref += sz * self.strand
            elif op == "I":
                pos_ref += sz * self.strand
            elif op == "D":
                pos_chr += sz
            else:
                raise AldyException("Invalid CIGAR string")

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
            for name, coord in yml["structure"]["regions"][self.genome].items():
                if name[0] == "e" and name[1:].isdigit():  # this is exon
                    num_exons += 1
                regions[name] = GRange(self.chr, coord[i * 2] - 1, coord[i * 2 + 1] - 1)
            for e in range(1, num_exons):  # fill introns
                r = [regions[f"e{e}"], regions[f"e{e + 1}"]][:: self.strand]
                regions[f"i{e}"] = GRange(self.chr, r[0][2], r[1][1])
            self.regions[i] = dict(
                sorted(regions.items(), key=lambda x: x[1])[:: self.strand]
            )

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

        # allele ID of the deletion allele (i.e. whole gene is missing).
        deletion_allele = None

        self.mutation_info: Dict[
            Tuple[int, str], Tuple[Optional[str], str, int, str]
        ] = {}
        for name, allele in yml["alleles"].items():
            name = allele_name(name)
            mutations: List[Mutation] = []
            if [self.name, "deletion"] in allele["mutations"]:
                deletion_allele = name
            else:
                for pos, op, *info in allele["mutations"]:
                    if pos in self.pseudogenes:
                        assert pos == self.pseudogenes[0]  # TODO: relax later
                        if op[-1] == "-":
                            fusions_left[name] = op[:-1]
                        else:
                            if op[-1] == "+":
                                op = op[:-1]
                            fusions_right[name] = op
                    else:
                        orig_op = op
                        rsid, function = info[0], info[1] if len(info) > 1 else None
                        if self.strand < 0:
                            if ">" in op:
                                l, r = op.split(">")
                                op = f"{rev_comp(l)}>{rev_comp(r)}"
                                pos = pos + len(l) - 1
                            elif op[:3] == "ins":
                                ins = op[3:]
                                while self.seq[pos - len(ins) : pos] == ins:
                                    pos -= len(ins)
                                pos += 1
                                op = f"ins{rev_comp(op[3:])}"
                            elif op[:3] == "del":
                                pos = pos + len(op) - 4
                                op = f"del{rev_comp(op[3:])}"
                        pos -= 1  # Cast to 0-based index
                        if pos not in self.ref_to_chr:
                            log.warn(
                                f"Ignoring {pos}.{op} in {name} (not in {self.refseq})"
                            )
                        elif self.region_at(self.ref_to_chr[pos])[1] == "-":
                            log.warn(
                                f"Ignoring {pos}.{op} in {name} (not in named region)"
                            )
                        else:
                            mutations.append(Mutation(self.ref_to_chr[pos], op))
                            self.mutation_info.setdefault(
                                (self.ref_to_chr[pos], op),
                                (function, rsid, pos, orig_op),
                            )
            alt_name = allele.get("label", None)
            if alt_name:
                alt_name = allele_name(alt_name)
            alleles[name] = MinorAllele(
                name,
                alt_name,
                mutations,
                allele.get("activity", None),
                allele.get("evidence", None),
                allele.get("pharmvar", None),
            )

        # TODO:
        self.common_tandems: List[tuple] = []
        if "common_tandems" in yml:
            self.common_tandems = [tuple(y) for y in yml["common_tandems"]]

        # Set copy number configurations
        # TODO: currently fusions only cover the space between the main gene and
        # the first pseudogene. Multi-pseudogene fusions are not supported.

        self.do_copy_number = deletion_allele or len(fusions_left) or len(fusions_right)
        self.cn_configs: Dict[str, CNConfig] = dict()

        for i, _ in enumerate(self.pseudogenes):
            if list(self.regions[i + 1].keys()) != list(self.regions[0].keys()):
                raise AldyException("Invalid database structure")
        rank = {a: ai for ai, a in enumerate(self.regions[0])}
        inverse_cn: Dict[tuple, str] = dict()
        # Left fusions are PSEUDOGENE + GENE fusions
        for a, brk in fusions_left.items():
            cn = [
                {r: int(rank[r] >= rank[brk]) for r in self.regions[0]},
                {r: int(rank[r] < rank[brk]) for r in self.regions[1]},
            ]
            key = freezekey(cn)
            if key not in inverse_cn:
                self.cn_configs[a] = CNConfig(
                    cn,
                    CNConfig.CNConfigType.LEFT_FUSION,
                    {a},
                    f"{self.pseudogenes[0]} fusion until {brk}",
                )
                inverse_cn[key] = a
            else:
                self.cn_configs[inverse_cn[key]].alleles.add(a)
        # Deletion is a special kind of left fusion
        if deletion_allele is not None:
            cn = [{r: 0 for r in self.regions[0]}]
            if len(self.pseudogenes) > 0:
                cn.append({r: 1 for r in self.regions[1]})
            self.cn_configs[deletion_allele] = CNConfig(
                cn,
                CNConfig.CNConfigType.DELETION,
                {deletion_allele},
                f"Whole {self.name} deletion",
            )
        # Right fusions GENE + PSEUDOGENE + whole copy of PSEUDOGENE fusions
        for a, brk in fusions_right.items():
            cn = [
                {r: int(rank[r] < rank[brk]) for r in self.regions[0]},
                {r: 1 + int(rank[r] >= rank[brk]) for r in self.regions[1]},
            ]
            key = freezekey(cn)
            if key not in inverse_cn:
                self.cn_configs[a] = CNConfig(
                    cn,
                    CNConfig.CNConfigType.RIGHT_FUSION,
                    {a},
                    f"{self.pseudogenes[0]} conservation after {brk}",
                )
                inverse_cn[key] = a
            else:
                self.cn_configs[inverse_cn[key]].alleles.add(a)

        # Normal CN case
        used_alleles = {
            a for _, (_, _, alleles, _) in self.cn_configs.items() for a in alleles
        }
        has_pseudogenes = len(self.pseudogenes) > 0
        default_cn = [
            {r: 1 for r in self.regions[g]} for g in range(1 + has_pseudogenes)
        ]
        self.cn_configs = {min(v.alleles): v for k, v in self.cn_configs.items()}
        self.cn_configs["1"] = CNConfig(
            default_cn,
            CNConfig.CNConfigType.DEFAULT,
            set(alleles.keys()) - used_alleles,
            "Standard copy-number configuration",
        )

        # Set CN of "empty" regions to zero (e.g. CYP2D6.pce)
        for _, conf in self.cn_configs.items():
            for g, r in enumerate(conf.cn):
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
            alleles_inverse[cn_config, tuple(fn_muts)].add(a)

        # Ensure that each distinct major allele has unique prefix
        used_names: Dict[str, int] = {}
        new_names = {}
        changed_cn_configs = {}
        for key, minors in alleles_inverse.items():
            an = min(minors)
            name = an.split(".")[0]  # remove ".0XX" suffix
            if name in used_names:  # special case
                name = alleles[an].alt_name if alleles[an].alt_name else an
            if name in used_names:  # Append letter
                used_names[name] += 1
                log.info(f"Renaming {name} -> {name}:{used_names[name]}")
                name += f":{used_names[name]}"
                used_names[name] = 1
            else:
                used_names[name] = 1
            if an in self.cn_configs and an != name:
                self.cn_configs[name] = self.cn_configs[an]
                changed_cn_configs[an] = name
                del self.cn_configs[an]
            new_names[key] = name
        self.alleles: Dict[str, MajorAllele] = dict()
        for key, minors in alleles_inverse.items():
            self.alleles[new_names[key]] = MajorAllele(
                name=new_names[key],
                cn_config=changed_cn_configs.get(key[0], key[0]),
                func_muts=key[1],
                minors={
                    sa: MinorAllele(
                        sa,
                        alleles[sa].alt_name,
                        set(alleles[sa].neutral_muts) - set(key[1]),
                        alleles[sa].activity,
                        alleles[sa].evidence,
                        alleles[sa].pharmvar,
                    )
                    for sa in minors
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
                new_name = "{}#{}".format(f, an)
                new_muts = preserved_mutations(f, a.func_muts)
                key = sorted_tuple(new_muts)
                if key in add:
                    add[key].minors.update(
                        {
                            san: MinorAllele(
                                san,
                                [],
                                preserved_mutations(f, sa.neutral_muts),
                                sa.activity,
                                sa.evidence,
                                sa.pharmvar,
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
                                san,
                                [],
                                preserved_mutations(f, sa.neutral_muts),
                                None,
                                None,
                                None,
                            )
                            for san, sa in a.minors.items()
                        },
                    )

            # Remove fusion (will be replaced at least by allele "1#{f}")
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
                        a.minors[min(sa)].alt_name,
                        set(nm),
                        a.minors[min(sa)].activity,
                        a.minors[min(sa)].evidence,
                        a.minors[min(sa)].pharmvar,
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

        self.exons = yml["reference"]["exons"]
        self.aminoacid = seq_to_amino("".join(self.seq[s:e] for [s, e] in self.exons))
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

    def get_functional(self, mut, infer=True) -> Optional[str]:
        """
        Returns:
            bool: ``True`` if a mutation is functional
            (i.e. does it affect the underlying aminoacid or not).
        """
        pos, op = mut
        if (pos, op) in self.mutation_info:
            return self.mutation_info[pos, op][0]

        # Calculate based on aminoacid change
        pos = self.chr_to_ref[pos]
        if infer and any(s <= pos < e for s, e in self.exons):
            if ">" not in op:
                return "indel"
            if self.strand < 0:
                op = self._reverse_op(op)
            seq = "".join(
                self.seq[s:pos] + op[2] + self.seq[pos + 1 : e]
                if s <= pos < e
                else self.seq[s:e]
                for s, e in self.exons
            )
            amino = seq_to_amino(seq)
            if amino != self.aminoacid:
                pos = next(i for i, a in enumerate(amino) if a != self.aminoacid[i])
                return f"{self.aminoacid[pos]}{pos}{amino[pos]}"
        return None

    def is_functional(self, mut, infer=True) -> bool:
        """
        Returns:
            bool: ``True`` if a mutation is functional
            (i.e. does it affect the underlying aminoacid or not).
        """
        return self.get_functional(mut, infer) is not None

    def _reverse_op(self, op: str) -> str:
        if ">" in op:
            l, r = op.split(">")
            return f"{rev_comp(l)}>{rev_comp(r)}"
        elif op[:3] == "ins":
            return f"ins{rev_comp(op[3:])}"
        elif op[:3] == "del":
            return f"del{rev_comp(op[3:])}"
        return op

    def get_dbsnp(self, *args) -> str:
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

    def get_refseq(self, *args, from_atg=False) -> str:
        """
        Returns:
            bool: ``True`` if a mutation is functional
            (i.e. does it affect the underlying aminoacid or not).
        """
        assert 1 <= len(args) <= 2
        atg_start = self.exons[0][0]
        if len(args) == 2:
            pos, op = args
        else:
            pos, op = args[0]
        if (pos, op) in self.mutation_info:
            pos, op = self.mutation_info[pos, op][2:4]
            if from_atg:
                if pos >= atg_start:
                    pos = pos - atg_start + 1
                else:
                    pos = pos - atg_start
            return f"{pos}{op}"
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
        return GRange(self.chr, mi, ma)

    def __contains__(self, i: int):
        return i in self.chr_to_ref

    def __getitem__(self, i: int):
        if i in self.chr_to_ref:
            s = self.seq[self.chr_to_ref[i]]
            return rev_comp(s) if self.strand < 0 else s
        return "N"

    def _print_mutation(self, m):
        fields = [
            self.get_dbsnp(m),
            self.get_refseq(m, from_atg=True),
            f"{self.refseq}:" + self.get_refseq(m),
        ]
        if m in self.mutation_info and self.mutation_info[m][0]:
            fields.append(self.mutation_info[m][0])
        if fields[0] == "-":
            fields = fields[1:]
        return f"{self.chr}:{m} ({', '.join(fields)})"

    def print_summary(self, query, full=True):
        minors = {m: a for a, al in self.alleles.items() for m in al.minors}
        alts = {
            mal.alt_name: (a, m)
            for a, al in self.alleles.items()
            for m, mal in al.minors.items()
            if mal.alt_name
        }
        if query:
            if query.lower() in map(str.lower, self.cn_configs):
                query = next(q for q in self.cn_configs if query.lower() == q.lower())
                self.print_cn(query)
            elif query.lower() in map(str.lower, self.alleles):
                query = next(q for q in self.alleles if query.lower() == q.lower())
                self.print_majors(query)
            elif query.lower() in map(str.lower, minors):
                query = next(q for q in minors if query.lower() == q.lower())
                self.print_minors(minors[query], query)
            elif query.lower() in map(str.lower, alts):
                query = next(q for q in alts if query.lower() == q.lower())
                self.print_minors(*alts[query])
            else:
                log.error(f"Cannot parse query '{query}'")
                sys.exit(1)
            return

        name = [f"Gene {self.name}"]
        if self.ensembl:
            name.append(f"ENSEMBL: {self.ensembl}")
        if self.pharmvar:
            name.append(f"PharmVar ID: {self.pharmvar}")
        log.info("-" * 80)
        log.info("{} {}", name[0], "" if len(name) == 1 else f"({', '.join(name[1:])})")
        log.info("-" * 80)

        st, ed = self.ref_to_chr[0], self.ref_to_chr[len(self.seq) - 1]
        log.info(f"{self.genome} genome locus: {self.chr}:{st}-{ed} ({self.refseq})")
        log.info("Strand: " + "3' (reverse)" if self.strand < 0 else "5'")
        pseudo = ", ".join(f"{p} (ID {i + 1})" for i, p in enumerate(self.pseudogenes))
        log.info(f"Pseudogenes: {pseudo if pseudo else 'none'}")
        amino = "\n  ".join(textwrap.wrap(self.aminoacid, width=78))
        log.info(f"Aminoacid:\n  {amino}")

        log.info("Regions of interest:")
        s = [f"  {'Region:':>10}  {self.name:>25}"]
        if self.pseudogenes:
            s[-1] += f" {self.pseudogenes[0]:>25}"
        for r in self.regions[0]:
            lg = str(self.regions[0].get(r, "-"))
            lp = str(self.regions.get(1, {}).get(r, ""))
            s.append(f"  {r:>10}: {lg:>25} {lp:>25}")
        log.info("\n".join(s))

        log.info(f"Structural alleles (deletions, conservations and fusions):")
        for c, cn in self.cn_configs.items():
            log.info(f"  {'*'+c+':':>8} {cn.description}")

        log.info("Major star-alleles:")
        for a, allele in natsorted(self.alleles.items()):
            if "#" in a:  # avoid fusions here
                continue
            log.info(f"  *{a}:")
            m = ",\n                   ".join(
                self._print_mutation(m) for m in sorted(allele.func_muts)
            )
            log.info(f"    Key mutations: {m}")

            minors = ",\n                   ".join(
                f"*{a}" + (f" (*{mi.alt_name})" if mi.alt_name else "")
                for a, mi in natsorted(allele.minors.items())
            )
            log.info(f"    Minor alleles: {minors}")

    def print_cn(self, major, full=False):
        if major not in self.cn_configs:
            raise AldyException(f"Structural allele {major} not found in {self.name}")

        config = self.cn_configs[major]
        log.info(f"Gene {self.name}, structural allele {self.name}*{major}:")

        lg = " ".join(f"{i:<8}" for i in [self.name] + self.pseudogenes)
        log.info(f"  Structure: {lg}")
        for r in self.regions[0]:
            lg = " ".join(f"{config.cn[i][r]:<8}" for i, _ in enumerate(self.regions))
            log.info(f"      {r:>5}: {lg}")

        if major in self.alleles:
            self.print_majors(major)

    def print_majors(self, major, full=False):
        if major not in self.alleles:
            raise AldyException(f"Major star-allele {major} not found in {self.name}")

        allele = self.alleles[major]
        log.info(f"Gene {self.name}, major star-allele {self.name}*{allele.name}:")
        log.info(f"  Structure: {allele.cn_config}")

        activities = {m.activity for m in allele.minors.values() if m.activity}
        if len(activities) == 1:
            log.info(f"  Activity: {activities.pop()}")
        elif len(activities) > 1:
            log.info(
                f"  Activity: {' or '.join(activities)} (please check minor alleles)"
            )

        ms = ",\n                 ".join(
            self._print_mutation(m) for m in sorted(allele.func_muts)
        )
        log.info(f"  Key mutations: {ms}")

        log.info("  Minor star-alleles:")
        for a, minor in natsorted(allele.minors.items()):
            log.info(f"    *{a}:")
            m = ",\n                        ".join(
                self._print_mutation(m) for m in sorted(minor.neutral_muts)
            )
            if minor.alt_name:
                log.info(f"      Legacy name: *{minor.alt_name}")
            log.info(f"      Silent mutations: {m}")

    def print_minors(self, major, minor):
        allele = self.alleles[major]
        if minor not in allele.minors:
            raise AldyException(f"Minor star-allele {minor} not found in {self.name}")

        log.info(f"\nGene {self.name}, minor star-allele {self.name}*{minor}:")
        log.info(f"  Structure: {allele.cn_config}")
        log.info(f"  Major star-allele: {self.name}*{allele.name}")
        minor = allele.minors[minor]
        if minor.activity:
            log.info(f"  Activity: {minor.activity}")
        if minor.evidence:
            evidence = {"L": "limited", "D": "definitive"}.get(
                minor.evidence, minor.evidence
            )
            log.info(f"  Evidence: {evidence}")
        if minor.pharmvar:
            log.info(f"  PharmVar details: {minor.pharmvar}")
        if minor.alt_name:
            log.info(f"  Legacy name: *{minor.alt_name}")

        ms = ",\n                 ".join(
            self._print_mutation(m) for m in allele.func_muts
        )
        log.info(f"  Key mutations: {ms}")

        m = ",\n                    ".join(
            self._print_mutation(m) for m in sorted(minor.neutral_muts)
        )
        log.info(f"  Silent mutations: {m}")
