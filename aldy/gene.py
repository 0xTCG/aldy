# 786
# Aldy source: gene.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Tuple, Dict, List, Optional, Set, NamedTuple
from dataclasses import dataclass, field
from enum import Enum

import os
import yaml
import collections
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


class Mutation(NamedTuple):
    """
    Mutation description. Immutable.

    :param pos: Reference genome position (0-based).
    :param op: Variant operation in HGVS-like format:
    - ``X>Y``:  a SNP from X to Y
    - ``insX``:  an insertion of X
    - ``delX``:  a deletion of X (X, Y are both of format ``[ACGT]+``).

    .. note:: Has custom printer (``__str__``).
    """

    pos: int
    op: str

    def __str__(self):
        return f"{self.pos + 1}.{self.op}"


@dataclass
class MinorAllele:
    """
    Minor allele description.

    :param name: Minor allale name.
    :param alt_name: List of alternative names.
    :param neutral_muts: Netural mutations that describe the minor allele.
    :param activity: Activity indicator (see PharmVar for details).
    :param evidence: Evidence indicator (see PharmVar for details).
    :param pharmvar: PharmVar ID.

    .. note:: Has custom printer (``__str__``).
    """

    name: str
    alt_name: Optional[str] = None
    neutral_muts: Set[Mutation] = field(default_factory=set)
    activity: Optional[str] = None
    evidence: Optional[str] = None
    pharmvar: Optional[str] = None

    def __str__(self):
        muts = ", ".join(map(str, self.neutral_muts))
        return f"Minor({self.name}; [{muts}])"


@dataclass
class MajorAllele:
    """
    Major allele description.

    :param name: Major allele name.
    :param cn_config: Copy-number configuration.
    :param func_muts: Functional mutations that describe the major allele.
    :param minors: Minor alleles (suballeles) that are derived from this major allele.
    """

    name: str
    cn_config: str = "1"
    func_muts: Set[Mutation] = field(default_factory=set)
    minors: Dict[str, MinorAllele] = field(default_factory=dict)

    def get_minor_mutations(self, minor: str):
        """Sequence of mutations that describe a minor allele."""
        for m in self.func_muts:
            yield m
        for m in self.minors[minor].neutral_muts:
            yield m


class CNConfigType(Enum):
    """Enumeration describing the type of a copy-number configuration"""

    DEFAULT = 0
    LEFT_FUSION = 1
    RIGHT_FUSION = 2
    DELETION = 3


@dataclass
class CNConfig:
    """
    Copy-number (CN) configuration description. Immutable.

    :param cn:
        Value of the expected region copy number in each gene.
        For example, ``cn[0]["e1"] == 1`` means that the exon 1 of the main gene (ID 0)
        has one copy (and thus should be present) in the configuration `cn`.
    :param kind: Type of the copy-number configuration.
    :param alleles: Allele IDs that have this CN configuration.
    :param description: Human-readable description (e.g. "deletion").

    .. note:: Has custom printer (``__str__``)
    """

    cn: List[Dict[str, int]]
    kind: CNConfigType
    alleles: Set[str]
    description: str = ""

    def __str__(self):
        cov = "|".join("".join(str(g[r]) for r in self.cn[0]) for g in self.cn)
        alleles = " ".join(natsorted(self.alleles))
        return f"CNConfig({str(self.kind)[13:]}; vector={cov}; alleles=[{alleles}])"


@dataclass
class Gene:
    r"""
    Gene (and associated pseudogenes) description relative to a reference genome.

    :param name: Gene name (e.g. _CYP2D6_).
    :param version: Database version.
    :param ensembl: ENSEMBL gene ID.
    :param pharmvar: PharmVar ID.
    :param refseq: RefSeq sequence ID.
    :param seq: RefSeq sequence (typically \*1 allele).
    :param genome: Reference genome version (e.g. hg19).
    :param chr: Reference genome chromosome.
    :param strand: Reference genome strand of RefSeq sequence.
    :param chr_to_ref: Location mapping from the reference to the RefSeq sequence.
    :param ref_to_chr: Location mapping from the RefSeq sequence to the reference.
    :param pseudogenes: Pseudogene names. A pseudogene has ID greater than zero.
    :param regions:
        Collection of regions (names and ranges) for each gene in the database.
        Maps a gene ID to a dictionary that maps a gene region name
        (e.g. "e9" for exon 9) to its region in the reference genome.
        Gene 0 is the main gene.
    :param exons: List of RefSeq exon coordinates.
    :param aminoacid: Aminoacid sequence of the main gene.

    :param mutations:
        Maps ``(position, mutation_type)`` to a corresponding :obj:`Mutation`.
    :param random_mutations:
        Set of (silent) mutations that can occur in any allele.
    :param mutation_info:
        Attributes of each mutation.
        These attributes currently consist:
        - function (None if the mutation is not functional),
        - rsID,
        - RefSeq position, and
        - RefSeq opcode.
    :param do_copy_number: Indicates if the gene is subject to structural variations.
    :param cn_configs:
        Copy-number configurations associated with the gene.
        `1` (akin to \*1) is the default CN config (no structural variations present).
    :param unique_regions: List of genic regions used for copy-number calling.
    :param alleles: Major star-alleles in the gene.
    :param common_tandems:
        List of common allele tandems. Used in diplotype assignment heuristics.
        For example, the fact that \*13 is always followed by \*1
        (encoded as ``('13', '1')``) will be used to group \*13 and \*1 together
        within the same haplotype (e.g. \*13+\*1).
    """

    # Basic information
    name: str
    version: str
    pharmvar: Optional[str]
    ensembl: Optional[str]

    # Structure
    refseq: str
    seq: str
    genome: str
    chr: str
    strand: int
    chr_to_ref: Dict[int, int]
    ref_to_chr: Dict[int, int]
    pseudogenes: List[str]
    regions: List[Dict[str, GRange]]
    exons: List[Tuple[int, int]]
    aminoacid: str

    # Mutations
    mutations: Dict[Tuple[int, str], Tuple[Optional[str], str, int, str]]
    random_mutations: Set[Tuple[Optional[str], str, int, str]]

    # Copy number configurations
    do_copy_number: bool
    cn_configs: Dict[str, CNConfig]
    unique_regions: List[str]

    # Alleles
    alleles: Dict[str, MajorAllele]
    common_tandems: List[Tuple]

    # _lookup_range / _lookup_seq / _region_at

    def __init__(
        self,
        path: Optional[str],
        name: Optional[str] = None,
        yml: Optional[str] = None,
        genome: Optional[str] = None,
    ) -> None:
        """
        Initialize the Gene class with the database description
        specified in YML file ``path``.

        :param path: Location of YML file
        :param name: Gene name
        :param yml: YML file contents (used if ``path`` is ``None``).
        :param genome: Reference genome version
        :raise :obj:`aldy.common.AldyException`:: if a path is not set.
        """

        if path and os.path.exists(path):
            with open(path) as file:
                yml = file.read()
            if not name:
                name = os.path.split(path)[-1].split(".")[0].upper()

        if not yml or not name:
            raise AldyException("Either a path or a name should be given")

        yml = yaml.safe_load(yml)
        if genome:
            self.genome = genome
        else:
            self.genome = list(yml["reference"]["mappings"].keys())[0]  # type: ignore
        self._init_basic(yml)
        self._init_regions(yml)
        self._init_alleles(yml)
        self._init_partials()

    def region_at(self, pos: int) -> Optional[Tuple[int, str]]:
        """
        :param pos: Position in the reference genome.
        :return: Gene ID and a region that covers the position.
        """
        return self._region_at.get(pos, None)

    def get_functional(self, mut, infer=True) -> Optional[str]:
        """
        :return: String describing the mutation effect if a mutation is functional;
                 otherwise None.
        """
        pos, op = mut
        if (pos, op) in self.mutations:
            return self.mutations[pos, op][0]

        def reverse_op(op: str) -> str:
            if ">" in op:
                l, r = op.split(">")
                return f"{rev_comp(l)}>{rev_comp(r)}"
            elif op[:3] == "ins":
                return f"ins{rev_comp(op[3:])}"
            elif op[:3] == "del":
                return f"del{rev_comp(op[3:])}"
            return op

        # Calculate based on aminoacid change
        pos = self.chr_to_ref[pos]
        if infer and any(s <= pos < e for s, e in self.exons):
            if ">" not in op:
                return "indel"
            if self.strand < 0:
                op = reverse_op(op)
            if op[2] == "N":
                return None
            seq = "".join(
                self.seq[s:pos] + op[2] + self.seq[pos + 1 : e]
                if s <= pos < e
                else self.seq[s:e]
                for s, e in self.exons
            )
            amino = seq_to_amino(seq)
            if amino != self.aminoacid:
                pos = next(i for i, a in enumerate(amino) if a != self.aminoacid[i])
                return f"{self.aminoacid[pos]}{pos + 1}{amino[pos]}"
        return None

    def is_functional(self, mut, infer=True) -> bool:
        """
        :return: ``True`` if a mutation is functional.
        """
        return self.get_functional(mut, infer) is not None

    def get_rsid(self, *args, default=True) -> str:
        """
        :return: rsID if a mutation has it; otherwise "-".
            bool: ``True`` if a mutation is functional
            (i.e. does it affect the underlying aminoacid or not).
        """
        assert 1 <= len(args) <= 2
        if len(args) == 2:
            pos, op = args
        else:
            pos, op = args[0]
        res = "-"
        if (pos, op) in self.mutations:
            res = self.mutations[pos, op][1]
        return res if res != "-" or not default else f"{pos + 1}.{op}"

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
        if (pos, op) in self.mutations:
            pos, op = self.mutations[pos, op][2:4]
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
                if cn.kind == CNConfigType.DELETION
            )
        except StopIteration:
            return None

    def has_coverage(self, a: str, pos: int):
        """
        :return: ``True`` if a major allele `a` covers the mutation `m`.
        """
        m = self.region_at(pos)
        if m:
            return self.cn_configs[self.alleles[a].cn_config].cn[m[0]][m[1]] > 0
        return False

    def get_wide_region(self):
        mi = min(r.start for g in self.regions for r in g.values())
        ma = max(r.end for g in self.regions for r in g.values())
        return GRange(self.chr, mi, ma)

    def __contains__(self, i: int):
        return i in self.chr_to_ref

    def __getitem__(self, i):
        s, e = self._lookup_range
        if isinstance(i, slice):
            i, j = i.start, i.stop
            if j <= s or i >= e:
                return "N" * (j - i)
            loff = max(0, s - i)
            roff = max(0, j - e)
            if loff:
                i = s
            if roff:
                j = e
            return ("N" * loff) + self._lookup_seq[i - s : j - s] + ("N" * roff)
        else:
            assert isinstance(i, int)
            if not s <= i < e:
                return "N"
            return self._lookup_seq[i - s]

    def __str__(self):
        return f"Gene({self.name})"

    def __repr__(self):
        return f"Gene({self.name})"

    def _init_basic(self, yml) -> None:
        """
        Read basic gene properties (``name``, ``seq`` and ``region``).

        All database YAML coordinates are indexed starting from 1.
        Aldy internally uses 0-based indexing.
        """
        self.name = yml["name"]
        self.version = f"{yml.get('version', 'v?')} ({yml.get('generated', '-')})"
        self.refseq = yml["reference"]["name"]
        self.pharmvar = yml.get("pharmvar", None)
        self.ensembl = yml.get("ensembl", None)

        self.seq = yml["reference"]["seq"].replace("\n", "")
        if "patches" in yml["reference"]:
            seq = list(self.seq)
            for [pos, nuc] in yml["reference"]["patches"]:
                seq[pos - 1] = nuc
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

        self._lookup_range = (start - 1, end - 1)
        self._lookup_seq = "".join(
            (
                rev_comp(self.seq[self.chr_to_ref[i]])
                if self.strand < 0
                else self.seq[self.chr_to_ref[i]]
            )
            if i in self.chr_to_ref
            else "N"
            for i in range(*self._lookup_range)
        )

        self.exons = sorted((s - 1, e - 1) for [s, e] in yml["reference"]["exons"])
        self.aminoacid = seq_to_amino("".join(self.seq[s:e] for [s, e] in self.exons))

    def _init_regions(self, yml) -> None:
        """
        Calculate the genic regions and pseudogenes
        (``regions``, ``unique_regions`` and ``pseudogenes``).
        Prepare ``region_at()`` call.
        """

        self.regions = []
        # Each pseudogene is associated with an index > 0
        self.pseudogenes = []

        for i, g in enumerate(yml["structure"]["genes"]):
            if i > 0:
                self.pseudogenes.append(g)
            regions = {}
            num_exons = 0
            for name, coord in yml["structure"]["regions"][self.genome].items():
                if name[0] == "e" and name[1:].isdigit():  # this is exon
                    num_exons += 1
                if not coord[i * 2] <= coord[i * 2 + 1]:
                    raise AldyException(
                        f"Malformed YML file {self.name} (structure:regions:{name})"
                    )
                regions[name] = GRange(self.chr, coord[i * 2] - 1, coord[i * 2 + 1] - 1)
            for e in range(1, num_exons):  # fill introns
                r = [regions[f"e{e}"], regions[f"e{e + 1}"]][:: self.strand]
                regions[f"i{e}"] = GRange(self.chr, r[0][2], r[1][1])
            self.regions.append(
                dict(sorted(regions.items(), key=lambda x: x[1])[:: self.strand])
            )

        #: dict[int, (int, str)]:
        #: reverse lookup (gene, region) of gene regions given a location
        #: within the reference genome.
        for d in self.regions:
            regs = sorted((r.start, r.end) for r in d.values())
            for ri in range(1, len(regs)):
                if regs[ri - 1][1] < regs[ri][0]:
                    raise AldyException(
                        "Region {}-{} is not annotated", regs[ri - 1][1], regs[ri][0]
                    )
        self._region_at = {
            i: (g, r)
            for g, d in enumerate(self.regions)
            for r, rng in d.items()
            for i in range(rng.start, rng.end)
        }
        for i in self.chr_to_ref:
            if i not in self._region_at:
                raise AldyException(
                    f"Position {i} not within a named region of {self.name}"
                )
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

        def process_mutation(pos, op, info):
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
                if info == []:
                    info = ["-"]
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
                    log.warn(f"Ignoring {pos}.{op} in {name} (not in {self.refseq})")
                elif self.region_at(self.ref_to_chr[pos]) is None:
                    log.warn(f"Ignoring {pos}.{op} in {name} (not in named region)")
                else:
                    self.mutations.setdefault(
                        (self.ref_to_chr[pos], op),
                        (function, rsid, pos, orig_op),
                    )
                    yield Mutation(self.ref_to_chr[pos], op)

        self.mutations = {}
        self.random_mutations = set()
        for name, allele in yml["alleles"].items():
            if name == "random":
                for pos, op, *info in allele:
                    for m in process_mutation(pos, op, info):
                        self.random_mutations.add(m)
                continue
            name = allele_name(name)
            mutations: Set[Mutation] = set()
            if [self.name, "deletion"] in allele["mutations"]:
                deletion_allele = name
            else:
                for pos, op, *info in allele["mutations"]:
                    for m in process_mutation(pos, op, info):
                        mutations.add(m)

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
        if "tandems" in yml["structure"]:
            self.common_tandems = [tuple(y) for y in yml["structure"]["tandems"]]

        # Set copy number configurations
        # TODO: currently fusions only cover the space between the main gene and
        # the first pseudogene. Multi-pseudogene fusions are not supported.

        self.do_copy_number = bool(
            deletion_allele or len(fusions_left) or len(fusions_right)
        )
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
                    CNConfigType.LEFT_FUSION,
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
                CNConfigType.DELETION,
                {deletion_allele},
                f"{self.name} deletion",
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
                    CNConfigType.RIGHT_FUSION,
                    {a},
                    f"{self.pseudogenes[0]} conservation after {brk}",
                )
                inverse_cn[key] = a
            else:
                self.cn_configs[inverse_cn[key]].alleles.add(a)

        # Normal CN case
        used_alleles = {a for cn in self.cn_configs.values() for a in cn.alleles}
        has_pseudogenes = len(self.pseudogenes) > 0
        default_cn = [
            {r: 1 for r in self.regions[g]} for g in range(1 + has_pseudogenes)
        ]
        self.cn_configs = {min(v.alleles): v for v in self.cn_configs.values()}
        self.cn_configs["1"] = CNConfig(
            default_cn,
            CNConfigType.DEFAULT,
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
                log.debug(f"[gene] renaming {name} -> {name}:{used_names[name]}")
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
                func_muts=set(key[1]),
                minors={
                    sa: MinorAllele(
                        sa,
                        alleles[sa].alt_name,
                        set(alleles[sa].neutral_muts) - set(key[1]),
                        alleles[sa].activity,
                        alleles[sa].evidence,
                        alleles[sa].pharmvar,
                    )
                    for sa in natsorted(minors)
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
                m = self.region_at(m.pos)
                if m:
                    return self.cn_configs[f].cn[m[0]][m[1]] > 0
                return False

            return set(filter(filter_f, m))

        for f in filter(
            lambda x: self.cn_configs[x].kind == CNConfigType.LEFT_FUSION,
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
                new_name = f"{f}#{an}"
                new_muts = preserved_mutations(f, a.func_muts)
                key = sorted_tuple(new_muts)
                if key in add:
                    add[key].minors.update(
                        {
                            f"{f}#{san}": MinorAllele(
                                f"{f}#{san}",
                                None,
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
                            f"{f}#{san}": MinorAllele(
                                f"{f}#{san}",
                                neutral_muts=preserved_mutations(f, sa.neutral_muts),
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
