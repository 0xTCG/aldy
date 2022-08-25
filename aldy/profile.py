# 786
# Aldy source: profile.py
#   This file is subject to the terms and conditions defined in
#   file 'LICENSE', which is part of this source code package.


from typing import Dict, List, Tuple, Optional
import pysam
import os
import os.path
import yaml

from collections import defaultdict
from natsort import natsorted
from .common import log, GRange, AldyException, script_path, chr_prefix
from .gene import Gene


class Profile:
    """Profile and model parameter information."""

    def __init__(self, name, cn_region=None, data=None, **kwargs):
        self.name = name
        """Name of the profile."""

        self.cn_region = cn_region
        """Location of the copy-number neutral region."""

        self.data = data
        """Profile coverage data."""

        self.gap = 0.0
        """
        Optimality gap. Non-zero values enable non-optimal solutions.
        Default: 0 (only optimal solutions)
        """

        self.cn_solution = None
        """
        User-specified copy-number configuration.
        Default: `None` (uses CN solver in :py:mod:`aldy.cn` for detection)
        """

        self.neutral_value = 0.0
        """
        Joint coverage of the copy-number neutral region.
        Default: 0 (typically specified in the profile's YAML file)
        """

        self.threshold = 0.5
        """
        Single-copy variant threshold.
        Its value indicate the fraction of total reads that contain the given variant
        in a single gene copy. For example, if two copies are given (maternal and
        paternal), and if the single copy coverage is 10, threshold of 0.5 will ensure
        that all variants with coverage less than 5 (i.e., 0.5 * 10) are filtered out.
        Default: 0.5
        """

        self.min_coverage = 5.0 if name == "illumina" else 2.0
        """
        Minimum coverage needed to call a variant.
        Default: 2 (5 for illumina/wgs)
        """

        self.min_quality = 10
        """
        Minimum base quality for a read base to be considered.
        Default: 10
        """

        self.min_mapq = 10
        """
        Minimum mapping quality for a read to be considered.
        Default: 10
        """

        self.phase = True
        """
        Set if the phasing model in :py:mod:`aldy.minor` is to be used.
        Default: `True`
        """

        self.sam_long_reads = False
        """
        Set if long reads should be split-mapped. Should be used when dealing
        with long PacBio or Nanopore reads
        Default: `False` (typically specified in the profile's YAML file)
        """

        self.sam_mappy_preset = "map-hifi"
        """
        Mappy preset to use for split-mapping.
        Default: map-hifi
        """

        self.cn_max = 20
        """
        Maximum possible copy number of a gene.
        Default: 20
        """

        self.cn_pce_penalty = 2.0
        """
        Error penalty applied to the PCE region during CN calling (1 for no penalty).
        Default: 2.0
        """

        self.cn_diff = 10.0
        """
        The first CN objective term (coverage fit) coefficient.
        Default: 10.0
        """

        self.cn_fit = 1.0
        """
        The second CN objective term (gene fit) coefficient.
        Default: 1.0
        """

        self.cn_parsimony = 0.5
        """
        The third CN objective term (max. parsimony) coefficient.
        Default: 0.5
        """

        self.cn_fusion_left = 0.5
        """
        Extra penalty for the left fusions.
        Default: 0.5.
        """

        self.cn_fusion_right = 0.25
        """
        Extra penalty for the right fusions.
        Default: 0.25.
        """

        self.major_novel = 21.0
        """
        Penalty for novel functional mutation (0 for no penalty).
        Should be large enough to avoid calling novel mutations unless really necessary.
        Default: 21.0 (i.e., `max_cn + 1`)
        """

        self.minor_miss = 1.5
        """
        Penalty for missed minor mutations (0 for no penalty).
        Ideally larger than `minor_add` as additions should be cheaper.
        Default: 1.5
        """

        self.minor_add = 1.0
        """
        Penalty for novel minor mutations (0 for no penalty).
        Zero penalty ensures that extra mutations are preferred over the coverage errors
        if the normalized variant slack coverage is >= 50%. Penalty of 1.0 prefers
        additions only if the variant slack coverage is >= 75%.
        Default: 1.0
        """

        self.minor_phase = 0.4
        """
        The minor star-allele calling model's phasing term coefficient.
        Default: 0.4
        """

        self.minor_phase_vars = 3_000
        """
        Number of variables to use during the phasing. Use lower number if the model
        takes too long to complete.
        Default: 3,000
        """

        self.male = False
        """
        Set if the sample is male (i.e., has two X chromosomes). Used for calling
        sex chromosome genes (e.g., G6PD) when the CN calling is disabled.
        Default: False
        """

        self.max_minor_solutions = 1
        """
        Maximum number of minor solutions to report for each major solution.
        Default: 1
        """

        self.update(kwargs)

    def update(self, kwargs):
        params = {}
        for n, v in kwargs.items():
            if v is not None and n in self.__dict__:
                if n == "cn_solution":
                    self.__dict__[n] = v
                else:
                    try:
                        typ = type(self.__dict__[n])
                        self.__dict__[n] = typ(v)
                    except (ValueError, TypeError):
                        raise AldyException(f"Invalid parameter {n}: {v}") from None
                params[n] = self.__dict__[n]
        ps = "; ".join(f"{n}={v}" for n, v in params.items())
        log.debug(f"[params] {ps}")
        return params

    @staticmethod
    def load(gene, profile, cn_region=None, **params):
        """
        Load the copy number profile and parameters from a profile file.

        :param gene: Gene instance.
        :param profile: A profile YAML or a SAM/BAM/CRAM file that contains profile
            data.
        :param cn_region: Copy-number neutral region.
        """
        prof = None
        is_yml = True
        if os.path.exists(profile) and os.path.isfile(profile):
            ext = os.path.splitext(profile)
            if ext[-1] in [".bam", ".sam", ".cram"]:
                assert gene and cn_region
                regions = {
                    (gene.name, r, gi): rng
                    for gi, gr in enumerate(gene.regions)
                    for r, rng in gr.items()
                }
                prof = Profile.get_sam_profile_data(
                    profile,
                    regions=regions,
                    genome=gene.genome,
                    cn_region=cn_region,
                )
                is_yml = False
            else:
                with open(profile) as f:
                    prof = yaml.safe_load(f)
        else:
            profile_path = script_path(
                "aldy.resources.profiles/{}.yml".format(profile.lower())
            )
            with open(profile_path) as f:
                prof = yaml.safe_load(f)
        if not prof:
            raise AldyException(f"Could not load profile from {profile}")
        if is_yml and profile != "illumina" and cn_region:
            raise AldyException("-n is only valid with illumina or BAM profile")
        if is_yml and profile == "illumina" and cn_region:
            prof["neutral"]["value"] = cn_region.end - cn_region.start
        if "neutral" not in prof or "value" not in prof["neutral"]:
            raise AldyException("Profile missing neutral region")
        if gene.name not in prof:
            raise AldyException(f"Profile missing {gene.name}")
        return Profile(
            profile,
            GRange(*prof["neutral"][gene.genome]),
            prof,
            neutral_value=prof["neutral"].get("value"),
            **dict(prof.get("options", {}), **params),
        )

    @staticmethod
    def get_sam_profile_data(
        sam_path: str,
        ref_path: Optional[str] = None,
        regions: Dict[Tuple[str, str, int], GRange] = dict(),
        cn_region: Optional[GRange] = None,
        genome: Optional[str] = "hg19",
        params: Dict = dict(),
    ) -> Dict[str, Dict[str, List[float]]]:
        """
        Load the profile information from a SAM/BAM/CRAM file.

        :param regions: List of regions to be extracted.
        :param cn_region: Copy-number neutral region.

        :return: list of tuples `(gene_name, chromosome, loci, coverage)`.


        .. note::
            Profile samples used in the original Aldy paper:

                1. PGRNseq-v1/v3: NA17642
                2. PGRNseq-v2: NA19789
                3. Illumina: by definition contains all ones (uniform coverage).
        """

        if not genome:
            genome = "hg19"
        if len(regions) == 0:
            import pkg_resources

            gene_regions = {}
            for g in sorted(pkg_resources.resource_listdir("aldy.resources", "genes")):
                if not g.endswith(".yml"):
                    continue
                log.debug("Loading {}...", g)
                gg = Gene(script_path(f"aldy.resources.genes/{g}"), genome=genome)
                for gi, gr in enumerate(gg.regions):
                    for r, rng in gr.items():
                        gene_regions[gg.name, r, gi] = rng
        else:
            gene_regions = regions

        default_cn_neutral_region = {
            "hg19": GRange("22", 42547463, 42548249),
            "hg38": GRange("22", 42151472, 42152258),
        }
        gene_regions["neutral", "value", 0] = (
            cn_region if cn_region else default_cn_neutral_region[genome]
        )

        chr_regions: Dict[str, Tuple[int, int]] = {}
        for c, s, e in gene_regions.values():
            if c not in chr_regions:
                chr_regions[c] = (s, e)
            else:
                chr_regions[c] = (min(s, chr_regions[c][0]), max(e, chr_regions[c][1]))

        cov: dict = defaultdict(lambda: defaultdict(int))
        for c, (s, e) in natsorted(chr_regions.items()):
            if sam_path == "<illumina>":
                continue
            with pysam.AlignmentFile(  # type: ignore
                sam_path, reference_filename=ref_path
            ) as sam:
                region = GRange(c, s, e).samtools(
                    pad_left=1000,
                    pad_right=1000,
                    prefix=chr_prefix(c, [x["SN"] for x in sam.header["SQ"]]),
                )
                log.info("Scanning {}...", region)
                try:
                    for read in sam.fetch(region=region):
                        start, s_start = read.reference_start, 0
                        if not read.cigartuples:
                            continue

                        for op, size in read.cigartuples:
                            if op == 2:
                                for i in range(size):
                                    cov[c][start + i] += 1
                                start += size
                            elif op == 1:
                                s_start += size
                            elif op == 4:
                                s_start += size
                            elif op in [0, 7, 8]:
                                for i in range(size):
                                    cov[c][start + i] += 1
                                start += size
                                s_start += size
                except ValueError:
                    log.warn("Cannot fetch {}", region)

        d: Dict = {}
        for (g, r, ri), (c, s, e) in gene_regions.items():
            if g not in d:
                d[g] = {}
            if r not in d[g]:
                d[g][r] = [0]
            if ri >= len(d[g][r]):
                d[g][r].append(0)
            if sam_path == "<illumina>":
                d[g][r][ri] = e - s
            else:
                d[g][r][ri] = sum(cov[c][i] for i in range(s, e))
        d["neutral"]["value"] = d["neutral"]["value"][0]
        d["neutral"][genome] = [*gene_regions["neutral", "value", 0]]
        if params:
            d["options"] = {}
            for k, v in Profile("")._parse_params(params).items():
                d["options"][k] = v
        return d
