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
    def __init__(self, name, cn_region=None, data=None, **params):
        self.name = name
        self.cn_region = cn_region
        self.data = data

        params = {k: v for k, v in params.items() if v is not None}
        self.cn_solution = params.get("cn_solution")
        self.neutral_value = float(params.get("neutral_value", 0))
        self.threshold = float(params.get("threshold", 0.5))
        self.min_coverage = float(
            params.get("min_coverage", 5.0 if name == "illumina" else 2.0)
        )
        self.min_quality = float(params.get("min_quality", 10.0))
        self.min_mapq = float(params.get("min_mapq", 10.0))
        self.phase = bool(params.get("phase", True))

        # Copy-number parameters
        self.cn_max = int(params.get("cn_max", 20))
        self.cn_pce_penalty = float(params.get("cn_pce_penalty", 2.0))
        """Error penalty applied to the PCE region (1 for no penalty)."""

        self.cn_diff = float(params.get("cn_diff", 10.0))
        self.cn_fit = float(params.get("cn_fit", 1.0))
        self.cn_parsimony = float(params.get("cn_parsimony", 0.5))
        self.cn_fusion_left = float(params.get("cn_fusion_left", 0.5))
        self.cn_fusion_right = float(params.get("cn_fusion_right", 0.25))

        self.major_novel = float(params.get("major_novel", 21.0))  # MAX_CN + 1
        """Penalty for each novel functional mutation (0 for no penalty).
            Should be large enough to prevent novel mutations unless really necessary."""

        self.minor_miss = float(params.get("minor_miss", 1.5))
        """ Penalty for each missed minor mutation (0 for no penalty).
            Ideally larger than `ADD_PENALTY_FACTOR` as additions should be cheaper.
        """
        self.minor_add = float(params.get("minor_add", 1.0))
        """ Penalty for each novel minor mutation (0 for no penalty).
            Zero penalty always prefers mutation additions over coverage errors if the
            normalized SNP slack coverage is >= 50%.
            Penalty of 1.0 prefers additions if the SNP slack coverage is >= 75%.
        """
        self.minor_phase = float(params.get("minor_phase", 0.4))

        log.debug(
            f"[params] "
            f"neutral={self.neutral_value}; "
            f"threshold={self.threshold}; "
            f"min_coverage={self.min_coverage}; "
            f"min_quality={self.min_quality}; "
            f"min_mapq={self.min_mapq}; "
            f"phase={self.phase}"
        )

    @staticmethod
    def load(gene, profile, cn_region=None, **params):
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
        regions: Dict[Tuple[str, str, int], GRange] = dict(),
        cn_region: Optional[GRange] = None,
        genome: Optional[str] = "hg19",
    ) -> Dict[str, Dict[str, List[float]]]:
        """
        Load the profile information from a SAM/BAM file.

        Returns:
            list[str, str, int, float]: list of tuples
            ``(gene_name, chromosome, loci, coverage)``.

        Params:
            factor (float):
                Scaling factor. Default is 2.0 (for two copies).
            regions (list[:obj:`GRange`], optional):
                List of regions to be extracted.

        Notes:
            Profiles that were used in Aldy paper:

                1. PGRNseq-v1/v3: NA17642 was used for all genes
                    (n.b. PGXT147 with rescale 2.52444127771 was used for CYP2B6 beta).
                2. PGRNseq-v2: NA19789.bam was used for all genes.
                3. Illumina: by definition contains all ones (uniform coverage profile).
        """

        if not genome:
            genome = "hg19"
        if len(regions) == 0:
            import pkg_resources

            gene_regions = {}
            for g in sorted(pkg_resources.resource_listdir("aldy.resources", "genes")):
                if g[-4:] != ".yml":
                    continue
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
            with pysam.AlignmentFile(
                sam_path, reference_filename="../data/ref/hg19.fa"
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

        d = {}
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
        return d
