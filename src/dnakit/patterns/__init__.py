"""Deterministic sequence-pattern annotation APIs."""

from dnakit.patterns.coding import scan_codon_sites, scan_orfs
from dnakit.patterns.crispr import (
    BUILTIN_PAM_RULES,
    PAM_CATALOG_VERSION,
    PAMRule,
    scan_pam_candidates,
)
from dnakit.patterns.motif import (
    PROMOTER_CATALOG_VERSION,
    PROMOTER_MOTIFS,
    PWM,
    scan_motif,
    scan_promoter_motifs,
    scan_pwm,
    scan_tf_pwm,
)
from dnakit.patterns.regions import find_cpg_islands, find_low_complexity_regions
from dnakit.patterns.repeats import (
    find_inverted_repeats,
    find_microsatellites,
    find_reverse_complement_palindromes,
    find_tandem_repeats,
)
from dnakit.patterns.restriction import (
    BUILTIN_RESTRICTION_ENZYMES,
    RESTRICTION_CATALOG_VERSION,
    RestrictionEnzyme,
    scan_restriction_sites,
)
from dnakit.patterns.results import (
    CodonSite,
    GuideCandidate,
    InvertedRepeatHit,
    LowComplexityResult,
    MotifHit,
    ORFHit,
    PalindromeHit,
    PatternResult,
    RegionHit,
    RestrictionSiteHit,
    TandemRepeatHit,
)

__all__ = [
    "BUILTIN_PAM_RULES",
    "BUILTIN_RESTRICTION_ENZYMES",
    "PAM_CATALOG_VERSION",
    "PROMOTER_CATALOG_VERSION",
    "PROMOTER_MOTIFS",
    "PWM",
    "RESTRICTION_CATALOG_VERSION",
    "CodonSite",
    "GuideCandidate",
    "InvertedRepeatHit",
    "LowComplexityResult",
    "MotifHit",
    "ORFHit",
    "PAMRule",
    "PalindromeHit",
    "PatternResult",
    "RegionHit",
    "RestrictionEnzyme",
    "RestrictionSiteHit",
    "TandemRepeatHit",
    "find_cpg_islands",
    "find_inverted_repeats",
    "find_low_complexity_regions",
    "find_microsatellites",
    "find_reverse_complement_palindromes",
    "find_tandem_repeats",
    "scan_codon_sites",
    "scan_motif",
    "scan_orfs",
    "scan_pam_candidates",
    "scan_promoter_motifs",
    "scan_pwm",
    "scan_restriction_sites",
    "scan_tf_pwm",
]
