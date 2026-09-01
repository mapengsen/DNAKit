"""Variant annotation and population-frequency APIs."""

from dnakit.annotation.variants import (
    annotate_rsid_vep,
    annotate_variant_vep,
    get_clinvar_significance,
    get_clinvar_variant,
    get_dbsnp_frequencies,
    get_dbsnp_variant,
    get_gnomad_population_frequencies,
    get_gnomad_variant,
    recode_variant,
    search_clinvar_variants,
    search_gnomad_variants,
)
from dnakit.core import ProviderResult

__all__ = [
    "ProviderResult",
    "annotate_rsid_vep",
    "annotate_variant_vep",
    "get_clinvar_significance",
    "get_clinvar_variant",
    "get_dbsnp_frequencies",
    "get_dbsnp_variant",
    "get_gnomad_population_frequencies",
    "get_gnomad_variant",
    "recode_variant",
    "search_clinvar_variants",
    "search_gnomad_variants",
]
