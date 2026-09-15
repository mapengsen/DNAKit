"""Reference-based gene annotation, variant annotation, and population frequencies."""

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.gene_models import (
    GeneAnnotationConfig,
    GeneAnnotationResult,
    GeneRecord,
    GeneSegment,
)
from dnakit.annotation.genes import annotate_genes
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
    "GeneAnnotationConfig",
    "GeneAnnotationResult",
    "GeneRecord",
    "GeneSegment",
    "ProviderResult",
    "annotate_genes",
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

_register_public_verbose(globals())
