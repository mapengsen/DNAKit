"""Explicit variant-annotation functions backed by optional public providers."""

from __future__ import annotations

from collections.abc import Sequence

from dnakit.backends.scientific import _run_scientific_function
from dnakit.core import ProviderResult
from dnakit.exceptions import ConfigurationError

_GNOMAD_DATASETS = frozenset(
    {
        "exac",
        "gnomad_r2_1",
        "gnomad_r2_1_controls",
        "gnomad_r2_1_non_cancer",
        "gnomad_r2_1_non_neuro",
        "gnomad_r2_1_non_topmed",
        "gnomad_r3",
        "gnomad_r3_controls_and_biobanks",
        "gnomad_r3_non_cancer",
        "gnomad_r3_non_neuro",
        "gnomad_r3_non_topmed",
        "gnomad_r3_non_v2",
        "gnomad_r4",
        "gnomad_r4_non_ukb",
    }
)


def _text(value: object, name: str, *, maximum: int = 2_048) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            f"{name} must be non-empty text.",
            code="INVALID_ANNOTATION_QUERY",
        )
    resolved = value.strip()
    if len(resolved) > maximum or "\x00" in resolved:
        raise ConfigurationError(
            f"{name} exceeds its text limit or contains NUL.",
            code="INVALID_ANNOTATION_QUERY",
        )
    return resolved


def _optional_text(value: str | None, name: str) -> str | None:
    return None if value is None else _text(value, name)


def _species(value: str) -> str:
    return _text(value, "species", maximum=256)


def _dataset(value: str) -> str:
    resolved = _text(value, "dataset", maximum=128)
    if resolved not in _GNOMAD_DATASETS:
        raise ConfigurationError(
            "Unsupported gnomAD dataset.",
            code="INVALID_GNOMAD_DATASET",
            context={"dataset": resolved},
        )
    return resolved


def _variant_names(value: str | Sequence[str] | None) -> str | tuple[str, ...] | None:
    if value is None or isinstance(value, str):
        return _optional_text(value, "variant_name")
    if not 1 <= len(value) <= 20:
        raise ConfigurationError(
            "variant_name must contain between 1 and 20 values.",
            code="INVALID_ANNOTATION_QUERY",
        )
    resolved = tuple(_text(item, "variant_name") for item in value)
    return resolved


def annotate_variant_vep(
    hgvs_notation: str,
    *,
    species: str = "human",
) -> ProviderResult:
    """Annotate one HGVS expression with the Ensembl VEP consequence service."""

    notation = _text(hgvs_notation, "hgvs_notation")
    organism = _species(species)
    arguments = {"hgvs_notation": notation, "species": organism}
    return _run_scientific_function(
        "annotate_variant_vep",
        arguments,
        parameters=arguments,
    )


def annotate_rsid_vep(
    rsid: str,
    *,
    species: str = "human",
) -> ProviderResult:
    """Annotate one rsID with Ensembl VEP."""

    variant_id = _text(rsid, "rsid", maximum=256)
    organism = _species(species)
    arguments = {"variant_id": variant_id, "species": organism}
    return _run_scientific_function(
        "annotate_rsid_vep",
        arguments,
        parameters=arguments,
    )


def recode_variant(
    variant_id: str,
    *,
    species: str = "human",
) -> ProviderResult:
    """Convert a variant identifier into equivalent Ensembl notations."""

    identifier = _text(variant_id, "variant_id")
    organism = _species(species)
    arguments = {"variant_id": identifier, "species": organism}
    return _run_scientific_function("recode_variant", arguments, parameters=arguments)


def search_clinvar_variants(
    *,
    gene: str | None = None,
    condition: str | None = None,
    variant_id: str | None = None,
    variant_name: str | Sequence[str] | None = None,
    clinical_significance: str | None = None,
    limit: int = 20,
) -> ProviderResult:
    """Search ClinVar using explicit gene, condition, identifier, or HGVS fields."""

    if isinstance(limit, bool) or not isinstance(limit, int) or not 1 <= limit <= 100:
        raise ConfigurationError(
            "limit must be an integer in [1, 100].",
            code="INVALID_ANNOTATION_QUERY",
        )
    names = _variant_names(variant_name)
    arguments: dict[str, object] = {
        "gene": _optional_text(gene, "gene"),
        "condition": _optional_text(condition, "condition"),
        "variant_id": _optional_text(variant_id, "variant_id"),
        "variant_name": list(names) if isinstance(names, tuple) else names,
        "clinical_significance": _optional_text(
            clinical_significance,
            "clinical_significance",
        ),
        "max_results": limit,
    }
    arguments = {key: value for key, value in arguments.items() if value is not None}
    if not any(key in arguments for key in ("gene", "condition", "variant_id", "variant_name")):
        raise ConfigurationError(
            "Provide at least one ClinVar search field.",
            code="CLINVAR_QUERY_REQUIRED",
        )
    return _run_scientific_function(
        "search_clinvar_variants",
        arguments,
        parameters=arguments,
    )


def get_clinvar_variant(variant_id: str) -> ProviderResult:
    """Return summary details for one ClinVar variation identifier."""

    arguments = {"variant_id": _text(variant_id, "variant_id", maximum=256)}
    return _run_scientific_function(
        "get_clinvar_variant",
        arguments,
        parameters=arguments,
    )


def get_clinvar_significance(variant_id: str) -> ProviderResult:
    """Return ClinVar clinical-significance assertions for one variation ID."""

    arguments = {"variant_id": _text(variant_id, "variant_id", maximum=256)}
    return _run_scientific_function(
        "get_clinvar_significance",
        arguments,
        parameters=arguments,
    )


def get_dbsnp_variant(rsid: str) -> ProviderResult:
    """Return coordinates and alleles for one dbSNP rsID."""

    arguments = {"rsid": _text(rsid, "rsid", maximum=256)}
    return _run_scientific_function("get_dbsnp_variant", arguments, parameters=arguments)


def get_dbsnp_frequencies(rsid: str) -> ProviderResult:
    """Return population allele-frequency records for one dbSNP rsID."""

    arguments = {"rsid": _text(rsid, "rsid", maximum=256)}
    return _run_scientific_function(
        "get_dbsnp_frequencies",
        arguments,
        parameters=arguments,
    )


def search_gnomad_variants(
    query: str,
    *,
    dataset: str = "gnomad_r4",
) -> ProviderResult:
    """Search gnomAD variants by rsID or supported free-text identifier."""

    arguments = {"query": _text(query, "query"), "dataset": _dataset(dataset)}
    return _run_scientific_function(
        "search_gnomad_variants",
        arguments,
        parameters=arguments,
    )


def get_gnomad_variant(
    variant_id: str,
    *,
    dataset: str = "gnomad_r4",
) -> ProviderResult:
    """Return aggregate gnomAD metadata and frequencies for one variant."""

    arguments = {
        "variant_id": _text(variant_id, "variant_id"),
        "dataset": _dataset(dataset),
    }
    return _run_scientific_function("get_gnomad_variant", arguments, parameters=arguments)


def get_gnomad_population_frequencies(
    variant_id: str,
    *,
    dataset: str = "gnomad_r4",
) -> ProviderResult:
    """Return population-stratified gnomAD frequencies for one variant."""

    arguments = {
        "variant_id": _text(variant_id, "variant_id"),
        "dataset": _dataset(dataset),
    }
    return _run_scientific_function(
        "get_gnomad_population_frequencies",
        arguments,
        parameters=arguments,
    )


__all__ = [
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
