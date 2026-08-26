"""Fixed-schema exact k-mer fingerprints."""

from __future__ import annotations

from itertools import product

from dnakit.core._json import FrozenDict
from dnakit.descriptors import canonical_kmer, kmer_statistics
from dnakit.exceptions import ConfigurationError
from dnakit.fingerprints._shared import (
    FingerprintAmbiguityPolicy,
    FingerprintRepresentation,
    KmerFingerprintMode,
    SequenceInput,
    coerce_enum,
    sequence_and_id,
    validate_bool,
    validate_positive_int,
)
from dnakit.fingerprints.results import FingerprintResult, FingerprintValues, Numeric

DEFAULT_MAX_DIMENSION = 1_000_000


def _coerce_mode(value: KmerFingerprintMode | str) -> KmerFingerprintMode:
    if value == "presence":
        return KmerFingerprintMode.BINARY
    return coerce_enum(value, KmerFingerprintMode, "k-mer fingerprint mode")


def _schema_dimension(k: int, canonical: bool, max_dimension: int) -> int:
    raw_limit = max_dimension * 2 if canonical else max_dimension
    raw_dimension = 1
    for _ in range(k):
        raw_dimension *= 4
        if raw_dimension > raw_limit:
            return max_dimension + 1
    if not canonical:
        return raw_dimension
    fixed_reverse_complements = 0 if k % 2 else 4 ** (k // 2)
    return (raw_dimension + fixed_reverse_complements) // 2


def _feature_names(k: int, canonical: bool) -> tuple[str, ...]:
    words = ("".join(symbols) for symbols in product("ACGT", repeat=k))
    if not canonical:
        return tuple(words)
    return tuple(sorted({canonical_kmer(word) for word in words}))


def _count_for(counts: FrozenDict, feature: str) -> int:
    value = counts.get(feature, 0)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(
            "The descriptor k-mer count result contained a non-integer value.",
            code="INVALID_INTERNAL_KMER_COUNT",
            context={"feature": feature, "value": value},
        )
    return value


def kmer(
    value: SequenceInput,
    *,
    k: int,
    canonical: bool = False,
    mode: KmerFingerprintMode | str = KmerFingerprintMode.COUNT,
    representation: FingerprintRepresentation | str = FingerprintRepresentation.DENSE,
    ambiguity_policy: FingerprintAmbiguityPolicy | str = FingerprintAmbiguityPolicy.ERROR,
    overlapping: bool = True,
    cross_gaps: bool = False,
    max_dimension: int = DEFAULT_MAX_DIMENSION,
) -> FingerprintResult:
    """Build a deterministic fixed-schema exact k-mer fingerprint.

    The feature order is lexicographic over all A/C/G/T words. Canonical mode
    collapses each word with its reverse complement and sorts the unique
    representatives. ``presence`` is accepted as an alias for ``binary``.
    """

    validate_positive_int(k, "k")
    validate_positive_int(max_dimension, "max_dimension")
    validate_bool(canonical, "canonical")
    validate_bool(overlapping, "overlapping")
    validate_bool(cross_gaps, "cross_gaps")
    resolved_mode = _coerce_mode(mode)
    resolved_representation = coerce_enum(
        representation,
        FingerprintRepresentation,
        "fingerprint representation",
    )
    resolved_ambiguity = coerce_enum(
        ambiguity_policy,
        FingerprintAmbiguityPolicy,
        "fingerprint ambiguity policy",
    )
    dimension = _schema_dimension(k, canonical, max_dimension)
    if dimension > max_dimension:
        raise ConfigurationError(
            "The requested exact k-mer schema exceeds max_dimension.",
            code="FINGERPRINT_DIMENSION_LIMIT",
            context={
                "k": k,
                "canonical": canonical,
                "dimension_lower_bound": dimension,
                "max_dimension": max_dimension,
            },
            hint="Reduce k, use canonical=True, or explicitly raise max_dimension.",
        )

    sequence, sequence_id = sequence_and_id(value)
    features = _feature_names(k, canonical)
    statistics = kmer_statistics(
        sequence,
        k,
        overlapping=overlapping,
        canonical=canonical,
        ambiguity_policy=resolved_ambiguity.value,
        cross_gaps=cross_gaps,
    )
    counts = statistics.counts
    count_values = tuple(_count_for(counts, feature) for feature in features)
    if resolved_mode is KmerFingerprintMode.COUNT:
        ordered_values: tuple[Numeric, ...] = count_values
    elif resolved_mode is KmerFingerprintMode.BINARY:
        ordered_values = tuple(int(count != 0) for count in count_values)
    else:
        denominator = statistics.denominator
        ordered_values = tuple(
            count / denominator if denominator else 0.0 for count in count_values
        )

    materialized: FingerprintValues
    if resolved_representation is FingerprintRepresentation.DENSE:
        materialized = ordered_values
    else:
        materialized = FrozenDict(
            {
                feature: numeric
                for feature, numeric in zip(features, ordered_values, strict=True)
                if numeric != 0
            }
        )
    return FingerprintResult(
        name="kmer",
        method="exact_fixed_schema_kmer",
        schema_version="dnakit.kmer.acgt.v1",
        sequence_id=sequence_id,
        symbol_length=sequence.symbol_length,
        gap_count=statistics.gap_count,
        unknown_gap_count=statistics.unknown_gap_count,
        k=k,
        canonical=canonical,
        mode=resolved_mode,
        representation=resolved_representation,
        ambiguity_policy=resolved_ambiguity,
        overlapping=overlapping,
        cross_gaps=cross_gaps,
        feature_names=features,
        values=materialized,
        observation_count=statistics.denominator,
        ignored_ambiguity_count=statistics.ignored_ambiguity_count,
        max_dimension=max_dimension,
    )


kmer_fingerprint = kmer


__all__ = ["DEFAULT_MAX_DIMENSION", "kmer", "kmer_fingerprint"]
