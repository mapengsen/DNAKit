"""Exact k-mer and numeric fingerprint comparisons."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Literal, TypeAlias, TypeGuard

from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.fingerprints import BitFingerprintResult, FingerprintResult
from dnakit.similarity._shared import (
    NamedVector,
    NumericVector,
    SequenceInput,
    sequence_and_id,
    validate_bool,
    validate_positive_int,
)
from dnakit.similarity.results import SimilarityResult

VectorMetric: TypeAlias = Literal["tanimoto", "jaccard", "cosine", "euclidean", "manhattan"]
KmerMetric: TypeAlias = Literal["jaccard", "containment", "cosine"]
KmerMode: TypeAlias = Literal["set", "count"]
ZeroVectorPolicy: TypeAlias = Literal["identity", "error"]
FingerprintVectorResult: TypeAlias = FingerprintResult | BitFingerprintResult
VectorInput: TypeAlias = FingerprintVectorResult | NumericVector | NamedVector
WeightInput: TypeAlias = NumericVector | NamedVector
KmerVectorInput: TypeAlias = FingerprintResult | NamedVector

_IUPAC_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)


def _validate_zero_policy(value: ZeroVectorPolicy) -> ZeroVectorPolicy:
    if value not in ("identity", "error"):
        raise ConfigurationError(
            "zero_vector_policy must be 'identity' or 'error'.",
            code="INVALID_ZERO_VECTOR_POLICY",
            context={"zero_vector_policy": value},
        )
    return value


def _finite_number(value: object, *, role: str, index: int | str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            "Vector values must be finite numbers.",
            code="INVALID_VECTOR_VALUE",
            context={"role": role, "index": index, "value": value},
        )
    return float(value)


def _is_fingerprint_result(value: object) -> TypeGuard[FingerprintVectorResult]:
    return isinstance(value, (FingerprintResult, BitFingerprintResult))


def _fingerprint_values(value: FingerprintVectorResult, *, role: str) -> tuple[float, ...]:
    return tuple(
        _finite_number(item, role=role, index=index)
        for index, item in enumerate(value.dense_values())
    )


def _mapping_values(
    left: Mapping[str, int | float],
    right: Mapping[str, int | float],
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[str, ...]]:
    if any(not isinstance(key, str) or not key for key in (*left.keys(), *right.keys())):
        raise ConfigurationError(
            "Named-vector feature names must be non-empty strings.",
            code="INVALID_VECTOR_FEATURE",
        )
    features = tuple(sorted(set(left) | set(right)))
    left_values = tuple(
        _finite_number(left.get(feature, 0), role="left", index=feature) for feature in features
    )
    right_values = tuple(
        _finite_number(right.get(feature, 0), role="right", index=feature) for feature in features
    )
    return left_values, right_values, features


def _sequence_values(value: object, *, role: str) -> tuple[float, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(
            "Numeric vectors must be non-text sequences of finite numbers.",
            code="INVALID_VECTOR_INPUT",
            context={"role": role, "input_type": type(value).__name__},
        )
    return tuple(_finite_number(item, role=role, index=index) for index, item in enumerate(value))


def _align_vectors(
    left: VectorInput,
    right: VectorInput,
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[str, ...] | None, str | None, str | None]:
    if _is_fingerprint_result(left) and _is_fingerprint_result(right):
        left_parameters = left.parameters if isinstance(left, BitFingerprintResult) else None
        right_parameters = right.parameters if isinstance(right, BitFingerprintResult) else None
        if (
            left.schema_version != right.schema_version
            or left.feature_names != right.feature_names
            or left_parameters != right_parameters
        ):
            raise ConfigurationError(
                "Fingerprint schemas must match exactly before comparison.",
                code="FINGERPRINT_SCHEMA_MISMATCH",
                context={
                    "left_schema": left.schema_version,
                    "right_schema": right.schema_version,
                },
            )
        return (
            _fingerprint_values(left, role="left"),
            _fingerprint_values(right, role="right"),
            left.feature_names,
            left.sequence_id,
            right.sequence_id,
        )
    if _is_fingerprint_result(left) or _is_fingerprint_result(right):
        raise ConfigurationError(
            "A versioned fingerprint cannot be mixed with an unversioned vector.",
            code="MIXED_VECTOR_TYPES",
        )
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        left_values, right_values, features = _mapping_values(left, right)
        return left_values, right_values, features, None, None
    if isinstance(left, Mapping) or isinstance(right, Mapping):
        raise ConfigurationError(
            "Named and positional vectors cannot be mixed.",
            code="MIXED_VECTOR_TYPES",
        )
    left_values = _sequence_values(left, role="left")
    right_values = _sequence_values(right, role="right")
    if len(left_values) != len(right_values):
        raise ConfigurationError(
            "Positional vectors must have equal dimensions.",
            code="VECTOR_DIMENSION_MISMATCH",
            context={"left_dimension": len(left_values), "right_dimension": len(right_values)},
        )
    return left_values, right_values, None, None, None


def _resolve_weights(
    weights: WeightInput | None,
    *,
    dimension: int,
    features: tuple[str, ...] | None,
) -> tuple[float, ...]:
    if weights is None:
        return (1.0,) * dimension
    if isinstance(weights, Mapping):
        if features is None:
            raise ConfigurationError(
                "Mapping weights require named vector features.",
                code="INVALID_VECTOR_WEIGHTS",
            )
        if any(not isinstance(feature, str) or not feature for feature in weights):
            raise ConfigurationError(
                "Weight feature names must be non-empty strings.",
                code="INVALID_VECTOR_WEIGHTS",
            )
        unknown = set(weights) - set(features)
        if unknown:
            raise ConfigurationError(
                "Weights contain unknown feature names.",
                code="UNKNOWN_WEIGHT_FEATURE",
                context={"features": sorted(unknown)},
            )
        resolved = tuple(
            _finite_number(weights.get(feature, 1), role="weight", index=feature)
            for feature in features
        )
    else:
        resolved = _sequence_values(weights, role="weight")
        if len(resolved) != dimension:
            raise ConfigurationError(
                "Weight dimension must match vector dimension.",
                code="WEIGHT_DIMENSION_MISMATCH",
                context={"weight_dimension": len(resolved), "dimension": dimension},
            )
    if any(weight < 0 for weight in resolved):
        raise ConfigurationError(
            "Feature weights must be non-negative.",
            code="NEGATIVE_VECTOR_WEIGHT",
        )
    return resolved


def _resolved_weight_parameter(
    weights: WeightInput | None,
    *,
    resolved: tuple[float, ...],
    features: tuple[str, ...] | None,
) -> tuple[float, ...] | FrozenDict | None:
    """Return the exact applied weights in a deterministic audit form."""

    if weights is None:
        return None
    if isinstance(weights, Mapping):
        if features is None:  # Guarded by ``_resolve_weights``.
            raise ConfigurationError(
                "Mapping weights require named vector features.",
                code="INVALID_VECTOR_WEIGHTS",
            )
        return FrozenDict(
            {feature: weight for feature, weight in zip(features, resolved, strict=True)}
        )
    return resolved


def _undefined_zero_similarity(
    *,
    left_zero: bool,
    right_zero: bool,
    policy: ZeroVectorPolicy,
    metric: str,
) -> float:
    if policy == "error":
        raise ConfigurationError(
            f"{metric} is undefined for an effective zero vector.",
            code="ZERO_VECTOR_UNDEFINED",
            context={"left_zero": left_zero, "right_zero": right_zero, "metric": metric},
        )
    return 1.0 if left_zero and right_zero else 0.0


def _compute_vector_metric(
    left: tuple[float, ...],
    right: tuple[float, ...],
    weights: tuple[float, ...],
    *,
    metric: VectorMetric,
    zero_policy: ZeroVectorPolicy,
) -> tuple[float, Literal["similarity", "distance"], FrozenDict]:
    if metric in ("tanimoto", "jaccard") and any(value < 0 for value in (*left, *right)):
        raise ConfigurationError(
            f"{metric} requires non-negative vector values.",
            code="NEGATIVE_VECTOR_NOT_ALLOWED",
            context={"metric": metric},
        )
    if metric == "jaccard":
        intersection = math.fsum(
            weight
            for left_value, right_value, weight in zip(left, right, weights, strict=True)
            if left_value > 0 and right_value > 0
        )
        union = math.fsum(
            weight
            for left_value, right_value, weight in zip(left, right, weights, strict=True)
            if left_value > 0 or right_value > 0
        )
        if union == 0:
            value = _undefined_zero_similarity(
                left_zero=True,
                right_zero=True,
                policy=zero_policy,
                metric=metric,
            )
        else:
            value = intersection / union
        return (
            value,
            "similarity",
            FrozenDict({"intersection_weight": intersection, "union_weight": union}),
        )

    dot = math.fsum(
        weight * left_value * right_value
        for left_value, right_value, weight in zip(left, right, weights, strict=True)
    )
    left_squared = math.fsum(
        weight * value * value for value, weight in zip(left, weights, strict=True)
    )
    right_squared = math.fsum(
        weight * value * value for value, weight in zip(right, weights, strict=True)
    )
    if metric == "tanimoto":
        denominator = math.fsum((left_squared, right_squared, -dot))
        if math.isclose(denominator, 0.0):
            value = _undefined_zero_similarity(
                left_zero=math.isclose(left_squared, 0.0),
                right_zero=math.isclose(right_squared, 0.0),
                policy=zero_policy,
                metric=metric,
            )
        else:
            value = dot / denominator
        return (
            value,
            "similarity",
            FrozenDict(
                {
                    "dot_product": dot,
                    "left_squared_norm": left_squared,
                    "right_squared_norm": right_squared,
                    "denominator": denominator,
                }
            ),
        )
    if metric == "cosine":
        denominator = math.sqrt(left_squared * right_squared)
        if math.isclose(denominator, 0.0):
            value = _undefined_zero_similarity(
                left_zero=math.isclose(left_squared, 0.0),
                right_zero=math.isclose(right_squared, 0.0),
                policy=zero_policy,
                metric=metric,
            )
        else:
            value = dot / denominator
            value = min(1.0, max(-1.0, value))
        return (
            value,
            "similarity",
            FrozenDict(
                {
                    "dot_product": dot,
                    "left_squared_norm": left_squared,
                    "right_squared_norm": right_squared,
                    "denominator": denominator,
                }
            ),
        )
    if metric == "euclidean":
        squared_distance = math.fsum(
            weight * (left_value - right_value) ** 2
            for left_value, right_value, weight in zip(left, right, weights, strict=True)
        )
        return (
            math.sqrt(squared_distance),
            "distance",
            FrozenDict({"squared_distance": squared_distance}),
        )
    if metric == "manhattan":
        distance = math.fsum(
            weight * abs(left_value - right_value)
            for left_value, right_value, weight in zip(left, right, weights, strict=True)
        )
        return distance, "distance", FrozenDict({"absolute_distance": distance})
    raise ConfigurationError(
        "Unknown fingerprint similarity metric.",
        code="INVALID_VECTOR_METRIC",
        context={"metric": metric},
    )


def fingerprint_similarity(
    left: VectorInput,
    right: VectorInput,
    *,
    metric: VectorMetric = "tanimoto",
    weights: WeightInput | None = None,
    zero_vector_policy: ZeroVectorPolicy = "identity",
) -> SimilarityResult:
    """Compare versioned fingerprints or homogeneous numeric vectors.

    ``jaccard`` compares non-zero feature presence. ``tanimoto`` uses the
    generalized non-negative dot-product formula. ``cosine`` allows signed
    vectors. Euclidean and Manhattan methods return distances. Under the
    default zero policy, two effective zero vectors have similarity 1 while
    exactly one zero vector has similarity 0. A ``FingerprintResult`` is
    compared as its already-materialized vector; its recorded upstream Gap
    traversal policy is not reinterpreted here.
    """

    if metric not in ("tanimoto", "jaccard", "cosine", "euclidean", "manhattan"):
        raise ConfigurationError(
            "Unknown fingerprint similarity metric.",
            code="INVALID_VECTOR_METRIC",
            context={"metric": metric},
        )
    zero_policy = _validate_zero_policy(zero_vector_policy)
    left_values, right_values, features, left_id, right_id = _align_vectors(left, right)
    resolved_weights = _resolve_weights(
        weights,
        dimension=len(left_values),
        features=features,
    )
    audited_weights = _resolved_weight_parameter(
        weights,
        resolved=resolved_weights,
        features=features,
    )
    value, value_kind, components = _compute_vector_metric(
        left_values,
        right_values,
        resolved_weights,
        metric=metric,
        zero_policy=zero_policy,
    )
    return SimilarityResult(
        name="fingerprint_similarity",
        method=metric,
        value=value,
        value_kind=value_kind,
        left_id=left_id,
        right_id=right_id,
        left_dimension=len(left_values),
        right_dimension=len(right_values),
        parameters=FrozenDict(
            {
                "weighted": weights is not None,
                "weights": audited_weights,
                "feature_schema": features,
            }
        ),
        components=components,
        zero_vector_policy=zero_policy if value_kind == "similarity" else None,
        iupac_matching=None,
    )


def _canonical_word(word: str) -> str:
    reverse_complement = word.translate(_IUPAC_COMPLEMENT)[::-1]
    return min(word, reverse_complement)


def _kmer_counts(
    text: str,
    *,
    k: int,
    overlapping: bool,
    canonical: bool,
) -> Counter[str]:
    step = 1 if overlapping else k
    words = (text[start : start + k] for start in range(0, len(text) - k + 1, step))
    return Counter(_canonical_word(word) if canonical else word for word in words)


def _kmer_components(
    left: Mapping[str, int | float],
    right: Mapping[str, int | float],
    *,
    mode: KmerMode,
    metric: KmerMetric,
    zero_policy: ZeroVectorPolicy,
) -> tuple[float, FrozenDict]:
    features = tuple(sorted(set(left) | set(right)))
    if mode == "set":
        left_values = {feature: float(left[feature] > 0) for feature in features}
        right_values = {feature: float(right[feature] > 0) for feature in features}
    else:
        left_values = {feature: float(left[feature]) for feature in features}
        right_values = {feature: float(right[feature]) for feature in features}
    shared = math.fsum(min(left_values[feature], right_values[feature]) for feature in features)
    left_total = math.fsum(left_values[feature] for feature in features)
    right_total = math.fsum(right_values[feature] for feature in features)
    union = math.fsum(max(left_values[feature], right_values[feature]) for feature in features)
    if metric == "jaccard":
        if union == 0:
            value = _undefined_zero_similarity(
                left_zero=True,
                right_zero=True,
                policy=zero_policy,
                metric="kmer_jaccard",
            )
        else:
            value = shared / union
        denominator = union
    elif metric == "containment":
        if left_total == 0:
            value = _undefined_zero_similarity(
                left_zero=True,
                right_zero=right_total == 0,
                policy=zero_policy,
                metric="kmer_containment",
            )
        else:
            value = shared / left_total
        denominator = left_total
    elif metric == "cosine":
        dot = math.fsum(left_values[feature] * right_values[feature] for feature in features)
        left_squared = math.fsum(left_values[feature] ** 2 for feature in features)
        right_squared = math.fsum(right_values[feature] ** 2 for feature in features)
        denominator = math.sqrt(left_squared * right_squared)
        if denominator == 0:
            value = _undefined_zero_similarity(
                left_zero=left_squared == 0,
                right_zero=right_squared == 0,
                policy=zero_policy,
                metric="kmer_cosine",
            )
        else:
            value = dot / denominator
    else:
        raise ConfigurationError(
            "Unknown k-mer similarity metric.",
            code="INVALID_KMER_METRIC",
            context={"metric": metric},
        )
    return value, FrozenDict(
        {
            "shared_weight": shared,
            "left_weight": left_total,
            "right_weight": right_total,
            "denominator": denominator,
        }
    )


def kmer_vector_similarity(
    left: KmerVectorInput,
    right: KmerVectorInput,
    *,
    metric: KmerMetric = "jaccard",
    mode: KmerMode = "set",
    zero_vector_policy: ZeroVectorPolicy = "identity",
) -> SimilarityResult:
    """Compare named exact k-mer vectors, including versioned fingerprints."""

    if metric not in ("jaccard", "containment", "cosine"):
        raise ConfigurationError(
            "Unknown k-mer similarity metric.",
            code="INVALID_KMER_METRIC",
            context={"metric": metric},
        )
    if mode not in ("set", "count"):
        raise ConfigurationError(
            "k-mer mode must be 'set' or 'count'.",
            code="INVALID_KMER_MODE",
            context={"mode": mode},
        )
    zero_policy = _validate_zero_policy(zero_vector_policy)
    left_values, right_values, features, left_id, right_id = _align_vectors(left, right)
    if features is None:
        raise ConfigurationError(
            "k-mer vector comparison requires named features.",
            code="UNNAMED_KMER_VECTOR",
        )
    if any(value < 0 for value in (*left_values, *right_values)):
        raise ConfigurationError(
            "k-mer vector values must be non-negative.",
            code="NEGATIVE_VECTOR_NOT_ALLOWED",
        )
    left_mapping = dict(zip(features, left_values, strict=True))
    right_mapping = dict(zip(features, right_values, strict=True))
    value, components = _kmer_components(
        left_mapping,
        right_mapping,
        mode=mode,
        metric=metric,
        zero_policy=zero_policy,
    )
    return SimilarityResult(
        name="kmer_vector_similarity",
        method=f"{mode}_{metric}",
        value=value,
        value_kind="similarity",
        left_id=left_id,
        right_id=right_id,
        left_dimension=len(features),
        right_dimension=len(features),
        parameters=FrozenDict(
            {
                "metric": metric,
                "mode": mode,
                "feature_schema": features,
                "source": "named_kmer_vector",
            }
        ),
        components=components,
        zero_vector_policy=zero_policy,
        iupac_matching=None,
    )


def kmer_similarity(
    left: SequenceInput,
    right: SequenceInput,
    *,
    k: int,
    metric: KmerMetric = "jaccard",
    mode: KmerMode = "set",
    canonical: bool = False,
    overlapping: bool = True,
    zero_vector_policy: ZeroVectorPolicy = "identity",
) -> SimilarityResult:
    """Compare exact literal-IUPAC k-mer sets or count multisets."""

    validate_positive_int(k, "k")
    validate_bool(canonical, "canonical")
    validate_bool(overlapping, "overlapping")
    if metric not in ("jaccard", "containment", "cosine"):
        raise ConfigurationError(
            "Unknown k-mer similarity metric.",
            code="INVALID_KMER_METRIC",
            context={"metric": metric},
        )
    if mode not in ("set", "count"):
        raise ConfigurationError(
            "k-mer mode must be 'set' or 'count'.",
            code="INVALID_KMER_MODE",
            context={"mode": mode},
        )
    zero_policy = _validate_zero_policy(zero_vector_policy)
    left_sequence, left_id = sequence_and_id(left, role="left")
    right_sequence, right_id = sequence_and_id(right, role="right")
    left_counts = _kmer_counts(
        left_sequence.symbols,
        k=k,
        overlapping=overlapping,
        canonical=canonical,
    )
    right_counts = _kmer_counts(
        right_sequence.symbols,
        k=k,
        overlapping=overlapping,
        canonical=canonical,
    )
    value, components = _kmer_components(
        left_counts,
        right_counts,
        mode=mode,
        metric=metric,
        zero_policy=zero_policy,
    )
    return SimilarityResult(
        name="kmer_similarity",
        method=f"{mode}_{metric}",
        value=value,
        value_kind="similarity",
        left_id=left_id,
        right_id=right_id,
        left_dimension=len(left_counts),
        right_dimension=len(right_counts),
        parameters=FrozenDict(
            {
                "k": k,
                "metric": metric,
                "mode": mode,
                "canonical": canonical,
                "overlapping": overlapping,
            }
        ),
        components=components,
        zero_vector_policy=zero_policy,
        iupac_matching="literal",
    )


__all__ = ["fingerprint_similarity", "kmer_similarity", "kmer_vector_similarity"]
