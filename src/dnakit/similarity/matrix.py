"""Bounded dense pairwise similarity and distance matrices."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from itertools import islice
from typing import TypeAlias

from dnakit.core import DNA, DNARecord, DNASequence
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.fingerprints import BitFingerprintResult, FingerprintResult
from dnakit.similarity._shared import validate_positive_int
from dnakit.similarity.compare import ComparisonInput, compare
from dnakit.similarity.results import DistanceResult, SimilarityMatrixResult, SimilarityResult
from dnakit.similarity.vector import KmerMode, WeightInput, ZeroVectorPolicy

DEFAULT_MAX_MATRIX_ITEMS = 1_000
_SUPPORTED_METHODS = frozenset(
    {
        "exact",
        "hamming",
        "edit",
        "kmer_jaccard",
        "kmer_containment",
        "kmer_cosine",
        "tanimoto",
        "jaccard",
        "cosine",
        "euclidean",
        "manhattan",
    }
)
_WEIGHTED_METHODS = frozenset({"tanimoto", "jaccard", "cosine", "euclidean", "manhattan"})
_ResolvedMatrixWeights: TypeAlias = tuple[float, ...] | dict[str, float]


def _label(value: ComparisonInput, index: int) -> str:
    if isinstance(value, DNA):
        return value.id
    if isinstance(value, DNARecord):
        return value.id
    if (
        isinstance(value, (FingerprintResult, BitFingerprintResult))
        and value.sequence_id is not None
    ):
        return value.sequence_id
    return f"item_{index}"


def _matrix_kind(method: str) -> str:
    return "distance" if method in ("hamming", "edit", "euclidean", "manhattan") else "similarity"


def _merge_resolved_weights(
    current: _ResolvedMatrixWeights | None,
    candidate: object,
) -> _ResolvedMatrixWeights:
    """Merge one scalar result's applied weights into matrix-level audit data."""

    if isinstance(candidate, Mapping):
        if current is not None and not isinstance(current, dict):
            raise ConfigurationError(
                "Matrix comparisons resolved incompatible weight forms.",
                code="INCONSISTENT_MATRIX_WEIGHTS",
            )
        merged = {} if current is None else dict(current)
        for feature, value in candidate.items():
            if (
                not isinstance(feature, str)
                or isinstance(value, bool)
                or not isinstance(value, (int, float))
            ):
                raise ConfigurationError(
                    "Matrix comparison returned invalid resolved weights.",
                    code="INVALID_RESOLVED_MATRIX_WEIGHTS",
                )
            number = float(value)
            if feature in merged and merged[feature] != number:
                raise ConfigurationError(
                    "Matrix comparisons resolved conflicting feature weights.",
                    code="INCONSISTENT_MATRIX_WEIGHTS",
                    context={"feature": feature},
                )
            merged[feature] = number
        return merged
    if isinstance(candidate, tuple) and all(
        not isinstance(value, bool) and isinstance(value, (int, float)) for value in candidate
    ):
        positional = tuple(float(value) for value in candidate)
        if current is not None and current != positional:
            raise ConfigurationError(
                "Matrix comparisons resolved conflicting positional weights.",
                code="INCONSISTENT_MATRIX_WEIGHTS",
            )
        return positional
    raise ConfigurationError(
        "Matrix comparison returned invalid resolved weights.",
        code="INVALID_RESOLVED_MATRIX_WEIGHTS",
    )


def _finalize_resolved_weights(
    weights: _ResolvedMatrixWeights | None,
) -> tuple[float, ...] | FrozenDict | None:
    if isinstance(weights, dict):
        return FrozenDict({feature: weights[feature] for feature in sorted(weights)})
    return weights


def similarity_matrix(
    items: Iterable[ComparisonInput],
    *,
    method: str,
    k: int | None = None,
    mode: KmerMode = "set",
    canonical: bool = False,
    overlapping: bool = True,
    reverse_complement: bool = False,
    weights: WeightInput | None = None,
    zero_vector_policy: ZeroVectorPolicy = "identity",
    max_items: int = DEFAULT_MAX_MATRIX_ITEMS,
) -> SimilarityMatrixResult:
    """Materialize a bounded full pairwise matrix in input order.

    The function intentionally refuses more than ``max_items`` inputs because
    the output itself is quadratic. K-mer containment is directional and
    therefore returns an asymmetric matrix; all other MVP methods are treated
    as symmetric.
    """

    validate_positive_int(max_items, "max_items")
    if method not in _SUPPORTED_METHODS:
        raise ConfigurationError(
            "Unknown pairwise matrix method.",
            code="INVALID_COMPARISON_METHOD",
            context={"method": method},
        )
    if weights is not None and method not in _WEIGHTED_METHODS:
        raise ConfigurationError(
            "weights are supported only for numeric fingerprint/vector matrix methods.",
            code="UNSUPPORTED_MATRIX_WEIGHTS",
            context={"method": method},
        )
    if isinstance(
        items,
        (str, bytes, DNASequence, DNARecord, FingerprintResult, BitFingerprintResult),
    ):
        raise ConfigurationError(
            "similarity_matrix items must be an iterable of comparison inputs.",
            code="INVALID_MATRIX_ITEMS",
        )
    try:
        iterator = iter(items)
    except TypeError as exc:
        raise ConfigurationError(
            "similarity_matrix items must be iterable.",
            code="INVALID_MATRIX_ITEMS",
        ) from exc
    values = tuple(islice(iterator, max_items + 1))
    if len(values) > max_items:
        raise ConfigurationError(
            "Pairwise matrix input exceeds max_items.",
            code="SIMILARITY_MATRIX_SIZE_LIMIT",
            context={
                "item_count": len(values),
                "item_count_is_lower_bound": True,
                "max_items": max_items,
            },
            hint="Reduce the input or use a future Top-k/indexed search API.",
        )
    if any(
        not isinstance(
            value,
            (
                DNASequence,
                DNARecord,
                FingerprintResult,
                BitFingerprintResult,
                tuple,
                list,
                dict,
            ),
        )
        for value in values
    ):
        raise ConfigurationError(
            "Matrix contains an unsupported comparison input.",
            code="INVALID_MATRIX_ITEM",
        )
    if weights is not None and not values:
        raise ConfigurationError(
            "A weighted matrix requires at least one item to resolve weight dimensions or keys.",
            code="EMPTY_WEIGHTED_MATRIX",
        )

    symmetric = method != "kmer_containment"
    size = len(values)
    rows = [[0.0] * size for _ in range(size)]
    resolved_matrix_weights: _ResolvedMatrixWeights | None = None
    for row in range(size):
        column_start = row if symmetric else 0
        for column in range(column_start, size):
            result = compare(
                values[row],
                values[column],
                method=method,
                k=k,
                mode=mode,
                canonical=canonical,
                overlapping=overlapping,
                reverse_complement=reverse_complement,
                weights=weights,
                zero_vector_policy=zero_vector_policy,
            )
            scalar = result.distance if isinstance(result, DistanceResult) else result.value
            if weights is not None and isinstance(result, SimilarityResult):
                candidate_weights = result.parameters.get("weights")
                if candidate_weights is not None:
                    resolved_matrix_weights = _merge_resolved_weights(
                        resolved_matrix_weights,
                        candidate_weights,
                    )
            rows[row][column] = scalar
            if symmetric:
                rows[column][row] = scalar
    kind = _matrix_kind(method)
    return SimilarityMatrixResult(
        name="similarity_matrix" if kind == "similarity" else "distance_matrix",
        method=method,
        value_kind=kind,  # type: ignore[arg-type]
        labels=tuple(_label(value, index) for index, value in enumerate(values)),
        values=tuple(tuple(row) for row in rows),
        symmetric=symmetric,
        max_items=max_items,
        parameters=FrozenDict(
            {
                "k": k,
                "mode": mode,
                "canonical": canonical,
                "overlapping": overlapping,
                "reverse_complement": reverse_complement,
                "weighted": weights is not None,
                "weights": _finalize_resolved_weights(resolved_matrix_weights),
                "zero_vector_policy": zero_vector_policy,
            }
        ),
    )


pairwise_matrix = similarity_matrix


__all__ = ["DEFAULT_MAX_MATRIX_ITEMS", "pairwise_matrix", "similarity_matrix"]
