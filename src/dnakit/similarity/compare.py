"""Unified scalar comparison dispatch for sequence and fingerprint inputs."""

from __future__ import annotations

from typing import TypeAlias, cast

from dnakit.core import DNA, DNARecord, DNASequence
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.similarity._shared import SequenceInput, sequence_and_id, validate_bool
from dnakit.similarity.distance import edit_distance, hamming_distance
from dnakit.similarity.results import DistanceResult, SimilarityResult
from dnakit.similarity.vector import (
    KmerMetric,
    KmerMode,
    KmerVectorInput,
    VectorInput,
    VectorMetric,
    WeightInput,
    ZeroVectorPolicy,
    fingerprint_similarity,
    kmer_similarity,
    kmer_vector_similarity,
)

ComparisonInput: TypeAlias = SequenceInput | VectorInput
ComparisonResult: TypeAlias = SimilarityResult | DistanceResult


def _as_sequence(value: ComparisonInput, *, role: str) -> SequenceInput:
    if not isinstance(value, (DNA, DNASequence, DNARecord)):
        raise ConfigurationError(
            f"{role} must be DNASequence or DNARecord for this method.",
            code="INVALID_COMPARISON_INPUT",
            context={"role": role, "input_type": type(value).__name__},
        )
    return value


def exact_similarity(
    left: SequenceInput,
    right: SequenceInput,
    *,
    reverse_complement: bool = False,
) -> SimilarityResult:
    """Return 1 for literal full-sequence equality, otherwise 0."""

    validate_bool(reverse_complement, "reverse_complement")
    left_sequence, left_id = sequence_and_id(left, role="left")
    right_sequence, right_id = sequence_and_id(right, role="right")
    forward_equal = left_sequence.symbols == right_sequence.symbols
    reverse_equal = (
        reverse_complement and left_sequence.reverse_complement().symbols == right_sequence.symbols
    )
    return SimilarityResult(
        name="exact_similarity",
        method="literal_full_sequence_equality",
        value=float(forward_equal or reverse_equal),
        value_kind="similarity",
        left_id=left_id,
        right_id=right_id,
        left_dimension=left_sequence.symbol_length,
        right_dimension=right_sequence.symbol_length,
        parameters=FrozenDict({"reverse_complement": reverse_complement}),
        components=FrozenDict(
            {"forward_equal": forward_equal, "reverse_complement_equal": reverse_equal}
        ),
        zero_vector_policy=None,
        iupac_matching="literal",
    )


def compare(
    left: ComparisonInput,
    right: ComparisonInput,
    *,
    method: str,
    k: int | None = None,
    mode: KmerMode = "set",
    canonical: bool = False,
    overlapping: bool = True,
    reverse_complement: bool = False,
    weights: WeightInput | None = None,
    zero_vector_policy: ZeroVectorPolicy = "identity",
    return_path: bool = False,
) -> ComparisonResult:
    """Dispatch one explicit native comparison method.

    Supported names are ``exact``, ``hamming``, ``edit``,
    ``kmer_jaccard``, ``kmer_containment``, ``kmer_cosine``, and the
    fingerprint/vector methods ``tanimoto``, ``jaccard``, ``cosine``,
    ``euclidean``, and ``manhattan``.
    """

    if not isinstance(method, str) or not method:
        raise ConfigurationError("Comparison method must be a non-empty string.")
    if method == "exact":
        return exact_similarity(
            _as_sequence(left, role="left"),
            _as_sequence(right, role="right"),
            reverse_complement=reverse_complement,
        )
    if method == "hamming":
        return hamming_distance(
            _as_sequence(left, role="left"),
            _as_sequence(right, role="right"),
        )
    if method == "edit":
        return edit_distance(
            _as_sequence(left, role="left"),
            _as_sequence(right, role="right"),
            return_path=return_path,
        )
    if method in ("kmer_jaccard", "kmer_containment", "kmer_cosine"):
        metric = cast(KmerMetric, method.removeprefix("kmer_"))
        left_is_sequence = isinstance(left, (DNA, DNASequence, DNARecord))
        right_is_sequence = isinstance(right, (DNA, DNASequence, DNARecord))
        if left_is_sequence != right_is_sequence:
            raise ConfigurationError(
                "k-mer comparison cannot mix a DNA sequence with a k-mer vector.",
                code="MIXED_COMPARISON_INPUTS",
            )
        if left_is_sequence and k is None:
            raise ConfigurationError(
                "A positive k is required for k-mer comparison.",
                code="MISSING_KMER_SIZE",
            )
        if left_is_sequence:
            assert k is not None
            return kmer_similarity(
                _as_sequence(left, role="left"),
                _as_sequence(right, role="right"),
                k=k,
                metric=metric,
                mode=mode,
                canonical=canonical,
                overlapping=overlapping,
                zero_vector_policy=zero_vector_policy,
            )
        return kmer_vector_similarity(
            cast(KmerVectorInput, left),
            cast(KmerVectorInput, right),
            metric=metric,
            mode=mode,
            zero_vector_policy=zero_vector_policy,
        )
    if method in ("tanimoto", "jaccard", "cosine", "euclidean", "manhattan"):
        return fingerprint_similarity(
            cast(VectorInput, left),
            cast(VectorInput, right),
            metric=cast(VectorMetric, method),
            weights=weights,
            zero_vector_policy=zero_vector_policy,
        )
    raise ConfigurationError(
        "Unknown comparison method.",
        code="INVALID_COMPARISON_METHOD",
        context={"method": method},
    )


__all__ = ["compare", "exact_similarity"]
