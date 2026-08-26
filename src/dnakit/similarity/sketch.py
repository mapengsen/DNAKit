"""Similarity metrics for compatible DNAKit sketches."""

from __future__ import annotations

from typing import Literal

from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError
from dnakit.fingerprints import SketchResult

from .results import SketchSimilarityResult

SketchMetric = Literal["jaccard", "containment"]


def _validate_compatible(left: SketchResult, right: SketchResult) -> None:
    if not isinstance(left, SketchResult) or not isinstance(right, SketchResult):
        raise ConfigurationError(
            "Sketch comparison requires two SketchResult objects.",
            code="INVALID_SKETCH_INPUT",
        )
    fields = ("schema_version", "k", "canonical", "seed", "selection")
    mismatches = tuple(name for name in fields if getattr(left, name) != getattr(right, name))
    if left.selection == "bottom_k" and left.num_hashes != right.num_hashes:
        mismatches += ("num_hashes",)
    if left.selection == "scaled" and left.scaled != right.scaled:
        mismatches += ("scaled",)
    if mismatches:
        raise ConfigurationError(
            "Sketch schemas are incompatible.",
            code="SKETCH_SCHEMA_MISMATCH",
            context={"fields": mismatches},
        )


def sketch_similarity(
    left: SketchResult,
    right: SketchResult,
    *,
    metric: SketchMetric = "jaccard",
    min_shared_hashes: int = 0,
) -> SketchSimilarityResult:
    """Compare compatible sketches as sets, with explicit empty-set semantics."""

    _validate_compatible(left, right)
    if metric not in {"jaccard", "containment"}:
        raise ConfigurationError("Unknown sketch metric.", code="UNKNOWN_SKETCH_METRIC")
    if (
        isinstance(min_shared_hashes, bool)
        or not isinstance(min_shared_hashes, int)
        or min_shared_hashes < 0
    ):
        raise ConfigurationError(
            "min_shared_hashes must be a non-negative integer.",
            code="INVALID_MIN_SHARED_HASHES",
        )
    left_set = set(left.hashes)
    right_set = set(right.hashes)
    shared = len(left_set & right_set)
    union = len(left_set | right_set)
    if metric == "jaccard":
        value = shared / union if union else 1.0
        denominator = union
    else:
        denominator = len(left_set)
        value = shared / denominator if denominator else (1.0 if not right_set else 0.0)
    passed = shared >= min_shared_hashes
    return SketchSimilarityResult(
        name="sketch_similarity",
        method=f"set-{metric}",
        value=value,
        metric=metric,
        left_id=left.sequence_id,
        right_id=right.sequence_id,
        left_hash_count=len(left_set),
        right_hash_count=len(right_set),
        shared_hash_count=shared,
        union_hash_count=union,
        denominator=denominator,
        passed_min_shared_hashes=passed,
        parameters=FrozenDict(
            {
                "schema_version": left.schema_version,
                "k": left.k,
                "canonical": left.canonical,
                "seed": left.seed,
                "selection": left.selection,
                "num_hashes": left.num_hashes,
                "scaled": left.scaled,
                "min_shared_hashes": min_shared_hashes,
            }
        ),
    )


__all__ = ["SketchMetric", "sketch_similarity"]
