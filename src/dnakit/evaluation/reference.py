"""Local versioned reference libraries, nearest references, novelty, and memorization."""

from __future__ import annotations

import heapq
from collections.abc import Mapping

from dnakit.core import DNA, DNASet
from dnakit.core._json import FrozenDict
from dnakit.exceptions import ConfigurationError

from ._shared import (
    EvaluationInput,
    digest_json,
    enforce_pair_limit,
    materialize_input,
    mean,
    pair_count,
    pair_similarity,
    record_for,
    report,
    sequence_digest_payload,
)
from .config import ReferenceSearchConfig
from .results import EvaluationEntry, EvaluationReport, ReferenceLibrary

_ReferenceHit = tuple[float, float | None, int, str, bool]
_RankedReferenceHit = tuple[float, int, _ReferenceHit]


def create_reference_library(
    records: DNA | DNASet,
    *,
    name: str,
    version: str,
    source: str,
    date: str | None = None,
    filters: Mapping[str, object] | None = None,
    index_parameters: Mapping[str, object] | None = None,
) -> ReferenceLibrary:
    """Bind caller-supplied DNA to version and content provenance.

    This function never performs a network request or downloads a database.
    The digest covers ordered IDs, sequence parts, alphabet, topology, and
    strandedness; metadata and annotations are deliberately out of scope.
    """

    if not isinstance(records, (DNA, DNASet)):
        raise ConfigurationError("records must be DNA or DNASet.", code="INVALID_REFERENCE_RECORDS")
    dataset = records.dataset if isinstance(records, DNA) else records
    digest_scope = "ordered IDs + sequence parts + alphabet + topology + strandedness"
    return ReferenceLibrary(
        records=dataset,
        name=name,
        version=version,
        source=source,
        digest=digest_json(sequence_digest_payload(dataset)),
        digest_scope=digest_scope,
        date=date,
        filters=FrozenDict(filters),
        index_parameters=FrozenDict(index_parameters),
    )


def _nearest_entries(
    queries: EvaluationInput,
    reference: ReferenceLibrary,
    config: ReferenceSearchConfig,
) -> tuple[EvaluationEntry, ...]:
    if not isinstance(reference, ReferenceLibrary):
        raise ConfigurationError(
            "reference must be ReferenceLibrary.", code="UNVERSIONED_REFERENCE"
        )
    query_items = materialize_input(queries, limits=config.limits)
    reference_items = materialize_input(reference.records, limits=config.limits)
    comparisons = pair_count(len(query_items), len(reference_items))
    enforce_pair_limit(comparisons, config.limits)
    entries: list[EvaluationEntry] = []
    for query_index, query_item in enumerate(query_items):
        query_record = record_for(query_item)
        ranked_hits: list[_RankedReferenceHit] = []
        for reference_index, reference_item in enumerate(reference_items):
            reference_record = record_for(reference_item)
            similarity, _ = pair_similarity(
                query_record,
                reference_record,
                method=config.method,
                k=config.k,
                canonical=config.canonical,
                max_alignment_cells=config.limits.max_alignment_cells,
            )
            exact = query_record.sequence.symbols == reference_record.sequence.symbols
            longest = max(
                query_record.sequence.symbol_length,
                reference_record.sequence.symbol_length,
            )
            coverage = (
                1.0
                if longest == 0
                else min(
                    query_record.sequence.symbol_length,
                    reference_record.sequence.symbol_length,
                )
                / longest
            )
            if similarity >= config.min_similarity and coverage >= config.min_coverage:
                hit: _ReferenceHit = (
                    similarity,
                    coverage,
                    reference_index,
                    reference_record.id,
                    exact,
                )
                ranked_hit: _RankedReferenceHit = (similarity, -reference_index, hit)
                if len(ranked_hits) < config.top_k:
                    heapq.heappush(ranked_hits, ranked_hit)
                elif (ranked_hit[0], ranked_hit[1]) > (
                    ranked_hits[0][0],
                    ranked_hits[0][1],
                ):
                    heapq.heapreplace(ranked_hits, ranked_hit)
        selected = sorted(
            (ranked_hit[2] for ranked_hit in ranked_hits),
            key=lambda hit: (-hit[0], hit[2]),
        )
        nearest = selected[0] if selected else None
        entries.append(
            EvaluationEntry(
                query_item.subject_id,
                query_index,
                FrozenDict(
                    {
                        "nearest_reference_id": nearest[3] if nearest else None,
                        "nearest_reference_index": nearest[2] if nearest else None,
                        "nearest_similarity": nearest[0] if nearest else None,
                        "nearest_coverage": nearest[1] if nearest else None,
                        "exact_match": nearest[4] if nearest else False,
                        "hits": tuple(
                            {
                                "reference_id": hit[3],
                                "reference_index": hit[2],
                                "similarity": hit[0],
                                "similarities": {
                                    resolved_method: hit[0],
                                    "exact": float(hit[4]),
                                },
                                "coverage": hit[1],
                                "exact": hit[4],
                            }
                            for hit in selected
                            for resolved_method in (config.method,)
                        ),
                    }
                ),
            )
        )
    return tuple(entries)


def nearest_reference(
    queries: EvaluationInput,
    reference: ReferenceLibrary,
    *,
    config: ReferenceSearchConfig | None = None,
) -> EvaluationReport:
    """Return deterministic top-k references from a bounded exhaustive local scan."""

    resolved = ReferenceSearchConfig() if config is None else config
    if not isinstance(resolved, ReferenceSearchConfig):
        raise TypeError("config must be ReferenceSearchConfig or None.")
    entries = _nearest_entries(queries, reference, resolved)
    found = tuple(
        float(value)
        for entry in entries
        if (value := entry.metrics["nearest_similarity"]) is not None
        and isinstance(value, (int, float))
        and not isinstance(value, bool)
    )
    return report(
        name="nearest_reference",
        method=f"bounded-exhaustive-{resolved.method}-scan",
        version="eval-nearest-reference-v1",
        parameters={
            "similarity_method": resolved.method,
            "k": resolved.k,
            "canonical": resolved.canonical,
            "top_k": resolved.top_k,
            "min_similarity": resolved.min_similarity,
            "min_coverage": resolved.min_coverage,
            "coverage_definition": "min(query_length,reference_length)/max(lengths)",
            "short_k_policy": "literal equality when either sequence is shorter than k",
            "reported_similarities": "selected method plus exact-symbol equality",
            "topology_policy": "linearized at stored origin; no circular rotation",
            "tie_break": "similarity-desc,reference-input-index-asc",
            "limits": resolved.limits,
            "reference": reference.to_dict(),
        },
        metrics={
            "query_count": len(entries),
            "reference_count": len(reference.records),
            "matched_query_count": len(found),
            "mean_nearest_similarity": mean(found),
            "pairwise_comparison_count": len(entries) * len(reference.records),
        },
        entries=entries,
        reference=reference,
    )


def evaluate_novelty(
    queries: EvaluationInput,
    reference: ReferenceLibrary,
    *,
    config: ReferenceSearchConfig | None = None,
) -> EvaluationReport:
    """Score novelty as ``1 - nearest similarity`` relative to one versioned library."""

    resolved = ReferenceSearchConfig() if config is None else config
    if not isinstance(resolved, ReferenceSearchConfig):
        raise TypeError("config must be ReferenceSearchConfig or None.")
    _validate_copy_search_thresholds(resolved)
    nearest_entries = _nearest_entries(queries, reference, resolved)
    entries: list[EvaluationEntry] = []
    for entry in nearest_entries:
        raw_similarity = entry.metrics["nearest_similarity"]
        similarity = (
            float(raw_similarity)
            if isinstance(raw_similarity, (int, float)) and not isinstance(raw_similarity, bool)
            else None
        )
        novelty = 1.0 - similarity if similarity is not None else 1.0
        metrics = dict(entry.metrics)
        metrics.update(
            {
                "novelty": novelty,
                "is_novel": similarity is None or similarity < resolved.copy_threshold,
                "threshold": resolved.copy_threshold,
            }
        )
        entries.append(EvaluationEntry(entry.subject_id, entry.input_index, FrozenDict(metrics)))
    materialized = tuple(entries)
    novel_count = sum(entry.metrics["is_novel"] is True for entry in materialized)
    novelty_values = tuple(
        value for entry in materialized if isinstance((value := entry.metrics["novelty"]), float)
    )
    return report(
        name="novelty",
        method=f"one-minus-nearest-{resolved.method}-similarity",
        version="eval-novelty-v1",
        parameters={
            "copy_threshold": resolved.copy_threshold,
            "threshold_rule": "novel iff no hit or nearest_similarity < copy_threshold",
            "k": resolved.k,
            "canonical": resolved.canonical,
            "min_similarity": resolved.min_similarity,
            "min_coverage": resolved.min_coverage,
            "topology_policy": "linearized at stored origin; no circular rotation",
            "limits": resolved.limits,
            "reference": reference.to_dict(),
        },
        metrics={
            "score": mean(novelty_values),
            "query_count": len(materialized),
            "novel_count": novel_count,
            "novel_fraction": novel_count / len(materialized) if materialized else None,
            "mean_novelty": mean(novelty_values),
            "pairwise_comparison_count": len(materialized) * len(reference.records),
        },
        entries=materialized,
        reference=reference,
    )


def evaluate_memorization(
    generated: EvaluationInput,
    training_reference: ReferenceLibrary,
    *,
    config: ReferenceSearchConfig | None = None,
) -> EvaluationReport:
    """Flag exact or threshold-near copies of a versioned training reference."""

    resolved = ReferenceSearchConfig() if config is None else config
    if not isinstance(resolved, ReferenceSearchConfig):
        raise TypeError("config must be ReferenceSearchConfig or None.")
    _validate_copy_search_thresholds(resolved)
    nearest_entries = _nearest_entries(generated, training_reference, resolved)
    entries: list[EvaluationEntry] = []
    exact_count = 0
    near_count = 0
    for entry in nearest_entries:
        raw_similarity = entry.metrics["nearest_similarity"]
        similarity = (
            float(raw_similarity)
            if isinstance(raw_similarity, (int, float)) and not isinstance(raw_similarity, bool)
            else None
        )
        exact = entry.metrics["exact_match"] is True
        copied = exact or (similarity is not None and similarity >= resolved.copy_threshold)
        exact_count += exact
        near_count += copied and not exact
        metrics = dict(entry.metrics)
        metrics.update(
            {
                "memorized": copied,
                "memorization_class": "exact" if exact else ("near" if copied else "none"),
                "threshold": resolved.copy_threshold,
            }
        )
        entries.append(EvaluationEntry(entry.subject_id, entry.input_index, FrozenDict(metrics)))
    materialized = tuple(entries)
    copied_count = exact_count + near_count
    fraction = copied_count / len(materialized) if materialized else None
    return report(
        name="memorization",
        method=f"exact-or-nearest-{resolved.method}-threshold",
        version="eval-memorization-v1",
        parameters={
            "copy_threshold": resolved.copy_threshold,
            "definition": "exact symbol equality OR nearest similarity >= copy_threshold",
            "length_stratification": False,
            "task_model_used": False,
            "min_similarity": resolved.min_similarity,
            "min_coverage": resolved.min_coverage,
            "topology_policy": "linearized at stored origin; no circular rotation",
            "limits": resolved.limits,
            "reference": training_reference.to_dict(),
        },
        metrics={
            "score": fraction,
            "generated_count": len(materialized),
            "memorized_count": copied_count,
            "memorization_fraction": fraction,
            "exact_copy_count": exact_count,
            "near_copy_count": near_count,
            "pairwise_comparison_count": len(materialized) * len(training_reference.records),
        },
        entries=materialized,
        reference=training_reference,
    )


def evaluate_reference_similarity(
    queries: EvaluationInput,
    reference: ReferenceLibrary,
    *,
    config: ReferenceSearchConfig | None = None,
) -> EvaluationReport:
    """Summarize nearest-reference similarity for a query collection."""

    resolved = ReferenceSearchConfig() if config is None else config
    if not isinstance(resolved, ReferenceSearchConfig):
        raise TypeError("config must be ReferenceSearchConfig or None.")
    entries = _nearest_entries(queries, reference, resolved)
    values = tuple(
        float(value)
        for entry in entries
        if (value := entry.metrics["nearest_similarity"]) is not None
        and isinstance(value, (int, float))
        and not isinstance(value, bool)
    )
    return report(
        name="reference_similarity",
        method=f"nearest-{resolved.method}-similarity",
        version="eval-reference-similarity-v1",
        parameters={
            "aggregation": "arithmetic mean of per-query nearest similarities",
            "missing_hit_policy": "omitted from mean and counted separately",
            "k": resolved.k,
            "canonical": resolved.canonical,
            "min_similarity": resolved.min_similarity,
            "min_coverage": resolved.min_coverage,
            "coverage_definition": "min(query_length,reference_length)/max(lengths)",
            "topology_policy": "linearized at stored origin; no circular rotation",
            "limits": resolved.limits,
            "reference": reference.to_dict(),
        },
        metrics={
            "score": mean(values),
            "query_count": len(entries),
            "matched_query_count": len(values),
            "mean_nearest_similarity": mean(values),
            "minimum_nearest_similarity": min(values) if values else None,
            "maximum_nearest_similarity": max(values) if values else None,
            "pairwise_comparison_count": len(entries) * len(reference.records),
        },
        entries=entries,
        reference=reference,
    )


def _validate_copy_search_thresholds(config: ReferenceSearchConfig) -> None:
    if config.min_similarity > config.copy_threshold:
        raise ConfigurationError(
            "min_similarity cannot exceed copy_threshold for novelty or memorization.",
            code="REFERENCE_COPY_THRESHOLD_CONFLICT",
            context={
                "min_similarity": config.min_similarity,
                "copy_threshold": config.copy_threshold,
            },
        )


__all__ = [
    "create_reference_library",
    "evaluate_memorization",
    "evaluate_novelty",
    "evaluate_reference_similarity",
    "nearest_reference",
]
