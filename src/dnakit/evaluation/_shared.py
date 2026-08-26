"""Shared bounded primitives for DNA evaluation implementations."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import TypeAlias, cast

from dnakit.alignment import AlignmentConfig, align_pairwise
from dnakit.core import (
    DNA,
    DNARecord,
    DNASequence,
    DNASet,
    ImplementationInfo,
    ImplementationLabel,
    Issue,
    IssueSeverity,
    Provenance,
    ReferenceInfo,
)
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.similarity import edit_distance, exact_similarity, kmer_similarity

from .config import EvaluationLimits, PairSimilarityMethod
from .results import EvaluationEntry, EvaluationReport, ReferenceLibrary

EvaluationInput: TypeAlias = DNA | DNASequence | DNARecord | DNASet | Iterable[DNARecord]


@dataclass(frozen=True, slots=True)
class InputItem:
    sequence: DNASequence
    subject_id: str
    record: DNARecord | None


def provenance(*, reference: ReferenceLibrary | None = None) -> Provenance:
    reference_info = None
    if reference is not None:
        reference_info = ReferenceInfo(
            reference.name,
            version=reference.version,
            date=reference.date,
            checksum=reference.digest,
            filters=reference.filters,
        )
    return Provenance(
        implementation=ImplementationInfo(label=ImplementationLabel.NATIVE),
        reference=reference_info,
    )


def issue(code: str, severity: IssueSeverity, message: str, **details: object) -> Issue:
    return Issue(code, severity, message, details=details)


def materialize_input(value: EvaluationInput, *, limits: EvaluationLimits) -> tuple[InputItem, ...]:
    if isinstance(value, DNA) and value.is_single:
        record = value.record
        items: tuple[InputItem, ...] = (InputItem(record.sequence, record.id, record),)
    elif isinstance(value, DNASequence):
        items = (InputItem(value, "sequence", None),)
    elif isinstance(value, DNARecord):
        items = (InputItem(value.sequence, value.id, value),)
    else:
        source = value.records if isinstance(value, (DNA, DNASet)) else value
        try:
            iterator = iter(source)
        except TypeError as exc:
            raise ConfigurationError(
                "Evaluation input must be a DNASequence, DNARecord, DNASet, or record iterable.",
                code="INVALID_EVALUATION_INPUT",
            ) from exc
        collected: list[InputItem] = []
        for index, record in enumerate(iterator):
            if not isinstance(record, DNARecord):
                raise ConfigurationError(
                    "Evaluation collections must contain DNARecord objects.",
                    code="INVALID_EVALUATION_RECORD",
                    context={"input_index": index, "type": type(record).__name__},
                )
            if index >= limits.max_records:
                raise ConfigurationError(
                    "Evaluation input exceeds max_records.",
                    code="EVALUATION_RECORD_LIMIT",
                    context={
                        "record_count_lower_bound": index + 1,
                        "max_records": limits.max_records,
                    },
                )
            collected.append(InputItem(record.sequence, record.id, record))
        items = tuple(collected)
    if len(items) > limits.max_records:
        raise ConfigurationError(
            "Evaluation input exceeds max_records.",
            code="EVALUATION_RECORD_LIMIT",
            context={"record_count": len(items), "max_records": limits.max_records},
        )
    total_symbols = sum(item.sequence.symbol_length for item in items)
    if total_symbols > limits.max_total_symbols:
        raise ConfigurationError(
            "Evaluation input exceeds max_total_symbols.",
            code="EVALUATION_SYMBOL_LIMIT",
            context={"total_symbols": total_symbols, "max_total_symbols": limits.max_total_symbols},
        )
    return items


def require_nonempty(items: tuple[InputItem, ...], operation: str) -> None:
    if not items:
        raise ConfigurationError(
            f"{operation} requires at least one record.", code="EMPTY_EVALUATION_DATASET"
        )


def pair_count(left_size: int, right_size: int | None = None) -> int:
    return left_size * (left_size - 1) // 2 if right_size is None else left_size * right_size


def enforce_pair_limit(count: int, limits: EvaluationLimits) -> None:
    if count > limits.max_pairwise_comparisons:
        raise ConfigurationError(
            "Evaluation workload exceeds max_pairwise_comparisons.",
            code="EVALUATION_PAIRWISE_LIMIT",
            context={"pairwise_comparisons": count, "limit": limits.max_pairwise_comparisons},
        )


def require_pairwise_compatible(sequence: DNASequence, *, role: str) -> None:
    if sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "Pairwise evaluation does not silently omit or bridge Gap objects.",
            code="EVALUATION_GAP_UNSUPPORTED",
            context={"role": role},
        )


def pair_similarity(
    left: DNARecord,
    right: DNARecord,
    *,
    method: PairSimilarityMethod,
    k: int,
    canonical: bool,
    max_alignment_cells: int,
) -> tuple[float, float | None]:
    require_pairwise_compatible(left.sequence, role="left")
    require_pairwise_compatible(right.sequence, role="right")
    if method == "exact":
        return exact_similarity(left, right).value, 1.0
    if method == "identity":
        result = align_pairwise(
            left,
            right,
            config=AlignmentConfig(mode="global", max_cells=max_alignment_cells),
        )
        return (1.0 if result.identity is None else result.identity), min(
            result.query_coverage, result.target_coverage
        )
    if method == "edit":
        distance_result = edit_distance(left, right, max_cells=max_alignment_cells)
        denominator = max(left.sequence.symbol_length, right.sequence.symbol_length)
        return (1.0 if denominator == 0 else 1.0 - distance_result.distance / denominator), None
    if left.sequence.symbol_length < k or right.sequence.symbol_length < k:
        return (1.0 if left.sequence.symbols == right.sequence.symbols else 0.0), None
    return kmer_similarity(left, right, k=k, canonical=canonical).value, None


def aggregate_numeric(entries: tuple[EvaluationEntry, ...], metric: str) -> FrozenDict:
    values = tuple(
        float(value)
        for entry in entries
        if (value := entry.metrics.get(metric)) is not None
        and not isinstance(value, bool)
        and isinstance(value, (int, float))
    )
    return FrozenDict(
        {
            "count": len(values),
            "mean": math.fsum(values) / len(values) if values else None,
            "minimum": min(values) if values else None,
            "maximum": max(values) if values else None,
        }
    )


def report(
    *,
    name: str,
    method: str,
    version: str,
    parameters: Mapping[str, object],
    metrics: Mapping[str, object],
    entries: Iterable[EvaluationEntry] = (),
    issues: Iterable[Issue] = (),
    reference: ReferenceLibrary | None = None,
) -> EvaluationReport:
    return EvaluationReport(
        name,
        method,
        version,
        FrozenDict(cast(Mapping[str, object], to_json_compatible(parameters))),
        FrozenDict(cast(Mapping[str, object], to_json_compatible(metrics))),
        tuple(entries),
        provenance(reference=reference),
        tuple(issues),
    )


def sequence_digest_payload(records: DNASet) -> dict[str, object]:
    return {
        "schema": "dnakit.reference-library-sequences.v1",
        "records": [
            {
                "id": record.id,
                "parts": [
                    part
                    if isinstance(part, str)
                    else {
                        "gap": {
                            "length": part.length,
                            "kind": part.kind.value,
                            "crossable": part.crossable,
                            "evidence": part.evidence,
                            "metadata": part.metadata,
                        }
                    }
                    for part in record.sequence.parts
                ],
                "alphabet": record.sequence.alphabet.value,
                "topology": record.sequence.topology.value,
                "strandedness": record.sequence.strandedness.value,
            }
            for record in records
        ],
    }


def digest_json(value: object) -> str:
    payload = json.dumps(
        to_json_compatible(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def record_for(item: InputItem) -> DNARecord:
    return item.record or DNARecord(item.sequence, item.subject_id)


def mean(values: Iterable[float]) -> float | None:
    materialized = tuple(values)
    return math.fsum(materialized) / len(materialized) if materialized else None


def as_float(mapping: Mapping[str, object], key: str) -> float:
    value = mapping[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError("Internal evaluation metric is not numeric.")
    return float(value)


__all__ = [
    "EvaluationInput",
    "InputItem",
    "aggregate_numeric",
    "as_float",
    "digest_json",
    "enforce_pair_limit",
    "issue",
    "materialize_input",
    "mean",
    "pair_count",
    "pair_similarity",
    "provenance",
    "record_for",
    "report",
    "require_nonempty",
    "sequence_digest_payload",
]
