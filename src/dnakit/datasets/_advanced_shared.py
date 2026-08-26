"""Shared bounded pairwise primitives for advanced dataset algorithms."""

from __future__ import annotations

import math
from collections.abc import Iterable

from dnakit.alignment import AlignmentConfig, align_pairwise
from dnakit.core import DNARecord, Gap, Topology
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.fingerprints import kmer_fingerprint
from dnakit.similarity import edit_distance, fingerprint_similarity, kmer_similarity

from .config import ClusterMethod


class UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root

    def groups(self) -> tuple[tuple[int, ...], ...]:
        grouped: dict[int, list[int]] = {}
        for index in range(len(self.parent)):
            grouped.setdefault(self.find(index), []).append(index)
        return tuple(tuple(values) for values in grouped.values())


def validate_max_records(value: int, maximum: int) -> None:
    if value > maximum:
        raise ConfigurationError(
            "Advanced dataset input exceeds max_records.",
            code="ADVANCED_DATASET_SIZE_LIMIT",
            context={"record_count": value, "max_records": maximum},
        )


def materialize_limited(
    records: Iterable[DNARecord],
    *,
    max_records: int,
) -> tuple[DNARecord, ...]:
    try:
        iterator = iter(records)
    except TypeError as exc:
        raise ConfigurationError(
            "records must be iterable.", code="INVALID_DATASET_RECORDS"
        ) from exc
    values: list[DNARecord] = []
    for index, record in enumerate(iterator):
        if not isinstance(record, DNARecord):
            raise ConfigurationError(
                "Every dataset item must be a DNARecord.",
                code="INVALID_DATASET_RECORD",
                context={"input_index": index},
            )
        values.append(record)
        if len(values) > max_records:
            validate_max_records(len(values), max_records)
    return tuple(values)


def validate_pairwise_sequences(records: tuple[DNARecord, ...], *, operation: str) -> None:
    for index, record in enumerate(records):
        if record.sequence.topology is Topology.CIRCULAR:
            raise ConfigurationError(
                f"{operation} requires linear sequences.",
                code="CIRCULAR_ADVANCED_DATASET_UNSUPPORTED",
                context={"input_index": index, "record_id": record.id},
            )
        if any(isinstance(part, Gap) for part in record.sequence.parts):
            raise UnsupportedGapOperationError(
                f"{operation} does not silently omit or bridge Gap objects.",
                code="ADVANCED_DATASET_GAP_UNSUPPORTED",
                context={"input_index": index, "record_id": record.id},
            )


def pair_count(size: int) -> int:
    return size * (size - 1) // 2


def ensure_pair_limit(size: int, limit: int) -> int:
    count = pair_count(size)
    if count > limit:
        raise ConfigurationError(
            "Pairwise workload exceeds max_pairwise_comparisons.",
            code="ADVANCED_PAIRWISE_LIMIT",
            context={"record_count": size, "pairwise_comparisons": count, "limit": limit},
        )
    return count


def quality_score(record: DNARecord) -> float:
    values = record.letter_annotations.get("phred_quality")
    return float("-inf") if not values else math.fsum(values) / len(values)


def pair_similarity(
    left: DNARecord,
    right: DNARecord,
    *,
    method: ClusterMethod,
    k: int,
    canonical: bool,
    max_alignment_cells: int,
) -> float:
    if method in {"kmer", "fingerprint"} and (
        left.sequence.symbol_length < k or right.sequence.symbol_length < k
    ):
        raise ConfigurationError(
            "k exceeds at least one sequence length in a pairwise dataset comparison.",
            code="DATASET_K_EXCEEDS_SEQUENCE_LENGTH",
            context={
                "k": k,
                "left_length": left.sequence.symbol_length,
                "right_length": right.sequence.symbol_length,
            },
        )
    if method == "identity":
        alignment_result = align_pairwise(
            left,
            right,
            config=AlignmentConfig(mode="global", max_cells=max_alignment_cells),
        )
        return 1.0 if alignment_result.identity is None else alignment_result.identity
    if method == "edit":
        distance_result = edit_distance(
            left,
            right,
            max_cells=max_alignment_cells,
        )
        denominator = max(left.sequence.symbol_length, right.sequence.symbol_length)
        return 1.0 if denominator == 0 else 1.0 - distance_result.distance / denominator
    if method == "kmer":
        return kmer_similarity(
            left,
            right,
            k=k,
            metric="jaccard",
            canonical=canonical,
        ).value
    left_fingerprint = kmer_fingerprint(left, k=k, canonical=canonical, mode="binary")
    right_fingerprint = kmer_fingerprint(right, k=k, canonical=canonical, mode="binary")
    return fingerprint_similarity(left_fingerprint, right_fingerprint, metric="tanimoto").value


def similarity_matrix(
    records: tuple[DNARecord, ...],
    *,
    method: ClusterMethod,
    k: int,
    canonical: bool,
    max_alignment_cells: int,
) -> tuple[tuple[float, ...], ...]:
    rows = [
        [1.0 if left == right else 0.0 for right in range(len(records))]
        for left in range(len(records))
    ]
    for left in range(len(records)):
        for right in range(left + 1, len(records)):
            value = pair_similarity(
                records[left],
                records[right],
                method=method,
                k=k,
                canonical=canonical,
                max_alignment_cells=max_alignment_cells,
            )
            rows[left][right] = rows[right][left] = value
    return tuple(tuple(row) for row in rows)


__all__ = [
    "UnionFind",
    "ensure_pair_limit",
    "materialize_limited",
    "pair_count",
    "pair_similarity",
    "quality_score",
    "similarity_matrix",
    "validate_max_records",
    "validate_pairwise_sequences",
]
