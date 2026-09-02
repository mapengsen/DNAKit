"""Tests for circular, IUPAC-aware, and approximate deduplication."""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from dnakit.core import DNAAlphabet, DNARecord, DNASequence, Gap, Topology
from dnakit.datasets import (
    ClusterConfig,
    IUPACDeduplicationConfig,
    deduplicate,
    deduplicate_approximate,
    deduplicate_iupac,
)
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError


def _record(record_id: str, symbols: str, *, circular: bool = False) -> DNARecord:
    alphabet = DNAAlphabet.STRICT if set(symbols) <= set("ACGT") else DNAAlphabet.IUPAC
    return DNARecord(
        DNASequence(
            symbols,
            alphabet=alphabet,
            topology=Topology.CIRCULAR if circular else Topology.LINEAR,
        ),
        record_id,
    )


def test_circular_deduplication_collapses_rotations_with_stable_representative() -> None:
    result = deduplicate(
        [
            _record("origin-0", "AACG", circular=True),
            _record("origin-2", "C GAA".replace(" ", ""), circular=True),
        ],
        equivalence="circular",
    )

    assert result.records.ids == ("origin-0",)
    assert result.groups[0].member_ids == ("origin-0", "origin-2")
    assert result.groups[0].rotation_offsets == (0, 2)
    assert result.groups[0].rotation_offset_definition is not None
    assert result.equivalence == "circular"


def test_circular_reverse_complement_and_domain_errors_are_explicit() -> None:
    result = deduplicate(
        [_record("forward", "AAGC", circular=True), _record("reverse", "GCTT", circular=True)],
        equivalence="circular_reverse_complement",
    )
    assert result.output_count == 1

    with pytest.raises(ConfigurationError) as topology_error:
        deduplicate([_record("linear", "AAGC")], equivalence="circular")
    assert topology_error.value.code == "CIRCULAR_DEDUPLICATION_TOPOLOGY_REQUIRED"
    gapped = DNARecord(DNASequence(["AA", Gap(1), "CG"], topology=Topology.CIRCULAR), "gap")
    with pytest.raises(ConfigurationError) as gap_error:
        deduplicate([gapped], equivalence="circular")
    assert gap_error.value.code == "CIRCULAR_DEDUPLICATION_GAP_UNSUPPORTED"


def test_iupac_dedup_uses_complete_link_to_avoid_incompatible_transitive_merge() -> None:
    result = deduplicate_iupac([_record("a", "A"), _record("n", "N"), _record("g", "G")])

    assert [group.member_ids for group in result.groups] == [("a", "n"), ("g",)]
    assert result.groups[0].relation == "compatible"
    assert result.groups[1].relation == "singleton"
    assert result.grouping_strategy == "stable-greedy-complete-link"
    assert result.pairwise_comparison_count == 3
    assert (
        result.identical_pair_count,
        result.compatible_pair_count,
        result.conflict_pair_count,
    ) == (
        0,
        2,
        1,
    )
    assert [relation.relation for relation in result.pair_relations] == [
        "compatible",
        "conflict",
        "compatible",
    ]


def test_iupac_dedup_rejects_gaps_and_pair_limit_before_comparison() -> None:
    with pytest.raises(UnsupportedGapOperationError):
        deduplicate_iupac([DNARecord(DNASequence(["A", Gap(1)]), "gap")])
    with pytest.raises(ConfigurationError) as limit_error:
        deduplicate_iupac(
            [_record("a", "A"), _record("b", "N"), _record("c", "G")],
            config=IUPACDeduplicationConfig(max_pairwise_comparisons=2),
        )
    assert limit_error.value.code == "ADVANCED_PAIRWISE_LIMIT"


@pytest.mark.parametrize("method", ["identity", "edit", "kmer", "fingerprint"])
def test_approximate_dedup_reuses_threshold_clustering_and_is_auditable(method: str) -> None:
    result = deduplicate_approximate(
        [_record("a", "AAAA"), _record("b", "AAAT"), _record("c", "CCCC")],
        config=ClusterConfig(method=method, threshold=0.5, k=2),  # type: ignore[arg-type]
    )

    assert result.labels[0] == result.labels[1]
    assert result.labels[2] != result.labels[0]
    assert result.representatives.ids == ("a", "c")
    assert result.pairwise_comparison_count == 3
    assert result.clustering_strategy == "exhaustive-threshold-graph-connected-components"


def test_approximate_dedup_consumes_only_limit_plus_one_records() -> None:
    consumed: list[int] = []

    def records() -> Iterator[DNARecord]:
        for index in range(20):
            consumed.append(index)
            yield _record(str(index), "AAAA")

    with pytest.raises(ConfigurationError) as error:
        deduplicate_approximate(records(), config=ClusterConfig(max_records=2))
    assert error.value.code == "ADVANCED_DATASET_SIZE_LIMIT"
    assert consumed == [0, 1, 2]
