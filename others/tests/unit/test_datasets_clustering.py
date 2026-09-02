"""Tests for clustering, hierarchy, and representative selection."""

from __future__ import annotations

import json

import pytest

from dnakit.core import DNARecord, DNASequence
from dnakit.datasets import (
    ClusterConfig,
    HierarchicalClusteringConfig,
    cluster_sequences,
    hierarchical_cluster,
    select_representatives,
)
from dnakit.exceptions import ConfigurationError


def _record(record_id: str, symbols: str, quality: int | None = None) -> DNARecord:
    annotations = {} if quality is None else {"phred_quality": (quality,) * len(symbols)}
    return DNARecord(DNASequence(symbols), record_id, letter_annotations=annotations)


def test_identity_cluster_connected_components_and_medoid_are_deterministic() -> None:
    records = [_record("a", "AAAA"), _record("b", "AAAT"), _record("c", "AATT")]
    config = ClusterConfig(method="identity", threshold=0.7, representative_policy="medoid")

    first = cluster_sequences(records, config=config)
    second = cluster_sequences(records, config=config)

    assert first == second
    assert first.labels == (0, 0, 0)
    assert first.representatives.ids == ("b",)
    assert json.loads(json.dumps(first.to_dict()))["seed"] == 0
    assert not first.seed_used


def test_representative_policies_preserve_cluster_label_order() -> None:
    records = [
        _record("long-low", "AAAA", 1),
        _record("short-high", "AA", 30),
        _record("other", "CCCC", 10),
    ]
    labels = [7, 7, 2]

    shortest = select_representatives(records, labels, policy="shortest")
    longest = select_representatives(records, labels, policy="longest")
    quality = select_representatives(records, labels, policy="best_quality")

    assert shortest.representative_ids == ("other", "short-high")
    assert longest.representative_ids == ("other", "long-low")
    assert quality.representative_ids == ("other", "short-high")


@pytest.mark.parametrize("linkage", ["single", "complete", "average"])
def test_hierarchical_clustering_outputs_valid_stable_linkage(linkage: str) -> None:
    records = [_record("a", "AAAA"), _record("b", "AAAT"), _record("c", "CCCC")]
    result = hierarchical_cluster(
        records,
        config=HierarchicalClusteringConfig(method="identity", linkage=linkage),  # type: ignore[arg-type]
    )

    assert len(result.linkage) == 2
    assert result.linkage[0].left_node == 0
    assert result.linkage[0].right_node == 1
    assert result.linkage[-1].member_count == 3
    assert result.pairwise_comparison_count == 3
    assert result.linkage_distance_update_count == 1


def test_hierarchical_empty_singleton_and_resource_limit() -> None:
    assert hierarchical_cluster([]).linkage == ()
    singleton = hierarchical_cluster([_record("a", "A")])
    assert singleton.record_ids == ("a",)
    assert singleton.linkage == ()
    with pytest.raises(ConfigurationError) as error:
        hierarchical_cluster(
            [_record("a", "A"), _record("b", "C"), _record("c", "G")],
            config=HierarchicalClusteringConfig(max_pairwise_comparisons=2),
        )
    assert error.value.code == "ADVANCED_PAIRWISE_LIMIT"


def test_cluster_rejects_alignment_cell_limit() -> None:
    with pytest.raises(ConfigurationError) as error:
        cluster_sequences(
            [_record("a", "AAAA"), _record("b", "AAAA")],
            config=ClusterConfig(method="identity", max_alignment_cells=4),
        )
    assert error.value.code == "ALIGNMENT_CELL_LIMIT"


def test_representative_selection_validates_limits_even_without_medoid() -> None:
    with pytest.raises(ConfigurationError):
        select_representatives([], [], max_records=0)
