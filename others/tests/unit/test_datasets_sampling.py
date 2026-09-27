"""Exercise proportional sampling, geometric representatives, and lossless file export."""

from __future__ import annotations

import csv
import json
from collections import Counter
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from dnakit import DNARecord, DNASequence
from dnakit.datasets import DiversitySamplingConfig, diversity_sample
from dnakit.exceptions import ConfigurationError
from dnakit.representations import RepresentationConfig


class MatrixBackend:
    def __init__(self, matrix: Any) -> None:
        self.matrix = np.asarray(matrix, dtype=np.float32)
        self.calls = 0

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        assert len(sequences) == len(self.matrix)
        self.calls += 1
        return self.matrix


def record(index: int, label: object = 0) -> DNARecord:
    sequence = "".join("ACGT"[(index // 4**i) % 4] for i in range(8))
    return DNARecord(
        DNASequence(sequence),
        f"seq_{index}",
        description=f"original {index}",
        metadata={"label": label, "source": {"row": index}},
    )


def config(n_samples: int, **kwargs: Any) -> DiversitySamplingConfig:
    return DiversitySamplingConfig(
        n_samples=n_samples,
        representation=RepresentationConfig(model="dnabert2", show_progress=False),
        **kwargs,
    )


@pytest.mark.parametrize("diversity", ["high", "low"])
def test_840_to_400_preserves_proportions_and_original_records(diversity: str) -> None:
    rng = np.random.default_rng(123)
    labels = rng.permutation([0] * 173 + [1] * 667)
    records = [record(i, int(label)) for i, label in enumerate(labels)]
    backend = MatrixBackend(rng.normal(size=(840, 4)))
    settings = config(400, seed=17, n_init=1, diversity=diversity)

    result = diversity_sample(iter(records), config=settings, backend=backend, verbose=False)
    assert backend.calls == 1
    assert {a.label: (a.input_count, a.n_samples) for a in result.allocations} == {
        0: (173, 82),
        1: (667, 318),
    }
    assert len(result.records) == 400
    assert len(result.clusters) == (400 if diversity == "high" else 10)
    assert len(set(result.selected_indices)) == 400
    assert result.selected_indices == tuple(sorted(result.selected_indices))
    assert Counter(r.metadata["label"] for r in result.records) == {0: 82, 1: 318}
    for selected, index in zip(result.records, result.selected_indices, strict=True):
        assert selected is records[index]
    covered = []
    for cluster in result.clusters:
        covered.extend(cluster.member_indices)
        assert len({records[i].metadata["label"] for i in cluster.member_indices}) == 1
        assert cluster.representative_index in cluster.member_indices
        assert all(
            result.cluster_labels[i] == cluster.cluster_index for i in cluster.member_indices
        )
    assert sorted(covered) == list(range(840))

    again = diversity_sample(records, config=settings, backend=backend, verbose=False)
    assert again.selected_indices == result.selected_indices
    assert again.cluster_labels == result.cluster_labels
    payload = json.loads(json.dumps(result.to_dict()))
    assert payload["model_name"] == "dnabert2"
    assert payload["n_samples"] == 400
    assert payload["diversity"] == diversity
    assert payload["low_n_clusters"] == (None if diversity == "high" else 5)


@pytest.mark.parametrize("n_samples", [1, 3, 6])
@pytest.mark.parametrize("with_config", [False, True])
def test_direct_sample_count_controls_total_clusters(n_samples: int, with_config: bool) -> None:
    settings = config(2, seed=19) if with_config else None
    result = diversity_sample(
        [record(i, i % 2) for i in range(6)],
        n_samples=n_samples,
        config=settings,
        backend=MatrixBackend([[i, 1] for i in range(6)]),
        verbose=False,
    )
    assert result.n_samples == len(result.clusters) == len(result.records) == n_samples
    assert sum(a.n_samples for a in result.allocations) == n_samples
    assert result.seed == (19 if with_config else 0)
    if settings is not None:
        assert settings.n_samples == 2


@pytest.mark.parametrize("n_samples", [0, -1, True, 1.5, 3])
def test_invalid_direct_sample_count_fails_before_inference(n_samples: Any) -> None:
    backend = MatrixBackend([[1, 0], [0, 1]])
    with pytest.raises(ConfigurationError):
        diversity_sample(
            [record(0), record(1)], n_samples=n_samples, backend=backend, verbose=False
        )
    assert backend.calls == 0


def test_representatives_are_nearest_own_class_centroid_with_global_indices() -> None:
    records = [record(i, i % 2) for i in range(6)]
    result = diversity_sample(
        records,
        config=config(2, normalize=False),
        backend=MatrixBackend([[0, 0], [10, 0], [2, 0], [12, 0], [9, 0], [19, 0]]),
        verbose=False,
    )
    assert result.selected_indices == (2, 3)
    assert result.cluster_labels == (0, 1, 0, 1, 0, 1)
    assert result.clusters[0].member_indices == (0, 2, 4)
    assert result.centers[0] == pytest.approx((11 / 3, 0))
    assert result.centers[1] == pytest.approx((41 / 3, 0))


def test_centroid_distance_ties_select_first_input_record() -> None:
    result = diversity_sample(
        [record(5), record(1)],
        config=config(1, normalize=False),
        backend=MatrixBackend([[0, 0], [2, 0]]),
        verbose=False,
    )
    assert result.selected_indices == (0,)
    assert result.records.ids == ("seq_5",)


def test_largest_remainder_ties_and_zero_allocations_are_explicit() -> None:
    result = diversity_sample(
        [record(0, "b"), record(1, "a"), record(2, "c")],
        config=config(2),
        backend=MatrixBackend([[1, 0], [0, 1], [1, 1]]),
        verbose=False,
    )
    assert [(a.label, a.n_samples) for a in result.allocations] == [("b", 1), ("a", 1), ("c", 0)]
    assert result.selected_indices == (0, 1)
    assert result.cluster_labels == (0, 1, None)


def test_label_types_do_not_merge_and_sampling_all_records_works() -> None:
    result = diversity_sample(
        [record(i, label) for i, label in enumerate([1, "1", True, 1.0])],
        config=config(4),
        backend=MatrixBackend([[1, 1]] * 4),
        verbose=False,
    )
    assert len(result.allocations) == 4
    assert result.selected_indices == (0, 1, 2, 3)
    assert [type(a.label) for a in result.allocations] == [int, str, bool, float]


def test_unlabeled_sampling_and_custom_label_key() -> None:
    records = [DNARecord(DNASequence("ACGT"), "a", metadata={"target": "positive"})]
    for key, expected_label in ((None, None), ("target", "positive")):
        result = diversity_sample(
            records,
            config=config(1, label_key=key),
            backend=MatrixBackend([[1, 0]]),
            verbose=False,
        )
        assert result.records[0] is records[0]
        assert result.allocations[0].label == expected_label


@pytest.mark.parametrize("label", [None, [1, 2], {"value": 0}])
def test_invalid_labels_fail_before_model_inference(label: object) -> None:
    backend = MatrixBackend([[1, 0]])
    with pytest.raises(ConfigurationError) as error:
        diversity_sample([record(0, label)], config=config(1), backend=backend, verbose=False)
    assert error.value.code == "INVALID_DIVERSITY_SAMPLING_LABEL"
    assert backend.calls == 0


def test_missing_label_and_invalid_size_fail_before_model_inference() -> None:
    backend = MatrixBackend([[1, 0]])
    for records, settings, code in (
        ([DNARecord(DNASequence("ACGT"), "a")], config(1), "INVALID_DIVERSITY_SAMPLING_LABEL"),
        ([record(0)], config(2), "INVALID_DIVERSITY_SAMPLE_SIZE"),
        ([], config(1), "EMPTY_NEURAL_CLUSTER_INPUT"),
        ([record(0), record(1)], config(1, max_records=1), "NEURAL_CLUSTER_RECORD_LIMIT"),
        ([record(0)], object(), "INVALID_DIVERSITY_SAMPLING_CONFIG"),
    ):
        with pytest.raises(ConfigurationError) as error:
            diversity_sample(records, config=settings, backend=backend, verbose=False)
        assert error.value.code == code
    assert backend.calls == 0


@pytest.mark.parametrize(
    "options",
    [
        {"n_samples": 0},
        {"n_samples": True},
        {"n_samples": 1.5},
        {"label_key": ""},
        {"label_key": 5},
        {"seed": -1},
        {"seed": 2**32},
        {"seed": True},
        {"normalize": 1},
        {"n_init": 0},
        {"max_iter": 0},
        {"tolerance": float("nan")},
        {"max_records": 0},
        {"representation": {}},
        {"diversity": "bad"},
        {"diversity": None},
        {"low_n_clusters": 0},
        {"low_n_clusters": True},
        {"low_n_clusters": 1.5},
    ],
)
def test_invalid_sampling_configuration(options: dict[str, Any]) -> None:
    with pytest.raises(ConfigurationError):
        DiversitySamplingConfig(**options)


@pytest.mark.parametrize("normalize", [False, True])
def test_insufficient_distinct_embeddings_never_silently_underfill(normalize: bool) -> None:
    matrix = [[1, 0], [2 if normalize else 1, 0]]
    with pytest.raises(ConfigurationError) as error:
        diversity_sample(
            [record(0), record(1)],
            config=config(2, normalize=normalize),
            backend=MatrixBackend(matrix),
            verbose=False,
        )
    assert error.value.code == "INSUFFICIENT_DIVERSITY_EMBEDDINGS"
    assert error.value.context["n_samples"] == 2


def test_zero_vectors_require_disabling_normalization() -> None:
    with pytest.raises(ConfigurationError) as error:
        diversity_sample([record(0)], config=config(1), backend=MatrixBackend([[0, 0]]))
    assert error.value.code == "ZERO_REPRESENTATION_VECTOR"


def test_representation_limit_is_respected_before_model_inference() -> None:
    settings = config(1)
    settings = replace(settings, representation=replace(settings.representation, max_records=1))
    backend = MatrixBackend([[1, 0], [0, 1]])
    with pytest.raises(ConfigurationError):
        diversity_sample([record(0), record(1)], config=settings, backend=backend, verbose=False)
    assert backend.calls == 0


@pytest.mark.parametrize("diversity", ["high", "low"])
def test_csv_sampling_preserves_original_rows_and_writes_audit_report(
    tmp_path: Path, diversity: Any
) -> None:
    source, output = tmp_path / "candidates.csv", tmp_path / "selected"
    source.write_text(
        "id,sequence,label,score\na,AAAA,0,0.10\nb,AAAT,0,0.20\nc,CCCC,1,0.30\nd,CCCT,1,0.40\n",
        encoding="utf-8",
    )
    result = diversity_sample(
        source,
        n_samples=2,
        diversity=diversity,
        output_path=output,
        config=config(400, normalize=False),
        backend=MatrixBackend([[0, 0], [2, 0], [10, 0], [12, 0]]),
        show_progress=False,
        verbose=False,
    )
    with (output / "diversity_sample.csv").open() as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        assert reader.fieldnames == ["id", "sequence", "label", "score"]
    assert rows == [
        {"id": "a", "sequence": "AAAA", "label": "0", "score": "0.10"},
        {"id": "c", "sequence": "CCCC", "label": "1", "score": "0.30"},
    ]
    report = json.loads((output / "result.json").read_text())
    assert report["selected_indices"] == list(result.selected_indices) == [0, 2]
    assert [(a["label"], a["n_samples"]) for a in report["allocations"]] == [("0", 1), ("1", 1)]
    assert report["diversity"] == diversity


def test_low_diversity_selects_compact_region_instead_of_spreading_across_clusters() -> None:
    # Five well-separated regions; the first is much tighter than the others.
    positions = np.concatenate(
        [
            100 * group + np.arange(-6, 6) * spread
            for group, spread in enumerate([0.001, 0.1, 0.2, 0.3, 0.4])
        ]
    )
    matrix = np.column_stack((positions, np.ones(len(positions))))
    records = [record(i) for i in range(len(matrix))]
    settings = config(6, normalize=False, n_init=2)
    low = diversity_sample(
        records, diversity="low", config=settings, backend=MatrixBackend(matrix), verbose=False
    )
    high = diversity_sample(
        records, diversity="high", config=settings, backend=MatrixBackend(matrix), verbose=False
    )
    assert len(low.clusters) == 5
    assert len(low.records) == 6
    assert all(index < 12 for index in low.selected_indices)
    assert np.ptp(positions[list(low.selected_indices)]) < 0.01
    assert np.ptp(positions[list(high.selected_indices)]) > 300
    assert low.sampling_strategy == "minimum-dispersion-centroid-neighborhood"
    assert settings.diversity == "high"


def test_low_diversity_spans_cluster_boundaries_to_fill_large_quota() -> None:
    matrix = [[value, 1] for value in [0, 1, 10, 11, 20, 21, 30, 31, 40, 41]]
    result = diversity_sample(
        [record(i) for i in range(10)],
        diversity="low",
        config=config(7, normalize=False),
        backend=MatrixBackend(matrix),
        verbose=False,
    )
    assert len(result.clusters) == 5
    assert all(len(cluster.member_indices) < 7 for cluster in result.clusters)
    assert len(result.records) == len(set(result.selected_indices)) == 7
    assert all(result.cluster_labels[index] is not None for index in result.selected_indices)


@pytest.mark.parametrize("n_samples", [1, 3, 6])
def test_low_diversity_accepts_duplicate_embeddings_and_small_classes(n_samples: int) -> None:
    result = diversity_sample(
        [record(i) for i in range(6)],
        n_samples=n_samples,
        diversity="low",
        backend=MatrixBackend([[1, 0]] * 6),
        verbose=False,
    )
    assert len(result.clusters) == 1
    assert result.selected_indices == tuple(range(n_samples))
    assert result.low_n_clusters == 5  # Requested count; actual count is capped.


def test_low_diversity_arguments_override_configuration_and_nearest_ties_are_stable() -> None:
    settings = config(4, normalize=False, diversity="high", low_n_clusters=5)
    result = diversity_sample(
        [record(i) for i in [9, 8, 7, 6]],
        n_samples=2,
        diversity="low",
        low_n_clusters=1,
        config=settings,
        backend=MatrixBackend([[-2, 0], [-1, 0], [1, 0], [2, 0]]),
        verbose=False,
    )
    assert len(result.clusters) == 1
    assert result.selected_indices == (1, 2)
    assert result.low_n_clusters == 1
    assert (settings.diversity, settings.low_n_clusters, settings.n_samples) == ("high", 5, 4)

    tied = diversity_sample(
        [record(9), record(8)],
        n_samples=1,
        diversity="low",
        low_n_clusters=1,
        config=settings,
        backend=MatrixBackend([[-1, 0], [1, 0]]),
        verbose=False,
    )
    assert tied.selected_indices == (0,)


@pytest.mark.parametrize("options", [{"diversity": "bad"}, {"low_n_clusters": 0}])
def test_invalid_low_diversity_arguments_fail_before_inference(options: dict[str, Any]) -> None:
    backend = MatrixBackend([[1, 0]])
    with pytest.raises(ConfigurationError):
        diversity_sample([record(0)], n_samples=1, backend=backend, verbose=False, **options)
    assert backend.calls == 0
