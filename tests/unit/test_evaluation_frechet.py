"""Tests for the LucaOne-backed Fréchet DNA representation distance."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Any, cast

import numpy as np
import pytest

from dnakit.core import DNARecord, DNASequence, DNASet
from dnakit.evaluation import FrechetDistanceConfig, evaluate_frechet_distance
from dnakit.exceptions import ConfigurationError
from dnakit.representations import RepresentationConfig


def _set(*pairs: tuple[str, str]) -> DNASet:
    return DNASet(DNARecord(DNASequence(symbols), record_id) for record_id, symbols in pairs)


class _MappingBackend:
    def __init__(self, vectors: Mapping[str, tuple[float, ...]]) -> None:
        self.vectors = vectors
        self.calls: list[tuple[tuple[str, ...], bool]] = []

    def extract(self, sequences: Sequence[str], *, show_progress: bool) -> Any:
        materialized = tuple(sequences)
        self.calls.append((materialized, show_progress))
        return np.asarray([self.vectors[sequence] for sequence in materialized], dtype=np.float32)


def test_frechet_defaults_to_lucaone_and_data_027_normalization() -> None:
    config = FrechetDistanceConfig()

    assert config.representation.model == "lucaone"
    assert config.normalize is True
    assert config.representation.show_progress is True

    records = _set(("a", "AAAA"), ("c", "CCCC"))
    with pytest.raises(ConfigurationError) as remote_code_error:
        evaluate_frechet_distance(records, records)
    assert remote_code_error.value.code == "MODEL_REMOTE_CODE_NOT_ALLOWED"


def test_frechet_known_one_dimensional_distance_is_symmetric() -> None:
    left = _set(("left-a", "AAAA"), ("left-c", "CCCC"))
    right = _set(("right-g", "GGGG"), ("right-t", "TTTT"))
    backend = _MappingBackend(
        {
            "AAAA": (0.0,),
            "CCCC": (2.0,),
            "GGGG": (1.0,),
            "TTTT": (3.0,),
        }
    )
    config = FrechetDistanceConfig(
        representation=RepresentationConfig(show_progress=False),
        normalize=False,
    )

    forward = evaluate_frechet_distance(left, right, config=config, backend=backend)
    reverse = evaluate_frechet_distance(right, left, config=config, backend=backend)

    assert forward.metrics["frechet_distance"] == pytest.approx(1.0)
    assert forward.metrics["mean_component"] == pytest.approx(1.0)
    assert forward.metrics["covariance_component"] == pytest.approx(0.0, abs=1e-12)
    assert reverse.metrics["frechet_distance"] == pytest.approx(
        forward.metrics["frechet_distance"]
    )
    assert backend.calls[0] == (("AAAA", "CCCC", "GGGG", "TTTT"), False)
    assert forward.parameters["encoder"] == "lucaone"
    assert forward.parameters["lower_is_better"] is True
    payload = forward.to_dict()
    assert payload["parameters"]["representation_config"]["dtype"] == "auto"
    json.dumps(payload)


def test_identical_normalized_representation_distributions_have_zero_distance() -> None:
    records = _set(("a", "AAAA"), ("c", "CCCC"), ("g", "GGGG"))
    backend = _MappingBackend(
        {
            "AAAA": (1.0, 0.0),
            "CCCC": (0.0, 2.0),
            "GGGG": (1.0, 1.0),
        }
    )
    report = evaluate_frechet_distance(
        records,
        records,
        config=FrechetDistanceConfig(
            representation=RepresentationConfig(show_progress=False)
        ),
        backend=backend,
    )

    assert report.metrics["frechet_distance"] == pytest.approx(0.0, abs=1e-12)
    assert report.metrics["embedding_dimension"] == 2
    assert report.metrics["left_covariance_rank_upper_bound"] == 2


def test_sample_space_result_matches_dense_covariance_formula() -> None:
    left = _set(("a", "AAAA"), ("c", "CCCC"), ("g", "GGGG"), ("t", "TTTT"))
    right = _set(
        ("ac", "ACAC"),
        ("ag", "AGAG"),
        ("at", "ATAT"),
        ("cg", "CGCG"),
        ("ct", "CTCT"),
    )
    vectors = {
        "AAAA": (0.0, 1.0, 2.0),
        "CCCC": (1.0, 2.0, 0.0),
        "GGGG": (2.0, 0.0, 1.0),
        "TTTT": (3.0, 2.0, 1.0),
        "ACAC": (0.5, 1.0, 1.5),
        "AGAG": (1.5, 1.0, 0.5),
        "ATAT": (2.5, 0.5, 1.0),
        "CGCG": (1.0, 2.5, 0.5),
        "CTCT": (3.0, 1.5, 2.0),
    }
    report = evaluate_frechet_distance(
        left,
        right,
        config=FrechetDistanceConfig(
            representation=RepresentationConfig(show_progress=False),
            normalize=False,
        ),
        backend=_MappingBackend(vectors),
    )
    left_matrix = np.asarray([vectors[record.sequence.symbols] for record in left], dtype=float)
    right_matrix = np.asarray(
        [vectors[record.sequence.symbols] for record in right], dtype=float
    )
    left_covariance = np.cov(left_matrix, rowvar=False, ddof=1)
    right_covariance = np.cov(right_matrix, rowvar=False, ddof=1)
    left_values, left_vectors = np.linalg.eigh(left_covariance)
    left_sqrt = (left_vectors * np.sqrt(np.clip(left_values, 0.0, None))) @ left_vectors.T
    middle = left_sqrt @ right_covariance @ left_sqrt
    cross_trace = float(np.sqrt(np.clip(np.linalg.eigvalsh(middle), 0.0, None)).sum())
    mean_delta = left_matrix.mean(axis=0) - right_matrix.mean(axis=0)
    expected = float(
        mean_delta @ mean_delta
        + np.trace(left_covariance)
        + np.trace(right_covariance)
        - 2.0 * cross_trace
    )

    assert report.metrics["frechet_distance"] == pytest.approx(expected, abs=1e-12)


def test_frechet_validates_sample_size_budget_and_zero_vectors() -> None:
    records = _set(("a", "AAAA"), ("c", "CCCC"))
    backend = _MappingBackend({"AAAA": (0.0, 0.0), "CCCC": (1.0, 0.0)})

    with pytest.raises(ConfigurationError) as sample_error:
        evaluate_frechet_distance(
            _set(("only", "AAAA")),
            records,
            backend=backend,
        )
    assert sample_error.value.code == "FRECHET_MINIMUM_SAMPLE_SIZE"
    assert backend.calls == []

    with pytest.raises(ConfigurationError) as budget_error:
        evaluate_frechet_distance(
            records,
            records,
            config=FrechetDistanceConfig(max_cross_gram_elements=3),
            backend=backend,
        )
    assert budget_error.value.code == "FRECHET_CROSS_GRAM_LIMIT"
    assert backend.calls == []

    with pytest.raises(ConfigurationError) as zero_error:
        evaluate_frechet_distance(records, records, backend=backend)
    assert zero_error.value.code == "ZERO_REPRESENTATION_VECTOR"


def test_frechet_config_rejects_invalid_values() -> None:
    with pytest.raises(ConfigurationError) as normalize_error:
        FrechetDistanceConfig(normalize=cast(Any, 1))
    assert normalize_error.value.code == "INVALID_FRECHET_DISTANCE_CONFIG"

    with pytest.raises(ConfigurationError):
        FrechetDistanceConfig(max_cross_gram_elements=0)
