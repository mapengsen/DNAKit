"""Numerical, filtering and input-contract checks for diploid dosage PCA."""

from __future__ import annotations

import copy
import json
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.genomics import population
from dnakit.genomics.population import genotype_pca

np = pytest.importorskip("numpy")


def test_known_covariance_eigenvalues_and_reconstruction() -> None:
    # Means are both 1 and scales sqrt(1/2), giving covariance [[2, 1], [1, 2]].
    result = genotype_pca([[0, 0], [1, 2], [2, 1]], verbose=False)
    np.testing.assert_allclose(result["eigenvalues"], [3, 1])
    np.testing.assert_allclose(result["explained_variance_ratio"], [0.75, 0.25])
    scores = np.asarray(result["scores"])
    loadings = np.asarray(result["loadings"])
    expected = np.sqrt(2) * np.array([[-1, -1], [0, 1], [1, 0]])
    np.testing.assert_allclose(scores @ loadings.T, expected, atol=1e-12)
    np.testing.assert_allclose(scores.T @ scores / 2, np.diag([3, 1]), atol=1e-12)
    np.testing.assert_allclose(loadings.T @ loadings, np.eye(2), atol=1e-12)
    np.testing.assert_allclose(scores.mean(axis=0), 0, atol=1e-12)
    assert result["sample_ids"] == ["sample_1", "sample_2", "sample_3"]
    assert result["marker_means"] == [1, 1]
    assert result["allele_frequencies"] == [0.5, 0.5]
    assert result["total_variance"] == pytest.approx(4)
    assert result["parameters"]["eigenvalue_denominator"] == "n_samples - 1"
    assert result["numerical_rank"] == 2
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("n_markers", [3, 11])
def test_covariance_eigendecomposition_matches_svd(n_markers: int) -> None:
    data = np.random.default_rng(89).uniform(0, 2, (8, n_markers))
    result = genotype_pca(data, components=2, verbose=False)
    frequencies = data.mean(axis=0) / 2
    standardized = (data - 2 * frequencies) / np.sqrt(2 * frequencies * (1 - frequencies))
    covariance = standardized.T @ standardized / 7
    eigenvalues, vectors = np.linalg.eigh(covariance)
    eigenvalues, vectors = eigenvalues[::-1], vectors[:, ::-1]
    loadings = np.asarray(result["loadings"])
    scores = np.asarray(result["scores"])
    np.testing.assert_allclose(result["eigenvalues"], eigenvalues[:2], rtol=1e-12)
    np.testing.assert_allclose(
        result["explained_variance_ratio"], eigenvalues[:2] / np.trace(covariance), rtol=1e-12
    )
    # Spectral projectors do not depend on arbitrary eigenvector signs.
    np.testing.assert_allclose(loadings @ loadings.T, vectors[:, :2] @ vectors[:, :2].T)
    np.testing.assert_allclose(scores, standardized @ loadings, atol=1e-12)
    assert sum(result["explained_variance_ratio"]) < 1
    assert np.all(loadings[np.argmax(np.abs(loadings), axis=0), np.arange(2)] >= 0)


def test_missing_imputation_filter_mapping_and_no_list_mutation() -> None:
    data = [
        [0, None, 1, 2, 0, 0.5],
        [None, None, 1, 2, 1, 0.5],
        [2, None, 1, 2, 2, None],
        [1, None, 1, 2, None, 0.5],
    ]
    original = copy.deepcopy(data)
    labels = ["a", "b", "c", "d"]
    result = genotype_pca(data, sample_ids=labels, verbose=False)
    assert data == original
    assert result["sample_ids"] == labels
    assert result["sample_ids"] is not labels
    assert result["marker_indices"] == [0, 4]
    assert result["dropped_all_missing_indices"] == [1]
    assert result["dropped_invariant_indices"] == [2, 3, 5]
    assert result["marker_means"] == [1, None, 1, 2, 1, 0.5]
    assert result["allele_frequencies"] == [0.5, None, 0.5, 1, 0.5, 0.25]
    assert result["imputed_value_count"] == 2
    assert result["imputed_counts_per_marker"] == [1, 1]
    assert result["n_markers"] == 6
    assert result["n_markers_used"] == 2
    assert result["numerical_rank"] == 1
    expected = np.sqrt(2) * np.array([[-1, -1], [0, 0], [1, 1], [0, 0]])
    np.testing.assert_allclose(
        np.asarray(result["scores"]) @ np.asarray(result["loadings"]).T, expected, atol=1e-12
    )
    json.dumps(result, allow_nan=False)


def test_numpy_missing_matches_null_and_input_is_read_only() -> None:
    data = np.array([[0, 0], [np.nan, 1], [2, 2], [1, np.nan]])
    original = data.copy()
    data.flags.writeable = False
    result = genotype_pca(data, verbose=False)
    from_json = genotype_pca([[0, 0], [None, 1], [2, 2], [1, None]], verbose=False)
    assert result == from_json
    np.testing.assert_equal(data, original)
    json.dumps(result, allow_nan=False)


def test_sample_and_marker_reordering_preserve_sample_geometry() -> None:
    data = np.array([[0, 1, 0.3], [1, 2, 1.2], [2, 0, 2], [0.5, 1.1, 0.4]])
    result = genotype_pca(data, components=3, verbose=False)
    permutation = [3, 0, 2, 1]
    reordered = genotype_pca(data[permutation][:, ::-1], components=3, verbose=False)
    scores = np.asarray(result["scores"])
    reordered_scores = np.asarray(reordered["scores"])
    np.testing.assert_allclose(result["eigenvalues"], reordered["eigenvalues"], atol=1e-12)
    np.testing.assert_allclose(
        scores[permutation] @ scores[permutation].T,
        reordered_scores @ reordered_scores.T,
        atol=1e-12,
    )


def test_allele_recoding_preserves_sample_geometry_and_variance() -> None:
    data = np.array([[0, 1, 0.3], [1, 2, 1.2], [2, 0, 2], [0.5, np.nan, 0.4]])
    original = genotype_pca(data, verbose=False)
    recoded_data = data.copy()
    recoded_data[:, 1:] = 2 - recoded_data[:, 1:]
    recoded = genotype_pca(recoded_data, verbose=False)
    original_scores, recoded_scores = np.asarray(original["scores"]), np.asarray(recoded["scores"])
    np.testing.assert_allclose(original["eigenvalues"], recoded["eigenvalues"])
    np.testing.assert_allclose(
        original_scores @ original_scores.T, recoded_scores @ recoded_scores.T, atol=1e-12
    )


def test_single_marker_has_known_scores_and_one_component() -> None:
    result = genotype_pca([[0], [1], [2]], components=1, verbose=False)
    np.testing.assert_allclose(result["scores"], [[-np.sqrt(2)], [0], [np.sqrt(2)]])
    np.testing.assert_allclose(result["loadings"], [[1]])
    np.testing.assert_allclose(result["eigenvalues"], [2])
    np.testing.assert_allclose(result["explained_variance_ratio"], [1])


@pytest.mark.parametrize(
    "data",
    [[], [0, 1], [[0]], [[], []], [[0], [1, 2]], [[[0]], [[2]]], None, "012"],
)
def test_rejects_nonmatrices_and_invalid_dimensions(data: Any) -> None:
    with pytest.raises(ConfigurationError):
        genotype_pca(data, components=1, verbose=False)


@pytest.mark.parametrize(
    "value", [True, False, np.bool_(True), "1", "NaN", "null", 1j, float("inf"), -1, 3]
)
def test_rejects_invalid_dosages(value: Any) -> None:
    with pytest.raises(ConfigurationError):
        genotype_pca([[0], [value], [2]], components=1, verbose=False)


@pytest.mark.parametrize("dtype", [bool, complex, str])
def test_rejects_invalid_numpy_matrix_types(dtype: Any) -> None:
    with pytest.raises(ConfigurationError):
        genotype_pca(np.array([[0], [1]], dtype=dtype), components=1, verbose=False)


@pytest.mark.parametrize(
    "data", [[[None], [None]], [[0], [0]], [[1], [1]], [[2], [2]], [[0.2], [None], [0.2]]]
)
def test_rejects_all_missing_or_invariant_markers(data: Any) -> None:
    with pytest.raises(ConfigurationError, match="No variable markers"):
        genotype_pca(data, components=1, verbose=False)


@pytest.mark.parametrize("components", [0, -1, True, 1.0, "1", None, 3])
def test_rejects_invalid_component_counts(components: Any) -> None:
    with pytest.raises(ConfigurationError, match="components"):
        genotype_pca([[0, 0], [1, 2], [2, 1]], components=components, verbose=False)


def test_component_limit_applies_after_filtering_and_centering() -> None:
    with pytest.raises(ConfigurationError, match="<= 1"):
        genotype_pca([[0, 1], [1, 1], [2, 1]], components=2, verbose=False)
    with pytest.raises(ConfigurationError, match="<= 1"):
        genotype_pca([[0, 1], [2, 0]], components=2, verbose=False)


@pytest.mark.parametrize(
    "labels", [["a"], ["a", "a", "c"], ["a", "", "c"], ["a", " ", "c"], [1, 2, 3], "abc"]
)
def test_rejects_invalid_sample_names(labels: Any) -> None:
    with pytest.raises(ConfigurationError, match="sample_ids"):
        genotype_pca([[0], [1], [2]], components=1, sample_ids=labels, verbose=False)


def test_rejects_nonboolean_verbosity_and_is_quiet_when_disabled(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(ConfigurationError, match="verbose"):
        genotype_pca([[0], [2]], components=1, verbose="false")  # type: ignore[arg-type]
    genotype_pca([[0], [2]], components=1, verbose=False)
    output = capsys.readouterr()
    assert output.out == output.err == ""


def test_failed_svd_is_reported_as_configuration_error(monkeypatch: pytest.MonkeyPatch) -> None:
    def failed_svd(*args: Any, **kwargs: Any) -> Any:
        raise np.linalg.LinAlgError("injected failure")

    monkeypatch.setattr(np.linalg, "svd", failed_svd)
    with pytest.raises(ConfigurationError, match="SVD did not converge"):
        genotype_pca([[0], [2]], components=1, verbose=False)


def test_numpy_dependency_uses_existing_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str, str]] = []

    def recorded_require(module: str, extra: str) -> Any:
        calls.append((module, extra))
        return np

    monkeypatch.setattr(population, "require", recorded_require)
    genotype_pca([[0], [2]], components=1, verbose=False)
    assert calls == [("numpy", "bio-analysis")]
