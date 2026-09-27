"""Known-value and validation checks for the independent breeding baselines."""

from __future__ import annotations

import json
from typing import Any

import pytest

from dnakit.exceptions import ConfigurationError
from dnakit.genomics.breeding import genomic_prediction, genomic_relationship_matrix

np = pytest.importorskip("numpy")


def test_vanraden_known_values_and_invariants() -> None:
    result = genomic_relationship_matrix([[0, 0], [1, 2], [2, 1]], verbose=False)
    expected = [[2, -1, -1], [-1, 1, 0], [-1, 0, 1]]
    np.testing.assert_allclose(result["relationship_matrix"], expected)
    matrix = np.asarray(result["relationship_matrix"])
    np.testing.assert_allclose(matrix.sum(axis=0), 0, atol=1e-12)
    assert np.linalg.eigvalsh(matrix).min() >= -1e-12
    assert result["denominator"] == 1
    assert result["allele_frequencies"] == [0.5, 0.5]
    json.dumps(result, allow_nan=False)


def test_grm_missing_values_filtering_and_no_input_mutation() -> None:
    data = np.array([[0, np.nan, 2], [np.nan, np.nan, 2], [2, np.nan, 2]])
    original = data.copy()
    result = genomic_relationship_matrix(data, verbose=False)
    np.testing.assert_allclose(result["relationship_matrix"], [[2, 0, -2], [0, 0, 0], [-2, 0, 2]])
    assert result["marker_indices"] == [0]
    assert result["dropped_all_missing_indices"] == [1]
    assert result["dropped_monomorphic_indices"] == [2]
    assert result["imputed_value_count"] == 1
    assert result["allele_frequencies"] == [0.5, None, 1.0]
    np.testing.assert_equal(data, original)
    json.dumps(result, allow_nan=False)


def test_all_heterozygous_column_is_not_allele_monomorphic() -> None:
    result = genomic_relationship_matrix([[0, 1], [2, 1]], verbose=False)
    np.testing.assert_allclose(result["relationship_matrix"], [[1, -1], [-1, 1]])
    assert result["n_markers_used"] == 2


def test_fixed_lambda_single_marker_analytical_solution() -> None:
    result = genomic_prediction(
        [[0], [1], [2]], [1, 2, 4], [[0], [1], [2]], ridge_lambda=2, verbose=False
    )
    # Centered covariance is 3, genotype sum of squares is 2: b = 3/(2+2).
    np.testing.assert_allclose(result["marker_effects"], [0.75])
    np.testing.assert_allclose(result["predictions"], [19 / 12, 28 / 12, 37 / 12])
    assert result["intercept"] == pytest.approx(7 / 3)
    assert result["raw_intercept"] == pytest.approx(19 / 12)
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("n_markers", [2, 8])
def test_ridge_matches_independent_augmented_least_squares(n_markers: int) -> None:
    rng = np.random.default_rng(742)
    genotypes = rng.uniform(0, 2, size=(5, n_markers))
    phenotypes = rng.normal(size=5)
    prediction = rng.uniform(0, 2, size=(3, n_markers))
    penalty = 1.75
    # Independent formulation: fit raw intercept and marker effects with an
    # augmented least-squares system, penalizing only the marker coefficients.
    design = np.column_stack((np.ones(5), genotypes))
    pseudo_rows = np.column_stack((np.zeros(n_markers), np.sqrt(penalty) * np.eye(n_markers)))
    coefficients = np.linalg.lstsq(
        np.vstack((design, pseudo_rows)), np.r_[phenotypes, np.zeros(n_markers)], rcond=None
    )[0]
    expected = np.column_stack((np.ones(3), prediction)) @ coefficients
    result = genomic_prediction(
        genotypes, phenotypes, prediction, ridge_lambda=penalty, verbose=False
    )
    np.testing.assert_allclose(result["predictions"], expected, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(result["marker_effects"], coefficients[1:], atol=1e-10)


def test_prediction_uses_only_training_statistics_and_does_not_mutate() -> None:
    training = np.array([[0, 0, np.nan, 2], [np.nan, 1, np.nan, 2], [2, 2, np.nan, 2]])
    original = training.copy()
    kwargs: dict[str, Any] = {"ridge_lambda": 2, "verbose": False}
    alone = genomic_prediction(training, [1, 2, 4], [[np.nan, 2, 1, 0]], **kwargs)
    batch = genomic_prediction(
        training, [1, 2, 4], [[np.nan, 2, 1, 0], [2, 2, 2, 2], [2, 0, 0, 1]], **kwargs
    )
    assert batch["predictions"][0] == pytest.approx(alone["predictions"][0])
    np.testing.assert_allclose(batch["marker_effects"], alone["marker_effects"])
    assert alone["allele_frequencies"] == [0.5, 0.5, None, 1]
    assert alone["marker_effects"][2:] == [0, 0]
    assert alone["imputed_training_value_count"] == 1
    assert alone["imputed_prediction_value_count"] == 1
    np.testing.assert_equal(training, original)


@pytest.mark.parametrize("ploidy", [1, 3, 4, True, 2.0])
def test_rejects_unsupported_ploidy(ploidy: Any) -> None:
    with pytest.raises(ConfigurationError, match="diploid"):
        genomic_relationship_matrix([[0], [2]], ploidy=ploidy, verbose=False)


@pytest.mark.parametrize(
    "values",
    [[], [[0]], [[0], [3]], [[0], [-1]], [[0], [float("inf")]], [[0], [1, 2]], [[0], [1j]]],
)
def test_rejects_invalid_genotype_matrix(values: object) -> None:
    with pytest.raises(ConfigurationError):
        genomic_relationship_matrix(values, verbose=False)


@pytest.mark.parametrize("values", [[[0], [0]], [[2], [2]], [[None], [None]]])
def test_rejects_all_unusable_markers(values: object) -> None:
    with pytest.raises(ConfigurationError, match="No observed polymorphic"):
        genomic_relationship_matrix(values, verbose=False)


@pytest.mark.parametrize("penalty", [0, -1, float("nan"), float("inf"), True, "bad"])
def test_rejects_invalid_penalty(penalty: Any) -> None:
    with pytest.raises(ConfigurationError, match="ridge_lambda"):
        genomic_prediction([[0], [2]], [1, 2], [[1]], ridge_lambda=penalty, verbose=False)


def test_missing_policy_phenotypes_and_prediction_columns() -> None:
    with pytest.raises(ConfigurationError, match="missing"):
        genomic_relationship_matrix([[0], [None]], missing="error", verbose=False)
    with pytest.raises(ConfigurationError, match="missing"):
        genomic_relationship_matrix([[0], [2]], missing="invalid", verbose=False)
    for response in ([1], [1, None], [[1], [2]], [1, float("inf")]):
        with pytest.raises(ConfigurationError, match="phenotype"):
            genomic_prediction([[0], [2]], response, [[1]], ridge_lambda=1, verbose=False)
    with pytest.raises(ConfigurationError, match="same columns"):
        genomic_prediction([[0], [2]], [1, 2], [[1, 2]], ridge_lambda=1, verbose=False)
    with pytest.raises(ConfigurationError, match="missing"):
        genomic_prediction(
            [[0], [2]], [1, 2], [[None]], ridge_lambda=1, missing="error", verbose=False
        )


def test_verbose_false_leaves_output_clean(capsys: pytest.CaptureFixture[str]) -> None:
    genomic_relationship_matrix([[0], [2]], verbose=False)
    captured = capsys.readouterr()
    assert captured.out == captured.err == ""


def test_json_null_dosages_match_numpy_missing_values_and_remain_json_safe() -> None:
    training = json.loads("[[0, null, 2], [null, null, 2], [2, null, 2]]")
    prediction = json.loads("[[null, 0, 0], [2, 2, 2]]")
    from_json = genomic_relationship_matrix(training, verbose=False)
    from_numpy = genomic_relationship_matrix(np.asarray(training, dtype=float), verbose=False)
    assert from_json == from_numpy
    predicted = genomic_prediction(training, [0, 1, 2], prediction, ridge_lambda=1, verbose=False)
    assert predicted["imputed_training_value_count"] == 1
    assert predicted["imputed_prediction_value_count"] == 1
    assert predicted["predictions"][0] == pytest.approx(1)
    assert predicted["allele_frequencies"] == [0.5, None, 1.0]
    json.dumps(predicted, allow_nan=False)


def test_allele_recoding_preserves_relationships_and_predictions() -> None:
    training = np.array([[0, 2, 0.2], [None, 1, 0.8], [2, 0, 1.4], [1, 1, 2]], dtype=float)
    prediction = np.array([[0, 0, 1.3], [2, None, 0.7]], dtype=float)
    flipped_training, flipped_prediction = training.copy(), prediction.copy()
    flipped_training[:, 1:] = 2 - flipped_training[:, 1:]
    flipped_prediction[:, 1:] = 2 - flipped_prediction[:, 1:]
    original_grm = genomic_relationship_matrix(training, verbose=False)
    flipped_grm = genomic_relationship_matrix(flipped_training, verbose=False)
    np.testing.assert_allclose(
        original_grm["relationship_matrix"], flipped_grm["relationship_matrix"], atol=1e-12
    )
    original = genomic_prediction(
        training, [2, 1, 3, 4], prediction, ridge_lambda=0.25, verbose=False
    )
    flipped = genomic_prediction(
        flipped_training, [2, 1, 3, 4], flipped_prediction, ridge_lambda=0.25, verbose=False
    )
    np.testing.assert_allclose(original["predictions"], flipped["predictions"], atol=1e-12)
    np.testing.assert_allclose(
        np.asarray(original["marker_effects"]) * [1, -1, -1], flipped["marker_effects"], atol=1e-12
    )


@pytest.mark.parametrize("value", [True, False, "0", "NaN", "null", "1.2"])
def test_json_text_and_boolean_are_never_coerced_into_dosages(value: object) -> None:
    with pytest.raises(ConfigurationError, match="numeric"):
        genomic_relationship_matrix([[0], [value], [2]], verbose=False)
    with pytest.raises(ConfigurationError, match="numeric"):
        genomic_prediction([[0], [2]], [1, 2], [[value]], ridge_lambda=1, verbose=False)


@pytest.mark.parametrize("value", [True, "1", "NaN", 1j])
def test_nonnumeric_phenotypes_are_rejected(value: object) -> None:
    with pytest.raises(ConfigurationError, match="phenotypes"):
        genomic_prediction([[0], [2]], [1, value], [[1]], ridge_lambda=1, verbose=False)


@pytest.mark.parametrize("value", ["1", "1e-3", np.bool_(True)])
def test_nonnumeric_penalty_is_not_coerced(value: object) -> None:
    with pytest.raises(ConfigurationError, match="ridge_lambda"):
        genomic_prediction([[0], [2]], [1, 2], [[1]], ridge_lambda=value, verbose=False)
