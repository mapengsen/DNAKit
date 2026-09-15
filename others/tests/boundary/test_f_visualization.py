"""F5: probability semantics, missing values, bins and CSV/CLI failures."""

from __future__ import annotations

import math
from typing import Any, cast

import pytest
from typer.testing import CliRunner

from dnakit.cli.app import app
from dnakit.exceptions import InputFormatError
from dnakit.visualization import DistributionPlotConfig, plot_csv_distributions, plot_distribution


@pytest.mark.parametrize("values", [[], [None, "", " "]])
@pytest.mark.parametrize("task", ["classification", "regression", "probability"])
def test_no_valid_samples_has_clear_error(values: Any, task: Any) -> None:
    with pytest.raises(InputFormatError, match="No valid samples"):
        plot_distribution(values, task=task, verbose=False)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf, "NaN", "Inf"])
def test_nonfinite_probability_is_rejected(value: Any) -> None:
    with pytest.raises(InputFormatError, match="finite"):
        plot_distribution([value], task="probability", verbose=False)


def test_missing_probability_denominator() -> None:
    plot = plot_distribution(
        [0.2, None, 0.8],
        task="probability",
        config=DistributionPlotConfig(bin_edges=[0, 0.5, 1]),
        verbose=False,
    )
    assert plot.metadata["valid_count"] == 2
    assert plot.metadata["missing_count"] == 1
    assert plot.metadata["proportions"] == (0.5, 0.5)


def test_probability_endpoints_and_inner_boundary() -> None:
    plot = plot_distribution(
        [0, 0.5, 1],
        task="probability",
        config=DistributionPlotConfig(bin_edges=[0, 0.5, 1]),
        verbose=False,
    )
    assert plot.metadata["counts"] == (1, 2)
    assert sum(cast(tuple[int, ...], plot.metadata["counts"])) == 3


@pytest.mark.parametrize(
    "values,task",
    [([0.5] * 100, "probability"), ([0.5], "probability"), (["1"] * 100, "classification")],
)
def test_constant_or_single_sample(values: Any, task: Any) -> None:
    plot = plot_distribution(values, task=task, verbose=False)
    assert sum(cast(tuple[int, ...], plot.metadata["counts"])) == len(values)
    assert sum(cast(tuple[float, ...], plot.metadata["proportions"])) == pytest.approx(1)
    assert "nan" not in plot.svg.lower()


@pytest.mark.parametrize("values", [[-0.1], [1.1], [-4, 2]])
def test_logits_or_out_of_range_probabilities_rejected(values: Any) -> None:
    with pytest.raises(InputFormatError, match="INVALID_PROBABILITY"):
        plot_distribution(values, task="probability", verbose=False)
    assert plot_distribution(values, task="regression", verbose=False).metadata[
        "valid_count"
    ] == len(values)


def test_csv_missing_column_and_out_of_range_probability(tmp_path: Any) -> None:
    source = tmp_path / "probabilities.csv"
    source.write_text("probability_1,logits\n0.2,-5\n0.8,4\n")
    with pytest.raises(InputFormatError, match="not found"):
        plot_csv_distributions(source, columns=["missing"], task="probability", verbose=False)
    with pytest.raises(InputFormatError, match="logits"):
        plot_csv_distributions(source, columns=["logits"], task="probability", verbose=False)
    assert (
        plot_csv_distributions(
            source, columns=["probability_1"], task="probability", verbose=False
        )["probability_1"].metadata["valid_count"]
        == 2
    )


def test_cli_invalid_probability_preserves_output(tmp_path: Any) -> None:
    source, target = tmp_path / "input.csv", tmp_path / "plot.svg"
    source.write_text("probability\n1.1\n")
    target.write_text("previous image")
    result = CliRunner().invoke(
        app,
        [
            "plot-distribution",
            "--input",
            str(source),
            "--column",
            "probability",
            "--task",
            "probability",
            "--output",
            str(target),
            "--overwrite",
        ],
    )
    assert result.exit_code != 0
    assert "[0, 1]" in result.output
    assert target.read_text() == "previous image"
