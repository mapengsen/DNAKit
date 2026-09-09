"""Distribution semantics, CSV validation, and CLI image export."""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import cast

import pytest
from typer.testing import CliRunner

from dnakit.cli.app import app
from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.visualization import (
    DistributionPlotConfig,
    SVGTheme,
    plot_csv_distributions,
    plot_distribution,
)


def test_classification_counts_and_missing_denominator() -> None:
    plot = plot_distribution(
        ["0", "1", "1", "2", "", None, " "],
        task="classification",
        config=DistributionPlotConfig(categories=("0", "1", "2", "3")),
    )
    assert plot.metadata["counts"] == (1, 2, 1, 0)
    assert plot.metadata["proportions"] == (0.25, 0.5, 0.25, 0.0)
    assert plot.metadata["valid_count"] == 4
    assert plot.metadata["total_count"] == 7
    assert plot.metadata["missing_count"] == 3
    assert plot.issues[0].code == "VIZ_MISSING_VALUES"


def test_regression_boundary_inclusion_and_unequal_bins() -> None:
    plot = plot_distribution(
        [-1, 0, 0.5, 1, 4, None],
        task="regression",
        config=DistributionPlotConfig(bin_edges=(-1, 0, 1, 4)),
    )
    assert plot.metadata["counts"] == (1, 2, 2)
    assert plot.metadata["proportions"] == (0.2, 0.4, 0.4)
    # Different bin widths still use sample proportions, not density.
    bars = [
        node for node in ET.fromstring(plot.svg).iter() if node.get("class") == "distribution-bar"
    ]
    assert bars[1].get("height") == bars[2].get("height")
    assert float(bars[2].attrib["width"]) == pytest.approx(float(bars[1].attrib["width"]) * 3)


def test_regression_equal_bins_and_constant_column() -> None:
    plot = plot_distribution(
        [0, 1, 2, 3, 4], task="regression", config=DistributionPlotConfig(bins=4)
    )
    assert plot.metadata["bin_edges"] == (0, 1, 2, 3, 4)
    assert plot.metadata["counts"] == (1, 1, 1, 2)
    constant = plot_distribution([2, 2, 2], task="regression")
    assert constant.metadata["counts"] == (3,)
    assert constant.metadata["proportions"] == (1.0,)


def test_histogram_matches_numpy_on_many_samples() -> None:
    np = pytest.importorskip("numpy")
    values = np.random.default_rng(123).normal(size=2000).tolist()
    counts, edges = np.histogram(values, bins=37)
    plot = plot_distribution(values, task="regression", config=DistributionPlotConfig(bins=37))
    assert plot.metadata["counts"] == tuple(counts.tolist())
    assert plot.metadata["bin_edges"] == pytest.approx(edges.tolist())
    assert sum(cast(tuple[float, ...], plot.metadata["proportions"])) == pytest.approx(1.0)


@pytest.mark.parametrize("values", [[], [None, " "], ["oops"], [float("nan")], [float("inf")]])
def test_invalid_regression_values(values: list[float | str | None]) -> None:
    with pytest.raises(InputFormatError):
        plot_distribution(values, task="regression")


def test_configuration_and_coverage_errors() -> None:
    with pytest.raises(ConfigurationError):
        DistributionPlotConfig(bin_edges=(0, 1, 1))
    with pytest.raises(ConfigurationError):
        DistributionPlotConfig(bin_edges=(0, math.inf))
    with pytest.raises(ConfigurationError):
        DistributionPlotConfig(bar_width=0)
    with pytest.raises(ConfigurationError):
        DistributionPlotConfig(categories=("a", "a"))
    with pytest.raises(ConfigurationError):
        plot_distribution([1, 2], task="regression", config=DistributionPlotConfig(max_values=1))
    with pytest.raises(InputFormatError, match="cover every"):
        plot_distribution(
            [0, 2], task="regression", config=DistributionPlotConfig(bin_edges=(0, 1))
        )
    with pytest.raises(InputFormatError, match="every observed"):
        plot_distribution(
            ["a", "b"], task="classification", config=DistributionPlotConfig(categories=("a",))
        )


def test_csv_multiple_columns_with_bom_and_quoted_label(tmp_path: Path) -> None:
    source = tmp_path / "data.csv"
    source.write_text('a,b\n"x,y",1\n,2\nz,2\n', encoding="utf-8-sig")
    plots = plot_csv_distributions(source, columns=["a", "b"], task="classification")
    assert plots["a"].metadata["labels"] == ("x,y", "z")
    assert plots["a"].metadata["valid_count"] == 2
    assert plots["b"].metadata["valid_count"] == 3
    assert plots["b"].metadata["proportions"] == (1 / 3, 2 / 3)


@pytest.mark.parametrize(
    "content,columns",
    [
        ("a,a\n1,2\n", ["a"]),
        ("a,b\n1\n", ["a"]),
        ("a,b\n1,2,3\n", ["a"]),
        ("a,b\n1,2\n", ["absent"]),
        ("a\n", ["a"]),
        ('a\n"unclosed', ["a"]),
    ],
)
def test_csv_errors(tmp_path: Path, content: str, columns: list[str]) -> None:
    source = tmp_path / "bad.csv"
    source.write_text(content)
    with pytest.raises(InputFormatError):
        plot_csv_distributions(source, columns=columns, task="classification")


def test_style_xml_escaping_and_many_categories() -> None:
    config = DistributionPlotConfig(
        size=900,
        font_size=20,
        tick_font_size=15,
        show_values=True,
        label_rotation=45,
        theme=SVGTheme(accent="#ff0000", font_family="serif"),
    )
    plot = plot_distribution(
        ["<label&>", *map(str, range(35))], task="classification", config=config
    )
    nodes = list(ET.fromstring(plot.svg).iter())
    bars = [node for node in nodes if node.get("class") == "distribution-bar"]
    assert len(bars) == 36
    assert bars[0].get("data-label") == "<label&>"
    assert bars[0].get("fill") == "#ff0000"
    assert plot.width == plot.height == 900
    assert len([node for node in nodes if node.get("class") == "category-tick"]) <= 12
    assert any(node.get("font-size") == "20" for node in nodes)
    assert not [node for node in nodes if node.get("class") == "plot-title"]


def test_cli_multiple_outputs_and_overwrite(tmp_path: Path) -> None:
    source = tmp_path / "data.csv"
    source.write_text("a,b\n0,1\n1,2\n1,3\n")
    first, second = tmp_path / "plots/a.svg", tmp_path / "plots/b.svg"
    args = [
        "plot-distribution",
        "--input",
        str(source),
        "--task",
        "classification",
        "--column",
        "a",
        "--column",
        "b",
        "--output",
        str(first),
        "--output",
        str(second),
    ]
    runner = CliRunner()
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert first.is_file() and second.is_file()
    ET.parse(first)
    original = first.read_bytes()
    result = runner.invoke(app, args)
    assert result.exit_code == 2
    assert first.read_bytes() == original
    assert runner.invoke(app, [*args, "--overwrite"]).exit_code == 0


def test_cli_reports_bad_column_before_writing(tmp_path: Path) -> None:
    source = tmp_path / "data.csv"
    source.write_text("a,b\n1,oops\n")
    first, second = tmp_path / "a.svg", tmp_path / "b.svg"
    result = CliRunner().invoke(
        app,
        [
            "plot-distribution",
            "--input",
            str(source),
            "--task",
            "regression",
            "--column",
            "a",
            "--column",
            "b",
            "--output",
            str(first),
            "--output",
            str(second),
        ],
    )
    assert result.exit_code == 2
    assert "Column 'b'" in result.output
    assert not first.exists() and not second.exists()


def test_cli_regression_export_png(tmp_path: Path) -> None:
    pytest.importorskip("cairosvg")
    image = pytest.importorskip("PIL.Image")
    source = tmp_path / "data.csv"
    source.write_text("value\n-1\n0\n1\n2\n")
    target = tmp_path / "histogram"
    result = CliRunner().invoke(
        app,
        [
            "plot-distribution",
            "--input",
            str(source),
            "--task",
            "regression",
            "--column",
            "value",
            "--output",
            str(target),
            "--bin-edges=-1,0,2",
            "--dpi",
            "96",
            "--size",
            "500",
        ],
    )
    assert result.exit_code == 0, result.output
    with image.open(target.with_suffix(".png")) as output:
        assert output.size == (500, 500)
