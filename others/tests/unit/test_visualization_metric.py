"""CSV metric comparison and optional fold/test inputs."""

from __future__ import annotations

from pathlib import Path

import pytest

from dnakit.exceptions import InputFormatError
from dnakit.visualization import plot_metric_comparison_csv


def test_metric_comparison_needs_only_one_csv_and_saves_image(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    source = tmp_path / "summary.csv"
    source.write_text(
        "model,mean,se,source\nA,0.51,0.02,Reference\nB,0.64,0.01,Reference\n",
        encoding="utf-8",
    )
    output = tmp_path / "plots" / "comparison.svg"
    assert plot_metric_comparison_csv(source, output=output, metric_label="MCC") == output
    assert output.read_text(encoding="utf-8").lstrip().startswith("<?xml")
    with pytest.raises(FileExistsError):
        plot_metric_comparison_csv(source, output=output)
    assert plot_metric_comparison_csv(source, output=output, overwrite=True) == output
    png = tmp_path / "plots" / "comparison.png"
    assert plot_metric_comparison_csv(source, output=png) == png
    assert png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_metric_comparison_computes_fold_se_and_filters_test_csv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    axes_module = pytest.importorskip("matplotlib.axes")
    reference = tmp_path / "reference.csv"
    reference.write_text(
        "dataset,source,model,mean,se\nd1,Table S3,Reference model,0.6,0.01\n",
        encoding="utf-8",
    )
    folds = tmp_path / "folds.csv"
    folds.write_text("fold,mcc\n1,0.4\n2,0.6\n", encoding="utf-8")
    tests = tmp_path / "tests.csv"
    tests.write_text("dataset,method,mcc\nd2,full,0.2\nd1,full,0.7\n", encoding="utf-8")
    seen: list[tuple[float, float | None]] = []
    original = axes_module.Axes.errorbar

    def capture(self, x, y, **kwargs):
        seen.append((y, kwargs["yerr"]))
        return original(self, x, y, **kwargs)

    monkeypatch.setattr(axes_module.Axes, "errorbar", capture)
    plot_metric_comparison_csv(
        reference,
        fold_csvs={"Our model": folds},
        test_csv=tests,
        test_model_names={"full": "Our model"},
        dataset="d1",
        output=tmp_path / "comparison.svg",
    )
    assert seen[0] == (0.6, 0.01)
    assert seen[1][0] == pytest.approx(0.5)
    assert seen[1][1] == pytest.approx(0.1)
    assert seen[2] == (0.7, None)


@pytest.mark.parametrize(
    "content",
    [
        "model,se\nA,0.1\n",
        "model,mean,se\nA,nan,0.1\n",
        "model,mean,se\nA,0.5,-0.1\n",
    ],
)
def test_metric_comparison_rejects_invalid_summary(tmp_path: Path, content: str) -> None:
    pytest.importorskip("matplotlib")
    source = tmp_path / "invalid.csv"
    source.write_text(content, encoding="utf-8")
    output = tmp_path / "plot.svg"
    with pytest.raises(InputFormatError):
        plot_metric_comparison_csv(source, output=output)
    assert not output.exists()


def test_metric_comparison_rejects_missing_reference_dataset(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    source = tmp_path / "summary.csv"
    source.write_text("dataset,model,mean\nd1,A,0.5\n", encoding="utf-8")
    with pytest.raises(InputFormatError, match="No reference metric rows"):
        plot_metric_comparison_csv(source, dataset="d2", output=tmp_path / "plot.svg")
