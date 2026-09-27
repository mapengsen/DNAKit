"""One-CSV learning curves and fold aggregation."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from dnakit.exceptions import ConfigurationError, InputFormatError
from dnakit.visualization import plot_learning_curves_csv


def test_learning_curves_group_csv_and_average_folds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    axes_module = pytest.importorskip("matplotlib.axes")
    source = tmp_path / "epochs.csv"
    source.write_text(
        "dataset,setting,method,fold,epoch,train_loss,validation_loss,train_mcc,validation_mcc\n"
        "enhancers,baseline,full,0,1,0.4,0.5,0.6,0.4\n"
        "enhancers,baseline,full,1,1,0.6,0.7,0.8,0.6\n"
        "enhancers,baseline,linear_probe,0,1,0.8,0.9,0.2,0.1\n"
        "enhancers,selected,full,0,1,0.3,0.4,0.9,0.8\n"
        "splice_sites_donors,baseline,full,0,1,0.7,0.8,0.5,0.4\n",
        encoding="utf-8",
    )
    plotted: list[tuple[list[int], list[float], str]] = []
    bands: list[tuple[list[int], list[float], list[float]]] = []
    original_plot = axes_module.Axes.plot
    original_fill = axes_module.Axes.fill_between

    def capture_plot(axis, x, y, **kwargs):
        plotted.append((list(x), list(y), kwargs["label"]))
        return original_plot(axis, x, y, **kwargs)

    def capture_fill(axis, x, lower, upper, **kwargs):
        bands.append((list(x), list(lower), list(upper)))
        return original_fill(axis, x, lower, upper, **kwargs)

    monkeypatch.setattr(axes_module.Axes, "plot", capture_plot)
    monkeypatch.setattr(axes_module.Axes, "fill_between", capture_fill)
    output = tmp_path / "figures"
    paths = plot_learning_curves_csv(source, output=output)

    assert set(paths) == {
        ("enhancers", "baseline"),
        ("enhancers", "selected"),
        ("splice_sites_donors", "baseline"),
    }
    assert {path.name for path in paths.values()} == {
        "enhancers_baseline.svg",
        "enhancers_selected.svg",
        "splice_sites_donors_baseline.svg",
    }
    assert all(
        path.read_text(encoding="utf-8").lstrip().startswith("<?xml") for path in paths.values()
    )
    assert "<text" in paths[("enhancers", "baseline")].read_text(encoding="utf-8")
    assert plotted[0] == ([1], [0.5], "full/train")
    assert bands[0][0] == [1]
    assert bands[0][1] == pytest.approx([0.5 - math.sqrt(0.02)])
    assert bands[0][2] == pytest.approx([0.5 + math.sqrt(0.02)])
    assert ([1], [0.7], "full/train") in plotted


def test_learning_curves_simple_csv_saves_one_image(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    source = tmp_path / "history.csv"
    source.write_text(
        "epoch,train_loss,validation_loss,train_mcc,validation_mcc\n"
        "1,0.6,0.7,0.4,0.3\n"
        "2,0.5,0.6,,0.4\n",
        encoding="utf-8",
    )
    target = tmp_path / "plot.png"

    assert plot_learning_curves_csv(source, output=target) == {("history", ""): target}
    assert target.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    with pytest.raises(FileExistsError):
        plot_learning_curves_csv(source, output=target)
    assert plot_learning_curves_csv(source, output=target, overwrite=True) == {
        ("history", ""): target
    }


@pytest.mark.parametrize(
    "content",
    [
        "epoch,train_loss,validation_loss,train_mcc\n1,0.3,0.4,0.5\n",
        "epoch,train_loss,validation_loss,train_mcc,validation_mcc\n",
        "epoch,train_loss,validation_loss,train_mcc,validation_mcc\nx,0.3,0.4,0.5,0.6\n",
        "epoch,train_loss,validation_loss,train_mcc,validation_mcc\n1,nan,0.4,0.5,0.6\n",
        "epoch,train_loss,validation_loss,train_mcc,validation_mcc\n1,,,,\n",
        "method,fold,epoch,train_loss,validation_loss,train_mcc,validation_mcc\n"
        "full,0,1,0.3,0.4,0.5,0.6\nfull,0,1,0.2,0.4,0.5,0.6\n",
    ],
)
def test_learning_curves_reject_invalid_csv(tmp_path: Path, content: str) -> None:
    source = tmp_path / "bad.csv"
    source.write_text(content, encoding="utf-8")

    with pytest.raises(InputFormatError):
        plot_learning_curves_csv(source, output=tmp_path / "plot.svg")
    assert not (tmp_path / "plot.svg").exists()


def test_learning_curves_need_directory_for_multiple_groups(tmp_path: Path) -> None:
    source = tmp_path / "epochs.csv"
    source.write_text(
        "dataset,epoch,train_loss,validation_loss,train_mcc,validation_mcc\n"
        "first,1,0.3,0.4,0.5,0.6\n"
        "second,1,0.2,0.3,0.7,0.8\n",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="output directory"):
        plot_learning_curves_csv(source, output=tmp_path / "one.svg")
    assert not (tmp_path / "one.svg").exists()
