"""Optional PNG output for uniqueness and diversity evaluation."""

from __future__ import annotations

import pytest

from dnakit import DNA
from dnakit.core import DNASet
from dnakit.evaluation import evaluate_diversity, evaluate_uniqueness


def test_collection_plots_are_opt_in_and_keep_scores(tmp_path, monkeypatch) -> None:
    pytest.importorskip("matplotlib.pyplot")
    monkeypatch.chdir(tmp_path)
    records = DNA(["AAAA", "AAAA", "AAAT", "CCCC"])

    uniqueness = evaluate_uniqueness(records, verbose=False)
    diversity = evaluate_diversity(records, verbose=False)
    assert list(tmp_path.iterdir()) == []

    plotted_uniqueness = evaluate_uniqueness(records, plot=True, verbose=False)
    plotted_diversity = evaluate_diversity(records, plot=True, verbose=False)
    assert plotted_uniqueness.metrics == uniqueness.metrics
    assert plotted_diversity.metrics == diversity.metrics
    assert {path.name for path in tmp_path.iterdir()} == {
        "uniqueness_ratio.png",
        "diversity_mean_bar.png",
        "diversity_pairwise_distance_histogram.png",
    }
    for path in tmp_path.iterdir():
        assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_diversity_plot_handles_an_undefined_singleton_score(tmp_path, monkeypatch) -> None:
    pytest.importorskip("matplotlib.pyplot")
    monkeypatch.chdir(tmp_path)

    result = evaluate_diversity(DNA(["AAAA"]), plot=True, verbose=False)

    assert result.score is None
    assert (tmp_path / "diversity_mean_bar.png").is_file()
    assert (tmp_path / "diversity_pairwise_distance_histogram.png").is_file()


def test_uniqueness_plot_handles_an_empty_collection(tmp_path, monkeypatch) -> None:
    pytest.importorskip("matplotlib.pyplot")
    monkeypatch.chdir(tmp_path)

    result = evaluate_uniqueness(DNASet([]), plot=True, verbose=False)

    assert result.score == 1.0
    assert (tmp_path / "uniqueness_ratio.png").is_file()


def test_diversity_histogram_receives_every_unordered_pair(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def capture(mean_distance: float | None, pair_distances: list[float]) -> None:
        captured["mean"] = mean_distance
        captured["distances"] = sorted(pair_distances)

    monkeypatch.setattr("dnakit.evaluation._plots.save_diversity_charts", capture)

    evaluate_diversity(DNA(["AAAA", "AAAT", "CCCC"]), plot=True, verbose=False)

    assert captured["mean"] == pytest.approx(3.0)
    assert captured["distances"] == [1.0, 4.0, 4.0]
