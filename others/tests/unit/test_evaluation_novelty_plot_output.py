"""Optional PNG output for novelty evaluation."""

from __future__ import annotations

import pytest

from dnakit import DNA
from dnakit.core import DNASet
from dnakit.evaluation import (
    ReferenceSearchConfig,
    create_reference_library,
    evaluate_novelty,
)


def _reference():
    return create_reference_library(
        DNA([{"id": "ref_a", "sequence": "AAAA"}, {"id": "ref_c", "sequence": "CCCC"}]),
        name="plot-reference",
        version="1",
        source="local:test",
        verbose=False,
    )


def test_novelty_plot_is_opt_in_and_saves_three_pngs(tmp_path, monkeypatch) -> None:
    pytest.importorskip("matplotlib.pyplot")
    monkeypatch.chdir(tmp_path)
    queries = DNA(
        [
            {"id": "copy", "sequence": "AAAA", "metadata": {"batch": "Copies"}},
            {"id": "near", "sequence": "AAAT", "metadata": {"batch": "Copies"}},
            {"id": "far", "sequence": "GGGG", "metadata": {"batch": "Other"}},
        ]
    )
    reference = _reference()

    plain = evaluate_novelty(
        queries, reference, config=ReferenceSearchConfig(method="identity"), verbose=False
    )
    assert list(tmp_path.iterdir()) == []

    plotted = evaluate_novelty(
        queries,
        reference,
        config=ReferenceSearchConfig(method="identity"),
        plot=True,
        verbose=False,
    )
    assert plotted.metrics == plain.metrics
    assert {path.name for path in tmp_path.iterdir()} == {
        "novelty_batch_summary.png",
        "novelty_distribution.png",
        "novelty_nearest_levenshtein.png",
    }
    for path in tmp_path.iterdir():
        assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_similarity_plot_receives_nearest_reference_edit_distances(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def capture(result, nearest_distances, *, batch_labels) -> None:
        captured["score"] = result.score
        captured["distances"] = nearest_distances
        captured["labels"] = batch_labels

    monkeypatch.setattr("dnakit.evaluation._plots.save_novelty_charts", capture)
    queries = DNA(
        [
            {"id": "copy", "sequence": "AAAA", "metadata": {"batch": "A"}},
            {"id": "near", "sequence": "AAAT", "metadata": {"batch": "A"}},
            {"id": "far", "sequence": "GGGG", "metadata": {"batch": "B"}},
        ]
    )

    evaluate_novelty(
        queries,
        _reference(),
        config=ReferenceSearchConfig(method="identity"),
        plot=True,
        verbose=False,
    )

    assert captured["score"] == pytest.approx(5 / 12)
    assert captured["distances"] == [0.0, 1.0, 4.0]
    assert captured["labels"] == ("A", "A", "B")


def test_levenshtein_plot_reuses_its_reported_distances(monkeypatch) -> None:
    captured: dict[str, object] = {}

    def capture(result, nearest_distances, *, batch_labels) -> None:
        captured["score"] = result.score
        captured["distances"] = nearest_distances
        captured["labels"] = batch_labels

    monkeypatch.setattr("dnakit.evaluation._plots.save_novelty_charts", capture)
    queries = DNA(["AAAA", "AAAT", "GGGG"])

    evaluate_novelty(
        queries,
        _reference(),
        config=ReferenceSearchConfig(novelty_calculation="levenshtein"),
        plot=True,
        verbose=False,
    )

    assert captured["score"] == pytest.approx(5 / 3)
    assert captured["distances"] == [0.0, 1.0, 4.0]
    assert captured["labels"] == ("Evaluated set",) * 3


def test_levenshtein_mode_saves_three_pngs(tmp_path, monkeypatch) -> None:
    pytest.importorskip("matplotlib.pyplot")
    monkeypatch.chdir(tmp_path)

    result = evaluate_novelty(
        DNA(["AAAA", "AAAT", "GGGG"]),
        _reference(),
        config=ReferenceSearchConfig(novelty_calculation="levenshtein"),
        plot=True,
        verbose=False,
    )

    assert result.score == pytest.approx(5 / 3)
    assert len(list(tmp_path.glob("novelty_*.png"))) == 3


def test_novelty_plot_handles_empty_queries_and_validates_flag(tmp_path, monkeypatch) -> None:
    pytest.importorskip("matplotlib.pyplot")
    monkeypatch.chdir(tmp_path)
    reference = _reference()

    with pytest.raises(TypeError, match="plot must be a boolean"):
        evaluate_novelty(DNA(["AAAA"]), reference, plot=1, verbose=False)  # type: ignore[arg-type]

    result = evaluate_novelty(DNASet([]), reference, plot=True, verbose=False)

    assert result.score is None
    assert len(list(tmp_path.glob("novelty_*.png"))) == 3
