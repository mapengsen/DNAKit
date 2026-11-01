"""Tests for tskit text-tree drawing (batch-0 coverage-gap review)."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

pytest.importorskip("tskit")

import tskit

from dnakit.exceptions import ConfigurationError
from dnakit.population import tskit_batch0_review_functions as review

TABLES = {
    "sequence_length": 100.0,
    "nodes": [
        {"time": 0.0, "sample": True},
        {"time": 0.0, "sample": True},
        {"time": 0.0, "sample": True},
        {"time": 1.0},
        {"time": 2.0},
    ],
    "edges": [
        {"left": 0.0, "right": 100.0, "parent": 3, "child": 0},
        {"left": 0.0, "right": 100.0, "parent": 3, "child": 1},
        {"left": 0.0, "right": 100.0, "parent": 4, "child": 2},
        {"left": 0.0, "right": 100.0, "parent": 4, "child": 3},
    ],
    "sites": [{"position": 10.0, "ancestral_state": "A"}],
    "mutations": [{"site": 0, "node": 0, "derived_state": "G"}],
}


def _reference_drawing(text_format: str) -> str:
    collection = tskit.TableCollection(100.0)
    for time, sample in [(0.0, True), (0.0, True), (0.0, True), (1.0, False), (2.0, False)]:
        collection.nodes.add_row(time=time, flags=tskit.NODE_IS_SAMPLE if sample else 0)
    for child, parent in [(0, 3), (1, 3), (2, 4), (3, 4)]:
        collection.edges.add_row(0.0, 100.0, parent, child)
    collection.sites.add_row(10.0, "A")
    collection.mutations.add_row(0, 0, "G")
    return str(collection.tree_sequence().at_index(0).draw(format=text_format))


def test_public_functions_are_verbose_keyword_only() -> None:
    for name in review.__all__:
        function = getattr(review, name)
        verbose = inspect.signature(function).parameters["verbose"]
        assert verbose.kind is inspect.Parameter.KEYWORD_ONLY
        assert verbose.default is True
        assert function.__doc__


@pytest.mark.parametrize("text_format", ["unicode", "ascii"])
def test_text_drawing_matches_tskit_draw(tmp_path: Path, text_format: str) -> None:
    target = tmp_path / f"tree_{text_format}.txt"
    result = review.draw_tree_text(TABLES, 0, target, text_format=text_format, verbose=False)
    assert result["status"] == "completed"
    assert result["operation"] == "draw_tree_text"
    assert result["image_path"] == str(target.resolve())
    assert target.read_text(encoding="utf-8") == _reference_drawing(text_format)
    assert result["text"] == _reference_drawing(text_format)
    assert result["bytes"] == target.stat().st_size
    assert result["format"] == text_format
    assert result["interval"] == [0.0, 100.0]
    assert result["num_roots"] == 1
    assert result["num_sites"] == 1
    assert result["newick"].startswith("(3:")
    assert result["provenance"]["package"].startswith("tskit ")


def test_default_format_is_unicode(tmp_path: Path) -> None:
    result = review.draw_tree_text(TABLES, 0, tmp_path / "default.txt", verbose=False)
    assert result["format"] == "unicode"
    assert "┏" in result["text"]


def test_checks_format_index_and_destination(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="text_format"):
        review.draw_tree_text(TABLES, 0, tmp_path / "svg.txt", text_format="svg", verbose=False)
    with pytest.raises(ConfigurationError, match="tree_index"):
        review.draw_tree_text(TABLES, 1, tmp_path / "index.txt", verbose=False)
    existing = tmp_path / "exists.txt"
    existing.write_text("keep", encoding="utf-8")
    with pytest.raises(ConfigurationError):
        review.draw_tree_text(TABLES, 0, existing, verbose=False)
    assert existing.read_text(encoding="utf-8") == "keep"
