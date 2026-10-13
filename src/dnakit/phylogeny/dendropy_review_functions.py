"""Review wrappers for DendroPy bipartition differences and tree-shape vectors on taxa.

- ``list_missing_bipartitions()`` lists the splits of one tree that the other tree lacks, using
  ``dendropy.calculate.treecompare.find_missing_bipartitions``.
- ``calculate_n_bar_index()``, ``calculate_node_ages()``, ``calculate_node_depths()``,
  ``calculate_coalescence_ages()`` and ``calculate_divergence_times()`` use
  ``dendropy.calculate.treemeasure``. The last four need a branch length on every edge.

Trees are Newick text or a path to a one-tree Newick file. Nothing is drawn, and no network or
external program is used. Packages are imported only when a function runs; a missing package
raises ``BackendUnavailableError`` with the ``pip install "dnakit[bio-dendropy]"`` hint.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.phylogeny.dendropy_functions import (
    _TREE_CONVENTION,
    _flag,
    _number,
    _read_tree,
    _tip_labels,
    _tree_pair,
    _tree_result,
)

_EXTRA = "bio-dendropy"
_MEASURE_CONVENTION = f"{_TREE_CONVENTION} Values are dendropy's; vectors are sorted ascending."


def _call(operation: str, function: Any) -> Any:
    """Run a DendroPy call and report its failures as input errors."""

    try:
        return function()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"dendropy could not {operation}: {exc}") from exc


def _require_lengths(tree: Any) -> None:
    """Every edge except the root edge needs a branch length for time-based measures."""

    for node in tree.preorder_node_iter():
        if node is not tree.seed_node and node.edge.length is None:
            raise ConfigurationError("This measure needs a branch length on every edge.")


def _tree(dendropy: Any, source: object, max_tips: int, max_file_bytes: int) -> Any:
    return _read_tree(dendropy, source, max_tips=max_tips, max_file_bytes=max_file_bytes)


def _values(vector: Sequence[object]) -> list[float | None]:
    return sorted((_number(value) for value in vector), key=lambda v: (v is None, v or 0.0))


def _split_tips(leafset: Sequence[str], all_tips: Sequence[str]) -> list[str]:
    """Report one side of a split: the smaller side, ties broken away from the first tip."""

    inside = sorted(leafset)
    outside = sorted(set(all_tips) - set(inside))
    first = min(all_tips)
    if len(inside) != len(outside):
        return inside if len(inside) < len(outside) else outside
    return outside if first in inside else inside


def list_missing_bipartitions(
    reference_tree: str | Path,
    comparison_tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """List the splits of the reference tree that the comparison tree does not contain."""

    dendropy = require("dendropy", _EXTRA)
    treecompare = require("dendropy.calculate.treecompare", _EXTRA)
    reference, comparison = _tree_pair(
        dendropy, reference_tree, comparison_tree, max_tips=max_tips, max_file_bytes=max_file_bytes
    )
    namespace = reference.taxon_namespace
    all_tips = _tip_labels(reference)
    missing = _call(
        "find missing bipartitions",
        lambda: treecompare.find_missing_bipartitions(reference, comparison),
    )
    splits = sorted(
        (
            _split_tips([taxon.label for taxon in item.leafset_taxa(namespace)], all_tips)
            for item in missing
        ),
        key=lambda tips: (len(tips), tips),
    )
    return _tree_result(
        "list_missing_bipartitions",
        reference,
        {"bipartitions": splits, "missing_count": len(splits)},
        convention=(
            "Each split is listed as the smaller of its two tip sets (ties: the side not "
            "containing the alphabetically first tip). Both trees must have the same tip labels."
        ),
    )


def calculate_n_bar_index(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Average number of nodes above a tip (the N-bar tree-shape statistic)."""

    dendropy = require("dendropy", _EXTRA)
    treemeasure = require("dendropy.calculate.treemeasure", _EXTRA)
    source = _tree(dendropy, tree, max_tips, max_file_bytes)
    value = _call("compute N-bar", lambda: treemeasure.N_bar(source))
    return _tree_result(
        "calculate_n_bar_index",
        source,
        {"n_bar": _number(value)},
        convention=_MEASURE_CONVENTION,
    )


def calculate_node_ages(
    tree: str | Path,
    *,
    internal_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Node ages in backward time (distance from each node to the deepest tip)."""

    dendropy = require("dendropy", _EXTRA)
    treemeasure = require("dendropy.calculate.treemeasure", _EXTRA)
    flag = _flag(internal_only, "internal_only")
    source = _tree(dendropy, tree, max_tips, max_file_bytes)
    _require_lengths(source)
    values = _call(
        "compute node ages",
        lambda: treemeasure.node_ages(source, is_internal_only=flag),
    )
    return _tree_result(
        "calculate_node_ages",
        source,
        {"internal_only": flag, "values": _values(values)},
        convention=_MEASURE_CONVENTION,
    )


def calculate_node_depths(
    tree: str | Path,
    *,
    internal_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Node depths in forward time (distance from the root to each node)."""

    dendropy = require("dendropy", _EXTRA)
    treemeasure = require("dendropy.calculate.treemeasure", _EXTRA)
    flag = _flag(internal_only, "internal_only")
    source = _tree(dendropy, tree, max_tips, max_file_bytes)
    _require_lengths(source)
    values = _call(
        "compute node depths",
        lambda: treemeasure.node_depths(source, is_internal_only=flag),
    )
    return _tree_result(
        "calculate_node_depths",
        source,
        {"internal_only": flag, "values": _values(values)},
        convention=_MEASURE_CONVENTION,
    )


def calculate_coalescence_ages(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Ages of the internal nodes (coalescence events) in backward time."""

    dendropy = require("dendropy", _EXTRA)
    treemeasure = require("dendropy.calculate.treemeasure", _EXTRA)
    source = _tree(dendropy, tree, max_tips, max_file_bytes)
    _require_lengths(source)
    values = _call("compute coalescence ages", lambda: treemeasure.coalescence_ages(source))
    return _tree_result(
        "calculate_coalescence_ages",
        source,
        {"values": _values(values)},
        convention=_MEASURE_CONVENTION,
    )


def calculate_divergence_times(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Times of the internal nodes (splits) in forward time from the root."""

    dendropy = require("dendropy", _EXTRA)
    treemeasure = require("dendropy.calculate.treemeasure", _EXTRA)
    source = _tree(dendropy, tree, max_tips, max_file_bytes)
    _require_lengths(source)
    values = _call("compute divergence times", lambda: treemeasure.divergence_times(source))
    return _tree_result(
        "calculate_divergence_times",
        source,
        {"values": _values(values)},
        convention=_MEASURE_CONVENTION,
    )


__all__ = [
    "calculate_coalescence_ages",
    "calculate_divergence_times",
    "calculate_n_bar_index",
    "calculate_node_ages",
    "calculate_node_depths",
    "list_missing_bipartitions",
]

_register_public_verbose(globals(), tuple(__all__))
