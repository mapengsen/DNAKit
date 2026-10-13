"""Optional toytree wrappers for tree rooting, editing, distances and taxon-level matrices.

- Tree editing with ``toytree.mod``: ``root_tree_on_midpoint()``,
  ``root_tree_on_balanced_midpoint()``, ``root_tree_on_minimal_ancestor_deviation()``,
  ``unroot_tree()``, ``ladderize_tree()``, ``extract_subtree()``, ``collapse_weak_nodes()``,
  ``remove_unary_nodes()`` and ``scale_tree_to_root_height()`` return the edited tree as Newick.
- Tree comparison with ``toytree.distance``: ``calculate_information_rf_distance()``,
  ``calculate_matching_split_distance()``,
  ``calculate_matching_split_information_distance()``,
  ``calculate_shared_phylogenetic_information_distance()`` and
  ``calculate_mutual_clustering_information_distance()`` compare two trees on the same tips.
- Taxon-level matrices and tables: ``calculate_tree_vcv_matrix()``,
  ``calculate_tree_correlation_matrix()``, ``list_tree_bipartitions()``,
  ``calculate_relative_evolutionary_divergence()`` and ``calculate_equal_splits()``.

Trees are given as Newick text or as a path to a one-tree Newick file. Nothing is drawn.
Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-toytree]"`` hint.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-toytree"
_MAX_FILE_BYTES_LIMIT = 200_000_000
_TREE_CONVENTION = (
    "Trees are read with toytree from Newick text or a one-tree Newick file. Tip labels must be "
    "unique and non-empty. The Newick output is toytree's write() with its default precision."
)


def _int(value: object, name: str, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be an integer in [{low}, {high}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _positive(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a positive number.")
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ConfigurationError(f"{name} must be a positive finite number.")
    return number


def _non_negative(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a non-negative number.")
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ConfigurationError(f"{name} must be a non-negative finite number.")
    return number


def _name(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 200 or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string of at most 200 characters.")
    return value


def _package() -> str:
    return f"toytree=={metadata.version('toytree')}"


def _number(value: object) -> float | None:
    number = float(value)
    return round(number, 12) if math.isfinite(number) else None


def _tree_source(value: object, max_file_bytes: int) -> tuple[str, Any]:
    """Return ("text", Newick string) or ("path", Path) after size checks."""

    max_bytes = _int(max_file_bytes, "max_file_bytes", 1, _MAX_FILE_BYTES_LIMIT)
    if isinstance(value, str) and ";" in value and "\x00" not in value:
        size = len(value.encode("utf-8"))
        if size > max_bytes:
            raise ConfigurationError(
                f"Newick text has {size} bytes; max_file_bytes is {max_bytes}."
            )
        return "text", value.strip()
    if not isinstance(value, (str, Path)) or "\x00" in str(value):
        raise ConfigurationError("tree must be Newick text containing ';' or a path to a file.")
    path = input_path(value)
    size = path.stat().st_size
    if size > max_bytes:
        raise ConfigurationError(f"Tree file has {size} bytes; max_file_bytes is {max_bytes}.")
    return "path", path


def _tip_labels(tree: Any) -> list[str]:
    labels = [str(label) for label in tree.get_tip_labels()]
    if any(not label for label in labels):
        raise ConfigurationError("Every tip needs a non-empty label.")
    if len(set(labels)) != len(labels):
        raise ConfigurationError("Tip labels must be unique.")
    return sorted(labels)


def _read_tree(toytree: Any, source: object, *, max_tips: int, max_file_bytes: int) -> Any:
    _int(max_tips, "max_tips", 2, 100_000)
    kind, value = _tree_source(source, max_file_bytes)
    try:
        tree = toytree.tree(str(value)) if kind == "path" else toytree.tree(value)
    except Exception as exc:
        raise ConfigurationError(f"toytree could not parse the Newick tree: {exc}") from exc
    labels = _tip_labels(tree)
    if len(labels) < 2:
        raise ConfigurationError("A tree needs at least two tips.")
    if len(labels) > max_tips:
        raise ConfigurationError(f"Tree has {len(labels)} tips; max_tips is {max_tips}.")
    return tree


def _tree_pair(
    toytree: Any, source_a: object, source_b: object, *, max_tips: int, max_file_bytes: int
) -> tuple[Any, Any]:
    tree_a = _read_tree(toytree, source_a, max_tips=max_tips, max_file_bytes=max_file_bytes)
    tree_b = _read_tree(toytree, source_b, max_tips=max_tips, max_file_bytes=max_file_bytes)
    if _tip_labels(tree_a) != _tip_labels(tree_b):
        raise ConfigurationError("Both trees must have the same tip labels.")
    return tree_a, tree_b


def _newick_result(
    operation: str, tree: Any, fields: dict[str, object] | None = None
) -> dict[str, object]:
    result: dict[str, object] = {
        "newick": str(tree.write()),
        "tip_count": len(_tip_labels(tree)),
        "tips": _tip_labels(tree),
    }
    result.update(fields or {})
    return {
        "status": "completed",
        "operation": operation,
        "result": result,
        "provenance": {"package": _package(), "convention": _TREE_CONVENTION},
    }


def _distance_result(
    operation: str, tree: Any, value: object, convention: str, fields: dict[str, object]
) -> dict[str, object]:
    return {
        "status": "completed",
        "operation": operation,
        "result": {"tip_count": len(_tip_labels(tree)), "distance": _number(value), **fields},
        "provenance": {"package": _package(), "convention": convention},
    }


def root_tree_on_midpoint(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Root a tree at the midpoint of its longest tip-to-tip path, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    return _newick_result("root_tree_on_midpoint", mod.root_on_midpoint(source))


def root_tree_on_balanced_midpoint(
    tree: str | Path,
    *,
    tolerance: float = 1e-6,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Root a tree at the balanced midpoint, which equalizes the depth of both sides."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    limit = _positive(tolerance, "tolerance")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    rooted = mod.root_on_balanced_midpoint(source, tolerance=limit)
    return _newick_result("root_tree_on_balanced_midpoint", rooted, {"tolerance": limit})


def root_tree_on_minimal_ancestor_deviation(
    tree: str | Path,
    *,
    min_dist: float = 1e-12,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Root a tree on the edge that minimizes ancestor deviation (MAD rooting)."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    limit = _positive(min_dist, "min_dist")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    rooted = mod.root_on_minimal_ancestor_deviation(source, min_dist=limit)
    return _newick_result("root_tree_on_minimal_ancestor_deviation", rooted, {"min_dist": limit})


def unroot_tree(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Remove the root split of a rooted tree and return the unrooted tree as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    return _newick_result("unroot_tree", mod.unroot(source))


def ladderize_tree(
    tree: str | Path,
    *,
    direction: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Order each node's children by clade size; direction=True reverses the order."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    flag = _flag(direction, "direction")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    return _newick_result(
        "ladderize_tree", mod.ladderize(source, direction=flag), {"direction": flag}
    )


def extract_subtree(
    tree: str | Path,
    tips: Sequence[str],
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Return the clade spanned by the given tips as Newick, keeping its stem length."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    if isinstance(tips, str) or not isinstance(tips, Sequence) or len(tips) < 2:
        raise ConfigurationError("tips must be a sequence of at least two tip labels.")
    names = [_name(label, "tip label") for label in tips]
    if len(set(names)) != len(names):
        raise ConfigurationError("tips must not repeat a label.")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    missing = sorted(set(names) - set(_tip_labels(source)))
    if missing:
        raise ConfigurationError(f"Tips are not in the tree: {', '.join(missing)}")
    clade = mod.extract_subtree(source, *names)
    return _newick_result("extract_subtree", clade, {"requested_tips": sorted(names)})


def collapse_weak_nodes(
    tree: str | Path,
    *,
    min_support: float,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Collapse internal nodes whose Newick support label is below min_support."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    support = _non_negative(min_support, "min_support")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    collapsed = mod.collapse_nodes(source, min_support=support)
    return _newick_result("collapse_weak_nodes", collapsed, {"min_support": support})


def remove_unary_nodes(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Remove single-child internal nodes and merge their branch lengths."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    return _newick_result("remove_unary_nodes", mod.remove_unary_nodes(source))


def scale_tree_to_root_height(
    tree: str | Path,
    *,
    height: float = 1.0,
    include_stem: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Rescale all branch lengths so the root-to-tip height equals height."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    target = _positive(height, "height")
    flag = _flag(include_stem, "include_stem")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    scaled = mod.edges_scale_to_root_height(source, treeheight=target, include_stem=flag)
    return _newick_result(
        "scale_tree_to_root_height", scaled, {"height": target, "include_stem": flag}
    )


def _square(frame: Any, labels: list[str]) -> list[list[float | None]]:
    return [[_number(frame.loc[row, column]) for column in labels] for row in labels]


def calculate_tree_vcv_matrix(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Tip variance-covariance matrix: shared root-to-node path lengths, rows sorted by tip."""

    toytree = require("toytree", _EXTRA)
    pcm = require("toytree.pcm", _EXTRA)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    labels = _tip_labels(source)
    frame = pcm.get_vcv_matrix_from_tree(source, df=True)
    return {
        "status": "completed",
        "operation": "calculate_tree_vcv_matrix",
        "result": {
            "tip_count": len(labels),
            "tips": labels,
            "matrix": _square(frame, labels),
        },
        "provenance": {
            "package": _package(),
            "convention": (
                "Diagonal entries are root-to-tip distances; off-diagonal entries are the length "
                "of the path shared by two tips from the root."
            ),
        },
    }


def calculate_tree_correlation_matrix(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Tip correlation matrix derived from the variance-covariance matrix, rows sorted by tip."""

    toytree = require("toytree", _EXTRA)
    pcm = require("toytree.pcm", _EXTRA)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    labels = _tip_labels(source)
    frame = pcm.get_corr_matrix_from_tree(source, df=True)
    return {
        "status": "completed",
        "operation": "calculate_tree_correlation_matrix",
        "result": {
            "tip_count": len(labels),
            "tips": labels,
            "matrix": _square(frame, labels),
        },
        "provenance": {
            "package": _package(),
            "convention": "Each VCV entry divided by the square root of its two tip variances.",
        },
    }


def _split_distance(
    operation: str,
    metric: Any,
    tree_a: str | Path,
    tree_b: str | Path,
    *,
    normalize: bool,
    max_tips: int,
    max_file_bytes: int,
    convention: str,
) -> dict[str, object]:
    toytree = require("toytree", _EXTRA)
    flag = _flag(normalize, "normalize")
    first, second = _tree_pair(
        toytree, tree_a, tree_b, max_tips=max_tips, max_file_bytes=max_file_bytes
    )
    try:
        value = metric(first, second, normalize=flag)
    except ValueError as exc:
        raise ConfigurationError(
            f"toytree could not compute {operation} for these trees: {exc}. Identical split "
            "sets can make the matching-based metric undefined."
        ) from exc
    return _distance_result(operation, first, value, convention, {"normalize": flag})


def calculate_information_rf_distance(
    tree_a: str | Path,
    tree_b: str | Path,
    *,
    normalize: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Information-corrected Robinson-Foulds distance (rfi) between two trees on the same tips."""

    distance = require("toytree.distance", _EXTRA)
    return _split_distance(
        "calculate_information_rf_distance",
        distance.get_treedist_rfi,
        tree_a,
        tree_b,
        normalize=normalize,
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
        convention="toytree get_treedist_rfi: split information corrected Robinson-Foulds.",
    )


def calculate_matching_split_distance(
    tree_a: str | Path,
    tree_b: str | Path,
    *,
    normalize: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Generalized Robinson-Foulds matching split distance between two trees on the same tips."""

    distance = require("toytree.distance", _EXTRA)
    return _split_distance(
        "calculate_matching_split_distance",
        distance.get_treedist_rfg_ms,
        tree_a,
        tree_b,
        normalize=normalize,
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
        convention=(
            "toytree get_treedist_rfg_ms (Bogdanowicz and Giaro matching split). toytree warns "
            "that this metric has no normalization method."
        ),
    )


def calculate_matching_split_information_distance(
    tree_a: str | Path,
    tree_b: str | Path,
    *,
    normalize: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Matching split information distance between two trees on the same tips."""

    distance = require("toytree.distance", _EXTRA)
    return _split_distance(
        "calculate_matching_split_information_distance",
        distance.get_treedist_rfg_msi,
        tree_a,
        tree_b,
        normalize=normalize,
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
        convention="toytree get_treedist_rfg_msi: matching split distance weighted by information.",
    )


def calculate_shared_phylogenetic_information_distance(
    tree_a: str | Path,
    tree_b: str | Path,
    *,
    normalize: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Generalized Robinson-Foulds shared phylogenetic information distance on the same tips."""

    distance = require("toytree.distance", _EXTRA)
    return _split_distance(
        "calculate_shared_phylogenetic_information_distance",
        distance.get_treedist_rfg_spi,
        tree_a,
        tree_b,
        normalize=normalize,
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
        convention="toytree get_treedist_rfg_spi: generalized RF, shared phylogenetic information.",
    )


def calculate_mutual_clustering_information_distance(
    tree_a: str | Path,
    tree_b: str | Path,
    *,
    normalize: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Generalized Robinson-Foulds mutual clustering information distance on the same tips."""

    distance = require("toytree.distance", _EXTRA)
    return _split_distance(
        "calculate_mutual_clustering_information_distance",
        distance.get_treedist_rfg_mci,
        tree_a,
        tree_b,
        normalize=normalize,
        max_tips=max_tips,
        max_file_bytes=max_file_bytes,
        convention="toytree get_treedist_rfg_mci: generalized RF, mutual clustering information.",
    )


def list_tree_bipartitions(
    tree: str | Path,
    *,
    include_singletons: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """List the splits of a tree as pairs of sorted tip sets, in a deterministic order."""

    toytree = require("toytree", _EXTRA)
    enum = require("toytree.enum", _EXTRA)
    flag = _flag(include_singletons, "include_singletons")
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    splits = []
    for side_a, side_b in enum.iter_bipartitions(source, include_singleton_partitions=flag):
        first = sorted(str(label) for label in side_a)
        second = sorted(str(label) for label in side_b)
        splits.append(sorted([first, second]))
    splits.sort()
    return {
        "status": "completed",
        "operation": "list_tree_bipartitions",
        "result": {
            "tip_count": len(_tip_labels(source)),
            "split_count": len(splits),
            "include_singletons": flag,
            "splits": splits,
        },
        "provenance": {
            "package": _package(),
            "convention": (
                "Each split is two sorted tip lists; the side with the smaller first tip comes "
                "first, and splits are sorted."
            ),
        },
    }


def calculate_relative_evolutionary_divergence(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Relative evolutionary divergence (RED) of every node, sorted by node index."""

    toytree = require("toytree", _EXTRA)
    pcm = require("toytree.pcm", _EXTRA)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    names = {int(node.idx): node.name or None for node in source.treenode.traverse()}
    values = pcm.get_relative_evolutionary_divergence(source)
    rows = [
        {
            "node_index": int(index),
            "name": names.get(int(index)),
            "red": _number(values[index]),
        }
        for index in sorted(values, key=int)
    ]
    return {
        "status": "completed",
        "operation": "calculate_relative_evolutionary_divergence",
        "result": {"tip_count": len(_tip_labels(source)), "nodes": rows},
        "provenance": {
            "package": _package(),
            "convention": "toytree get_relative_evolutionary_divergence, keyed by node index.",
        },
    }


def calculate_equal_splits(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Equal splits (ES) metric for every tip, sorted by tip label."""

    toytree = require("toytree", _EXTRA)
    pcm = require("toytree.pcm", _EXTRA)
    source = _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)
    tips = set(_tip_labels(source))
    series = pcm.get_equal_splits(source)
    rows = sorted(
        (
            {"tip": str(label), "equal_splits": _number(value)}
            for label, value in series.items()
            if str(label) in tips
        ),
        key=lambda row: str(row["tip"]),
    )
    return {
        "status": "completed",
        "operation": "calculate_equal_splits",
        "result": {"tip_count": len(rows), "tips": rows},
        "provenance": {
            "package": _package(),
            "convention": "toytree get_equal_splits (Redding and Mooers 2006), one value per tip.",
        },
    }


__all__ = [
    "calculate_equal_splits",
    "calculate_information_rf_distance",
    "calculate_matching_split_distance",
    "calculate_matching_split_information_distance",
    "calculate_mutual_clustering_information_distance",
    "calculate_relative_evolutionary_divergence",
    "calculate_shared_phylogenetic_information_distance",
    "calculate_tree_correlation_matrix",
    "calculate_tree_vcv_matrix",
    "collapse_weak_nodes",
    "extract_subtree",
    "ladderize_tree",
    "list_tree_bipartitions",
    "remove_unary_nodes",
    "root_tree_on_balanced_midpoint",
    "root_tree_on_midpoint",
    "root_tree_on_minimal_ancestor_deviation",
    "scale_tree_to_root_height",
    "unroot_tree",
]

_register_public_verbose(globals(), tuple(__all__))
