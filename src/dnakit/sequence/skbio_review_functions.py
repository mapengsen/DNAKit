"""scikit-bio tree functions added by review group G18o (sequence area).

Tree builders and tree-improvement steps take an aligned DNA mapping (name to aligned sequence,
gaps allowed). They compute a p-distance matrix themselves: mismatches divided by the number of
sites where both sequences carry an unambiguous base. Tree comparison functions take Newick
strings. Each function imports scikit-bio only when it runs.
"""

from __future__ import annotations

import importlib.metadata
import io
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-skibio"
_PACKAGE = "scikit-bio"
_UNAMBIGUOUS = frozenset("ACGT")
_ALIGNED_SYMBOLS = frozenset("ACGTRYSWKMBDHVN-")
_MIN_TAXA = 3
_MAX_TAXA = 500
_MAX_ALIGNMENT_LENGTH = 100_000
_MAX_TREES = 200
_MAX_TREE_TEXT = 1_000_000
_TREE_METRICS = ("cityblock", "euclidean", "sqeuclidean")


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _alignment(sequences: object) -> dict[str, str]:
    """Validate an aligned DNA mapping: unique names, equal lengths, 3 to 500 taxa."""

    if not isinstance(sequences, Mapping):
        raise ConfigurationError("sequences must be a mapping of taxon name to aligned DNA.")
    if not _MIN_TAXA <= len(sequences) <= _MAX_TAXA:
        raise ConfigurationError(f"sequences must hold {_MIN_TAXA} to {_MAX_TAXA} taxa.")
    cleaned: dict[str, str] = {}
    for name, value in sequences.items():
        if not isinstance(name, str) or not name or any(c.isspace() for c in name):
            raise ConfigurationError("taxon names must be non-empty strings without spaces.")
        if not isinstance(value, str) or not value:
            raise ConfigurationError(f"sequence for {name!r} must be a non-empty string.")
        text = value.upper()
        if len(text) > _MAX_ALIGNMENT_LENGTH:
            raise ConfigurationError(
                f"aligned sequences may be at most {_MAX_ALIGNMENT_LENGTH} long."
            )
        if set(text) - _ALIGNED_SYMBOLS:
            raise ConfigurationError(f"sequence for {name!r} contains symbols outside DNA and '-'.")
        cleaned[name] = text
    if len({len(text) for text in cleaned.values()}) != 1:
        raise ConfigurationError("aligned sequences must all have the same length.")
    return cleaned


def _p_distance_matrix(sequences: Mapping[str, str]) -> Any:
    """Return the p-distance matrix of aligned DNA sequences, taxa sorted by name."""

    skbio = require("skbio", _EXTRA)
    names = sorted(sequences)
    size = len(names)
    data = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i + 1, size):
            left, right = sequences[names[i]], sequences[names[j]]
            compared = 0
            differences = 0
            for a, b in zip(left, right, strict=True):
                if a in _UNAMBIGUOUS and b in _UNAMBIGUOUS:
                    compared += 1
                    differences += a != b
            if compared == 0:
                raise ConfigurationError(
                    f"{names[i]!r} and {names[j]!r} share no unambiguous aligned base."
                )
            value = differences / compared
            data[i][j] = data[j][i] = value
    return skbio.DistanceMatrix(data, names)


def _newick(tree: Any) -> str:
    """Newick text with branch lengths rounded to 10 decimals (removes floating-point noise)."""

    for node in tree.postorder(include_self=True):
        if node.length is not None:
            node.length = round(float(node.length), 10)
    return str(tree).strip()


def _tree_from_newick(newick: object, name: str) -> Any:
    tree_module = require("skbio", _EXTRA).TreeNode
    if not isinstance(newick, str) or not newick.strip():
        raise ConfigurationError(f"{name} must be a non-empty Newick string.")
    if len(newick) > _MAX_TREE_TEXT:
        raise ConfigurationError(f"{name} is longer than {_MAX_TREE_TEXT} characters.")
    try:
        return tree_module.read(io.StringIO(newick), format="newick")
    except Exception as exc:
        raise ConfigurationError(f"{name} is not a valid Newick tree: {exc}") from exc


def _tree_list(newicks: object) -> list[Any]:
    if not isinstance(newicks, Sequence) or isinstance(newicks, str):
        raise ConfigurationError("trees must be a sequence of Newick strings.")
    if not 2 <= len(newicks) <= _MAX_TREES:
        raise ConfigurationError(f"trees must contain 2 to {_MAX_TREES} Newick strings.")
    return [_tree_from_newick(text, f"trees[{i}]") for i, text in enumerate(newicks)]


def _ids(values: object, count: int) -> list[str] | None:
    if values is None:
        return None
    if not isinstance(values, Sequence) or isinstance(values, str) or len(values) != count:
        raise ConfigurationError(f"ids must list one label per tree ({count} labels).")
    labels = [str(item) for item in values]
    if any(not label or any(c.isspace() for c in label) for label in labels):
        raise ConfigurationError("ids must be non-empty labels without spaces.")
    if len(set(labels)) != len(labels):
        raise ConfigurationError("ids must be unique.")
    return labels


def _matrix_result(operation: str, matrix: Any, metric: str, method: str) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": operation,
        "metric": metric,
        "ids": [str(item) for item in matrix.ids],
        "matrix": [[float(value) for value in row] for row in matrix.data.tolist()],
        "provenance": _provenance(method=method),
    }


def build_bme_tree(
    sequences: Mapping[str, str],
    *,
    neg_as_zero: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a balanced minimum evolution tree from the p-distances of aligned DNA.

    ``neg_as_zero`` clamps negative branch lengths to zero, as in scikit-bio. Taxa are sorted by
    name. The Newick string is unrooted, with a multifurcating root as scikit-bio writes it.
    """

    alignment = _alignment(sequences)
    flag = _flag(neg_as_zero, "neg_as_zero")
    tree_module = require("skbio.tree", _EXTRA)
    matrix = _p_distance_matrix(alignment)
    tree = tree_module.bme(matrix, neg_as_zero=flag)
    return {
        "status": "completed",
        "operation": "build_bme_tree",
        "taxa": sorted(alignment),
        "distance": "p-distance",
        "neg_as_zero": flag,
        "newick": _newick(tree),
        "provenance": _provenance(method="skbio.tree.bme on p-distances"),
    }


def build_gme_tree(
    sequences: Mapping[str, str],
    *,
    neg_as_zero: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a greedy minimum evolution tree from the p-distances of aligned DNA."""

    alignment = _alignment(sequences)
    flag = _flag(neg_as_zero, "neg_as_zero")
    tree_module = require("skbio.tree", _EXTRA)
    matrix = _p_distance_matrix(alignment)
    tree = tree_module.gme(matrix, neg_as_zero=flag)
    return {
        "status": "completed",
        "operation": "build_gme_tree",
        "taxa": sorted(alignment),
        "distance": "p-distance",
        "neg_as_zero": flag,
        "newick": _newick(tree),
        "provenance": _provenance(method="skbio.tree.gme on p-distances"),
    }


def improve_tree_nni(
    newick: str,
    sequences: Mapping[str, str],
    *,
    balanced: bool = True,
    neg_as_zero: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Apply nearest-neighbour interchanges to a Newick tree that improve its fit to p-distances.

    The tip names of the tree must match the taxa of the aligned DNA exactly. ``balanced`` picks
    the balanced minimum-evolution criterion of scikit-bio instead of the ordinary one.
    """

    alignment = _alignment(sequences)
    balance = _flag(balanced, "balanced")
    flag = _flag(neg_as_zero, "neg_as_zero")
    tree = _tree_from_newick(newick, "newick")
    tips = sorted(tip.name for tip in tree.tips())
    if tips != sorted(alignment):
        raise ConfigurationError("tree tip names must match the aligned sequence names exactly.")
    if len(tips) < 4:
        raise ConfigurationError(
            "NNI needs at least four tips: a three-tip tree has no internal edge."
        )
    tree_module = require("skbio.tree", _EXTRA)
    matrix = _p_distance_matrix(alignment)
    try:
        improved = tree_module.nni(tree, matrix, balanced=balance, neg_as_zero=flag)
    except ValueError as exc:
        raise ConfigurationError(f"skbio could not improve the tree: {exc}") from exc
    return {
        "status": "completed",
        "operation": "improve_tree_nni",
        "tip_count": len(tips),
        "balanced": balance,
        "neg_as_zero": flag,
        "distance": "p-distance",
        "newick": _newick(improved),
        "provenance": _provenance(method="skbio.tree.nni on p-distances"),
    }


def tree_path_distances(
    newicks: Sequence[str],
    *,
    ids: Sequence[str] | None = None,
    metric: str = "euclidean",
    use_length: bool = True,
    shared_by_all: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare trees by their tip-to-tip path lengths (scikit-bio ``path_dists``).

    For each pair of trees, the metric compares the path lengths between every pair of shared
    tips. ``use_length`` counts branch lengths; otherwise each edge counts as 1. Sampling is not
    offered, so the result is exact and deterministic.
    """

    trees = _tree_list(newicks)
    labels = _ids(ids, len(trees))
    chosen = _choice(metric, _TREE_METRICS, "metric")
    use = _flag(use_length, "use_length")
    shared = _flag(shared_by_all, "shared_by_all")
    tree_module = require("skbio.tree", _EXTRA)
    matrix = tree_module.path_dists(
        trees, ids=labels, shared_by_all=shared, metric=chosen, use_length=use
    )
    result = _matrix_result(
        "tree_path_distances", matrix, chosen, "skbio.tree.path_dists (exact, no sampling)"
    )
    result["use_length"] = use
    result["shared_by_all"] = shared
    return result


def weighted_robinson_foulds_distances(
    newicks: Sequence[str],
    *,
    ids: Sequence[str] | None = None,
    metric: str = "cityblock",
    rooted: bool = False,
    include_tips: bool = True,
    shared_by_all: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare trees by weighted Robinson-Foulds distances (scikit-bio ``wrf_dists``).

    Each split carries its branch length. ``metric`` combines the length differences of shared
    splits with the lengths of unshared ones. ``include_tips`` counts the terminal branches.
    """

    trees = _tree_list(newicks)
    labels = _ids(ids, len(trees))
    chosen = _choice(metric, _TREE_METRICS, "metric")
    root_flag = _flag(rooted, "rooted")
    tips = _flag(include_tips, "include_tips")
    shared = _flag(shared_by_all, "shared_by_all")
    tree_module = require("skbio.tree", _EXTRA)
    matrix = tree_module.wrf_dists(
        trees,
        ids=labels,
        shared_by_all=shared,
        metric=chosen,
        rooted=root_flag,
        include_tips=tips,
    )
    result = _matrix_result(
        "weighted_robinson_foulds_distances",
        matrix,
        chosen,
        "skbio.tree.wrf_dists (weighted Robinson-Foulds)",
    )
    result["rooted"] = root_flag
    result["include_tips"] = tips
    result["shared_by_all"] = shared
    return result


__all__ = [
    "build_bme_tree",
    "build_gme_tree",
    "improve_tree_nni",
    "tree_path_distances",
    "weighted_robinson_foulds_distances",
]

_register_public_verbose(globals(), tuple(__all__))
