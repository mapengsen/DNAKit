"""Review wrappers for toytree node queries, tree enumeration, tree editing, consensus and tables.

- Node queries with ``toytree.distance``: ``list_tree_node_path()``,
  ``calculate_tree_node_distance()``, ``calculate_all_node_distance_matrix()``,
  ``calculate_internal_node_distance_matrix()``, ``calculate_descendant_distances()``,
  ``find_tree_farthest_node()`` and ``calculate_farthest_node_distance()``.
- Enumeration with ``toytree.enum``: ``list_tree_quadripartitions()`` and ``list_tree_quartets()``.
- Consensus with ``toytree.infer``: ``map_consensus_tree_features()`` maps split support and
  distance statistics from a list of Newick trees onto a target tree.
- Editing with ``toytree.mod``: ``extend_tree_tips_to_align()``, ``set_tree_node_heights()``,
  ``rotate_tree_node()``, ``drop_tree_tips()``, ``bisect_tree()``, ``resolve_tree_polytomies()``,
  ``add_tree_internal_node()``, ``add_tree_child_node()``, ``add_tree_sister_node()``,
  ``add_tree_internal_node_and_child()``, ``add_tree_internal_node_and_subtree()``,
  ``remove_tree_nodes()`` and ``merge_tree_nodes()`` return the edited tree as Newick.
- Networks and tables: ``parse_tree_network()`` uses ``toytree.network``;
  ``calculate_tip_distance_matrix_from_vcv()`` and ``calculate_tip_diversification_rates()`` use
  ``toytree.pcm``.

Trees are Newick text or a path to a one-tree Newick file. A node query is a node name or an
integer node index. Nothing is drawn. Packages are imported only when a function runs; a missing
package raises ``BackendUnavailableError`` with the ``pip install "dnakit[bio-toytree]"`` hint.
"""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.phylogeny.toytree_functions import (
    _TREE_CONVENTION,
    _flag,
    _int,
    _name,
    _newick_result,
    _non_negative,
    _number,
    _package,
    _positive,
    _read_tree,
    _tip_labels,
    _tree_source,
)

_EXTRA = "bio-toytree"
_CONVENTION = (
    f"{_TREE_CONVENTION} Node queries are node names or integer node indexes; "
    "matrices are ordered by node index."
)
_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]{0,49}\Z")


def _query(value: object, name: str) -> str | int:
    """Return a node name or a non-negative integer node index."""

    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise ConfigurationError(f"{name} must be a node name or an integer node index.")
    if isinstance(value, int):
        if value < 0:
            raise ConfigurationError(f"{name} index must be non-negative.")
        return value
    if not value.strip() or len(value) > 200 or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty node name of at most 200 characters.")
    return value


def _optional_query(value: object, name: str) -> str | int | None:
    return None if value is None else _query(value, name)


def _optional_name(value: object, name: str) -> str | None:
    return None if value is None else _name(value, name)


def _optional_distance(value: object, name: str) -> float | None:
    return None if value is None else _non_negative(value, name)


def _call(operation: str, function: Callable[[], Any]) -> Any:
    """Run a toytree call and report its failures as input errors."""

    try:
        return function()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise ConfigurationError(f"toytree could not {operation}: {exc}") from exc


def _result(
    operation: str, fields: Mapping[str, object], convention: str = _CONVENTION
) -> dict[str, object]:
    return {
        "status": "completed",
        "operation": operation,
        "result": dict(fields),
        "provenance": {"package": _package(), "convention": convention},
    }


_EDIT_CONVENTION = (
    f"{_TREE_CONVENTION} toytree's write() omits internal node names, so a name given to an "
    "internal node is echoed in the result fields but not in the Newick text."
)


def _edited(
    operation: str, tree: Any, fields: Mapping[str, object] | None = None
) -> dict[str, object]:
    """Newick result for an edited tree, with the write() naming limitation in its provenance."""

    result = _newick_result(operation, tree, dict(fields or {}))
    result["provenance"] = {"package": _package(), "convention": _EDIT_CONVENTION}
    return result


def _node_record(node: Any) -> dict[str, object]:
    return {
        "idx": int(node.idx),
        "name": str(node.name) if node.name else None,
        "is_tip": bool(node.is_leaf()),
    }


def _labels(nodes: Sequence[Any]) -> list[str]:
    """One label per node: its name, or node_<idx> when an internal node has no name."""

    labels = [str(node.name) if node.name else f"node_{node.idx}" for node in nodes]
    if len(set(labels)) != len(labels):
        raise ConfigurationError("Node labels must be unique; rename duplicate internal nodes.")
    return labels


def _matrix(values: Any) -> list[list[float | None]]:
    return [[_number(value) for value in row] for row in values]


def _tree(toytree: Any, tree: object, max_tips: int, max_file_bytes: int) -> Any:
    return _read_tree(toytree, tree, max_tips=max_tips, max_file_bytes=max_file_bytes)


def _tip_set(tree: Any) -> set[str]:
    return set(_tip_labels(tree))


def _require_known_tips(tree: Any, tips: Sequence[str]) -> None:
    known = _tip_set(tree)
    unknown = sorted(set(tips) - known)
    if unknown:
        raise ConfigurationError(f"Unknown tip labels: {unknown}.")


def _check_remaining_tips(tree: Any, removed: int) -> None:
    remaining = len(_tip_labels(tree)) - removed
    if remaining < 2:
        raise ConfigurationError("The edited tree must keep at least two tips.")


def list_tree_node_path(
    tree: str | Path,
    start: str | int,
    end: str | int,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """List the nodes on the path between two nodes, including both ends."""

    toytree = require("toytree", _EXTRA)
    distance = require("toytree.distance", _EXTRA)
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    first, last = _query(start, "start"), _query(end, "end")
    path = _call("find the node path", lambda: distance.get_node_path(source, first, last))
    return _result(
        "list_tree_node_path",
        {
            "start": first,
            "end": last,
            "nodes": [_node_record(node) for node in path],
            "edge_count": len(path) - 1,
        },
    )


def calculate_tree_node_distance(
    tree: str | Path,
    start: str | int,
    end: str | int,
    *,
    topology_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Patristic (branch length) or topological distance between two nodes of a tree."""

    toytree = require("toytree", _EXTRA)
    distance = require("toytree.distance", _EXTRA)
    flag = _flag(topology_only, "topology_only")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    first, last = _query(start, "start"), _query(end, "end")
    value = _call(
        "compute the node distance",
        lambda: distance.get_node_distance(source, first, last, topology_only=flag),
    )
    return _result(
        "calculate_tree_node_distance",
        {"start": first, "end": last, "topology_only": flag, "distance": _number(value)},
    )


def calculate_all_node_distance_matrix(
    tree: str | Path,
    *,
    topology_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Pairwise distances between all tips and internal nodes, ordered by node index."""

    toytree = require("toytree", _EXTRA)
    distance = require("toytree.distance", _EXTRA)
    flag = _flag(topology_only, "topology_only")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    nodes = list(source.get_nodes())
    labels = _labels(nodes)
    values = _call(
        "compute the node distance matrix",
        lambda: distance.get_node_distance_matrix(source, topology_only=flag),
    )
    return _result(
        "calculate_all_node_distance_matrix",
        {"labels": labels, "matrix": _matrix(values), "topology_only": flag},
    )


def calculate_internal_node_distance_matrix(
    tree: str | Path,
    *,
    topology_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Pairwise distances between the internal nodes (including the root), by node index."""

    toytree = require("toytree", _EXTRA)
    distance = require("toytree.distance", _EXTRA)
    flag = _flag(topology_only, "topology_only")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    nodes = [node for node in source.get_nodes() if not node.is_leaf()]
    labels = _labels(nodes)
    values = _call(
        "compute the internal node distance matrix",
        lambda: distance.get_internal_node_distance_matrix(source, topology_only=flag),
    )
    return _result(
        "calculate_internal_node_distance_matrix",
        {"labels": labels, "matrix": _matrix(values), "topology_only": flag},
    )


def calculate_descendant_distances(
    tree: str | Path,
    node: str | int | None = None,
    *,
    topology_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Distance from a node (default: the root) to each of its descendants, sorted by index."""

    toytree = require("toytree", _EXTRA)
    distance = require("toytree.distance", _EXTRA)
    flag = _flag(topology_only, "topology_only")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _optional_query(node, "node")
    mapping = _call(
        "compute descendant distances",
        lambda: distance.get_descendant_dists(source, query, topology_only=flag),
    )
    rows = [
        {**_node_record(item), "distance": _number(value)}
        for item, value in sorted(mapping.items(), key=lambda pair: int(pair[0].idx))
    ]
    return _result(
        "calculate_descendant_distances",
        {"node": query, "descendants": rows, "topology_only": flag},
    )


def find_tree_farthest_node(
    tree: str | Path,
    node: str | int,
    *,
    topology_only: bool = False,
    descendants_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Find the node farthest from a given node, with its distance."""

    toytree = require("toytree", _EXTRA)
    distance = require("toytree.distance", _EXTRA)
    topo = _flag(topology_only, "topology_only")
    below = _flag(descendants_only, "descendants_only")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _query(node, "node")
    far = _call(
        "find the farthest node",
        lambda: distance.get_farthest_node(
            source, query, topology_only=topo, descendants_only=below
        ),
    )
    value = _call(
        "compute the farthest distance",
        lambda: distance.get_farthest_node_distance(
            source, query, topology_only=topo, descendants_only=below
        ),
    )
    return _result(
        "find_tree_farthest_node",
        {
            "node": query,
            "farthest": _node_record(far),
            "distance": _number(value),
            "topology_only": topo,
            "descendants_only": below,
        },
    )


def calculate_farthest_node_distance(
    tree: str | Path,
    node: str | int | None = None,
    *,
    topology_only: bool = False,
    descendants_only: bool = False,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Distance from a node (default: the root) to the farthest node."""

    toytree = require("toytree", _EXTRA)
    distance = require("toytree.distance", _EXTRA)
    topo = _flag(topology_only, "topology_only")
    below = _flag(descendants_only, "descendants_only")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _optional_query(node, "node")
    value = _call(
        "compute the farthest distance",
        lambda: distance.get_farthest_node_distance(
            source, query, topology_only=topo, descendants_only=below
        ),
    )
    return _result(
        "calculate_farthest_node_distance",
        {
            "node": query,
            "distance": _number(value),
            "topology_only": topo,
            "descendants_only": below,
        },
    )


def list_tree_quadripartitions(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Quadripartitions: the four tip groups around each internal edge, sorted by content."""

    toytree = require("toytree", _EXTRA)
    enum = require("toytree.enum", _EXTRA)
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    groups = _call("enumerate quadripartitions", lambda: list(enum.iter_quadripartitions(source)))
    rows = sorted(
        ([sorted(group) for pair in item for group in pair] for item in groups),
        key=lambda row: [tuple(group) for group in row],
    )
    return _result("list_tree_quadripartitions", {"quadripartitions": rows, "count": len(rows)})


def list_tree_quartets(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Quartets induced by the internal edges: two tip groups on each side of the edge."""

    toytree = require("toytree", _EXTRA)
    enum = require("toytree.enum", _EXTRA)
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    items = _call("enumerate quartets", lambda: list(enum.iter_quartets(source)))
    rows = sorted(
        ({"left": sorted(left), "right": sorted(right)} for left, right in items),
        key=lambda row: (row["left"], row["right"]),
    )
    return _result("list_tree_quartets", {"quartets": rows, "count": len(rows)})


def map_consensus_tree_features(
    tree: str | Path,
    trees: Sequence[str | Path],
    *,
    features: Sequence[str] = ("support",),
    max_trees: int = 1_000,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Map split support and branch statistics from a tree list onto a target tree."""

    toytree = require("toytree", _EXTRA)
    infer = require("toytree.infer", _EXTRA)
    limit = _int(max_trees, "max_trees", 1, 100_000)
    if isinstance(trees, (str, bytes)) or not isinstance(trees, Sequence):
        raise ConfigurationError("trees must be a list of Newick trees or tree file paths.")
    if not 1 <= len(trees) <= limit:
        raise ConfigurationError(f"trees must contain between 1 and {limit} trees.")
    if isinstance(features, (str, bytes)) or not features:
        raise ConfigurationError("features must be a non-empty list of feature names.")
    names = [_name(item, "feature") for item in features]
    if len(set(names)) != len(names):
        raise ConfigurationError("features must be unique.")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    tips = _tip_set(source)
    members = []
    for item in trees:
        member = _tree(toytree, item, max_tips, max_file_bytes)
        if _tip_set(member) != tips:
            raise ConfigurationError(
                "Every tree in the list must have the same tips as the target."
            )
        members.append(member)
    multitree = _call("build the tree list", lambda: toytree.mtree(members))
    mapped = _call(
        "map consensus features",
        lambda: infer.get_consensus_features(source, multitree, features=names),
    )
    data = mapped.get_node_data()
    for feature in names:
        if feature not in data.columns:
            raise ConfigurationError(f"Unknown consensus feature: {feature!r}.")
    values = {
        int(index): {feature: _number(row[feature]) for feature in names}
        for index, row in zip(data["idx"], data.to_dict("records"), strict=True)
    }
    nodes = [
        {
            **_node_record(node),
            "tips": sorted(node.get_leaf_names()),
            "values": values[int(node.idx)],
        }
        for node in mapped.get_nodes()
    ]
    return _edited(
        "map_consensus_tree_features",
        mapped,
        {"features": names, "tree_count": len(members), "nodes": nodes},
    )


def extend_tree_tips_to_align(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Extend every tip branch so that all tips end at the same depth, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    aligned = _call("extend tips", lambda: mod.edges_extend_tips_to_align(source))
    return _edited("extend_tree_tips_to_align", aligned)


def set_tree_node_heights(
    tree: str | Path,
    heights: Mapping[str, float],
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Set the height (distance from the root) of named or indexed nodes, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    if not isinstance(heights, Mapping) or not heights or len(heights) > 10_000:
        raise ConfigurationError("heights must be a non-empty mapping of node to height.")
    data = {
        _query(key, "height node"): _non_negative(value, f"height for {key!r}")
        for key, value in heights.items()
    }
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    edited = _call("set node heights", lambda: mod.edges_set_node_heights(source, data))
    rows = [
        {"node": key, "height": value}
        for key, value in sorted(data.items(), key=lambda pair: str(pair[0]))
    ]
    return _edited("set_tree_node_heights", edited, {"heights": rows})


def rotate_tree_node(
    tree: str | Path,
    node: str | int,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Reverse the child order of one node, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _query(node, "node")
    rotated = _call("rotate the node", lambda: mod.rotate_node(source, query))
    return _edited("rotate_tree_node", rotated, {"node": query})


def drop_tree_tips(
    tree: str | Path,
    tips: Sequence[str],
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Remove the named tips from a tree and return the remaining tree as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    if isinstance(tips, (str, bytes)) or not isinstance(tips, Sequence) or not tips:
        raise ConfigurationError("tips must be a non-empty list of tip labels.")
    names = [_name(item, "tip") for item in tips]
    if len(set(names)) != len(names):
        raise ConfigurationError("tips must be unique.")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    _require_known_tips(source, names)
    _check_remaining_tips(source, len(names))
    dropped = _call("drop the tips", lambda: mod.drop_tips(source, *names))
    return _edited("drop_tree_tips", dropped, {"dropped": sorted(names)})


def bisect_tree(
    tree: str | Path,
    node: str | int,
    *,
    dist_partition: float = 1.0,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Split a tree on the edge above a node into the clade and the remaining tree."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    limit = _positive(dist_partition, "dist_partition")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _query(node, "node")
    clade, remainder = _call(
        "bisect the tree",
        lambda: mod.bisect(source, query, reroot=False, dist_partition=limit),
    )
    return _result(
        "bisect_tree",
        {
            "node": query,
            "dist_partition": limit,
            "subtree_newick": str(clade.write()),
            "subtree_tips": _tip_labels(clade),
            "remainder_newick": str(remainder.write()),
            "remainder_tips": _tip_labels(remainder),
        },
    )


def resolve_tree_polytomies(
    tree: str | Path,
    *,
    node: str | int | None = None,
    seed: int,
    dist: float = 1e-6,
    recursive: bool = True,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Resolve multifurcations at random with a fixed seed, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**32:
        raise ConfigurationError("seed must be an integer in [0, 2**32).")
    length = _positive(dist, "dist")
    flag = _flag(recursive, "recursive")
    query = _optional_query(node, "node")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    args = () if query is None else (query,)
    resolved = _call(
        "resolve polytomies",
        lambda: mod.resolve_polytomies(source, *args, dist=length, recursive=flag, seed=seed),
    )
    return _edited(
        "resolve_tree_polytomies",
        resolved,
        {"node": query, "seed": seed, "dist": length, "recursive": flag},
    )


def add_tree_internal_node(
    tree: str | Path,
    node: str | int,
    *,
    name: str | None = None,
    dist: float | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Insert an internal node on the edge above a node, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    label = _optional_name(name, "name")
    length = _optional_distance(dist, "dist")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _query(node, "node")
    edited = _call(
        "add the internal node",
        lambda: mod.add_internal_node(source, query, dist=length, name=label),
    )
    return _edited("add_tree_internal_node", edited, {"node": query, "name": label, "dist": length})


def add_tree_child_node(
    tree: str | Path,
    node: str | int,
    *,
    name: str | None = None,
    dist: float | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Add a new tip as a child of a node, returned as Newick.

    Given a tip, the tip becomes an unnamed internal node with the new child as its only child,
    so the original tip label is no longer a tip (toytree behaviour).
    """

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    label = _optional_name(name, "name")
    length = _optional_distance(dist, "dist")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _query(node, "node")
    edited = _call(
        "add the child node",
        lambda: mod.add_child_node(source, query, name=label, dist=length),
    )
    return _edited("add_tree_child_node", edited, {"node": query, "name": label, "dist": length})


def add_tree_sister_node(
    tree: str | Path,
    node: str | int,
    *,
    name: str | None = None,
    dist: float | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Add a new tip as a sister of a node, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    label = _optional_name(name, "name")
    length = _optional_distance(dist, "dist")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _query(node, "node")
    edited = _call(
        "add the sister node",
        lambda: mod.add_sister_node(source, query, name=label, dist=length),
    )
    return _edited("add_tree_sister_node", edited, {"node": query, "name": label, "dist": length})


def add_tree_internal_node_and_child(
    tree: str | Path,
    node: str | int,
    *,
    name: str | None = None,
    dist: float | None = None,
    parent_name: str | None = None,
    parent_dist: float | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Split the edge above a node with a new internal node that gets a new tip child."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    label = _optional_name(name, "name")
    length = _optional_distance(dist, "dist")
    parent = _optional_name(parent_name, "parent_name")
    parent_length = _optional_distance(parent_dist, "parent_dist")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    query = _query(node, "node")
    edited = _call(
        "add the internal node and child",
        lambda: mod.add_internal_node_and_child(
            source,
            query,
            name=label,
            dist=length,
            parent_dist=parent_length,
            parent_name=parent,
        ),
    )
    return _edited(
        "add_tree_internal_node_and_child",
        edited,
        {
            "node": query,
            "name": label,
            "dist": length,
            "parent_name": parent,
            "parent_dist": parent_length,
        },
    )


def add_tree_internal_node_and_subtree(
    tree: str | Path,
    node: str | int,
    subtree: str | Path,
    *,
    subtree_stem_dist: float | None = None,
    subtree_rescale: bool = False,
    parent_name: str | None = None,
    parent_dist: float | None = None,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Graft a Newick subtree onto the edge above a node through a new internal node."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    stem = _optional_distance(subtree_stem_dist, "subtree_stem_dist")
    rescale = _flag(subtree_rescale, "subtree_rescale")
    parent = _optional_name(parent_name, "parent_name")
    parent_length = _optional_distance(parent_dist, "parent_dist")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    graft = _tree(toytree, subtree, max_tips, max_file_bytes)
    if _tip_set(source) & _tip_set(graft):
        raise ConfigurationError("The subtree tips must not repeat tips of the tree.")
    query = _query(node, "node")
    edited = _call(
        "add the subtree",
        lambda: mod.add_internal_node_and_subtree(
            source,
            query,
            subtree=graft,
            subtree_stem_dist=stem,
            subtree_rescale=rescale,
            parent_dist=parent_length,
            parent_name=parent,
        ),
    )
    return _edited(
        "add_tree_internal_node_and_subtree",
        edited,
        {
            "node": query,
            "subtree_tips": _tip_labels(graft),
            "subtree_stem_dist": stem,
            "subtree_rescale": rescale,
            "parent_name": parent,
            "parent_dist": parent_length,
        },
    )


def _removed_tip_names(tree: Any, queries: Sequence[str | int]) -> set[str]:
    """Names of the tips that a list of node queries would remove."""

    names = set()
    for query in queries:
        for node in tree.get_nodes(query):
            if node.is_leaf():
                names.add(str(node.name))
    return names


def remove_tree_nodes(
    tree: str | Path,
    nodes: Sequence[str | int],
    *,
    preserve_dists: bool = True,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Remove named or indexed nodes, keeping their children, returned as Newick."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    flag = _flag(preserve_dists, "preserve_dists")
    if isinstance(nodes, (str, bytes)) or not isinstance(nodes, Sequence) or not nodes:
        raise ConfigurationError("nodes must be a non-empty list of node queries.")
    queries = [_query(item, "node") for item in nodes]
    if len(queries) > 10_000:
        raise ConfigurationError("nodes may contain at most 10000 queries.")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    removed_tips = _call("find the nodes", lambda: _removed_tip_names(source, queries))
    _check_remaining_tips(source, len(removed_tips))
    edited = _call(
        "remove the nodes",
        lambda: mod.remove_nodes(source, *queries, preserve_dists=flag),
    )
    return _edited(
        "remove_tree_nodes",
        edited,
        {"removed": [str(item) for item in queries], "preserve_dists": flag},
    )


def merge_tree_nodes(
    tree: str | Path,
    feature: str,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Merge clades whose tips share one node feature, keeping one representative per clade."""

    toytree = require("toytree", _EXTRA)
    mod = require("toytree.mod", _EXTRA)
    if not isinstance(feature, str) or not _IDENTIFIER.fullmatch(feature):
        raise ConfigurationError("feature must be an identifier naming a node attribute.")
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    before = len(source.get_nodes())
    merged = _call(
        "merge the nodes",
        lambda: mod.merge_nodes(source, merge_method=feature),
    )
    return _edited(
        "merge_tree_nodes",
        merged,
        {"feature": feature, "removed_node_count": before - len(merged.get_nodes())},
    )


def parse_tree_network(
    network: str | Path,
    *,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Parse an extended Newick network with hybrid nodes into a major tree and admixture edges."""

    network_module = require("toytree.network", _EXTRA)
    kind, value = _tree_source(network, max_file_bytes)
    argument = value if kind == "text" else str(value)
    major, admix = _call(
        "parse the network",
        lambda: network_module.parse_network_to_tree_and_admix(argument),
    )
    _tip_labels(major)
    edges = []
    for destination, source_tips, _ignored, _extra, gamma_text in admix:
        try:
            gamma = float(gamma_text)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(f"toytree gave an unreadable gamma: {gamma_text!r}") from exc
        edges.append(
            {
                "destination_tips": sorted(destination),
                "source_tips": sorted(source_tips),
                "gamma": _number(gamma),
            }
        )
    return _edited(
        "parse_tree_network",
        major,
        {"admixture_edges": edges, "admixture_count": len(edges)},
    )


def calculate_tip_distance_matrix_from_vcv(
    vcv: Sequence[Sequence[float]],
    labels: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, object]:
    """Euclidean distances between tips computed from a variance-covariance matrix."""

    pcm = require("toytree.pcm", _EXTRA)
    numpy = require("numpy", _EXTRA)
    if isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence):
        raise ConfigurationError("labels must be a list of tip labels.")
    names = [_name(item, "label") for item in labels]
    if len(set(names)) != len(names) or len(names) < 2 or len(names) > 2_000:
        raise ConfigurationError("labels must be 2 to 2000 unique tip labels.")
    try:
        array = numpy.asarray(vcv, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("vcv must be a rectangular numeric matrix.") from exc
    size = len(names)
    if array.shape != (size, size) or not numpy.isfinite(array).all():
        raise ConfigurationError(f"vcv must be a finite {size} by {size} matrix.")
    if not numpy.allclose(array, array.T, rtol=0, atol=1e-9):
        raise ConfigurationError("vcv must be symmetric.")
    distances = _call(
        "compute distances from the VCV matrix",
        lambda: pcm.get_distance_matrix_from_vcv_matrix(array),
    )
    return _result(
        "calculate_tip_distance_matrix_from_vcv",
        {"labels": names, "matrix": _matrix(distances)},
        convention=(
            "Entry (i, j) is V[i,i] + V[j,j] - 2 V[i,j] exactly as toytree.pcm computes it. "
            "That is the squared Euclidean distance between tips: upstream does not take a square "
            "root, although its docstring says Euclidean."
        ),
    )


def calculate_tip_diversification_rates(
    tree: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Tip-level diversification rates sensu Jetz 2012 for every tip of a tree."""

    toytree = require("toytree", _EXTRA)
    pcm = require("toytree.pcm", _EXTRA)
    source = _tree(toytree, tree, max_tips, max_file_bytes)
    tips = _tip_labels(source)
    table = _call(
        "compute diversification rates", lambda: pcm.get_tip_level_diversification(source)
    )
    series = table["DR"] if hasattr(table, "columns") else table
    values = {str(label): value for label, value in zip(series.index, series, strict=True)}
    rates = [
        {"tip": tip, "diversification_rate": _number(values[tip])} for tip in tips if tip in values
    ]
    return _result(
        "calculate_tip_diversification_rates",
        {"rates": rates, "tip_count": len(rates)},
        convention=(
            "Rates are toytree's DR values (Jetz et al. 2012) for tips only; the unnamed root "
            "row that toytree also returns is omitted."
        ),
    )


__all__ = [
    "add_tree_child_node",
    "add_tree_internal_node",
    "add_tree_internal_node_and_child",
    "add_tree_internal_node_and_subtree",
    "add_tree_sister_node",
    "bisect_tree",
    "calculate_all_node_distance_matrix",
    "calculate_descendant_distances",
    "calculate_farthest_node_distance",
    "calculate_internal_node_distance_matrix",
    "calculate_tip_distance_matrix_from_vcv",
    "calculate_tip_diversification_rates",
    "calculate_tree_node_distance",
    "drop_tree_tips",
    "extend_tree_tips_to_align",
    "find_tree_farthest_node",
    "list_tree_node_path",
    "list_tree_quadripartitions",
    "list_tree_quartets",
    "map_consensus_tree_features",
    "merge_tree_nodes",
    "parse_tree_network",
    "remove_tree_nodes",
    "resolve_tree_polytomies",
    "rotate_tree_node",
    "set_tree_node_heights",
]

_register_public_verbose(globals(), tuple(__all__))
