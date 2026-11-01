"""scikit-bio tree and taxonomy-tree operations (pack P20, part 2).

Trees go in and out as Newick strings. A node is selected by its name (tip or named internal
node) or by a list of tip names, in which case their last common ancestor is used. Every
function imports scikit-bio only when it runs.
"""

from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _skbio_p20_common as common
from dnakit.exceptions import ConfigurationError

_MAX_TREE_TEXT = 5_000_000
_MAX_TIPS = 100_000
_REROOT_METHODS = ("root_at", "unrooted_copy", "unrooted_move")


def _parse(newick: object, name: str = "newick") -> Any:
    package = common.skbio()
    value = common.text(newick, name, maximum=_MAX_TREE_TEXT)
    try:
        tree = package.TreeNode.read(io.StringIO(value))
    except Exception as exc:
        raise ConfigurationError(f"{name} is not a valid Newick tree: {exc}") from exc
    if tree.count(tips=True) > _MAX_TIPS:
        raise ConfigurationError(f"{name} may have at most {_MAX_TIPS} tips.")
    return tree


def _write(tree: Any) -> str:
    return str(tree).strip()


def _find(tree: Any, node: str | Sequence[str], name: str = "node") -> Any:
    try:
        if isinstance(node, str):
            return tree.find(node)
        if isinstance(node, Sequence) and node and all(isinstance(i, str) for i in node):
            return tree.lca(list(node))
    except Exception as exc:
        raise ConfigurationError(f"{name} {node!r} was not found in the tree: {exc}") from exc
    raise ConfigurationError(f"{name} must be a node name or a list of tip names.")


def _label(node: Any) -> str:
    """Tip name, or the comma-joined sorted tip names of an internal node."""

    if node.is_tip():
        return str(node.name)
    return ",".join(sorted(str(tip.name) for tip in node.tips()))


def shuffle_tree_tips(
    newick: str,
    *,
    names: Sequence[str] | None = None,
    k: int | None = None,
    seed: int | None = None,
    n: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shuffle tip labels over a fixed topology (``TreeNode.shuffle``), for null models.

    ``names`` restricts the shuffle to those tips (default all tips), ``k`` limits how many
    labels move, ``seed`` makes it reproducible, ``n`` is the number of shuffled trees.
    """

    tree = _parse(newick)
    count = common.integer(n, "n", minimum=1, maximum=1000)
    kwargs: dict[str, Any] = {"n": count}
    if names is not None:
        kwargs["names"] = [str(v) for v in names]
    if k is not None:
        kwargs["k"] = common.integer(k, "k", minimum=1, maximum=_MAX_TIPS)
    if seed is not None:
        kwargs["shuffler"] = common.integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    try:
        trees = [_write(item) for item in tree.shuffle(**kwargs)]
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not shuffle the tree: {exc}") from exc
    return {
        "status": "completed",
        "operation": "shuffle_tree_tips",
        "trees": trees,
        "provenance": common.provenance(method="skbio.TreeNode.shuffle"),
    }


def reroot_tree(
    newick: str,
    node: str | Sequence[str],
    *,
    above: bool = False,
    method: str = "root_at",
    verbose: bool = True,
) -> dict[str, Any]:
    """Reroot a tree at a node.

    ``method`` is ``root_at`` (rerooting at the node, or on the branch above it when ``above`` is
    True, which places the root half-way along that branch), ``unrooted_copy`` (a rerooted copy
    that keeps branch support and length) or ``unrooted_move`` (the in-place variant of the
    same operation, applied to a copy here).
    """

    tree = _parse(newick)
    chosen = common.choice(method, _REROOT_METHODS, "method")
    above_flag = common.flag(above, "above")
    target = _find(tree, node)
    try:
        if chosen == "root_at":
            rerooted = tree.root_at(target, above=above_flag)
        elif chosen == "unrooted_copy":
            rerooted = target.unrooted_copy()
        else:
            target.unrooted_move()
            rerooted = target
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not reroot the tree: {exc}") from exc
    return {
        "status": "completed",
        "operation": "reroot_tree",
        "method": chosen,
        "newick": _write(rerooted),
        "tip_count": rerooted.count(tips=True),
        "provenance": common.provenance(method=f"skbio.TreeNode.{chosen}"),
    }


def tree_tip_subsets(
    newick: str,
    *,
    node: str | Sequence[str] | None = None,
    within: Sequence[str] | None = None,
    include_full: bool = False,
    include_tips: bool = False,
    map_to_length: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the tip sets (clades) of a tree.

    With ``node`` the result is the tip set below that node (``TreeNode.subset``); without it
    every internal clade is listed (``TreeNode.subsets``). ``within`` limits tips,
    ``include_full``/``include_tips`` add the whole-tree and single-tip sets, and
    ``map_to_length`` also reports the branch length of each clade.
    """

    tree = _parse(newick)
    if node is not None:
        target = _find(tree, node)
        return {
            "status": "completed",
            "operation": "tree_tip_subsets",
            "tips": sorted(target.subset(include_self=True)),
            "provenance": common.provenance(method="skbio.TreeNode.subset"),
        }
    try:
        found = tree.subsets(
            within=None if within is None else [str(v) for v in within],
            include_full=common.flag(include_full, "include_full"),
            include_tips=common.flag(include_tips, "include_tips"),
            map_to_length=common.flag(map_to_length, "map_to_length"),
        )
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not list clades: {exc}") from exc
    items = found.items() if isinstance(found, dict) else ((s, None) for s in found)
    clades = sorted(
        (
            {"tips": sorted(s), **({"length": length} if length is not None else {})}
            for s, length in items
        ),
        key=lambda c: (len(c["tips"]), c["tips"]),
    )
    return {
        "status": "completed",
        "operation": "tree_tip_subsets",
        "clades": clades,
        "count": len(clades),
        "provenance": common.provenance(method="skbio.TreeNode.subsets"),
    }


def extract_tree_supports(newick: str, *, verbose: bool = True) -> dict[str, Any]:
    """Move numeric internal-node labels into branch supports (``TreeNode.assign_supports``).

    Returns, for every internal node with a numeric label, its clade and support value, and the
    tree rewritten without those labels.
    """

    tree = _parse(newick)
    try:
        tree.assign_supports()
    except Exception as exc:
        raise ConfigurationError(f"node labels cannot be read as supports: {exc}") from exc
    supports = [
        {"clade": _label(n), "support": float(n.support)}
        for n in tree.non_tips(include_self=False)
        if n.support is not None
    ]
    return {
        "status": "completed",
        "operation": "extract_tree_supports",
        "supports": supports,
        "newick": _write(tree),
        "provenance": common.provenance(method="skbio.TreeNode.assign_supports"),
    }


def propagate_tip_counts(
    newick: str,
    tip_counts: Mapping[str, int],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sum tip counts into every ancestor (``TreeNode.observed_node_counts``).

    ``tip_counts`` maps tip names to observed counts; the result lists a count for each node
    that has at least one observation, named by tip or by the sorted tip names of the clade.
    """

    tree = _parse(newick)
    if not isinstance(tip_counts, Mapping) or not tip_counts:
        raise ConfigurationError("tip_counts must be a mapping of tip name to count.")
    counts = {
        str(k): common.integer(v, f"count of {k}", minimum=0, maximum=2**40)
        for k, v in tip_counts.items()
    }
    try:
        found = tree.observed_node_counts(counts)
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not propagate counts: {exc}") from exc
    nodes = sorted(
        ({"node": _label(n), "count": int(c)} for n, c in found.items()),
        key=lambda r: str(r["node"]),
    )
    return {
        "status": "completed",
        "operation": "propagate_tip_counts",
        "nodes": nodes,
        "provenance": common.provenance(method="skbio.TreeNode.observed_node_counts"),
    }


def tree_diameter(newick: str, *, use_length: bool = True, verbose: bool = True) -> dict[str, Any]:
    """Largest tip-to-tip distance and the two tips that realise it (``TreeNode.maxdist``).

    ``use_length=False`` counts branches instead of summing branch lengths.
    """

    tree = _parse(newick)
    distance, (first, second) = tree.maxdist(use_length=common.flag(use_length, "use_length"))
    return {
        "status": "completed",
        "operation": "tree_diameter",
        "distance": float(distance),
        "tips": sorted([str(first.name), str(second.name)]),
        "provenance": common.provenance(method="skbio.TreeNode.maxdist"),
    }


def tree_from_linkage_matrix(
    linkage: Sequence[Sequence[float]],
    names: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert a SciPy hierarchical-clustering linkage matrix into a Newick tree.

    ``linkage`` is the (n-1) x 4 matrix returned by ``scipy.cluster.hierarchy.linkage`` and
    ``names`` the n leaf labels in the order of the clustered observations.
    """

    package = common.skbio()
    array = common.matrix(linkage, "linkage")
    if array.shape[1] != 4:
        raise ConfigurationError("linkage must have four columns.")
    leaves = common.labels(names, array.shape[0] + 1, "names", "L")
    try:
        tree = package.TreeNode.from_linkage_matrix(array, leaves)
    except Exception as exc:
        raise ConfigurationError(f"linkage is not a valid linkage matrix: {exc}") from exc
    return {
        "status": "completed",
        "operation": "tree_from_linkage_matrix",
        "newick": _write(tree),
        "tip_count": len(leaves),
        "provenance": common.provenance(method="skbio.TreeNode.from_linkage_matrix"),
    }


def build_taxonomy_tree(
    lineages: Mapping[str, Sequence[str]],
    *,
    extract_rank: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a taxonomy tree from lineages (``TreeNode.from_taxonomy``).

    ``lineages`` maps an identifier (the tip) to its ranks from the root down, for example
    ``{"g1": ["k__Bacteria", "p__Firmicutes"]}``. ``extract_rank`` splits ``k__`` style prefixes
    into a rank attribute and strips them from the node names.
    """

    package = common.skbio()
    if not isinstance(lineages, Mapping) or not lineages or len(lineages) > _MAX_TIPS:
        raise ConfigurationError("lineages must be a non-empty mapping of id to a list of ranks.")
    items = []
    for key, ranks in lineages.items():
        if isinstance(ranks, (str, bytes)) or not isinstance(ranks, Sequence) or not ranks:
            raise ConfigurationError(f"lineage of {key!r} must be a non-empty list of strings.")
        items.append((str(key), [str(r) for r in ranks]))
    try:
        tree = package.TreeNode.from_taxonomy(
            items, extract_rank=common.flag(extract_rank, "extract_rank")
        )
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not build the taxonomy tree: {exc}") from exc
    return {
        "status": "completed",
        "operation": "build_taxonomy_tree",
        "newick": _write(tree),
        "tip_count": tree.count(tips=True),
        "node_count": tree.count(),
        "provenance": common.provenance(method="skbio.TreeNode.from_taxonomy"),
    }


def export_tree_lineages(
    newick: str,
    *,
    allow_empty: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the lineage (names from the root down) of every tip (``TreeNode.to_taxonomy``)."""

    tree = _parse(newick)
    try:
        rows = [
            {"id": str(tip.name), "lineage": [str(r) for r in lineage]}
            for tip, lineage in tree.to_taxonomy(
                allow_empty=common.flag(allow_empty, "allow_empty")
            )
        ]
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not export lineages: {exc}") from exc
    return {
        "status": "completed",
        "operation": "export_tree_lineages",
        "lineages": rows,
        "count": len(rows),
        "provenance": common.provenance(method="skbio.TreeNode.to_taxonomy"),
    }


def build_tree_from_taxdump(
    nodes_path: str,
    *,
    names_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build the NCBI taxonomy tree from ``nodes.dmp`` (and optionally ``names.dmp``).

    Reads the dump tables with scikit-bio's ``taxdump`` reader and calls
    ``TreeNode.from_taxdump``; nodes are named by scientific name when ``names_path`` is given
    and by taxonomy id otherwise.
    """

    import pandas as pd

    package = common.skbio()
    nodes_file = common.input_path(nodes_path, "nodes_path")
    names_file = None if names_path is None else common.input_path(names_path, "names_path")
    try:
        nodes = package.io.read(
            str(nodes_file), format="taxdump", into=pd.DataFrame, scheme="nodes"
        )
        names = (
            None
            if names_file is None
            else package.io.read(
                str(names_file), format="taxdump", into=pd.DataFrame, scheme="names"
            )
        )
        tree = package.TreeNode.from_taxdump(nodes, names)
    except Exception as exc:
        raise ConfigurationError(f"scikit-bio could not build the taxdump tree: {exc}") from exc
    return {
        "status": "completed",
        "operation": "build_tree_from_taxdump",
        "newick": _write(tree),
        "node_count": tree.count(),
        "tip_count": tree.count(tips=True),
        "provenance": common.provenance(method="skbio.TreeNode.from_taxdump"),
    }


__all__ = [
    "build_taxonomy_tree",
    "build_tree_from_taxdump",
    "export_tree_lineages",
    "extract_tree_supports",
    "propagate_tip_counts",
    "reroot_tree",
    "shuffle_tree_tips",
    "tree_diameter",
    "tree_from_linkage_matrix",
    "tree_tip_subsets",
]

_register_public_verbose(globals(), tuple(__all__))
