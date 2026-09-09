"""Native UPGMA/neighbor joining trees, UniFrac, and optional tree file codecs."""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import input_path, matrix, require
from dnakit.bioinformatics.statistics import distance_matrix
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True)
class TreeNode:
    """Branch length belongs to the edge from parent to this node."""

    name: str | None = None
    length: float = 0.0
    children: tuple[TreeNode, ...] = ()

    def __post_init__(self) -> None:
        if not math.isfinite(self.length):
            raise ConfigurationError("Tree branch lengths must be finite.")
        if not self.children and not self.name:
            raise ConfigurationError("Tree tips need names.")

    def walk(self) -> Iterator[TreeNode]:
        yield self
        for child in self.children:
            yield from child.walk()

    @property
    def tips(self) -> tuple[str, ...]:
        return tuple(
            node.name for node in self.walk() if not node.children and node.name is not None
        )


def build_tree(distances: Any, ids: Iterable[str], *, method: str = "nj") -> TreeNode:
    """Deterministic neighbor joining or size-weighted UPGMA (native NumPy).

    NJ retains negative branch estimates rather than silently altering distances.
    UPGMA assumes an ultrametric molecular clock. NJ output root is arbitrary.
    """
    np = require("numpy", "bio-analysis")
    data = distance_matrix(distances).copy()
    names = tuple(ids)
    if len(names) != len(data) or len(set(names)) != len(names) or any(not name for name in names):
        raise ConfigurationError("Unique tip IDs must match the distance matrix.")
    if method not in {"nj", "upgma"}:
        raise ConfigurationError("method must be nj or upgma.")
    nodes = [TreeNode(name) for name in names]
    sizes = [1] * len(nodes)
    heights = [0.0] * len(nodes)
    while len(nodes) > 2:
        n = len(nodes)
        criterion = data.copy()
        totals = data.sum(axis=1)
        if method == "nj":
            criterion = (n - 2) * data - totals[:, None] - totals[None, :]
        np.fill_diagonal(criterion, np.inf)
        left, right = np.unravel_index(np.argmin(criterion), criterion.shape)
        if left > right:
            left, right = right, left
        height = float(data[left, right] / 2)
        if method == "nj":
            left_length = height + float((totals[left] - totals[right]) / (2 * (n - 2)))
            right_length = float(data[left, right] - left_length)
        else:
            left_length, right_length = height - heights[left], height - heights[right]
        parent = TreeNode(
            children=(
                replace(nodes[left], length=left_length),
                replace(nodes[right], length=right_length),
            )
        )
        keep = [index for index in range(n) if index not in {left, right}]
        if method == "nj":
            distances_to_parent = (data[left, keep] + data[right, keep] - data[left, right]) / 2
        else:
            distances_to_parent = (
                sizes[left] * data[left, keep] + sizes[right] * data[right, keep]
            ) / (sizes[left] + sizes[right])
        updated = np.zeros((n - 1, n - 1))
        updated[:-1, :-1] = data[np.ix_(keep, keep)]
        updated[:-1, -1] = updated[-1, :-1] = distances_to_parent
        merged_size = sizes[left] + sizes[right]
        nodes = [nodes[index] for index in keep] + [parent]
        heights = [heights[index] for index in keep] + [height]
        sizes = [sizes[index] for index in keep] + [merged_size]
        data = updated
    half = float(data[0, 1] / 2)
    return TreeNode(
        children=(
            replace(nodes[0], length=half - (heights[0] if method == "upgma" else 0)),
            replace(nodes[1], length=half - (heights[1] if method == "upgma" else 0)),
        )
    )


def read_tree(path: str | Path, *, format: str = "newick") -> TreeNode:
    """Biopython handles Newick, Nexus, phyloXML, NeXML; tree algorithms remain native."""
    parsed = require("Bio.Phylo", "bio-formats").read(str(input_path(path)), format)

    def convert(node: Any) -> TreeNode:
        return TreeNode(
            node.name,
            float(node.branch_length or 0),
            tuple(convert(child) for child in node.clades),
        )

    tree = convert(parsed.root)
    if len(set(tree.tips)) != len(tree.tips):
        raise ConfigurationError("Tree tips must have unique names.")
    return tree


def write_tree(tree: TreeNode, path: str | Path, *, format: str = "newick") -> Path:
    base = require("Bio.Phylo.BaseTree", "bio-formats")

    def convert(node: TreeNode) -> Any:
        return base.Clade(
            branch_length=node.length,
            name=node.name,
            clades=[convert(child) for child in node.children],
        )

    target = Path(path)
    with target.open("x", encoding="utf-8") as handle:
        require("Bio.Phylo", "bio-formats").write(base.Tree(root=convert(tree)), handle, format)
    return target


def unifrac(
    counts: Any,
    taxa: Iterable[str],
    tree: TreeNode,
    *,
    weighted: bool = False,
    normalized: bool = True,
) -> Any:
    """Native rooted weighted/unweighted UniFrac over descendant branch abundance."""
    np = require("numpy", "bio-analysis")
    data = matrix(counts, min_rows=2)
    taxa = tuple(taxa)
    if (
        len(taxa) != data.shape[1]
        or len(set(taxa)) != len(taxa)
        or set(taxa) != set(tree.tips)
        or len(set(tree.tips)) != len(tree.tips)
        or (data < 0).any()
        or (data.sum(1) == 0).any()
    ):
        raise ConfigurationError(
            "Positive-total abundances and unique taxa must exactly match tree tips."
        )
    nodes = [node for node in tree.walk() if node is not tree]
    if any(node.length < 0 for node in nodes):
        raise ConfigurationError("UniFrac requires nonnegative branch lengths and a rooted tree.")
    lengths = np.asarray([node.length for node in nodes])
    lookup = {name: index for index, name in enumerate(taxa)}
    abundance = np.stack(
        [data[:, [lookup[name] for name in node.tips]].sum(axis=1) for node in nodes], axis=1
    )
    abundance = abundance / data.sum(axis=1, keepdims=True) if weighted else abundance > 0
    result = np.zeros((len(data), len(data)))
    for left in range(len(data)):
        for right in range(left):
            a, b = abundance[left], abundance[right]
            if weighted:
                numerator = float((lengths * np.abs(a - b)).sum())
                denominator = float((lengths * (a + b)).sum()) if normalized else 1.0
            else:
                numerator = float(lengths[np.logical_xor(a, b)].sum())
                denominator = float(lengths[np.logical_or(a, b)].sum())
            result[left, right] = result[right, left] = (
                numerator / denominator if denominator else 0.0
            )
    return result


def plot_tree(tree: TreeNode, *, axes: Any = None) -> Any:
    """Draw branch-length phylogram; return a matplotlib Axes without displaying it."""
    plt = require("matplotlib.pyplot", "bio-plots")
    if axes is None:
        _, axes = plt.subplots()
    position = 0

    def draw(node: TreeNode, x: float) -> float:
        nonlocal position
        if not node.children:
            y = float(position)
            position += 1
            axes.text(x, y, " " + str(node.name), va="center")
            return y
        ys = []
        for child in node.children:
            y = draw(child, x + child.length)
            axes.plot([x, x + child.length], [y, y], color="black")
            ys.append(y)
        axes.plot([x, x], [min(ys), max(ys)], color="black")
        return sum(ys) / len(ys)

    draw(tree, 0.0)
    axes.set_xlabel("Branch length")
    axes.set_yticks([])
    return axes


def _tip_set(tree: TreeNode) -> frozenset[str]:
    tips = tree.tips
    if len(tips) != len(set(tips)):
        raise ConfigurationError("Tree tip names must be unique.")
    return frozenset(tips)


def robinson_foulds(
    left: TreeNode, right: TreeNode, *, rooted: bool = False, normalized: bool = False
) -> float:
    """RF symmetric difference of nontrivial clades/splits; identical tip sets required.

    Normalization divides by the total number of splits in the input trees.
    Branch lengths and internal labels do not contribute to this distance.
    """
    taxa = _tip_set(left)
    if taxa != _tip_set(right):
        raise ConfigurationError("RF distance requires identical tip sets.")

    def splits(tree: TreeNode) -> set[frozenset[str]]:
        output = set()
        for node in tree.walk():
            group = frozenset(node.tips)
            if rooted:
                if 1 < len(group) < len(taxa):
                    output.add(group)
            elif 1 < len(group) < len(taxa) - 1:
                complement = taxa - group
                canonical = min((group, complement), key=lambda part: (len(part), sorted(part)))
                output.add(canonical)
        return output

    a, b = splits(left), splits(right)
    difference = float(len(a ^ b))
    return difference / (len(a) + len(b)) if normalized and a | b else difference


def consensus_tree(trees: Iterable[TreeNode], *, threshold: float = 0.5) -> TreeNode:
    """Rooted strict-majority consensus, or strict consensus at threshold=1.

    All inputs must share taxa and root interpretation. Clades with frequency
    strictly above threshold survive (all trees at threshold=1). Internal names
    store support fractions; lengths are means over trees containing each clade.
    """
    from collections import defaultdict

    collection = tuple(trees)
    if not collection or not 0.5 <= threshold <= 1:
        raise ConfigurationError("Provide trees and a consensus threshold in [0.5,1].")
    taxa = _tip_set(collection[0])
    observed: dict[frozenset[str], list[float]] = defaultdict(list)
    for tree in collection:
        if _tip_set(tree) != taxa:
            raise ConfigurationError("Consensus trees must share all tip names.")
        local: dict[frozenset[str], float] = defaultdict(float)
        for node in tree.walk():
            local[frozenset(node.tips)] += node.length
        for group, length in local.items():
            observed[group].append(length)
    accepted = {
        group
        for group, lengths in observed.items()
        if len(lengths) == len(collection) or len(lengths) / len(collection) > threshold
    }

    def build(group: frozenset[str]) -> TreeNode:
        lengths = observed[group]
        if len(group) == 1:
            return TreeNode(next(iter(group)), sum(lengths) / len(lengths))
        candidates = [subset for subset in accepted if subset < group]
        children = [
            subset for subset in candidates if not any(subset < other for other in candidates)
        ]
        return TreeNode(
            f"{len(lengths) / len(collection):.6g}",
            0.0 if group == taxa else sum(lengths) / len(lengths),
            tuple(build(child) for child in sorted(children, key=lambda x: sorted(x))),
        )

    return build(taxa)


def prune_tree(tree: TreeNode, retain: Iterable[str]) -> TreeNode:
    """Retain named tips, collapsing unary branches while preserving path lengths."""
    selected = frozenset(retain)
    if not selected or not selected <= _tip_set(tree):
        raise ConfigurationError("Retain a nonempty subset of existing tree tips.")

    def prune(node: TreeNode) -> TreeNode | None:
        if not node.children:
            return node if node.name in selected else None
        children = tuple(result for child in node.children if (result := prune(child)) is not None)
        if not children:
            return None
        if len(children) == 1:
            return replace(children[0], length=children[0].length + node.length)
        return replace(node, children=children)

    result = prune(tree)
    assert result is not None
    return replace(result, length=0.0)
