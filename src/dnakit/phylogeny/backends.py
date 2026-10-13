"""Optional-backend tree summaries, layout coordinates and nucleotide simulation.

- ``read_tree_summary()``, ``calculate_patristic_distances()`` and
  ``calculate_branch_length_distance()`` use ``dendropy`` to read Newick or
  Nexus trees, report tips and branch lengths, and compare trees.
- ``calculate_node_coordinates()``, ``calculate_quartet_distance()`` and
  ``root_tree_on_outgroup()`` use ``toytree``. Layout coordinates come from
  toytree's layout classes, so nothing is drawn.
- ``simulate_nucleotide_alignment()``, ``simulate_nucleotide_sequences()`` and
  ``simulate_substitution_counts()`` use ``pyvolve`` with an explicit seed.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
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

_DENDROPY_EXTRA = "bio-dendropy"
_TOYTREE_EXTRA = "bio-toytree"
_PYVOLVE_EXTRA = "bio-pyvolve"
_TREE_SCHEMAS = ("newick", "nexus")
_LAYOUTS = ("r", "l", "u", "d", "c", "unrooted")
_NUCLEOTIDES = ("A", "C", "G", "T")
_TRANSITIONS = frozenset({"AG", "CT"})
_EXCHANGEABILITY_KEYS = ("AC", "AG", "AT", "CG", "CT", "GT")
_MAX_FILE_BYTES_LIMIT = 200_000_000
_MAX_SEED = 2**32 - 1
_SIMULATION_CONVENTION = (
    "Branch lengths are expected substitutions per site (pyvolve scales the rate matrix to a "
    "mean rate of 1). Transitions A<->G and C<->T are scaled by kappa, which is the HKY85 rate "
    "matrix. The seed fixes the random generator, so identical inputs give identical output."
)


def _bounded_int(value: object, name: str, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be an integer in [{low}, {high}].")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _label(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or len(value) > 200 or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string of at most 200 characters.")
    return value


def _package(name: str) -> str:
    return f"{name}=={metadata.version(name)}"


def _json_number(value: object) -> float | None:
    number = float(value)
    return number if math.isfinite(number) else None


def _tree_file(value: object, max_bytes: int) -> Path:
    if not isinstance(value, (str, Path)) or "\x00" in str(value):
        raise ConfigurationError("path must be a string or Path to an existing tree file.")
    _bounded_int(max_bytes, "max_file_bytes", 1, _MAX_FILE_BYTES_LIMIT)
    path = input_path(value)
    size = path.stat().st_size
    if size > max_bytes:
        raise ConfigurationError(f"Tree file has {size} bytes; max_file_bytes is {max_bytes}.")
    return path


def _fresh_output(value: object) -> Path:
    if not isinstance(value, (str, Path)) or "\x00" in str(value):
        raise ConfigurationError("output_path must be a string or Path.")
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Output already exists and will not be overwritten: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    return target


def _edge_length(value: object) -> float | None:
    """Finite branch length, or None when the file omitted it or wrote a non-finite value."""

    if value is None:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _require_lengths(lengths: Sequence[float | None], operation: str) -> None:
    missing = sum(1 for length in lengths if length is None)
    if missing:
        raise ConfigurationError(
            f"{operation} needs a branch length on every edge; {missing} edge(s) have none."
        )


def _load_dendropy_trees(
    dendropy: Any, path: Path, schema: str, taxon_namespace: Any = None
) -> Any:
    options: dict[str, Any] = {"path": str(path), "schema": schema, "preserve_underscores": True}
    if taxon_namespace is not None:
        options["taxon_namespace"] = taxon_namespace
    try:
        return dendropy.TreeList.get(**options)
    except Exception as exc:
        raise ConfigurationError(
            f"dendropy could not parse {path.name} as {schema}: {exc}"
        ) from exc


def _select_tree(trees: Any, tree_index: int) -> tuple[Any, int]:
    count = len(trees)
    if count == 0:
        raise ConfigurationError("The tree file contains no trees.")
    if tree_index >= count:
        raise ConfigurationError(f"tree_index {tree_index} is out of range for {count} tree(s).")
    return trees[tree_index], count


def _dendropy_leaves(tree: Any) -> list[tuple[str, Any, float | None]]:
    """Tips as (label, taxon, pendant length), sorted by label; labels must be unique."""

    leaves = []
    for leaf in tree.leaf_node_iter():
        taxon = leaf.taxon
        if taxon is None or not taxon.label:
            raise ConfigurationError("Every tip needs a non-empty label.")
        leaves.append((taxon.label, taxon, _edge_length(leaf.edge.length)))
    if len({label for label, _, _ in leaves}) != len(leaves):
        raise ConfigurationError("Tip labels must be unique.")
    return sorted(leaves, key=lambda item: item[0])


def _dendropy_edges(tree: Any) -> list[float | None]:
    """Branch lengths of every edge except the root's own edge, in preorder."""

    root = tree.seed_node
    return [
        _edge_length(node.edge.length) for node in tree.preorder_node_iter() if node is not root
    ]


def _rooting(tree: Any) -> str:
    if tree.is_rooted:
        return "rooted"
    if tree.is_unrooted:
        return "unrooted"
    return "undefined"


def _read_toytree(toytree: Any, path: Path) -> Any:
    try:
        return toytree.tree(path)
    except Exception as exc:
        raise ConfigurationError(f"toytree could not parse {path.name} as Newick: {exc}") from exc


def _outgroup_names(value: object) -> list[str]:
    items = [value] if isinstance(value, str) else value
    if not isinstance(items, (list, tuple)) or not 1 <= len(items) <= 1000:
        raise ConfigurationError(
            "outgroup must be a tip label or a non-empty list of tip labels (at most 1000)."
        )
    names = [_label(item, "outgroup tip") for item in items]
    if len(set(names)) != len(names):
        raise ConfigurationError("outgroup tips must be unique.")
    return names


def _nucleotide_inputs(
    length: object, kappa: object, state_freqs: object, seed: object, max_length: int
) -> tuple[int, float, list[float], int]:
    length = _bounded_int(length, "length", 1, max_length)
    if isinstance(kappa, bool) or not isinstance(kappa, (int, float)) or not math.isfinite(kappa):
        raise ConfigurationError("kappa must be a finite positive number.")
    if not 0 < kappa <= 1000:
        raise ConfigurationError("kappa must be in (0, 1000].")
    if not isinstance(state_freqs, (list, tuple)) or len(state_freqs) != len(_NUCLEOTIDES):
        raise ConfigurationError("state_freqs must be four frequencies for A, C, G and T.")
    freqs = []
    for value in state_freqs:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ConfigurationError("state_freqs must be finite positive numbers.")
        if not math.isfinite(value) or value <= 0:
            raise ConfigurationError("state_freqs must be finite positive numbers.")
        freqs.append(float(value))
    total = sum(freqs)
    if abs(total - 1.0) > 1e-6:
        raise ConfigurationError(f"state_freqs must sum to 1; they sum to {total}.")
    freqs = [value / total for value in freqs]
    seed = _bounded_int(seed, "seed", 0, _MAX_SEED)
    return length, float(kappa), freqs, seed


def _model_summary(kappa: float, freqs: Sequence[float]) -> dict[str, object]:
    return {
        "name": "HKY85-style nucleotide rate matrix",
        "kappa": kappa,
        "state_freqs": {base: freq for base, freq in zip(_NUCLEOTIDES, freqs, strict=True)},
    }


def _nucleotide_evolver(
    pyvolve: Any,
    tree_path: Path,
    *,
    length: int,
    kappa: float,
    freqs: Sequence[float],
    max_bases: int,
) -> tuple[Any, list[str]]:
    """Return a pyvolve Evolver for a nucleotide model on a Newick tree, and sorted tip names."""

    try:
        tree = pyvolve.read_tree(file=str(tree_path))
    except Exception as exc:
        raise ConfigurationError(
            f"pyvolve could not parse {tree_path.name} as Newick: {exc}"
        ) from exc
    names: list[str] = []
    stack = [tree]
    while stack:
        node = stack.pop()
        if node.children:
            stack.extend(node.children)
        else:
            names.append(node.name)
    if len(names) < 2:
        raise ConfigurationError("The tree needs at least two tips.")
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ConfigurationError("Tip names must be non-empty and unique.")
    if len(names) * length > max_bases:
        raise ConfigurationError(
            f"tips * length is {len(names) * length}; max_bases is {max_bases}."
        )
    mu = {key: (kappa if key in _TRANSITIONS else 1.0) for key in _EXCHANGEABILITY_KEYS}
    model = pyvolve.Model("nucleotide", {"mu": mu, "state_freqs": list(freqs)})
    partition = pyvolve.Partition(models=model, size=length)
    return pyvolve.Evolver(tree=tree, partitions=partition), sorted(names)


def read_tree_summary(
    path: str | Path,
    *,
    schema: str = "newick",
    tree_index: int = 0,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Read one tree from a Newick or Nexus file and report its tips and branch lengths."""

    dendropy = require("dendropy", _DENDROPY_EXTRA)
    schema = _choice(schema, _TREE_SCHEMAS, "schema")
    tree_index = _bounded_int(tree_index, "tree_index", 0, 100_000)
    max_tips = _bounded_int(max_tips, "max_tips", 1, 100_000)
    file_path = _tree_file(path, max_file_bytes)
    tree, count = _select_tree(_load_dendropy_trees(dendropy, file_path, schema), tree_index)
    leaves = _dendropy_leaves(tree)
    if len(leaves) > max_tips:
        raise ConfigurationError(f"Tree has {len(leaves)} tips; max_tips is {max_tips}.")
    lengths = _dendropy_edges(tree)
    known = [length for length in lengths if length is not None]
    return {
        "status": "completed",
        "operation": "read_tree_summary",
        "result": {
            "tree_count": count,
            "tree_index": tree_index,
            "rooting": _rooting(tree),
            "tip_count": len(leaves),
            "tips": [{"label": label, "branch_length": length} for label, _, length in leaves],
            "branch_count": len(lengths),
            "branch_lengths": lengths,
            "missing_branch_lengths": len(lengths) - len(known),
            "total_branch_length": float(sum(known)),
        },
        "provenance": {
            "package": _package("dendropy"),
            "convention": (
                "branch_lengths lists every edge except the root's own edge, in preorder; "
                "None marks a missing or non-finite length."
            ),
        },
    }


def calculate_patristic_distances(
    path: str | Path,
    *,
    schema: str = "newick",
    tree_index: int = 0,
    max_taxa: int = 200,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Tip-to-tip patristic distance matrix (summed branch lengths along the path) for one tree."""

    dendropy = require("dendropy", _DENDROPY_EXTRA)
    schema = _choice(schema, _TREE_SCHEMAS, "schema")
    tree_index = _bounded_int(tree_index, "tree_index", 0, 100_000)
    max_taxa = _bounded_int(max_taxa, "max_taxa", 2, 1000)
    file_path = _tree_file(path, max_file_bytes)
    tree, count = _select_tree(_load_dendropy_trees(dendropy, file_path, schema), tree_index)
    leaves = _dendropy_leaves(tree)
    if len(leaves) > max_taxa:
        raise ConfigurationError(f"Tree has {len(leaves)} tips; max_taxa is {max_taxa}.")
    _require_lengths(_dendropy_edges(tree), "Patristic distances")
    distances = tree.phylogenetic_distance_matrix()
    taxa = [taxon for _, taxon, _ in leaves]
    size = len(taxa)
    matrix: list[list[float | None]] = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row + 1, size):
            value = _json_number(distances.patristic_distance(taxa[row], taxa[column]))
            matrix[row][column] = value
            matrix[column][row] = value
    return {
        "status": "completed",
        "operation": "calculate_patristic_distances",
        "result": {
            "tree_count": count,
            "tree_index": tree_index,
            "taxa": [label for label, _, _ in leaves],
            "matrix": matrix,
        },
        "provenance": {
            "package": _package("dendropy"),
            "convention": (
                "matrix[i][j] is the path length between taxa i and j in the input's branch "
                "length units; the root's own edge lies on no tip-to-tip path."
            ),
        },
    }


def calculate_branch_length_distance(
    path_a: str | Path,
    path_b: str | Path,
    *,
    schema: str = "newick",
    tree_index: int = 0,
    max_taxa: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Felsenstein 2004 branch-length distance between two trees with identical tip labels."""

    dendropy = require("dendropy", _DENDROPY_EXTRA)
    treecompare = require("dendropy.calculate.treecompare", _DENDROPY_EXTRA)
    schema = _choice(schema, _TREE_SCHEMAS, "schema")
    tree_index = _bounded_int(tree_index, "tree_index", 0, 100_000)
    max_taxa = _bounded_int(max_taxa, "max_taxa", 2, 100_000)
    file_a = _tree_file(path_a, max_file_bytes)
    file_b = _tree_file(path_b, max_file_bytes)
    namespace = dendropy.TaxonNamespace()
    tree_a, _ = _select_tree(_load_dendropy_trees(dendropy, file_a, schema, namespace), tree_index)
    tree_b, _ = _select_tree(_load_dendropy_trees(dendropy, file_b, schema, namespace), tree_index)
    leaves_a = _dendropy_leaves(tree_a)
    labels = [label for label, _, _ in leaves_a]
    if labels != [label for label, _, _ in _dendropy_leaves(tree_b)]:
        raise ConfigurationError("Both trees must have the same tip labels.")
    if len(labels) > max_taxa:
        raise ConfigurationError(f"Trees have {len(labels)} tips; max_taxa is {max_taxa}.")
    _require_lengths(_dendropy_edges(tree_a), "Branch-length distance")
    _require_lengths(_dendropy_edges(tree_b), "Branch-length distance")
    distance = _json_number(treecompare.euclidean_distance(tree_a, tree_b))
    return {
        "status": "completed",
        "operation": "calculate_branch_length_distance",
        "result": {
            "tree_index": tree_index,
            "tip_count": len(labels),
            "distance": distance,
        },
        "provenance": {
            "package": _package("dendropy"),
            "convention": (
                "Euclidean branch-length distance on unrooted splits: the two root-adjacent "
                "edges are merged, so the placement of the root does not change the value."
            ),
        },
    }


def calculate_node_coordinates(
    path: str | Path,
    *,
    layout: str = "r",
    use_edge_lengths: bool = True,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Node and tip coordinates for a toytree layout, computed without drawing a figure."""

    toytree = require("toytree", _TOYTREE_EXTRA)
    layout_module = require("toytree.layout", _TOYTREE_EXTRA)
    style_module = require("toytree.style", _TOYTREE_EXTRA)
    layout_name = _choice(layout, _LAYOUTS, "layout")
    if not isinstance(use_edge_lengths, bool):
        raise ConfigurationError("use_edge_lengths must be a boolean.")
    max_tips = _bounded_int(max_tips, "max_tips", 1, 100_000)
    tree = _read_toytree(toytree, _tree_file(path, max_file_bytes))
    if tree.ntips > max_tips:
        raise ConfigurationError(f"Tree has {tree.ntips} tips; max_tips is {max_tips}.")
    style = style_module.TreeStyle(layout=layout_name, use_edge_lengths=use_edge_lengths)
    if layout_name == "c":
        layout_class = layout_module.CircularLayout
    elif layout_name == "unrooted":
        layout_class = layout_module.UnrootedLayout
    else:
        layout_class = layout_module.LinearLayout
    coords = layout_class(tree, style).coords
    nodes = sorted(tree.treenode.traverse(), key=lambda node: node.idx)
    if coords is None or len(coords) != len(nodes):
        raise ConfigurationError("toytree returned coordinates for an unexpected node count.")
    rows = []
    for node in nodes:
        x_value, y_value = coords[node.idx]
        rows.append(
            {
                "index": int(node.idx),
                "name": node.name or None,
                "is_tip": bool(node.is_leaf()),
                "parent": None if node.up is None else int(node.up.idx),
                "x": _json_number(x_value),
                "y": _json_number(y_value),
            }
        )
    return {
        "status": "completed",
        "operation": "calculate_node_coordinates",
        "result": {
            "layout": layout_name,
            "use_edge_lengths": use_edge_lengths,
            "tip_count": int(tree.ntips),
            "node_count": len(rows),
            "nodes": rows,
        },
        "provenance": {
            "package": _package("toytree"),
            "convention": (
                "index is toytree's node index; x and y are layout units, not pixels. "
                "No figure is drawn."
            ),
        },
    }


def calculate_quartet_distance(
    path_a: str | Path,
    path_b: str | Path,
    *,
    max_taxa: int = 40,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Quartet counts that agree or differ between two unrooted trees with identical tips."""

    toytree = require("toytree", _TOYTREE_EXTRA)
    distance = require("toytree.distance", _TOYTREE_EXTRA)
    max_taxa = _bounded_int(max_taxa, "max_taxa", 4, 200)
    tree_a = _read_toytree(toytree, _tree_file(path_a, max_file_bytes))
    tree_b = _read_toytree(toytree, _tree_file(path_b, max_file_bytes))
    labels = sorted(tree_a.get_tip_labels())
    if labels != sorted(tree_b.get_tip_labels()):
        raise ConfigurationError("Both trees must have the same tip labels.")
    size = len(labels)
    if size < 4:
        raise ConfigurationError("Quartet distance needs at least four tips.")
    if size > max_taxa:
        raise ConfigurationError(f"Trees have {size} tips; max_taxa is {max_taxa}.")
    metrics = distance.get_treedist_quartets(tree_a, tree_b)
    total = round(float(metrics["Q"]))
    different = round(float(metrics["D"]))
    return {
        "status": "completed",
        "operation": "calculate_quartet_distance",
        "result": {
            "tip_count": size,
            "total_quartets": total,
            "shared_quartets": round(float(metrics["S"])),
            "different_quartets": different,
            "different_fraction": different / total,
        },
        "provenance": {
            "package": _package("toytree"),
            "convention": (
                "Quartets are the 4-tip subsets. A quartet is shared when both trees induce the "
                "same resolved topology on it, and different when their topologies disagree."
            ),
        },
    }


def root_tree_on_outgroup(
    path: str | Path,
    outgroup: str | Sequence[str],
    output_path: str | Path,
    *,
    max_tips: int = 10_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Root a Newick tree on the edge above the outgroup and write the rooted Newick file."""

    toytree = require("toytree", _TOYTREE_EXTRA)
    names = _outgroup_names(outgroup)
    max_tips = _bounded_int(max_tips, "max_tips", 1, 100_000)
    tree = _read_toytree(toytree, _tree_file(path, max_file_bytes))
    if tree.ntips > max_tips:
        raise ConfigurationError(f"Tree has {tree.ntips} tips; max_tips is {max_tips}.")
    tips = set(tree.get_tip_labels())
    missing = [name for name in names if name not in tips]
    if missing:
        raise ConfigurationError(f"Outgroup tips are not in the tree: {', '.join(missing)}")
    if len(names) >= tree.ntips:
        raise ConfigurationError("The outgroup must leave at least one ingroup tip.")
    target = _fresh_output(output_path)
    try:
        rooted = tree.root(*names)
    except Exception as exc:
        raise ConfigurationError(f"toytree could not root on the outgroup: {exc}") from exc
    newick = rooted.write()
    if not newick.endswith("\n"):
        newick += "\n"
    with target.open("x", encoding="utf-8") as handle:
        handle.write(newick)
    return {
        "status": "completed",
        "operation": "root_tree_on_outgroup",
        "result": {
            "output_path": str(target),
            "outgroup": names,
            "tip_count": int(rooted.ntips),
            "newick": newick,
        },
        "provenance": {
            "package": _package("toytree"),
            "convention": (
                "The new root is placed on the edge above the most recent common ancestor of "
                "the outgroup tips; the root edge length is split between the two sides."
            ),
        },
    }


def simulate_nucleotide_alignment(
    tree_path: str | Path,
    output_path: str | Path,
    *,
    length: int = 1000,
    kappa: float = 1.0,
    state_freqs: Sequence[float] = (0.25, 0.25, 0.25, 0.25),
    seed: int = 0,
    max_length: int = 100_000,
    max_bases: int = 10_000_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Simulate a DNA alignment along a Newick tree with a fixed seed and write it as FASTA."""

    pyvolve = require("pyvolve", _PYVOLVE_EXTRA)
    length, kappa, freqs, seed = _nucleotide_inputs(length, kappa, state_freqs, seed, max_length)
    max_bases = _bounded_int(max_bases, "max_bases", 1, 1_000_000_000)
    tree_file = _tree_file(tree_path, max_file_bytes)
    target = _fresh_output(output_path)
    evolver, taxa = _nucleotide_evolver(
        pyvolve, tree_file, length=length, kappa=kappa, freqs=freqs, max_bases=max_bases
    )
    evolver(seed=seed, seqfile=str(target), seqfmt="fasta", ratefile=False, infofile=False)
    return {
        "status": "completed",
        "operation": "simulate_nucleotide_alignment",
        "result": {
            "output_path": str(target),
            "format": "fasta",
            "taxa": taxa,
            "sequence_count": len(taxa),
            "site_count": length,
            "model": _model_summary(kappa, freqs),
            "seed": seed,
        },
        "provenance": {"package": _package("pyvolve"), "convention": _SIMULATION_CONVENTION},
    }


def simulate_nucleotide_sequences(
    tree_path: str | Path,
    *,
    length: int = 1000,
    kappa: float = 1.0,
    state_freqs: Sequence[float] = (0.25, 0.25, 0.25, 0.25),
    seed: int = 0,
    max_length: int = 100_000,
    max_bases: int = 10_000_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Simulate DNA sequences at the tips of a Newick tree with a fixed seed, without a file."""

    pyvolve = require("pyvolve", _PYVOLVE_EXTRA)
    length, kappa, freqs, seed = _nucleotide_inputs(length, kappa, state_freqs, seed, max_length)
    max_bases = _bounded_int(max_bases, "max_bases", 1, 1_000_000_000)
    tree_file = _tree_file(tree_path, max_file_bytes)
    evolver, taxa = _nucleotide_evolver(
        pyvolve, tree_file, length=length, kappa=kappa, freqs=freqs, max_bases=max_bases
    )
    evolver(seed=seed, seqfile=False, ratefile=False, infofile=False)
    sequences = evolver.get_sequences()
    return {
        "status": "completed",
        "operation": "simulate_nucleotide_sequences",
        "result": {
            "taxa": taxa,
            "site_count": length,
            "sequences": {taxon: str(sequences[taxon]) for taxon in taxa},
            "model": _model_summary(kappa, freqs),
            "seed": seed,
        },
        "provenance": {"package": _package("pyvolve"), "convention": _SIMULATION_CONVENTION},
    }


def simulate_substitution_counts(
    tree_path: str | Path,
    *,
    length: int = 1000,
    kappa: float = 1.0,
    state_freqs: Sequence[float] = (0.25, 0.25, 0.25, 0.25),
    seed: int = 0,
    max_length: int = 100_000,
    max_bases: int = 10_000_000,
    max_file_bytes: int = 10_000_000,
    verbose: bool = True,
) -> dict[str, object]:
    """Simulate a nucleotide alignment with a fixed seed and count substitutions on each branch."""

    pyvolve = require("pyvolve", _PYVOLVE_EXTRA)
    length, kappa, freqs, seed = _nucleotide_inputs(length, kappa, state_freqs, seed, max_length)
    max_bases = _bounded_int(max_bases, "max_bases", 1, 1_000_000_000)
    tree_file = _tree_file(tree_path, max_file_bytes)
    evolver, taxa = _nucleotide_evolver(
        pyvolve, tree_file, length=length, kappa=kappa, freqs=freqs, max_bases=max_bases
    )
    # algorithm=1 uses the exact Gillespie jump chain; matrix exponentiation undercounts.
    evolver(seed=seed, algorithm=1, seqfile=False, ratefile=False, infofile=False)
    per_model = list(evolver.branch_substitution_counts.values())
    if len(per_model) != 1:
        raise ConfigurationError(
            "pyvolve returned substitution counts for an unexpected model count."
        )
    counts = {branch: int(count) for branch, count in sorted(per_model[0].items())}
    return {
        "status": "completed",
        "operation": "simulate_substitution_counts",
        "result": {
            "taxa": taxa,
            "site_count": length,
            "branch_substitution_counts": counts,
            "total_substitutions": sum(counts.values()),
            "model": _model_summary(kappa, freqs),
            "seed": seed,
        },
        "provenance": {
            "package": _package("pyvolve"),
            "convention": (
                _SIMULATION_CONVENTION + " Branch keys are tip names or pyvolve internalNode names;"
                " the root edge is not counted."
            ),
        },
    }


__all__ = [
    "calculate_branch_length_distance",
    "calculate_node_coordinates",
    "calculate_patristic_distances",
    "calculate_quartet_distance",
    "read_tree_summary",
    "root_tree_on_outgroup",
    "simulate_nucleotide_alignment",
    "simulate_nucleotide_sequences",
    "simulate_substitution_counts",
]

_register_public_verbose(globals(), tuple(__all__))
