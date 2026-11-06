"""Biopython GenePop editing and ``Bio.Cluster`` clustering operations added by integration pack P19.

``Bio.PopGen.GenePop`` operations read a GenePop file and write new GenePop files (a new output file or
directory is always created, nothing is overwritten). ``Bio.Cluster`` operations (the C Clustering
Library, Eisen's Cluster 3.0) take a data matrix as a list of rows; rows are the items clustered unless
``transpose=True``. Missing values are marked with ``mask`` (0 = missing, 1 = present). Distances are named
``euclidean``, ``cityblock``, ``pearson``, ``abs_pearson``, ``uncentered``, ``abs_uncentered``,
``spearman`` and ``kendall``. Biopython is imported only when a function runs.

The random starting points of k-means, k-medians, k-medoids and the self-organizing map come from the
C library and cannot be seeded: pass ``initial_clusters`` to make k-means/k-medians/k-medoids
deterministic.
"""

from __future__ import annotations

import importlib.metadata
import math
import re
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_MAX_CELLS = 5_000_000
_MAX_FILE_BYTES = 256 * 1024 * 1024
_DISTANCES = {
    "euclidean": "e",
    "cityblock": "b",
    "pearson": "c",
    "abs_pearson": "a",
    "uncentered": "u",
    "abs_uncentered": "x",
    "spearman": "s",
    "kendall": "k",
}
_LINKAGES = {"single": "s", "complete": "m", "average": "a", "centroid": "c"}
_CENTER_METHODS = {"mean": "a", "median": "m"}
_CLUSTER_DISTANCE_METHODS = {
    "mean": "a",
    "median": "m",
    "single": "s",
    "complete": "x",
    "average": "v",
}
_NAME = re.compile(r"[A-Za-z0-9_.\-]{1,100}")


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _bio(module: str) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return require(module, _EXTRA)


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _key(table: Mapping[str, str], value: object, name: str) -> str:
    if not isinstance(value, str) or value not in table:
        raise ConfigurationError(f"{name} must be one of {list(table)}; got {value!r}.")
    return table[value]


def _matrix(value: object, name: str, *, minimum_rows: int = 1, minimum_columns: int = 1) -> Any:
    """Validate a rectangular list of finite numbers and return a float numpy array."""

    np = require("numpy", _EXTRA)
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of rows.")
    width = None
    rows = []
    for row in value:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence):
            raise ConfigurationError(f"{name} must be a list of rows of numbers.")
        if width is None:
            width = len(row)
        if len(row) != width or width == 0:
            raise ConfigurationError(f"{name} rows must all have the same, non-zero length.")
        for item in row:
            if (
                isinstance(item, bool)
                or not isinstance(item, (int, float))
                or not math.isfinite(item)
            ):
                raise ConfigurationError(f"{name} must hold finite numbers only.")
        rows.append([float(item) for item in row])
    if len(rows) < minimum_rows or (width or 0) < minimum_columns:
        raise ConfigurationError(
            f"{name} needs at least {minimum_rows} rows and {minimum_columns} columns."
        )
    if len(rows) * (width or 0) > _MAX_CELLS:
        raise ConfigurationError(f"{name} may hold at most {_MAX_CELLS} values.")
    return np.array(rows, dtype=float)


def _mask(value: object, shape: tuple[int, int]) -> Any:
    if value is None:
        return None
    np = require("numpy", _EXTRA)
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError("mask must be a list of rows of 0/1 integers.")
    array = np.array(value, dtype=int) if value else np.zeros((0, 0), dtype=int)
    if array.shape != shape or ((array != 0) & (array != 1)).any():
        raise ConfigurationError(
            f"mask must have the shape of the data {shape} and hold only 0 and 1."
        )
    return array


def _weights(value: object, length: int) -> Any:
    if value is None:
        return None
    np = require("numpy", _EXTRA)
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != length:
        raise ConfigurationError(f"weights must be a list of {length} numbers.")
    if any(
        isinstance(w, bool) or not isinstance(w, (int, float)) or not math.isfinite(w) or w < 0
        for w in value
    ):
        raise ConfigurationError("weights must be finite non-negative numbers.")
    return np.array(value, dtype=float)


def _setup(
    data: object, mask: object, weights: object, transpose: object
) -> tuple[Any, Any, Any, bool]:
    array = _matrix(data, "data")
    flip = _flag(transpose, "transpose")
    return (
        array,
        _mask(mask, array.shape),
        _weights(weights, array.shape[0] if flip else array.shape[1]),
        flip,
    )


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return function(*args, **kwargs)
    except (ValueError, TypeError, RuntimeError, MemoryError) as exc:
        raise ConfigurationError(f"Bio.Cluster could not complete the operation: {exc}") from exc


def _list(array: Any) -> Any:
    return array.tolist() if hasattr(array, "tolist") else array


def _initial(value: object, count: int, n_clusters: int, name: str = "initial_clusters") -> Any:
    if value is None:
        return None
    np = require("numpy", _EXTRA)
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) != count:
        raise ConfigurationError(f"{name} must hold one cluster number per item ({count}).")
    for item in value:
        _integer(item, name, minimum=0, maximum=n_clusters - 1)
    return np.array(value, dtype=np.int32)


# --------------------------------------------------------------------------------------------
# GenePop
# --------------------------------------------------------------------------------------------


def _read_genepop(path: str | Path) -> tuple[Any, Path]:
    source = input_path(path)
    if source.stat().st_size > _MAX_FILE_BYTES:
        raise ConfigurationError(f"{source.name} is larger than {_MAX_FILE_BYTES} bytes.")
    genepop = _bio("Bio.PopGen.GenePop")
    with source.open(encoding="utf-8") as handle:
        try:
            return genepop.read(handle), source
        except (ValueError, IndexError, AttributeError) as exc:
            raise ConfigurationError(f"{source.name} is not a valid GenePop file: {exc}") from exc


def _safe_stem(value: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9_.\-]", "_", value.strip())[:100]
    return stem or "item"


def split_genepop_by_population(
    genepop_path: str | Path,
    output_dir: str | Path,
    *,
    population_names: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one GenePop file per population (``Bio.PopGen.GenePop.Record.split_in_pops``).

    Population names in GenePop files are not reliable, so by default the files are called
    ``pop1.gen``, ``pop2.gen``, ... in file order; pass ``population_names`` (one per population) to
    name them. Every output keeps all loci. ``output_dir`` must not exist yet.
    """

    record, source = _read_genepop(genepop_path)
    count = len(record.populations)
    if population_names is None:
        names = [f"pop{index}" for index in range(1, count + 1)]
    else:
        if isinstance(population_names, (str, bytes)) or len(population_names) != count:
            raise ConfigurationError(f"population_names must hold {count} names.")
        names = [str(item) for item in population_names]
        if len({_safe_stem(item) for item in names}) != count:
            raise ConfigurationError("population_names must be distinct.")
    parts = record.split_in_pops(names)
    directory = cli.new_directory(output_dir, "output_dir")
    files = {}
    for name, part in parts.items():
        target = directory / f"{_safe_stem(name)}.gen"
        target.write_text(str(part), encoding="utf-8")
        files[name] = str(target)
    return {
        "status": "completed",
        "operation": "split_genepop_by_population",
        "source": str(source),
        "population_count": count,
        "locus_count": len(record.loci_list),
        "files": files,
        "provenance": _provenance(parser="Bio.PopGen.GenePop.Record.split_in_pops"),
    }


def split_genepop_by_locus(
    genepop_path: str | Path, output_dir: str | Path, *, verbose: bool = True
) -> dict[str, Any]:
    """Write one GenePop file per locus (``Bio.PopGen.GenePop.Record.split_in_loci``).

    Each output ``<locus>.gen`` holds a single locus and all populations. ``output_dir`` must not exist
    yet.
    """

    record, source = _read_genepop(genepop_path)
    parts = record.split_in_loci(record)
    stems = [_safe_stem(name) for name in parts]
    if len(set(stems)) != len(stems):
        raise ConfigurationError("locus names are not distinct after removing unsafe characters.")
    directory = cli.new_directory(output_dir, "output_dir")
    files = {}
    for (name, part), stem in zip(parts.items(), stems, strict=False):
        target = directory / f"{stem}.gen"
        target.write_text(str(part), encoding="utf-8")
        files[name] = str(target)
    return {
        "status": "completed",
        "operation": "split_genepop_by_locus",
        "source": str(source),
        "population_count": len(record.populations),
        "locus_count": len(parts),
        "files": files,
        "provenance": _provenance(parser="Bio.PopGen.GenePop.Record.split_in_loci"),
    }


def remove_genepop_populations(
    genepop_path: str | Path,
    output_path: str | Path,
    positions: Sequence[int],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a copy of a GenePop file without the populations at the given 0-based positions.

    Uses ``Bio.PopGen.GenePop.Record.remove_population`` (the file-based ``FileRecord`` variant writes
    wrongly padded alleles in Biopython 1.88, so it is not used). ``output_path`` must not exist.
    """

    record, source = _read_genepop(genepop_path)
    count = len(record.populations)
    if isinstance(positions, (str, bytes)) or not isinstance(positions, Sequence) or not positions:
        raise ConfigurationError("positions must be a non-empty list of integers.")
    chosen = sorted(
        {_integer(item, "positions", minimum=0, maximum=count - 1) for item in positions}
    )
    if len(chosen) == count:
        raise ConfigurationError("positions would remove every population.")
    target = cli.new_file(output_path, "output_path")
    for position in reversed(chosen):
        record.remove_population(position)
    target.write_text(str(record), encoding="utf-8")
    return {
        "status": "completed",
        "operation": "remove_genepop_populations",
        "source": str(source),
        "output": str(target),
        "removed_positions": chosen,
        "population_count": count - len(chosen),
        "provenance": _provenance(parser="Bio.PopGen.GenePop.Record.remove_population"),
    }


def remove_genepop_loci(
    genepop_path: str | Path,
    output_path: str | Path,
    *,
    positions: Sequence[int] | None = None,
    names: Sequence[str] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a copy of a GenePop file without the given loci, selected by 0-based position or by name.

    Uses ``Bio.PopGen.GenePop.Record.remove_locus_by_position`` / ``remove_locus_by_name``. Exactly one of
    ``positions`` and ``names`` is required. ``output_path`` must not exist.
    """

    if (positions is None) == (names is None):
        raise ConfigurationError("Provide exactly one of positions and names.")
    record, source = _read_genepop(genepop_path)
    loci = list(record.loci_list)
    if positions is not None:
        if (
            isinstance(positions, (str, bytes))
            or not isinstance(positions, Sequence)
            or not positions
        ):
            raise ConfigurationError("positions must be a non-empty list of integers.")
        chosen = sorted(
            {_integer(item, "positions", minimum=0, maximum=len(loci) - 1) for item in positions},
            reverse=True,
        )
        removed = [loci[position] for position in reversed(chosen)]
        for position in chosen:
            record.remove_locus_by_position(position)
    else:
        if isinstance(names, (str, bytes)) or not isinstance(names, Sequence) or not names:
            raise ConfigurationError("names must be a non-empty list of locus names.")
        unknown = [item for item in names if item not in loci]
        if unknown:
            raise ConfigurationError(f"names holds loci that are not in the file: {unknown}.")
        removed = list(dict.fromkeys(names))
        for name in removed:
            record.remove_locus_by_name(name)
    if not record.loci_list:
        raise ConfigurationError("the selection would remove every locus.")
    target = cli.new_file(output_path, "output_path")
    target.write_text(str(record), encoding="utf-8")
    return {
        "status": "completed",
        "operation": "remove_genepop_loci",
        "source": str(source),
        "output": str(target),
        "removed_loci": removed,
        "remaining_loci": list(record.loci_list),
        "provenance": _provenance(
            parser="Bio.PopGen.GenePop.Record.remove_locus_by_position / remove_locus_by_name"
        ),
    }


# --------------------------------------------------------------------------------------------
# Bio.Cluster
# --------------------------------------------------------------------------------------------


def kmeans_cluster(
    data: Sequence[Sequence[float]],
    n_clusters: int,
    *,
    method: str = "mean",
    distance: str = "euclidean",
    n_passes: int = 1,
    initial_clusters: Sequence[int] | None = None,
    mask: Sequence[Sequence[int]] | None = None,
    weights: Sequence[float] | None = None,
    transpose: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """K-means (``method="mean"``) or k-medians (``method="median"``) clustering (``Bio.Cluster.kcluster``).

    Runs the EM algorithm ``n_passes`` times from random starts and keeps the best solution, or once from
    ``initial_clusters`` (one cluster number per item, fully deterministic). Returns the cluster number of
    every item, the within-cluster sum of distances (``error``) and how often the best solution was found.
    """

    array, mask_array, weight_array, flip = _setup(data, mask, weights, transpose)
    count = array.shape[1] if flip else array.shape[0]
    k = _integer(n_clusters, "n_clusters", minimum=1, maximum=max(count, 1))
    cluster = _bio("Bio.Cluster")
    start = _initial(initial_clusters, count, k)
    clusterid, error, found = _call(
        cluster.kcluster,
        array,
        nclusters=k,
        mask=mask_array,
        weight=weight_array,
        transpose=flip,
        npass=_integer(n_passes, "n_passes", minimum=1, maximum=100_000),
        method=_key(_CENTER_METHODS, method, "method"),
        dist=_key(_DISTANCES, distance, "distance"),
        initialid=start,
    )
    return {
        "status": "completed",
        "operation": "kmeans_cluster",
        "method": method,
        "distance": distance,
        "cluster_ids": [int(item) for item in clusterid],
        "error": float(error),
        "times_found": int(found),
        "provenance": _provenance(parser="Bio.Cluster.kcluster"),
    }


def kmedoids_cluster(
    distance_matrix: Sequence[Sequence[float]],
    n_clusters: int,
    *,
    n_passes: int = 1,
    initial_clusters: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """K-medoids clustering of items given their distance matrix (``Bio.Cluster.kmedoids``).

    ``distance_matrix`` is a square symmetric matrix (only the lower triangle is read). ``initial_clusters``
    gives a starting cluster number (0 to ``n_clusters`` - 1) per item. The cluster number returned for an
    item is the index of its medoid item. Returns the assignments, the within-cluster sum of
    distances and how often the best solution was found.
    """

    array = _matrix(distance_matrix, "distance_matrix")
    if array.shape[0] != array.shape[1]:
        raise ConfigurationError("distance_matrix must be square.")
    size = array.shape[0]
    k = _integer(n_clusters, "n_clusters", minimum=1, maximum=size)
    cluster = _bio("Bio.Cluster")
    clusterid, error, found = _call(
        cluster.kmedoids,
        array,
        nclusters=k,
        npass=_integer(n_passes, "n_passes", minimum=1, maximum=100_000),
        initialid=_initial(initial_clusters, size, k),
    )
    return {
        "status": "completed",
        "operation": "kmedoids_cluster",
        "cluster_ids": [int(item) for item in clusterid],
        "error": float(error),
        "times_found": int(found),
        "provenance": _provenance(parser="Bio.Cluster.kmedoids"),
    }


def _tree_nodes(tree: Any) -> list[dict[str, Any]]:
    return [
        {
            "left": int(tree[i].left),
            "right": int(tree[i].right),
            "distance": float(tree[i].distance),
        }
        for i in range(len(tree))
    ]


def _build_tree(nodes: object) -> Any:
    cluster = _bio("Bio.Cluster")
    if isinstance(nodes, (str, bytes)) or not isinstance(nodes, Sequence) or not nodes:
        raise ConfigurationError(
            "nodes must be a non-empty list of {left, right, distance} records."
        )
    built = []
    for item in nodes:
        if not isinstance(item, Mapping) or not {"left", "right"} <= set(item):
            raise ConfigurationError(
                "each node needs 'left' and 'right' (and optionally 'distance')."
            )
        distance = item.get("distance", 0.0)
        if (
            isinstance(distance, bool)
            or not isinstance(distance, (int, float))
            or not math.isfinite(distance)
        ):
            raise ConfigurationError("node distance must be a finite number.")
        built.append(
            cluster.Node(
                _integer(item["left"], "left", minimum=-len(nodes), maximum=len(nodes)),
                _integer(item["right"], "right", minimum=-len(nodes), maximum=len(nodes)),
                float(distance),
            )
        )
    return _call(cluster.Tree, built)


def hierarchical_cluster(
    data: Sequence[Sequence[float]] | None = None,
    *,
    distance_matrix: Sequence[Sequence[float]] | None = None,
    linkage: str = "complete",
    distance: str = "euclidean",
    mask: Sequence[Sequence[int]] | None = None,
    weights: Sequence[float] | None = None,
    transpose: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pairwise hierarchical clustering (``Bio.Cluster.treecluster``).

    ``linkage`` is ``single``, ``complete`` (the Biopython default), ``average`` or ``centroid``. Give
    either ``data`` (items by measurements; all four linkages) or a square ``distance_matrix``
    (single, complete and average only). The tree is returned as nodes in merge order: item ``i >= 0`` is
    a leaf, ``-k`` refers to the node created at step ``k`` (1-based), as in Eisen's Cluster 3.0.
    """

    method = _key(_LINKAGES, linkage, "linkage")
    if (data is None) == (distance_matrix is None):
        raise ConfigurationError("Provide exactly one of data and distance_matrix.")
    cluster = _bio("Bio.Cluster")
    if data is not None:
        array, mask_array, weight_array, flip = _setup(data, mask, weights, transpose)
        tree = _call(
            cluster.treecluster,
            array,
            mask=mask_array,
            weight=weight_array,
            transpose=flip,
            method=method,
            dist=_key(_DISTANCES, distance, "distance"),
        )
        size = array.shape[1] if flip else array.shape[0]
    else:
        if method == "c":
            raise ConfigurationError("centroid linkage needs data, not a distance matrix.")
        matrix = _matrix(distance_matrix, "distance_matrix", minimum_rows=2)
        if matrix.shape[0] != matrix.shape[1]:
            raise ConfigurationError("distance_matrix must be square.")
        tree = _call(cluster.treecluster, None, method=method, distancematrix=matrix)
        size = matrix.shape[0]
    return {
        "status": "completed",
        "operation": "hierarchical_cluster",
        "linkage": linkage,
        "item_count": size,
        "nodes": _tree_nodes(tree),
        "provenance": _provenance(parser="Bio.Cluster.treecluster"),
    }


def cut_cluster_tree(
    nodes: Sequence[Mapping[str, Any]], n_clusters: int, *, verbose: bool = True
) -> dict[str, Any]:
    """Cut a hierarchical clustering tree into ``n_clusters`` clusters (``Bio.Cluster.Tree.cut``).

    ``nodes`` is the node list returned by ``hierarchical_cluster``. Returns the cluster number of each
    item.
    """

    tree = _build_tree(nodes)
    k = _integer(n_clusters, "n_clusters", minimum=1, maximum=len(nodes) + 1)
    ids = _call(tree.cut, k)
    return {
        "status": "completed",
        "operation": "cut_cluster_tree",
        "n_clusters": k,
        "cluster_ids": [int(item) for item in ids],
        "provenance": _provenance(parser="Bio.Cluster.Tree.cut"),
    }


def sort_cluster_tree(
    nodes: Sequence[Mapping[str, Any]],
    *,
    order: Sequence[float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Order the leaves of a hierarchical clustering tree (``Bio.Cluster.Tree.sort``).

    Left and right subnodes are swapped so that the leaf order tends to follow ``order`` (one number per
    item; default is the item index). Returns the item indices in left-to-right order and the sorted node
    list.
    """

    tree = _build_tree(nodes)
    np = require("numpy", _EXTRA)
    values = None
    if order is not None:
        if isinstance(order, (str, bytes)) or len(order) != len(nodes) + 1:
            raise ConfigurationError(f"order must hold {len(nodes) + 1} numbers.")
        values = np.array([float(item) for item in order], dtype=float)
    indices = _call(tree.sort) if values is None else _call(tree.sort, values)
    return {
        "status": "completed",
        "operation": "sort_cluster_tree",
        "leaf_order": [int(item) for item in indices],
        "nodes": _tree_nodes(tree),
        "provenance": _provenance(parser="Bio.Cluster.Tree.sort"),
    }


def som_cluster(
    data: Sequence[Sequence[float]],
    *,
    nx_grid: int = 2,
    ny_grid: int = 1,
    init_tau: float = 0.02,
    n_iterations: int = 100,
    distance: str = "euclidean",
    mask: Sequence[Sequence[int]] | None = None,
    weights: Sequence[float] | None = None,
    transpose: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Self-organizing map on a rectangular grid (``Bio.Cluster.somcluster``).

    Returns the grid cell (x, y) of every item and the centroid vector of every cell
    (``cell_data[x][y]``). The map starts from random vectors; results vary between calls.
    """

    array, mask_array, weight_array, flip = _setup(data, mask, weights, transpose)
    cluster = _bio("Bio.Cluster")
    tau = float(init_tau)
    if not math.isfinite(tau) or tau <= 0:
        raise ConfigurationError("init_tau must be a positive number.")
    cell_ids, cell_data = _call(
        cluster.somcluster,
        array,
        mask=mask_array,
        weight=weight_array,
        transpose=flip,
        nxgrid=_integer(nx_grid, "nx_grid", minimum=1, maximum=1000),
        nygrid=_integer(ny_grid, "ny_grid", minimum=1, maximum=1000),
        inittau=tau,
        niter=_integer(n_iterations, "n_iterations", minimum=1, maximum=10_000_000),
        dist=_key(_DISTANCES, distance, "distance"),
    )
    return {
        "status": "completed",
        "operation": "som_cluster",
        "cell_ids": _list(cell_ids),
        "cell_data": _list(cell_data),
        "provenance": _provenance(parser="Bio.Cluster.somcluster"),
    }


def cluster_distance(
    data: Sequence[Sequence[float]],
    index1: Sequence[int],
    index2: Sequence[int],
    *,
    method: str = "mean",
    distance: str = "euclidean",
    mask: Sequence[Sequence[int]] | None = None,
    weights: Sequence[float] | None = None,
    transpose: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distance between two clusters of items (``Bio.Cluster.clusterdistance``).

    ``index1`` and ``index2`` list the member items of each cluster. ``method`` defines the cluster
    distance: ``mean`` (distance between the means), ``median`` (between the medians), ``single``
    (smallest pairwise distance), ``complete`` (largest pairwise distance) or ``average`` (mean of
    pairwise distances).
    """

    array, mask_array, weight_array, flip = _setup(data, mask, weights, transpose)
    count = array.shape[1] if flip else array.shape[0]
    np = require("numpy", _EXTRA)
    groups = []
    for name, value in (("index1", index1), ("index2", index2)):
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
            raise ConfigurationError(f"{name} must be a non-empty list of item indices.")
        groups.append(
            np.array(
                [_integer(item, name, minimum=0, maximum=count - 1) for item in value],
                dtype=np.int32,
            )
        )
    cluster = _bio("Bio.Cluster")
    value = _call(
        cluster.clusterdistance,
        array,
        mask=mask_array,
        weight=weight_array,
        index1=groups[0],
        index2=groups[1],
        method=_key(_CLUSTER_DISTANCE_METHODS, method, "method"),
        dist=_key(_DISTANCES, distance, "distance"),
        transpose=flip,
    )
    return {
        "status": "completed",
        "operation": "cluster_distance",
        "method": method,
        "distance_type": distance,
        "distance": float(value),
        "provenance": _provenance(parser="Bio.Cluster.clusterdistance"),
    }


def cluster_centroids(
    data: Sequence[Sequence[float]],
    cluster_ids: Sequence[int],
    *,
    method: str = "mean",
    mask: Sequence[Sequence[int]] | None = None,
    transpose: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Centroid of every cluster as the mean or median over its items (``Bio.Cluster.clustercentroids``).

    ``cluster_ids`` gives the non-negative cluster number of each item. Returns one centroid row per
    cluster and the matching 0/1 mask (0 = no data).
    """

    array = _matrix(data, "data")
    flip = _flag(transpose, "transpose")
    count = array.shape[1] if flip else array.shape[0]
    if (
        isinstance(cluster_ids, (str, bytes))
        or not isinstance(cluster_ids, Sequence)
        or len(cluster_ids) != count
    ):
        raise ConfigurationError(f"cluster_ids must hold one cluster number per item ({count}).")
    np = require("numpy", _EXTRA)
    ids = np.array(
        [_integer(item, "cluster_ids", minimum=0, maximum=count) for item in cluster_ids],
        dtype=np.int32,
    )
    cluster = _bio("Bio.Cluster")
    centroids, centroid_mask = _call(
        cluster.clustercentroids,
        array,
        mask=_mask(mask, array.shape),
        clusterid=ids,
        method=_key(_CENTER_METHODS, method, "method"),
        transpose=flip,
    )
    return {
        "status": "completed",
        "operation": "cluster_centroids",
        "method": method,
        "centroids": _list(centroids),
        "centroid_mask": _list(centroid_mask),
        "provenance": _provenance(parser="Bio.Cluster.clustercentroids"),
    }


def profile_distance_matrix(
    data: Sequence[Sequence[float]],
    *,
    distance: str = "euclidean",
    mask: Sequence[Sequence[int]] | None = None,
    weights: Sequence[float] | None = None,
    transpose: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Distance matrix between the expression/profile rows (or columns) of a data matrix.

    ``distance`` is one of ``euclidean``, ``cityblock``, ``pearson`` (1 - Pearson r), ``abs_pearson``
    (1 - |r|), ``uncentered`` (1 - uncentered correlation), ``abs_uncentered``, ``spearman`` (1 - rank
    correlation) or ``kendall`` (1 - Kendall tau); this is ``Bio.Cluster.distancematrix``. Returns the
    full symmetric matrix.
    """

    array, mask_array, weight_array, flip = _setup(data, mask, weights, transpose)
    cluster = _bio("Bio.Cluster")
    lower = _call(
        cluster.distancematrix,
        array,
        mask=mask_array,
        weight=weight_array,
        transpose=flip,
        dist=_key(_DISTANCES, distance, "distance"),
    )
    size = array.shape[1] if flip else array.shape[0]
    full = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row):
            full[row][column] = full[column][row] = float(lower[row][column])
    return {
        "status": "completed",
        "operation": "profile_distance_matrix",
        "distance": distance,
        "item_count": size,
        "matrix": full,
        "provenance": _provenance(parser="Bio.Cluster.distancematrix"),
    }


def cluster_pca(data: Sequence[Sequence[float]], *, verbose: bool = True) -> dict[str, Any]:
    """Principal component analysis of a data matrix (``Bio.Cluster.pca``).

    Returns the column means, the principal components (one per row), the coordinates of every row along
    the components and the eigenvalues, all sorted by decreasing eigenvalue. Adding the column means to
    ``coordinates @ components`` reproduces the data.
    """

    array = _matrix(data, "data", minimum_rows=2, minimum_columns=2)
    cluster = _bio("Bio.Cluster")
    means, coordinates, components, eigenvalues = _call(
        cluster.pca, array
    )  # Biopython's return order
    return {
        "status": "completed",
        "operation": "cluster_pca",
        "column_means": _list(means),
        "components": _list(components),
        "coordinates": _list(coordinates),
        "eigenvalues": _list(eigenvalues),
        "provenance": _provenance(parser="Bio.Cluster.pca"),
    }


__all__ = [
    "cluster_centroids",
    "cluster_distance",
    "cluster_pca",
    "cut_cluster_tree",
    "hierarchical_cluster",
    "kmeans_cluster",
    "kmedoids_cluster",
    "profile_distance_matrix",
    "remove_genepop_loci",
    "remove_genepop_populations",
    "som_cluster",
    "sort_cluster_tree",
    "split_genepop_by_locus",
    "split_genepop_by_population",
]

_register_public_verbose(globals(), tuple(__all__))
