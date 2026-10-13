"""Optional tskit and dadi adapters for tree-sequence and frequency-spectrum statistics.

The functions in this module are explicit adapters over optional packages:

- ``build_tree_sequence()``, ``read_tree_sequence()`` and ``write_tree_sequence()``
  use ``tskit`` to build a tree sequence from explicit tables, or to read and write
  a ``.trees`` file. No coalescent simulation is run.
- ``calculate_diversity()``, ``calculate_tajimas_d()``, ``calculate_fst()`` and
  ``calculate_site_frequency_spectrum()`` call tskit's own statistics on the tree
  sequence built from the tables.
- ``expected_neutral_sfs()``, ``fold_sfs()``, ``estimate_sfs_theta()`` and
  ``calculate_sfs_log_likelihood()`` use ``dadi`` for the constant-size neutral
  spectrum, spectrum folding and the Poisson likelihood of an observed spectrum.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Mapping, Sequence
from numbers import Integral, Real
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_TSKIT_EXTRA = "bio-tskit"
_DADI_EXTRA = "bio-dadi"
_DEFAULT_MAX_ROWS = 100_000
_MAX_ROWS_LIMIT = 1_000_000
_DEFAULT_MAX_BYTES = 100_000_000
_MAX_BYTES_LIMIT = 1_000_000_000
_DEFAULT_MAX_SAMPLE_SIZE = 1000
_MAX_SAMPLE_SIZE_LIMIT = 2000
_DEFAULT_GRID_POINTS = 1000
_MIN_GRID_POINTS = 50
_MAX_GRID_POINTS = 10_000
_TABLE_KEYS = frozenset({"sequence_length", "nodes", "edges", "sites", "mutations"})
_NODE_KEYS = frozenset({"time", "sample"})
_EDGE_KEYS = frozenset({"left", "right", "parent", "child"})
_SITE_KEYS = frozenset({"position", "ancestral_state"})
_MUTATION_KEYS = frozenset({"site", "node", "derived_state", "parent", "time"})


def _integer(value: object, name: str, *, minimum: int = 0, maximum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ConfigurationError(f"{name} must be an integer.")
    number = int(value)
    if number < minimum or (maximum is not None and number > maximum):
        bound = f"[{minimum}, {maximum}]" if maximum is not None else f">= {minimum}"
        raise ConfigurationError(f"{name} must be an integer {bound}; got {number}.")
    return number


def _number(value: object, name: str, *, positive: bool = False, minimum: float = 0.0) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ConfigurationError(f"{name} must be a finite number.")
    number = float(value)
    if not math.isfinite(number):
        raise ConfigurationError(f"{name} must be a finite number.")
    if positive and number <= 0:
        raise ConfigurationError(f"{name} must be greater than zero.")
    if number < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum:g}.")
    return number


def _text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string without NUL characters.")
    return value


def _limit(value: object, name: str, maximum: int) -> int:
    return _integer(value, name, minimum=1, maximum=maximum)


def _finite(value: object) -> float | None:
    number = float(value)
    return number if math.isfinite(number) else None


def _provenance(package: str, convention: str) -> dict[str, str]:
    return {
        "package": f"{package} {importlib.metadata.version(package)}",
        "convention": convention,
    }


_TABLE_CONVENTION = (
    "intervals are half-open [left, right) on 0 to sequence_length; tables are sorted before "
    "building; a sample set lists sample node IDs"
)


def _rows(value: object, table: str, max_rows: int) -> list[Mapping[str, Any]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{table} must be a list of row mappings.")
    if len(value) > max_rows:
        raise ConfigurationError(
            f"{table} has {len(value)} rows, above max_rows={max_rows}; "
            "increase the explicit limit."
        )
    rows = []
    for row in value:
        if not isinstance(row, Mapping):
            raise ConfigurationError(f"each {table} row must be a mapping.")
        rows.append(row)
    return rows


def _check_row(
    row: Mapping[str, Any], table: str, allowed: frozenset[str], required: frozenset[str]
) -> None:
    unknown = sorted(str(key) for key in row if key not in allowed)
    if unknown:
        raise ConfigurationError(f"{table} row has unknown keys: {unknown}.")
    missing = sorted(required - set(row))
    if missing:
        raise ConfigurationError(f"{table} row is missing keys: {missing}.")


def _build_tables(tskit: Any, data: object, max_rows: int) -> Any:
    """Turn explicit table rows into a sorted tskit TableCollection."""

    if not isinstance(data, Mapping):
        raise ConfigurationError(
            "tables must be a mapping with sequence_length, nodes, edges, sites and mutations."
        )
    unknown = sorted(str(key) for key in data if key not in _TABLE_KEYS)
    if unknown:
        raise ConfigurationError(f"tables has unknown keys: {unknown}.")
    if "sequence_length" not in data or "nodes" not in data:
        raise ConfigurationError("tables must include sequence_length and nodes.")
    length = _number(data["sequence_length"], "sequence_length", positive=True)
    limit = _limit(max_rows, "max_rows", _MAX_ROWS_LIMIT)
    nodes = _rows(data["nodes"], "nodes", limit)
    edges = _rows(data.get("edges", []), "edges", limit)
    sites = _rows(data.get("sites", []), "sites", limit)
    mutations = _rows(data.get("mutations", []), "mutations", limit)
    if not nodes:
        raise ConfigurationError("nodes must contain at least one row.")

    collection = tskit.TableCollection(sequence_length=length)
    for row in nodes:
        _check_row(row, "nodes", _NODE_KEYS, frozenset({"time"}))
        sample = row.get("sample", False)
        if not isinstance(sample, bool):
            raise ConfigurationError("node sample must be true or false.")
        collection.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE if sample else 0,
            time=_number(row["time"], "node time", minimum=0.0),
        )
    for row in edges:
        _check_row(row, "edges", _EDGE_KEYS, _EDGE_KEYS)
        collection.edges.add_row(
            left=_number(row["left"], "edge left", minimum=0.0),
            right=_number(row["right"], "edge right", minimum=0.0),
            parent=_integer(row["parent"], "edge parent"),
            child=_integer(row["child"], "edge child"),
        )
    for row in sites:
        _check_row(row, "sites", _SITE_KEYS, _SITE_KEYS)
        collection.sites.add_row(
            position=_number(row["position"], "site position", minimum=0.0),
            ancestral_state=_text(row["ancestral_state"], "ancestral_state"),
        )
    for row in mutations:
        _check_row(row, "mutations", _MUTATION_KEYS, frozenset({"site", "node", "derived_state"}))
        options: dict[str, float] = {}
        if "time" in row:
            options["time"] = _number(row["time"], "mutation time", minimum=0.0)
        collection.mutations.add_row(
            site=_integer(row["site"], "mutation site"),
            node=_integer(row["node"], "mutation node"),
            derived_state=_text(row["derived_state"], "derived_state"),
            parent=_integer(row.get("parent", -1), "mutation parent", minimum=-1),
            **options,
        )
    try:
        collection.sort()
    except (tskit.TskitException, ValueError) as exc:
        raise ConfigurationError(f"tskit could not sort the tables: {exc}") from exc
    return collection


def _tree_sequence(tskit: Any, data: object, max_rows: int) -> Any:
    tables = _build_tables(tskit, data, max_rows)
    try:
        return tables.tree_sequence()
    except (tskit.TskitException, ValueError) as exc:
        raise ConfigurationError(f"tskit rejected the tables: {exc}") from exc


def _summary(ts: Any) -> dict[str, Any]:
    return {
        "sequence_length": float(ts.sequence_length),
        "num_nodes": int(ts.num_nodes),
        "num_edges": int(ts.num_edges),
        "num_samples": int(ts.num_samples),
        "num_trees": int(ts.num_trees),
        "num_sites": int(ts.num_sites),
        "num_mutations": int(ts.num_mutations),
        "num_populations": int(ts.num_populations),
    }


def _sample_sets(ts: Any, sample_sets: object) -> dict[str, list[int]]:
    """Validate named sample sets of sample node IDs; the default is all samples."""

    samples = {int(node) for node in ts.samples()}
    if sample_sets is None:
        if len(samples) < 2:
            raise ConfigurationError("the tree sequence needs at least two sample nodes.")
        return {"all": sorted(samples)}
    if not isinstance(sample_sets, Mapping) or not sample_sets:
        raise ConfigurationError("sample_sets must map set names to lists of sample node IDs.")
    result: dict[str, list[int]] = {}
    for name, members in sample_sets.items():
        label = _text(name, "sample set name")
        if isinstance(members, (str, bytes)) or not isinstance(members, Sequence):
            raise ConfigurationError(f"sample set {label!r} must be a list of node IDs.")
        nodes = [_integer(node, f"node in sample set {label!r}") for node in members]
        if len(nodes) < 2:
            raise ConfigurationError(f"sample set {label!r} needs at least two sample nodes.")
        if len(set(nodes)) != len(nodes):
            raise ConfigurationError(f"sample set {label!r} lists a node more than once.")
        not_samples = [node for node in nodes if node not in samples]
        if not_samples:
            raise ConfigurationError(f"sample set {label!r} has non-sample nodes: {not_samples}.")
        result[label] = nodes
    return result


def build_tree_sequence(
    tables: Mapping[str, Any],
    *,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a tree sequence from explicit node, edge, site and mutation rows, using tskit.

    ``tables`` has ``sequence_length``, ``nodes``, ``edges``, ``sites`` and ``mutations``.
    Nodes take ``time`` and optional ``sample``; edges take ``left``, ``right``,
    ``parent`` and ``child``; sites take ``position`` and ``ancestral_state``;
    mutations take ``site``, ``node``, ``derived_state`` and optional ``parent``/``time``.
    The rows are sorted and validated by tskit. ``max_rows`` bounds each table.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    return {
        "status": "completed",
        "operation": "build_tree_sequence",
        **_summary(ts),
        "provenance": _provenance("tskit", _TABLE_CONVENTION),
    }


def read_tree_sequence(
    source: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a tskit ``.trees`` file and summarize it, using tskit.

    The file size is checked against ``max_bytes`` before tskit opens it.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    path = input_path(source)
    limit = _limit(max_bytes, "max_bytes", _MAX_BYTES_LIMIT)
    size = path.stat().st_size
    if size > limit:
        raise ConfigurationError(f"{path.name} has {size} bytes, above max_bytes={limit}.")
    try:
        ts = tskit.load(str(path))
    except (tskit.TskitException, tskit.FileFormatError, OSError, ValueError) as exc:
        raise ConfigurationError(f"tskit could not read {path.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_tree_sequence",
        "source": str(path),
        "bytes": size,
        **_summary(ts),
        "provenance": _provenance("tskit", _TABLE_CONVENTION),
    }


def write_tree_sequence(
    tables: Mapping[str, Any],
    destination: str | Path,
    *,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a tree sequence from tables and write it as a tskit ``.trees`` file.

    An existing file is never overwritten, and the parent directory must exist.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if not isinstance(destination, (str, Path)):
        raise ConfigurationError("destination must be a path.")
    target = Path(_text(str(destination), "destination")).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Output file already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Output directory does not exist: {target.parent}")
    ts = _tree_sequence(tskit, tables, max_rows)
    ts.dump(str(target))
    return {
        "status": "completed",
        "operation": "write_tree_sequence",
        "destination": str(target),
        "bytes": target.stat().st_size,
        **_summary(ts),
        "provenance": _provenance("tskit", _TABLE_CONVENTION),
    }


def calculate_diversity(
    tables: Mapping[str, Any],
    *,
    sample_sets: Mapping[str, Sequence[int]] | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide diversity per base for each sample set, using tskit.

    Values are tskit's site-mode mean pairwise differences divided by ``sequence_length``.
    With ``sample_sets`` omitted, all sample nodes form one set named ``"all"``.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    sets = _sample_sets(ts, sample_sets)
    values = ts.diversity(sample_sets=list(sets.values()), mode="site", span_normalise=True)
    return {
        "status": "completed",
        "operation": "calculate_diversity",
        "sequence_length": float(ts.sequence_length),
        "diversity_per_base": {
            name: _finite(value) for name, value in zip(sets, values, strict=True)
        },
        "provenance": _provenance("tskit", "site-mode diversity divided by sequence_length"),
    }


def calculate_tajimas_d(
    tables: Mapping[str, Any],
    *,
    sample_sets: Mapping[str, Sequence[int]] | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tajima's D and segregating-site counts for each sample set, using tskit.

    Values are None when tskit returns NaN, for example when a set has no segregating sites.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    sets = _sample_sets(ts, sample_sets)
    groups = list(sets.values())
    d_values = ts.Tajimas_D(sample_sets=groups, mode="site")
    counts = ts.segregating_sites(sample_sets=groups, mode="site", span_normalise=False)
    return {
        "status": "completed",
        "operation": "calculate_tajimas_d",
        "tajimas_d": {name: _finite(value) for name, value in zip(sets, d_values, strict=True)},
        "segregating_sites": {
            name: round(float(value)) for name, value in zip(sets, counts, strict=True)
        },
        "provenance": _provenance("tskit", "Tajima's D from tskit site-mode diversity and S"),
    }


def calculate_fst(
    tables: Mapping[str, Any],
    *,
    sample_sets: Mapping[str, Sequence[int]],
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fst between exactly two sample sets, using tskit's site-mode estimator.

    tskit computes ``1 - 2(d_X + d_Y) / (d_X + 2 d_XY + d_Y)``, where ``d_X`` and ``d_Y``
    are within-set diversities and ``d_XY`` is the between-set divergence. The value
    is None when the denominator is zero.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    ts = _tree_sequence(tskit, tables, max_rows)
    sets = _sample_sets(ts, sample_sets)
    if len(sets) != 2:
        raise ConfigurationError("calculate_fst needs exactly two sample sets.")
    (name_a, nodes_a), (name_b, nodes_b) = sets.items()
    value = ts.Fst(sample_sets=[nodes_a, nodes_b], mode="site")
    return {
        "status": "completed",
        "operation": "calculate_fst",
        "fst": _finite(value),
        "population_a": {"name": name_a, "samples": len(nodes_a)},
        "population_b": {"name": name_b, "samples": len(nodes_b)},
        "provenance": _provenance(
            "tskit", "site-mode Fst, 1 - 2(d_X + d_Y) / (d_X + 2 d_XY + d_Y)"
        ),
    }


def calculate_site_frequency_spectrum(
    tables: Mapping[str, Any],
    *,
    sample_sets: Mapping[str, Sequence[int]] | None = None,
    folded: bool = False,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Site frequency spectrum counts for each sample set, using tskit.

    ``spectrum[i]`` counts sites whose derived allele is carried by ``i`` of the set's
    sample nodes. With ``folded=True`` the minor-allele spectrum is returned, so counts
    of ``i`` and ``n - i`` are merged into index ``min(i, n - i)``.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if not isinstance(folded, bool):
        raise ConfigurationError("folded must be true or false.")
    ts = _tree_sequence(tskit, tables, max_rows)
    sets = _sample_sets(ts, sample_sets)
    spectra = {}
    for name, nodes in sets.items():
        values = ts.allele_frequency_spectrum(
            sample_sets=[nodes], mode="site", span_normalise=False, polarised=not folded
        )
        spectra[name] = [round(float(value)) for value in values]
    return {
        "status": "completed",
        "operation": "calculate_site_frequency_spectrum",
        "folded": folded,
        "spectra": spectra,
        "provenance": _provenance("tskit", "site counts from tskit allele_frequency_spectrum"),
    }


_DADI_CONVENTION = (
    "spectrum[i] is the count for derived-allele frequency i of n sampled chromosomes, "
    "i from 0 to n; entries 0 and n are monomorphic and masked (None) where dadi masks them"
)


def _grid(value: object) -> int:
    return _integer(value, "grid_points", minimum=_MIN_GRID_POINTS, maximum=_MAX_GRID_POINTS)


def _sample_size(value: object, limit: int) -> int:
    return _integer(value, "sample_size", minimum=2, maximum=limit)


def _limit_sample_size(value: object) -> int:
    return _limit(value, "max_sample_size", _MAX_SAMPLE_SIZE_LIMIT)


def _spectrum_values(np: Any, spectrum: Any) -> list[float | None]:
    """Plain list of floats, with masked entries as None."""

    data = np.ma.getdata(spectrum)
    mask = np.ma.getmaskarray(spectrum)
    return [None if masked else _finite(value) for value, masked in zip(data, mask, strict=True)]


def _polymorphic_counts(counts: list[float]) -> None:
    if sum(counts[1:-1]) <= 0:
        raise ConfigurationError("observed needs at least one count between indices 1 and n-1.")


def _observed_counts(value: object, limit: int) -> list[float]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) < 3:
        raise ConfigurationError(
            "observed must be a list of counts for derived-allele frequencies 0..n with n >= 2."
        )
    if len(value) - 1 > limit:
        raise ConfigurationError(
            f"observed has sample size {len(value) - 1}, above max_sample_size={limit}."
        )
    counts = []
    for item in value:
        if isinstance(item, bool) or not isinstance(item, Real) or not math.isfinite(float(item)):
            raise ConfigurationError("observed counts must be finite numbers.")
        if float(item) < 0:
            raise ConfigurationError("observed counts must be non-negative.")
        counts.append(float(item))
    return counts


def expected_neutral_sfs(
    sample_size: int,
    *,
    theta: float = 1.0,
    folded: bool = False,
    grid_points: int = _DEFAULT_GRID_POINTS,
    max_sample_size: int = _DEFAULT_MAX_SAMPLE_SIZE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected site frequency spectrum for a constant-size neutral population, using dadi.

    dadi's ``Demographics1D.snm`` integrates the standard neutral model on a grid of
    ``grid_points``, so values approximate ``theta / i``. The relative error falls as the
    grid grows. The result is multiplied by ``theta``; ``folded=True`` applies dadi's fold.
    """

    dadi = require("dadi", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    limit = _limit_sample_size(max_sample_size)
    n = _sample_size(sample_size, limit)
    scale = _number(theta, "theta", positive=True)
    if not isinstance(folded, bool):
        raise ConfigurationError("folded must be true or false.")
    points = _grid(grid_points)
    spectrum = dadi.Demographics1D.snm(None, (n,), points) * scale
    if folded:
        spectrum = spectrum.fold()
    return {
        "status": "completed",
        "operation": "expected_neutral_sfs",
        "sample_size": n,
        "theta": scale,
        "folded": folded,
        "spectrum": _spectrum_values(np, spectrum),
        "provenance": _provenance(
            "dadi", _DADI_CONVENTION + "; grid-integrated, E[xi_i] = theta / i"
        ),
    }


def fold_sfs(
    spectrum: Sequence[float],
    *,
    max_sample_size: int = _DEFAULT_MAX_SAMPLE_SIZE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fold an unfolded site frequency spectrum into the minor-allele spectrum, using dadi.

    Entries 0 and n are masked by dadi and do not contribute. Entries that dadi masks
    after folding are returned as None.
    """

    dadi = require("dadi", _DADI_EXTRA)
    np = require("numpy", _DADI_EXTRA)
    limit = _limit_sample_size(max_sample_size)
    counts = _observed_counts(spectrum, limit)
    n = len(counts) - 1
    folded = dadi.Spectrum(counts).fold()
    return {
        "status": "completed",
        "operation": "fold_sfs",
        "sample_size": n,
        "spectrum": _spectrum_values(np, folded),
        "provenance": _provenance("dadi", _DADI_CONVENTION + "; minor-allele folded"),
    }


def estimate_sfs_theta(
    observed: Sequence[float],
    *,
    grid_points: int = _DEFAULT_GRID_POINTS,
    max_sample_size: int = _DEFAULT_MAX_SAMPLE_SIZE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Maximum-likelihood theta for an observed spectrum under the neutral model, using dadi.

    The estimate is dadi's ``optimal_sfs_scaling``: the observed total over the model
    total on the polymorphic classes 1 to n-1.
    """

    dadi = require("dadi", _DADI_EXTRA)
    limit = _limit_sample_size(max_sample_size)
    counts = _observed_counts(observed, limit)
    _polymorphic_counts(counts)
    n = len(counts) - 1
    data = dadi.Spectrum(counts)
    model = dadi.Demographics1D.snm(None, (n,), _grid(grid_points))
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    return {
        "status": "completed",
        "operation": "estimate_sfs_theta",
        "sample_size": n,
        "theta": _finite(theta),
        "provenance": _provenance(
            "dadi", "maximum-likelihood theta under the constant-size neutral model"
        ),
    }


def calculate_sfs_log_likelihood(
    observed: Sequence[float],
    *,
    theta: float | None = None,
    grid_points: int = _DEFAULT_GRID_POINTS,
    max_sample_size: int = _DEFAULT_MAX_SAMPLE_SIZE,
    verbose: bool = True,
) -> dict[str, Any]:
    """Poisson log-likelihood of an observed spectrum under the neutral model, using dadi.

    With ``theta`` omitted, theta is set to its maximum-likelihood value first, as in
    dadi's ``ll_multinom``. With ``theta`` given, it is used as is. The log-likelihood is
    dadi's ``Inference.ll`` over the polymorphic classes 1 to n-1.
    """

    dadi = require("dadi", _DADI_EXTRA)
    limit = _limit_sample_size(max_sample_size)
    counts = _observed_counts(observed, limit)
    _polymorphic_counts(counts)
    n = len(counts) - 1
    data = dadi.Spectrum(counts)
    model = dadi.Demographics1D.snm(None, (n,), _grid(grid_points))
    if theta is None:
        scale = float(dadi.Inference.optimal_sfs_scaling(model, data))
        source = "maximum likelihood"
    else:
        scale = _number(theta, "theta", positive=True)
        source = "given"
    value = dadi.Inference.ll(model * scale, data)
    return {
        "status": "completed",
        "operation": "calculate_sfs_log_likelihood",
        "sample_size": n,
        "theta": scale,
        "theta_source": source,
        "log_likelihood": _finite(value),
        "provenance": _provenance(
            "dadi", "Poisson log-likelihood, dadi Inference.ll, classes 1 to n-1"
        ),
    }


__all__ = [
    "build_tree_sequence",
    "calculate_diversity",
    "calculate_fst",
    "calculate_sfs_log_likelihood",
    "calculate_site_frequency_spectrum",
    "calculate_tajimas_d",
    "estimate_sfs_theta",
    "expected_neutral_sfs",
    "fold_sfs",
    "read_tree_sequence",
    "write_tree_sequence",
]

_register_public_verbose(globals(), tuple(__all__))
