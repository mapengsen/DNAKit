"""tskit tree-sequence editing, statistics and export functions (integration pack P31).

Every function takes a tree sequence either as the explicit table rows of ``build_tree_sequence()``
(a mapping with ``sequence_length``, ``nodes``, ``edges``, ``sites``, ``mutations``) or as the path of
a tskit ``.trees`` file. Functions that edit a tree sequence return its summary and the edited
tables in the same row format, so the result can be passed on to the next function. Intervals are
half-open ``[left, right)``. tskit is imported only when a function runs; no network is used.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import (
    _DEFAULT_MAX_BYTES,
    _DEFAULT_MAX_ROWS,
    _MAX_BYTES_LIMIT,
    _TABLE_CONVENTION,
    _TSKIT_EXTRA,
    _build_tables,
    _integer,
    _limit,
    _number,
    _provenance,
    _summary,
    _text,
    _tree_sequence,
)

_MODES = ("site", "branch", "node")
_STAT_MODES = ("site", "branch")
_SIMPLE_STATS = ("Y1", "Y2", "Y3", "f4")
_SUMMARIES = ("identity", "square", "product")
_TRAIT_STATS = ("correlation", "covariance", "linear_model")
_PAIR_OUTPUTS = ("counts", "quantiles", "rates")
_LD_OUTPUTS = ("matrix", "pair", "array")
_RATE_QUANTITIES = ("rate", "cumulative_mass")
_TRIM_SIDES = ("both", "left", "right")
_MUTATION_TIME_METHODS = ("tables", "min")
_MAP_METHODS = ("parsimony",)
_TreeInput = Mapping[str, Any] | str | Path


def _tskit() -> Any:
    return require("tskit", _TSKIT_EXTRA)


def _guard(call: Any, what: str) -> Any:
    """Run ``call()`` and turn tskit failures into ``ConfigurationError``."""

    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:  # tskit raises TskitException, ValueError and _tskit.LibraryError
        raise ConfigurationError(f"tskit could not run {what}: {exc}") from exc


def _load(tskit: Any, data: object, max_rows: int, max_bytes: int = _DEFAULT_MAX_BYTES) -> Any:
    """Tree sequence from explicit table rows or from a ``.trees`` file."""

    if isinstance(data, (str, Path)):
        path = input_path(data)
        limit = _limit(max_bytes, "max_bytes", _MAX_BYTES_LIMIT)
        if path.stat().st_size > limit:
            raise ConfigurationError(f"{path.name} is larger than max_bytes={limit}.")
        return _guard(lambda: tskit.load(str(path)), f"reading {path.name}")
    return _tree_sequence(tskit, data, max_rows)


def _plain(value: Any) -> Any:
    """numpy values to plain JSON-serialisable Python values; NaN becomes None."""

    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _table_rows(tskit: Any, ts: Any) -> dict[str, Any]:
    """Tables of ``ts`` in the row format accepted by ``build_tree_sequence()``."""

    nodes = [
        {"time": float(n.time), "sample": bool(n.flags & tskit.NODE_IS_SAMPLE)} for n in ts.nodes()
    ]
    edges = [
        {"left": float(e.left), "right": float(e.right), "parent": int(e.parent), "child": int(e.child)}
        for e in ts.edges()
    ]
    sites = [
        {"position": float(s.position), "ancestral_state": s.ancestral_state} for s in ts.sites()
    ]
    mutations = []
    for m in ts.mutations():
        row: dict[str, Any] = {
            "site": int(m.site),
            "node": int(m.node),
            "derived_state": m.derived_state,
            "parent": int(m.parent),
        }
        if not tskit.is_unknown_time(m.time):
            row["time"] = float(m.time)
        mutations.append(row)
    return {
        "sequence_length": float(ts.sequence_length),
        "nodes": nodes,
        "edges": edges,
        "sites": sites,
        "mutations": mutations,
    }


def _edited(op: str, tskit: Any, ts: Any, convention: str, **extra: Any) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": op,
        **extra,
        "result_summary": _summary(ts),
        "tables": _table_rows(tskit, ts),
        "provenance": _provenance("tskit", f"{_TABLE_CONVENTION}; {convention}"),
    }


def _report(op: str, convention: str, **fields: Any) -> dict[str, Any]:
    return {
        "status": "completed",
        "operation": op,
        **fields,
        "provenance": _provenance("tskit", convention),
    }


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _choice(value: object, allowed: Sequence[str], name: str) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {list(allowed)}; got {value!r}.")
    return value


def _nodes(value: object, name: str, *, minimum: int = 0) -> list[int]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ConfigurationError(f"{name} must be a list of node IDs.")
    nodes = [_integer(v, f"node in {name}") for v in value]
    if len(nodes) < minimum:
        raise ConfigurationError(f"{name} needs at least {minimum} node(s).")
    return nodes


def _node_sets(ts: Any, sets: object, name: str = "sample_sets") -> list[list[int]]:
    """Lists of node IDs; the default is one set with every sample."""

    if sets is None:
        return [[int(n) for n in ts.samples()]]
    if isinstance(sets, (str, bytes)) or not isinstance(sets, Sequence) or not sets:
        raise ConfigurationError(f"{name} must be a non-empty list of lists of node IDs.")
    return [_nodes(s, f"{name} entry", minimum=1) for s in sets]


def _indexes(value: object, width: int, count: int) -> list[list[int]] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("indexes must be a non-empty list of index tuples.")
    rows = []
    for row in value:
        items = _nodes(row, "indexes entry")
        if len(items) != width:
            raise ConfigurationError(f"each indexes entry must have {width} sample-set indexes.")
        if any(i < 0 or i >= count for i in items):
            raise ConfigurationError("an indexes entry points outside sample_sets.")
        rows.append(items)
    return rows


def _windows(value: object) -> list[float] | None:
    if value is None:
        return None
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or len(value) < 2:
        raise ConfigurationError("windows must be an increasing list of at least two breakpoints.")
    return [_number(v, "window breakpoint", minimum=0.0) for v in value]


def _intervals(value: object, name: str = "intervals") -> list[list[float]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of [left, right] pairs.")
    result = []
    for pair in value:
        if isinstance(pair, (str, bytes)) or not isinstance(pair, Sequence) or len(pair) != 2:
            raise ConfigurationError(f"each entry of {name} must be a [left, right] pair.")
        result.append([_number(pair[0], "interval left", minimum=0.0),
                       _number(pair[1], "interval right", minimum=0.0)])
    return result


def _matrix(value: object, name: str, rows: int | None = None) -> list[list[float]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of numeric rows.")
    result = []
    for row in value:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or not row:
            raise ConfigurationError(f"{name} must be a list of non-empty numeric lists.")
        result.append([_number(v, f"value in {name}", minimum=-math.inf) for v in row])
    if len({len(r) for r in result}) != 1:
        raise ConfigurationError(f"all rows of {name} must have the same length.")
    if rows is not None and len(result) != rows:
        raise ConfigurationError(f"{name} needs one row per sample ({rows}); got {len(result)}.")
    return result


def _site_ids(ts: Any, ids: object) -> list[int]:
    items = _nodes(ids, "site_ids", minimum=1)
    bad = [i for i in items if i < 0 or i >= ts.num_sites]
    if bad:
        raise ConfigurationError(f"site IDs outside 0..{ts.num_sites - 1}: {bad}.")
    return items


def _stat_options(
    windows: object, mode: str, span_normalise: bool, modes: Sequence[str] = _STAT_MODES
) -> dict[str, Any]:
    return {
        "windows": _windows(windows),
        "mode": _choice(mode, modes, "mode"),
        "span_normalise": _flag(span_normalise, "span_normalise"),
    }


def _stat_result(op: str, ts: Any, values: Any, options: Mapping[str, Any], note: str, **extra: Any):
    windows = options["windows"]
    return _report(
        op,
        f"{note}; windows are [breakpoint_i, breakpoint_i+1); mode={options['mode']}; "
        f"span_normalise={options['span_normalise']}",
        mode=options["mode"],
        windows=windows,
        span_normalise=options["span_normalise"],
        values=_plain(values),
        **extra,
    )


# ---------------------------------------------------------------------------------------------
# genealogy and mutation tables
# ---------------------------------------------------------------------------------------------


def tskit_link_ancestors(
    tables: _TreeInput,
    *,
    samples: Sequence[int],
    ancestors: Sequence[int],
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Link the ancestral segments of ``ancestors`` to the ``samples`` that inherit them, using tskit.

    Runs ``TableCollection.link_ancestors()`` (also reachable as the tree-sequence method). The
    result is an edge list ``{left, right, parent, child}`` where ``parent`` is an ancestor node and
    ``child`` a sample node.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    sample_ids = _nodes(samples, "samples", minimum=1)
    ancestor_ids = _nodes(ancestors, "ancestors", minimum=1)

    def run() -> Any:
        return ts.dump_tables().link_ancestors(sample_ids, ancestor_ids)

    edges = _guard(run, "link_ancestors")
    rows = [
        {"left": float(e.left), "right": float(e.right), "parent": int(e.parent), "child": int(e.child)}
        for e in edges
    ]
    return _report(
        "tskit_link_ancestors",
        "edges map ancestor nodes to the sample nodes that inherit their genome; half-open intervals",
        num_edges=len(rows),
        edges=rows,
    )


def tskit_compute_mutation_parents(
    tables: _TreeInput, *, max_rows: int = _DEFAULT_MAX_ROWS, verbose: bool = True
) -> dict[str, Any]:
    """Compute the parent mutation of every mutation (``TableCollection.compute_mutation_parents``).

    The mutation parent is the closest mutation above it at the same site; ``-1`` means none. The
    result gives the parent array and the tables with the parents filled in.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)

    def run() -> Any:
        tc = ts.dump_tables()
        tc.compute_mutation_parents()
        return tc.tree_sequence()

    out = _guard(run, "compute_mutation_parents")
    return _edited(
        "tskit_compute_mutation_parents",
        tskit,
        out,
        "parent -1 means no parent mutation",
        mutation_parents=[int(m.parent) for m in out.mutations()],
    )


def tskit_compute_mutation_times(
    tables: _TreeInput,
    *,
    method: str = "tables",
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute or impute mutation times from node times, using tskit.

    ``method="tables"`` runs ``TableCollection.compute_mutation_times()``, which sets every mutation
    time (unknown times are placed between the node time and its parent's time and known times are
    overwritten). ``method="min"`` runs ``TreeSequence.impute_unknown_mutations_time()``: unknown
    times take the time of the node below the mutation and known times are kept; the result lists
    the times without changing the tables.
    """

    tskit = _tskit()
    _choice(method, _MUTATION_TIME_METHODS, "method")
    ts = _load(tskit, tables, max_rows)
    if method == "min":
        times = _guard(lambda: ts.impute_unknown_mutations_time(method="min"), "impute")
        return _report(
            "tskit_compute_mutation_times",
            "method=min: unknown mutation time = time of the node below the mutation",
            method=method,
            mutation_times=_plain(times),
        )

    def run() -> Any:
        tc = ts.dump_tables()
        tc.compute_mutation_times()
        return tc.tree_sequence()

    out = _guard(run, "compute_mutation_times")
    return _edited(
        "tskit_compute_mutation_times",
        tskit,
        out,
        "method=tables: times computed by tskit for all mutations",
        method=method,
        mutation_times=_plain([m.time for m in out.mutations()]),
    )


def tskit_deduplicate_sites(
    tables: _TreeInput, *, max_rows: int = _DEFAULT_MAX_ROWS, verbose: bool = True
) -> dict[str, Any]:
    """Merge sites at identical positions (``TableCollection.deduplicate_sites``).

    The first site of each position is kept and its mutations absorb those of the duplicates; mutation parents are recomputed afterwards.
    A tree sequence cannot hold duplicate sites, so the input is read as tables (explicit rows,
    which are sorted, or the tables of a ``.trees`` file) before deduplication.
    """

    tskit = _tskit()
    if isinstance(tables, (str, Path)):
        path = input_path(tables)
        tc = _guard(lambda: tskit.TableCollection.load(str(path)), f"reading {path.name}")
    else:
        tc = _build_tables(tskit, tables, max_rows)
    before = int(tc.sites.num_rows)

    def run() -> Any:
        tc.deduplicate_sites()
        tc.build_index()
        tc.compute_mutation_parents()  # merged sites need consistent mutation parents
        return tc.tree_sequence()

    out = _guard(run, "deduplicate_sites")
    return _edited(
        "tskit_deduplicate_sites",
        tskit,
        out,
        "sites at one position are merged",
        sites_before=before,
        sites_removed=before - int(out.num_sites),
    )


def tskit_delete_sites(
    tables: _TreeInput,
    *,
    site_ids: Sequence[int],
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Delete the sites with the given IDs and their mutations (``TreeSequence.delete_sites``)."""

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    ids = _site_ids(ts, site_ids)
    out = _guard(lambda: ts.delete_sites(ids, record_provenance=False), "delete_sites")
    return _edited("tskit_delete_sites", tskit, out, "site IDs renumbered after deletion", deleted=sorted(set(ids)))


def tskit_select_intervals(
    tables: _TreeInput,
    *,
    intervals: Sequence[Sequence[float]],
    action: str = "delete",
    simplify: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Delete or keep genomic intervals of a tree sequence, using tskit.

    ``action="delete"`` runs ``delete_intervals()`` (the genealogy in the intervals is removed and
    the remaining genome is closed up); ``action="keep"`` runs ``keep_intervals()`` (everything
    outside the intervals is removed, coordinates unchanged). ``simplify`` is tskit's flag.
    """

    tskit = _tskit()
    _choice(action, ("delete", "keep"), "action")
    flag = _flag(simplify, "simplify")
    pairs = _intervals(intervals)
    ts = _load(tskit, tables, max_rows)
    method = ts.delete_intervals if action == "delete" else ts.keep_intervals
    out = _guard(lambda: method(pairs, simplify=flag, record_provenance=False), f"{action}_intervals")
    return _edited("tskit_select_intervals", tskit, out, f"action={action}; half-open intervals", action=action)


def tskit_remove_older_genealogy(
    tables: _TreeInput,
    *,
    time: float,
    add_roots: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Remove the genealogy older than ``time``, using tskit.

    ``add_roots=True`` runs ``TreeSequence.decapitate()`` (the lineages crossing ``time`` get new
    root nodes at ``time``); ``add_roots=False`` runs ``TableCollection.delete_older()`` (edges and
    mutations at or above ``time`` are only deleted).
    """

    tskit = _tskit()
    moment = _number(time, "time", minimum=0.0)
    roots = _flag(add_roots, "add_roots")
    ts = _load(tskit, tables, max_rows)
    if roots:
        out = _guard(lambda: ts.decapitate(moment), "decapitate")
    else:
        def run() -> Any:
            tc = ts.dump_tables()
            tc.delete_older(moment)
            tc.sort()
            tc.build_index()
            return tc.tree_sequence()

        out = _guard(run, "delete_older")
    return _edited("tskit_remove_older_genealogy", tskit, out, f"time={moment}; add_roots={roots}", time=moment)


def tskit_trim_genome(
    tables: _TreeInput,
    *,
    side: str = "both",
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim empty flanks of the genome, using tskit.

    ``side="left"`` runs ``ltrim()`` (shift the first tree start to 0), ``side="right"`` runs
    ``rtrim()`` (set the sequence length to the end of the last tree) and ``side="both"`` runs
    ``trim()`` (both). Sites in the trimmed flanks are not allowed by tskit.
    """

    tskit = _tskit()
    _choice(side, _TRIM_SIDES, "side")
    ts = _load(tskit, tables, max_rows)
    method = {"left": ts.ltrim, "right": ts.rtrim, "both": ts.trim}[side]
    out = _guard(lambda: method(record_provenance=False), f"{side} trim")
    return _edited("tskit_trim_genome", tskit, out, f"side={side}", side=side)


def tskit_shift_coordinates(
    tables: _TreeInput,
    *,
    value: float,
    sequence_length: float | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shift all genomic coordinates by ``value`` (``TreeSequence.shift``).

    ``sequence_length`` sets the new length (default: old length plus ``value``).
    """

    tskit = _tskit()
    shift = _number(value, "value", minimum=-math.inf)
    length = None if sequence_length is None else _number(sequence_length, "sequence_length", positive=True)
    ts = _load(tskit, tables, max_rows)
    out = _guard(lambda: ts.shift(shift, sequence_length=length, record_provenance=False), "shift")
    return _edited("tskit_shift_coordinates", tskit, out, f"shift={shift}", value=shift)


def tskit_simplify(
    tables: _TreeInput,
    *,
    samples: Sequence[int] | None = None,
    reduce_to_site_topology: bool = False,
    keep_unary: bool = False,
    keep_input_roots: bool = False,
    filter_sites: bool = True,
    filter_nodes: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Simplify the genealogy to the history of ``samples`` (``TreeSequence.simplify``).

    Defaults are those of tskit. The result gives the node map (old ID to new ID, -1 if removed).
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    chosen = None if samples is None else _nodes(samples, "samples", minimum=1)
    options = {
        "reduce_to_site_topology": _flag(reduce_to_site_topology, "reduce_to_site_topology"),
        "keep_unary": _flag(keep_unary, "keep_unary"),
        "keep_input_roots": _flag(keep_input_roots, "keep_input_roots"),
        "filter_sites": _flag(filter_sites, "filter_sites"),
        "filter_nodes": _flag(filter_nodes, "filter_nodes"),
    }
    out, node_map = _guard(
        lambda: ts.simplify(chosen, map_nodes=True, record_provenance=False, **options), "simplify"
    )
    return _edited(
        "tskit_simplify", tskit, out, "node_map maps old node IDs to new IDs (-1 removed)",
        node_map=_plain(node_map), options=options,
    )


def tskit_subset(
    tables: _TreeInput,
    *,
    nodes: Sequence[int],
    reorder_populations: bool = True,
    remove_unreferenced: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep only the listed nodes and the edges, sites and mutations between them (``subset``).

    Node IDs are renumbered in the order given; ``nodes`` need not be sorted.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    ids = _nodes(nodes, "nodes", minimum=1)
    out = _guard(
        lambda: ts.subset(
            ids,
            record_provenance=False,
            reorder_populations=_flag(reorder_populations, "reorder_populations"),
            remove_unreferenced=_flag(remove_unreferenced, "remove_unreferenced"),
        ),
        "subset",
    )
    return _edited("tskit_subset", tskit, out, "nodes renumbered in the order given", kept_nodes=ids)


def tskit_union(
    tables: _TreeInput,
    other: _TreeInput,
    *,
    node_mapping: Sequence[int],
    all_edges: bool = False,
    all_mutations: bool = False,
    check_shared_equality: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add the non-shared part of ``other`` to ``tables`` (``TreeSequence.union``).

    ``node_mapping[i]`` is the node of ``tables`` equal to node ``i`` of ``other`` or ``-1`` if the
    node is new.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    second = _load(tskit, other, max_rows)
    mapping = [_integer(v, "node_mapping value", minimum=-1) for v in node_mapping]
    out = _guard(
        lambda: ts.union(
            second,
            mapping,
            check_shared_equality=_flag(check_shared_equality, "check_shared_equality"),
            record_provenance=False,
            all_edges=_flag(all_edges, "all_edges"),
            all_mutations=_flag(all_mutations, "all_mutations"),
        ),
        "union",
    )
    return _edited("tskit_union", tskit, out, "-1 in node_mapping marks nodes that are new")


def tskit_concatenate(
    tables: _TreeInput,
    others: Sequence[_TreeInput],
    *,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Join tree sequences end to end along the genome (``TreeSequence.concatenate``).

    The samples must be the same nodes in all inputs; lengths are added in the order given.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    if isinstance(others, (str, bytes)) or not isinstance(others, Sequence) or not others:
        raise ConfigurationError("others must be a non-empty list of tree sequences.")
    rest = [_load(tskit, item, max_rows) for item in others]
    out = _guard(lambda: ts.concatenate(*rest, record_provenance=False), "concatenate")
    return _edited("tskit_concatenate", tskit, out, "inputs joined in order along the genome")


def tskit_extend_haplotypes(
    tables: _TreeInput,
    *,
    max_iter: int = 10,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extend ancestral haplotype segments across unary regions (``extend_haplotypes``)."""

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    steps = _integer(max_iter, "max_iter", minimum=1, maximum=10_000)
    out = _guard(lambda: ts.extend_haplotypes(max_iter=steps), "extend_haplotypes")
    return _edited("tskit_extend_haplotypes", tskit, out, f"max_iter={steps}", max_iter=steps)


def tskit_split_edges(
    tables: _TreeInput,
    *,
    time: float,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split every edge that spans ``time`` by inserting a unary node there (``split_edges``)."""

    tskit = _tskit()
    moment = _number(time, "time", minimum=0.0)
    ts = _load(tskit, tables, max_rows)
    out = _guard(lambda: ts.split_edges(moment), "split_edges")
    return _edited("tskit_split_edges", tskit, out, f"time={moment}", time=moment)


def tskit_map_mutations(
    tables: _TreeInput,
    *,
    genotypes: Sequence[int],
    tree_index: int = 0,
    alleles: Sequence[str] | None = None,
    ancestral_state: str | int | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Place mutations on one tree that explain observed genotypes by parsimony (``Tree.map_mutations``).

    ``genotypes`` has one allele index per sample (-1 missing) in sample order, ``alleles`` names
    the indexes (default ``["0", "1", ...]``) and ``ancestral_state`` fixes the ancestral allele.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    index = _integer(tree_index, "tree_index", minimum=0, maximum=ts.num_trees - 1)
    values = [_integer(v, "genotype", minimum=-1) for v in genotypes]
    if len(values) != ts.num_samples:
        raise ConfigurationError(f"genotypes needs {ts.num_samples} values; got {len(values)}.")
    names = list(alleles) if alleles is not None else [str(i) for i in range(max(values + [1]) + 1)]
    names = [_text(a, "allele") for a in names]
    tree = ts.at_index(index)
    anc, muts = _guard(
        lambda: tree.map_mutations(values, names, ancestral_state=ancestral_state), "map_mutations"
    )
    rows = [
        {"node": int(m.node), "derived_state": m.derived_state, "parent": int(m.parent),
         "time": None if tskit.is_unknown_time(m.time) else float(m.time)}
        for m in muts
    ]
    return _report(
        "tskit_map_mutations",
        "Hartigan-style parsimony as implemented by tskit; parent indexes the returned mutations",
        tree_index=index,
        ancestral_state=anc,
        alleles=names,
        num_mutations=len(rows),
        mutations=rows,
    )


# ---------------------------------------------------------------------------------------------
# listings, genotypes, sequences
# ---------------------------------------------------------------------------------------------


def _check_listing_size(count: int, max_rows: int, what: str) -> None:
    if count > max_rows:
        raise ConfigurationError(f"{count} {what} is above max_rows={max_rows}; raise the limit.")


def tskit_list_mutations(
    tables: _TreeInput,
    *,
    tree_index: int | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List mutations of the whole tree sequence (``TreeSequence.mutations``) or one tree (``Tree.mutations``)."""

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    if tree_index is None:
        source = ts.mutations()
    else:
        source = ts.at_index(_integer(tree_index, "tree_index", minimum=0, maximum=ts.num_trees - 1)).mutations()
    rows = []
    for m in source:
        rows.append({
            "id": int(m.id), "site": int(m.site), "node": int(m.node),
            "derived_state": m.derived_state, "parent": int(m.parent),
            "time": None if tskit.is_unknown_time(m.time) else float(m.time),
            "position": float(ts.site(m.site).position),
        })
    _check_listing_size(len(rows), max_rows, "mutations")
    return _report("tskit_list_mutations", "position is the position of the mutation's site",
                   tree_index=tree_index, num_mutations=len(rows), mutations=rows)


def tskit_list_sites(
    tables: _TreeInput,
    *,
    tree_index: int | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """List variant sites of the whole tree sequence (``TreeSequence.sites``) or one tree (``Tree.sites``)."""

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    if tree_index is None:
        source = ts.sites()
    else:
        source = ts.at_index(_integer(tree_index, "tree_index", minimum=0, maximum=ts.num_trees - 1)).sites()
    rows = [
        {"id": int(s.id), "position": float(s.position), "ancestral_state": s.ancestral_state,
         "num_mutations": len(s.mutations)}
        for s in source
    ]
    _check_listing_size(len(rows), max_rows, "sites")
    return _report("tskit_list_sites", "positions are in genome coordinates",
                   tree_index=tree_index, num_sites=len(rows), sites=rows)


def tskit_genotype_matrix(
    tables: _TreeInput,
    *,
    samples: Sequence[int] | None = None,
    isolated_as_missing: bool = True,
    alleles: Sequence[str] | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotype matrix (sites by samples) of allele indexes (``TreeSequence.genotype_matrix``).

    ``-1`` marks missing data. ``alleles`` fixes the allele order for every site (default: the
    ancestral allele first, then derived alleles in order of appearance).
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    chosen = None if samples is None else _nodes(samples, "samples", minimum=1)
    fixed = None if alleles is None else tuple(_text(a, "allele") for a in alleles)
    matrix = _guard(
        lambda: ts.genotype_matrix(
            samples=chosen,
            isolated_as_missing=_flag(isolated_as_missing, "isolated_as_missing"),
            alleles=fixed,
        ),
        "genotype_matrix",
    )
    _check_listing_size(int(matrix.size), max_rows * 10, "genotype cells")
    return _report(
        "tskit_genotype_matrix", "rows are sites in order, columns are samples in sample order",
        num_sites=int(matrix.shape[0]), num_samples=int(matrix.shape[1]),
        positions=_plain(ts.tables.sites.position), matrix=_plain(matrix),
    )


def tskit_variants(
    tables: _TreeInput,
    *,
    samples: Sequence[int] | None = None,
    isolated_as_missing: bool = True,
    alleles: Sequence[str] | None = None,
    left: float | None = None,
    right: float | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-site genotype calls of a tree sequence (``TreeSequence.variants`` / ``Variant.decode``).

    Each record gives the site, position, allele list, the genotype index per sample, the allele
    state per sample (``Variant.states()``) and the allele counts (``Variant.counts()``).
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    chosen = None if samples is None else _nodes(samples, "samples", minimum=1)
    fixed = None if alleles is None else tuple(_text(a, "allele") for a in alleles)
    lo = None if left is None else _number(left, "left", minimum=0.0)
    hi = None if right is None else _number(right, "right", minimum=0.0)

    def run() -> list[dict[str, Any]]:
        out = []
        for v in ts.variants(
            samples=chosen, isolated_as_missing=_flag(isolated_as_missing, "isolated_as_missing"),
            alleles=fixed, left=lo, right=hi,
        ):
            out.append({
                "site": int(v.site.id), "position": float(v.site.position),
                "alleles": list(v.alleles), "genotypes": _plain(v.genotypes),
                "states": list(v.states()),
                "counts": {("missing" if k is None else k): int(c) for k, c in v.counts().items()},
            })
        return out

    records = _guard(run, "variants")
    _check_listing_size(len(records), max_rows, "variants")
    return _report(
        "tskit_variants", "genotype -1 and state null mean missing; counts use allele strings",
        num_variants=len(records), variants=records,
    )


def tskit_alignments(
    tables: _TreeInput,
    *,
    reference_sequence: str,
    missing_data_character: str = "N",
    isolated_as_missing: bool = True,
    samples: Sequence[int] | None = None,
    left: int | None = None,
    right: int | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Full-length sample sequence alignments over a reference (``TreeSequence.alignments``).

    ``reference_sequence`` supplies the bases at positions without a site; its length must equal
    the sequence length. Returns one string per sample (``left``/``right`` select a sub-region).
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    ref = _text(reference_sequence, "reference_sequence")
    miss = _text(missing_data_character, "missing_data_character")
    if len(miss) != 1:
        raise ConfigurationError("missing_data_character must be a single character.")
    chosen = None if samples is None else _nodes(samples, "samples", minimum=1)
    lo = None if left is None else _integer(left, "left")
    hi = None if right is None else _integer(right, "right")
    sequences = _guard(
        lambda: list(ts.alignments(
            reference_sequence=ref, missing_data_character=miss,
            isolated_as_missing=_flag(isolated_as_missing, "isolated_as_missing"),
            samples=chosen, left=lo, right=hi,
        )),
        "alignments",
    )
    ids = chosen if chosen is not None else [int(n) for n in ts.samples()]
    return _report(
        "tskit_alignments", "one string per sample over [left, right); positions without a site come from the reference",
        num_samples=len(sequences), sample_nodes=ids,
        alignments=[{"node": n, "sequence": s} for n, s in zip(ids, sequences, strict=True)],
    )


def tskit_haplotypes(
    tables: _TreeInput,
    *,
    missing_data_character: str = "-",
    isolated_as_missing: bool = True,
    samples: Sequence[int] | None = None,
    left: float | None = None,
    right: float | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Haplotype string per sample over the variable sites only (``TreeSequence.haplotypes``).

    Requires single-character alleles. Each string has one character per site in ``[left, right)``.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    miss = _text(missing_data_character, "missing_data_character")
    chosen = None if samples is None else _nodes(samples, "samples", minimum=1)
    lo = None if left is None else _number(left, "left", minimum=0.0)
    hi = None if right is None else _number(right, "right", minimum=0.0)
    out = _guard(
        lambda: list(ts.haplotypes(
            isolated_as_missing=_flag(isolated_as_missing, "isolated_as_missing"),
            missing_data_character=miss, samples=chosen, left=lo, right=hi,
        )),
        "haplotypes",
    )
    ids = chosen if chosen is not None else [int(n) for n in ts.samples()]
    return _report(
        "tskit_haplotypes", "one character per site, in site order",
        num_samples=len(out), haplotypes=[{"node": n, "haplotype": h} for n, h in zip(ids, out, strict=True)],
    )


def tskit_breakpoints(
    tables: _TreeInput, *, max_rows: int = _DEFAULT_MAX_ROWS, verbose: bool = True
) -> dict[str, Any]:
    """Genomic breakpoints between genealogies (``TreeSequence.breakpoints``), including 0 and the length."""

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    points = _plain(ts.breakpoints(as_array=True))
    return _report("tskit_breakpoints", "first value 0, last value the sequence length",
                   num_trees=int(ts.num_trees), breakpoints=points)


def tskit_map_to_vcf_model(
    tables: _TreeInput,
    *,
    ploidy: int | None = None,
    contig_id: str = "1",
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map sample nodes to VCF individuals (``TreeSequence.map_to_vcf_model``).

    Without individuals in the tree sequence ``ploidy`` is required and nodes are grouped in
    sample order. Returns the individual names, the node IDs of each individual (``individual_nodes``), the contig
    length and the contig ID.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    ploid = None if ploidy is None else _integer(ploidy, "ploidy", minimum=1, maximum=64)
    model = _guard(lambda: ts.map_to_vcf_model(ploidy=ploid, contig_id=_text(contig_id, "contig_id")), "map_to_vcf_model")
    return _report(
        "tskit_map_to_vcf_model", "VCF positions are integer genome coordinates",
        contig_id=model.contig_id, contig_length=int(model.contig_length),
        individual_names=list(model.individuals_name),
        individual_nodes=_plain(model.individuals_nodes),
    )


def tskit_sample_nodes_by_ploidy(
    tables: _TreeInput,
    *,
    ploidy: int,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group the sample nodes into rows of ``ploidy`` nodes (``TreeSequence.sample_nodes_by_ploidy``)."""

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    p = _integer(ploidy, "ploidy", minimum=1, maximum=64)
    groups = _guard(lambda: ts.sample_nodes_by_ploidy(p), "sample_nodes_by_ploidy")
    return _report("tskit_sample_nodes_by_ploidy", "each row is one individual", ploidy=p, groups=_plain(groups))


# ---------------------------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------------------------


def tskit_ld_r2(
    tables: _TreeInput,
    *,
    output: str = "matrix",
    site_a: int | None = None,
    site_b: int | None = None,
    focal_site: int | None = None,
    direction: int = 1,
    max_mutations: int | None = None,
    max_distance: float | None = None,
    max_sites: int | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Linkage disequilibrium r2 between sites, using ``tskit.LdCalculator``.

    ``output="pair"`` gives r2 of ``site_a`` and ``site_b`` (``r2``); ``output="array"`` gives r2
    between ``focal_site`` and the following (``direction=1``) or preceding (``-1``) sites, bounded
    by ``max_mutations``, ``max_distance`` and ``max_sites`` (``r2_array``); ``output="matrix"``
    gives r2 for all pairs of sites (``r2_matrix``).
    """

    tskit = _tskit()
    _choice(output, _LD_OUTPUTS, "output")
    ts = _load(tskit, tables, max_rows)
    calc = tskit.LdCalculator(ts)
    if output == "matrix":
        values = _guard(calc.r2_matrix, "r2_matrix")
        return _report("tskit_ld_r2", "r2 matrix over all sites, 1.0 on the diagonal",
                       output=output, num_sites=int(ts.num_sites), values=_plain(values))
    if output == "pair":
        if site_a is None or site_b is None:
            raise ConfigurationError("output='pair' needs site_a and site_b.")
        a = _integer(site_a, "site_a", maximum=max(ts.num_sites - 1, 0))
        b = _integer(site_b, "site_b", maximum=max(ts.num_sites - 1, 0))
        return _report("tskit_ld_r2", "r2 of two sites", output=output, site_a=a, site_b=b,
                       values=_plain(_guard(lambda: calc.r2(a, b), "r2")))
    if focal_site is None:
        raise ConfigurationError("output='array' needs focal_site.")
    focal = _integer(focal_site, "focal_site", maximum=max(ts.num_sites - 1, 0))
    if direction not in (1, -1):
        raise ConfigurationError("direction must be 1 or -1.")
    options = {
        "max_mutations": None if max_mutations is None else _integer(max_mutations, "max_mutations", minimum=1),
        "max_distance": None if max_distance is None else _number(max_distance, "max_distance", minimum=0.0),
        "max_sites": None if max_sites is None else _integer(max_sites, "max_sites", minimum=1),
    }
    values = _guard(lambda: calc.r2_array(focal, direction=direction, **options), "r2_array")
    return _report("tskit_ld_r2", "r2 between the focal site and the sites in the chosen direction",
                   output=output, focal_site=focal, direction=direction, values=_plain(values))


def tskit_sample_set_statistic(
    tables: _TreeInput,
    *,
    statistic: str,
    sample_sets: Sequence[Sequence[int]] | None = None,
    indexes: Sequence[Sequence[int]] | None = None,
    windows: Sequence[float] | None = None,
    mode: str = "site",
    span_normalise: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the tskit multi-sample-set statistic Y1, Y2, Y3 or f4.

    ``statistic`` selects ``Y1`` (one set), ``Y2`` (pairs of sets), ``Y3`` (triples) or ``f4`` (four
    sets). ``indexes`` lists tuples of ``sample_sets`` positions (not used by Y1; Y2 needs pairs,
    Y3 triples, f4 four-tuples; default: all of them when ``sample_sets`` has exactly that many).
    """

    tskit = _tskit()
    _choice(statistic, _SIMPLE_STATS, "statistic")
    ts = _load(tskit, tables, max_rows)
    sets = _node_sets(ts, sample_sets)
    options = _stat_options(windows, mode, span_normalise)
    width = {"Y2": 2, "Y3": 3, "f4": 4}.get(statistic)
    method = getattr(ts, statistic)
    if width is None:
        values = _guard(lambda: method(sets, **options), statistic)
        idx = None
    else:
        idx = _indexes(indexes, width, len(sets))
        values = _guard(lambda: method(sets, indexes=idx, **options), statistic)
    return _stat_result("tskit_sample_set_statistic", ts, values, options,
                        f"statistic={statistic}; one column per sample set or index tuple",
                        statistic=statistic, indexes=idx)


def tskit_segregating_sites(
    tables: _TreeInput,
    *,
    sample_sets: Sequence[Sequence[int]] | None = None,
    windows: Sequence[float] | None = None,
    mode: str = "site",
    span_normalise: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Segregating site density per window and sample set (``TreeSequence.segregating_sites``).

    ``span_normalise=True`` (tskit's default) divides by the window span; use ``False`` for counts.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    sets = _node_sets(ts, sample_sets)
    options = _stat_options(windows, mode, span_normalise)
    values = _guard(lambda: ts.segregating_sites(sets, **options), "segregating_sites")
    return _stat_result("tskit_segregating_sites", ts, values, options, "one value per window and sample set")


def tskit_general_stat(
    tables: _TreeInput,
    *,
    weights: Sequence[Sequence[float]] | None = None,
    sample_sets: Sequence[Sequence[int]] | None = None,
    summary: str = "identity",
    polarised: bool = False,
    strict: bool = True,
    windows: Sequence[float] | None = None,
    mode: str = "site",
    span_normalise: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Windowed statistic from sample weights or sample-set counts (``general_stat`` / ``sample_count_stat``).

    Pass ``weights`` (one row per sample) to run ``TreeSequence.general_stat()`` or ``sample_sets``
    to run ``TreeSequence.sample_count_stat()`` (counts of each set below a branch or allele).
    tskit takes a Python function as the summary; this function offers the named summaries
    ``identity`` f(x)=x, ``square`` f(x)=x*x and ``product`` f(x)=prod(x) (one output per summary). With ``strict=True`` (tskit's default) tskit
    requires the summary to be zero for zero weight and for the total weight, so centre the weights
    (columns summing to zero) or pass ``strict=False``.
    """

    tskit = _tskit()
    np = require("numpy", _TSKIT_EXTRA)
    _choice(summary, _SUMMARIES, "summary")
    if (weights is None) == (sample_sets is None):
        raise ConfigurationError("give exactly one of weights and sample_sets.")
    ts = _load(tskit, tables, max_rows)
    options = _stat_options(windows, mode, span_normalise)
    pol = _flag(polarised, "polarised")
    strict = _flag(strict, "strict")
    func = {
        "identity": lambda x: x,
        "square": lambda x: x * x,
        "product": lambda x: np.array([np.prod(x)]),
    }[summary]
    if weights is not None:
        w = np.array(_matrix(weights, "weights", rows=ts.num_samples), dtype=float)
        dim = 1 if summary == "product" else w.shape[1]
        values = _guard(lambda: ts.general_stat(w, func, dim, polarised=pol, strict=strict, **options), "general_stat")
        source = "weights"
    else:
        sets = _node_sets(ts, sample_sets)
        dim = 1 if summary == "product" else len(sets)
        values = _guard(lambda: ts.sample_count_stat(sets, func, dim, polarised=pol, strict=strict, **options), "sample_count_stat")
        source = "sample_sets"
    return _stat_result("tskit_general_stat", ts, values, options,
                        f"input={source}; summary={summary}; polarised={pol}", input=source, summary=summary)


def tskit_genetic_relatedness(
    tables: _TreeInput,
    *,
    sample_sets: Sequence[Sequence[int]] | None = None,
    indexes: Sequence[Sequence[int]] | None = None,
    windows: Sequence[float] | None = None,
    mode: str = "site",
    span_normalise: bool = True,
    centre: bool = True,
    proportion: bool = True,
    polarised: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genetic relatedness between sample sets (``genetic_relatedness`` / ``genetic_relatedness_matrix``).

    With ``indexes`` (pairs of sample-set positions) this runs ``genetic_relatedness()`` and returns
    one value per pair; without ``indexes`` it runs ``genetic_relatedness_matrix()`` and returns the
    full sample-set by sample-set matrix (``centre``, ``proportion`` and ``polarised`` then do not apply).
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    sets = _node_sets(ts, sample_sets) if indexes is not None or sample_sets is not None else None
    options = _stat_options(windows, mode, span_normalise)
    if indexes is None:
        values = _guard(lambda: ts.genetic_relatedness_matrix(sets, **options), "genetic_relatedness_matrix")
        return _stat_result("tskit_genetic_relatedness", ts, values, options, "matrix over sample sets",
                            output="matrix", indexes=None)
    idx = _indexes(indexes, 2, len(sets))
    values = _guard(
        lambda: ts.genetic_relatedness(
            sets, indexes=idx, centre=_flag(centre, "centre"),
            proportion=_flag(proportion, "proportion"), polarised=_flag(polarised, "polarised"), **options,
        ),
        "genetic_relatedness",
    )
    return _stat_result("tskit_genetic_relatedness", ts, values, options,
                        f"centre={centre}; proportion={proportion}; polarised={polarised}",
                        output="pairs", indexes=idx)


def tskit_genetic_relatedness_vector(
    tables: _TreeInput,
    *,
    weights: Sequence[Sequence[float]],
    windows: Sequence[float] | None = None,
    span_normalise: bool = True,
    centre: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Product of the genetic relatedness matrix with a weight matrix (``genetic_relatedness_vector``).

    ``weights`` has one row per sample; the result has the same shape per window. tskit supports
    only branch mode for this method, so branch mode is used.
    """

    tskit = _tskit()
    np = require("numpy", _TSKIT_EXTRA)
    ts = _load(tskit, tables, max_rows)
    w = np.array(_matrix(weights, "weights", rows=ts.num_samples), dtype=float)
    options = _stat_options(windows, "branch", span_normalise)
    values = _guard(lambda: ts.genetic_relatedness_vector(w, centre=_flag(centre, "centre"), **options),
                    "genetic_relatedness_vector")
    return _stat_result("tskit_genetic_relatedness_vector", ts, values, options, f"centre={centre}")


def tskit_genetic_relatedness_weighted(
    tables: _TreeInput,
    *,
    weights: Sequence[Sequence[float]],
    indexes: Sequence[Sequence[int]] | None = None,
    windows: Sequence[float] | None = None,
    mode: str = "site",
    span_normalise: bool = True,
    polarised: bool = False,
    centre: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Weighted genetic relatedness (``genetic_relatedness_weighted``).

    ``weights`` has one row per sample; ``indexes`` are pairs of weight-column positions
    (default: the pair (0, 1) when there are exactly two columns).
    """

    tskit = _tskit()
    np = require("numpy", _TSKIT_EXTRA)
    ts = _load(tskit, tables, max_rows)
    w = np.array(_matrix(weights, "weights", rows=ts.num_samples), dtype=float)
    idx = _indexes(indexes, 2, w.shape[1])
    options = _stat_options(windows, mode, span_normalise)
    values = _guard(
        lambda: ts.genetic_relatedness_weighted(
            w, indexes=idx, polarised=_flag(polarised, "polarised"), centre=_flag(centre, "centre"), **options),
        "genetic_relatedness_weighted",
    )
    return _stat_result("tskit_genetic_relatedness_weighted", ts, values, options,
                        f"polarised={polarised}; centre={centre}", indexes=idx)


def tskit_trait_statistic(
    tables: _TreeInput,
    *,
    statistic: str,
    weights: Sequence[Sequence[float]],
    covariates: Sequence[Sequence[float]] | None = None,
    windows: Sequence[float] | None = None,
    mode: str = "site",
    span_normalise: bool = True,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Association of genotype with traits: correlation, covariance or linear model.

    ``statistic`` is ``correlation`` (``trait_correlation``), ``covariance`` (``trait_covariance``)
    or ``linear_model`` (``trait_linear_model``, optionally with ``covariates``). ``weights`` is
    one row per sample and one column per trait.
    """

    tskit = _tskit()
    np = require("numpy", _TSKIT_EXTRA)
    _choice(statistic, _TRAIT_STATS, "statistic")
    ts = _load(tskit, tables, max_rows)
    w = np.array(_matrix(weights, "weights", rows=ts.num_samples), dtype=float)
    options = _stat_options(windows, mode, span_normalise)
    if statistic == "linear_model":
        z = None if covariates is None else np.array(_matrix(covariates, "covariates", rows=ts.num_samples), dtype=float)
        values = _guard(lambda: ts.trait_linear_model(w, Z=z, **options), "trait_linear_model")
    else:
        if covariates is not None:
            raise ConfigurationError("covariates are only used by statistic='linear_model'.")
        method = ts.trait_correlation if statistic == "correlation" else ts.trait_covariance
        values = _guard(lambda: method(w, **options), f"trait_{statistic}")
    return _stat_result("tskit_trait_statistic", ts, values, options, f"statistic={statistic}", statistic=statistic)


def tskit_mean_descendants(
    tables: _TreeInput,
    *,
    sample_sets: Sequence[Sequence[int]],
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean number of descendant samples of each sample set for every node (``mean_descendants``).

    The result has one row per node and one column per sample set; values are averaged over the genome.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    sets = _node_sets(ts, sample_sets)
    values = _guard(lambda: ts.mean_descendants(sets), "mean_descendants")
    return _report("tskit_mean_descendants", "rows are nodes, columns are sample sets, averaged over the genome",
                   num_nodes=int(ts.num_nodes), values=_plain(values))


def tskit_genealogical_nearest_neighbours(
    tables: _TreeInput,
    *,
    focal: Sequence[int],
    sample_sets: Sequence[Sequence[int]],
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Proportion of each focal node's genealogical nearest neighbours that lie in each sample set.

    Runs ``TreeSequence.genealogical_nearest_neighbours()``; rows follow ``focal``, columns follow
    ``sample_sets`` and the sample sets must not overlap.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    nodes = _nodes(focal, "focal", minimum=1)
    sets = _node_sets(ts, sample_sets)
    values = _guard(lambda: ts.genealogical_nearest_neighbours(nodes, sets, num_threads=0), "gnn")
    return _report("tskit_genealogical_nearest_neighbours", "rows sum to 1 where neighbours exist",
                   focal=nodes, values=_plain(values))


def tskit_pair_coalescence(
    tables: _TreeInput,
    *,
    output: str = "rates",
    time_windows: Sequence[float] | str = "nodes",
    quantiles: Sequence[float] | None = None,
    sample_sets: Sequence[Sequence[int]] | None = None,
    indexes: Sequence[Sequence[int]] | None = None,
    windows: Sequence[float] | None = None,
    span_normalise: bool = True,
    pair_normalise: bool = False,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pair coalescence counts, quantiles or rates over time (tskit pair_coalescence_* methods).

    ``output="counts"`` runs ``pair_coalescence_counts`` (``time_windows`` may be ``"nodes"`` or
    breakpoints; ``span_normalise``, ``pair_normalise``); ``output="rates"`` runs
    ``pair_coalescence_rates`` (``time_windows`` is an increasing list of times; tskit needs the last
    window to end at infinity, so infinity is appended as the final breakpoint);
    ``output="quantiles"`` runs ``pair_coalescence_quantiles`` and needs ``quantiles`` in [0, 1].
    ``indexes`` are pairs of sample-set positions (default: the one pair when two sets are given).
    """

    tskit = _tskit()
    np = require("numpy", _TSKIT_EXTRA)
    _choice(output, _PAIR_OUTPUTS, "output")
    ts = _load(tskit, tables, max_rows)
    sets = _node_sets(ts, sample_sets)
    idx = _indexes(indexes, 2, len(sets))
    win = _windows(windows)
    if output == "quantiles":
        if quantiles is None:
            raise ConfigurationError("output='quantiles' needs quantiles.")
        q = [_number(v, "quantile", minimum=0.0) for v in quantiles]
        values = _guard(lambda: ts.pair_coalescence_quantiles(np.array(q), sample_sets=sets, indexes=idx, windows=win),
                        "pair_coalescence_quantiles")
        return _report("tskit_pair_coalescence", "output=quantiles; times of the requested coalescence quantiles",
                       output=output, quantiles=q, values=_plain(values))
    if output == "rates":
        if isinstance(time_windows, str):
            raise ConfigurationError("output='rates' needs time_windows as an increasing list of times.")
        tw = [_number(v, "time window breakpoint", minimum=0.0) for v in time_windows]
        tw.append(math.inf)  # tskit requires the last time window to end at infinity
        values = _guard(lambda: ts.pair_coalescence_rates(np.array(tw), sample_sets=sets, indexes=idx, windows=win),
                        "pair_coalescence_rates")
        return _report("tskit_pair_coalescence", "output=rates; rate per time window, shape windows x time windows x pairs",
                       output=output, time_windows=tw[:-1] + ["inf"], values=_plain(values))
    tw2: Any = time_windows if isinstance(time_windows, str) else [_number(v, "time window breakpoint", minimum=0.0) for v in time_windows]
    if isinstance(tw2, str) and tw2 != "nodes":
        raise ConfigurationError("time_windows must be 'nodes' or a list of times.")
    values = _guard(
        lambda: ts.pair_coalescence_counts(
            sample_sets=sets, indexes=idx, windows=win, span_normalise=_flag(span_normalise, "span_normalise"),
            pair_normalise=_flag(pair_normalise, "pair_normalise"), time_windows=tw2 if isinstance(tw2, str) else np.array(tw2)),
        "pair_coalescence_counts",
    )
    return _report("tskit_pair_coalescence", "output=counts; coalescing pairs per time window",
                   output=output, time_windows=tw2, values=_plain(values))


def tskit_ibd_segments(
    tables: _TreeInput,
    *,
    within: Sequence[int] | None = None,
    between: Sequence[Sequence[int]] | None = None,
    max_time: float | None = None,
    min_span: float | None = None,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Identity-by-descent segments from the genealogy (``TreeSequence.ibd_segments``).

    ``within`` limits to pairs inside one list of samples; ``between`` is two or more lists of
    samples and reports pairs across lists. A segment is a region where a pair has the same most
    recent common ancestor; ``max_time`` and ``min_span`` filter them.
    """

    tskit = _tskit()
    ts = _load(tskit, tables, max_rows)
    if within is not None and between is not None:
        raise ConfigurationError("give within or between, not both.")
    inside = None if within is None else _nodes(within, "within", minimum=2)
    across = None if between is None else _node_sets(ts, between, "between")
    options = {
        "max_time": None if max_time is None else _number(max_time, "max_time", minimum=0.0),
        "min_span": None if min_span is None else _number(min_span, "min_span", minimum=0.0),
    }

    def run() -> list[dict[str, Any]]:
        result = ts.ibd_segments(within=inside, between=across, store_pairs=True, store_segments=True, **options)
        rows = []
        for pair, segs in result.items():
            for seg in segs:
                rows.append({"node_a": int(pair[0]), "node_b": int(pair[1]), "left": float(seg.left),
                             "right": float(seg.right), "mrca": int(seg.node)})
        return rows

    rows = _guard(run, "ibd_segments")
    _check_listing_size(len(rows), max_rows, "IBD segments")
    return _report("tskit_ibd_segments", "segments are half-open [left, right); mrca is the common ancestor node",
                   num_segments=len(rows), segments=rows)


# ---------------------------------------------------------------------------------------------
# recombination map
# ---------------------------------------------------------------------------------------------


def _rate_rows(rate_map: Any) -> list[dict[str, Any]]:
    return [
        {"left": float(l), "right": float(r), "rate": None if not math.isfinite(float(v)) else float(v)}
        for l, r, v in zip(rate_map.left, rate_map.right, rate_map.rate, strict=True)
    ]


def _hapmap(tskit: Any, source: object, sequence_length: float | None, has_header: bool,
            position_col: int | None, rate_col: int | None, map_col: int | None) -> Any:
    path = input_path(source)

    def run() -> Any:
        with open(path) as handle:
            return tskit.RateMap.read_hapmap(
                handle, sequence_length=sequence_length, has_header=has_header,
                position_col=position_col, rate_col=rate_col, map_col=map_col)

    return _guard(run, f"read_hapmap on {path.name}")


def tskit_read_hapmap(
    source: str | Path,
    *,
    sequence_length: float | None = None,
    has_header: bool = True,
    position_col: int | None = None,
    rate_col: int | None = None,
    map_col: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a HapMap-style recombination map file (``tskit.RateMap.read_hapmap``).

    The file is whitespace-delimited with a header (``has_header``); columns are found from the
    header or given as 0-based ``position_col``, ``rate_col`` (cM/Mb) and ``map_col`` (cM).
    Returns the intervals with their rate in per-base-pair units, the total and mean rate.
    """

    tskit = _tskit()
    length = None if sequence_length is None else _number(sequence_length, "sequence_length", positive=True)
    rm = _hapmap(tskit, source, length, _flag(has_header, "has_header"), position_col, rate_col, map_col)
    return _report(
        "tskit_read_hapmap", "rates are per base pair per generation as in tskit; NaN (missing) is null",
        sequence_length=float(rm.sequence_length), total_mass=float(rm.total_mass),
        mean_rate=float(rm.mean_rate), intervals=_rate_rows(rm),
    )


def tskit_recombination_map_lookup(
    positions: Sequence[float],
    *,
    quantity: str = "rate",
    hapmap_file: str | Path | None = None,
    map_positions: Sequence[float] | None = None,
    map_rates: Sequence[float] | None = None,
    sequence_length: float | None = None,
    has_header: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up a recombination map at genome positions, using ``tskit.RateMap``.

    ``quantity="rate"`` returns the rate at each position (``get_rate``); ``quantity="cumulative_mass"``
    returns the cumulative genetic map position (``get_cumulative_mass``). The map comes from a
    HapMap file (``hapmap_file``) or from breakpoints ``map_positions`` (n + 1 values) with
    ``map_rates`` (n values).
    """

    tskit = _tskit()
    np = require("numpy", _TSKIT_EXTRA)
    _choice(quantity, _RATE_QUANTITIES, "quantity")
    x = [_number(p, "position", minimum=0.0) for p in positions]
    if not x:
        raise ConfigurationError("positions must not be empty.")
    if (hapmap_file is None) == (map_positions is None):
        raise ConfigurationError("give exactly one of hapmap_file and map_positions/map_rates.")
    length = None if sequence_length is None else _number(sequence_length, "sequence_length", positive=True)
    if hapmap_file is not None:
        rm = _hapmap(tskit, hapmap_file, length, _flag(has_header, "has_header"), None, None, None)
    else:
        if map_rates is None:
            raise ConfigurationError("map_rates is required with map_positions.")
        pos = np.array([_number(p, "map position", minimum=0.0) for p in map_positions], dtype=float)
        rates = np.array([_number(r, "map rate", minimum=0.0) for r in map_rates], dtype=float)
        rm = _guard(lambda: tskit.RateMap(position=pos, rate=rates), "RateMap")
    getter = rm.get_rate if quantity == "rate" else rm.get_cumulative_mass
    values = _guard(lambda: getter(np.array(x, dtype=float)), f"get_{quantity}")
    return _report(
        "tskit_recombination_map_lookup",
        "rate per base pair per generation; cumulative_mass is the integral of rate from 0",
        quantity=quantity, positions=x, values=_plain(values),
    )


_register_public_verbose(globals())
