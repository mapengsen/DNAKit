"""tskit wrappers for family-variant rows that no existing dnakit wrapper reaches (group H1).

Two tskit functions have a DNA-level result that ``dnakit.population.tskit_functions`` does not
give: the topology counts of one chosen tree (``tskit.combinatorics.tree_count_topologies``) and
the NEXUS encoding that includes the tree blocks (``tskit.write_nexus`` with ``include_trees``).

Inputs are the explicit table mappings of ``build_tree_sequence()``. Output files are never
overwritten, no network is used, and tskit is imported only when a function runs.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.population.backends import (
    _DEFAULT_MAX_ROWS,
    _TSKIT_EXTRA,
    _provenance,
    _sample_sets,
    _text,
    _tree_sequence,
)
from dnakit.population.tskit_functions import _jsonable, _new_file, _precision, _tree_index


def count_tree_topologies(
    tables: Mapping[str, Any],
    sample_sets: Mapping[str, Sequence[int]],
    *,
    tree_index: int = 0,
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the embedded topologies of named sample sets in ONE tree, using tskit.

    ``sample_sets`` maps names to lists of sample node IDs; every set needs at least two
    samples. ``tree_index`` picks the tree (default the first). Counts are keyed by the set
    names joined with ``+`` and, within a key, by tskit's topology rank (``shape`` and ``label``),
    sorted by rank. The counts equal the entry of ``count_sample_topologies()`` for the same tree,
    but only this tree is evaluated.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if sample_sets is None:
        raise ConfigurationError("sample_sets is required, for example {'A': [0, 1], 'B': [2, 3]}.")
    ts = _tree_sequence(tskit, tables, max_rows)
    sets = _sample_sets(ts, sample_sets)
    index = _tree_index(tree_index, ts.num_trees)
    names = list(sets)
    groups = [sets[name] for name in names]
    counter = tskit.combinatorics.tree_count_topologies(ts.at_index(index), groups)
    counts: dict[str, list[dict[str, Any]]] = {}
    for key in sorted(counter.topologies):
        label = "+".join(names[position] for position in key)
        ordered = sorted(counter.topologies[key].items(), key=lambda p: (p[0].shape, p[0].label))
        counts[label] = [
            {"shape": _jsonable(rank.shape), "label": _jsonable(rank.label), "count": int(n)}
            for rank, n in ordered
        ]
    tree = ts.at_index(index)
    return {
        "status": "completed",
        "operation": "count_tree_topologies",
        "tree_index": index,
        "interval": [float(tree.interval.left), float(tree.interval.right)],
        "sample_sets": {name: len(sets[name]) for name in names},
        "counts": counts,
        "provenance": _provenance("tskit", "tskit tree_count_topologies for one tree"),
    }


def export_tree_sequence_nexus(
    tables: Mapping[str, Any],
    destination: str | Path,
    *,
    include_trees: bool = True,
    include_alignments: bool | None = None,
    precision: int | None = None,
    missing_data_character: str = "?",
    max_rows: int = _DEFAULT_MAX_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a tree sequence as a NEXUS file with its TREES block, using tskit's ``write_nexus()``.

    ``include_trees=True`` writes one Newick tree per marginal tree (named ``t<left>^<right>``);
    ``include_alignments`` adds the sample haplotypes (``None`` adds them when the tree sequence
    has at least one site). ``export_sample_sequences(file_format="nexus")`` writes the
    alignments only. Trees with several roots cannot be written. The destination must not exist.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if not isinstance(include_trees, bool):
        raise ConfigurationError("include_trees must be true or false.")
    if include_alignments is not None and not isinstance(include_alignments, bool):
        raise ConfigurationError("include_alignments must be true, false or None.")
    if not include_trees and include_alignments is False:
        raise ConfigurationError("include_trees and include_alignments cannot both be false.")
    digits = None if precision is None else _precision(precision)
    if len(_text(missing_data_character, "missing_data_character")) != 1:
        raise ConfigurationError("missing_data_character must be a single character.")
    ts = _tree_sequence(tskit, tables, max_rows)
    target = _new_file(destination, "destination")
    try:
        tskit.write_nexus(
            ts,
            str(target),
            precision=digits,
            include_trees=include_trees,
            include_alignments=include_alignments,
            reference_sequence=None,
            missing_data_character=missing_data_character,
            node_labels=None,
        )
    except (ValueError, tskit.TskitException) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"tskit could not write NEXUS: {exc}") from exc
    text = target.read_text(encoding="utf8")
    return {
        "status": "completed",
        "operation": "export_tree_sequence_nexus",
        "destination": str(target),
        "bytes": target.stat().st_size,
        "num_samples": int(ts.num_samples),
        "num_trees": int(ts.num_trees),
        "num_sites": int(ts.num_sites),
        "has_trees_block": "BEGIN TREES;" in text,
        "has_data_block": "BEGIN DATA;" in text,
        "provenance": _provenance(
            "tskit", "NEXUS; trees named t<left>^<right>; taxa n<node id>; missing data character"
        ),
    }


__all__ = [
    "count_tree_topologies",
    "export_tree_sequence_nexus",
]

_register_public_verbose(globals(), tuple(__all__))
