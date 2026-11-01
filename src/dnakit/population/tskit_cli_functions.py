"""Command-line facade for tskit's commands on a ``.trees`` file (wave 1d, coverage gap).

Each public function performs the operation of one ``python3 -m tskit`` subcommand in ``tskit.cli``:

- ``tree_file_info()``: ``info``, the summary of the tree sequence (the text of
  ``str(TreeSequence)``).
- ``tree_file_trees()``: ``trees``, the index, interval and site count of each tree, with an
  optional unicode drawing.
- ``export_tree_file_<table>_table()``: ``individuals``, ``nodes``, ``edges``, ``sites``,
  ``mutations``, ``populations``, ``migrations`` and ``provenances``, the tabular text that
  ``TreeSequence.dump_text()`` writes for one table.
- ``export_tree_file_fasta()`` and ``export_tree_file_vcf()``: ``run_fasta()`` and ``run_vcf()``,
  which call ``TreeSequence.write_fasta()`` and ``TreeSequence.write_vcf()``. tskit 1.0.3 has the
  ``fasta`` subcommand commented out, so only the Python-level operation is available there.

The file is loaded with ``tskit.load()`` after a size check. Output files are written to new paths
and are never overwritten. No network is used, and tskit is imported only when a function runs.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.population.backends import (
    _DEFAULT_MAX_BYTES,
    _MAX_BYTES_LIMIT,
    _TSKIT_EXTRA,
    _integer,
    _limit,
    _provenance,
    _text,
)
from dnakit.population.tskit_functions import _new_file

_INFO_CONVENTION = (
    "tskit summary of the tree sequence: the text is str(TreeSequence) (tskit's info table); "
    "counts are the table sizes of the file"
)
_TREES_CONVENTION = (
    "one record per tree in genomic order; interval is [left, right) rounded to precision decimal "
    "places, as in tskit's trees command"
)
_TABLE_CONVENTION = (
    "tabular text of one table from TreeSequence.dump_text(), with a header line; "
    "row_count excludes the header"
)
_PROVENANCE_CONVENTION = (
    "tabular text of the provenance table from TreeSequence.dump_text(); with human_readable each "
    "record is printed as in tskit's provenances command (id, timestamp, pretty-printed JSON)"
)
_FASTA_CONVENTION = (
    "FASTA of sample haplotypes from tskit alignments; wrap_width=0 disables wrapping"
)
_VCF_CONVENTION = "VCF of tskit mutations; samples grouped into individuals of the given ploidy"


def _load_tree_file(tskit: Any, source: str | Path, max_bytes: int) -> tuple[Path, Any]:
    """Check the file size against ``max_bytes`` and load the tree sequence with tskit."""

    path = input_path(source)
    limit = _limit(max_bytes, "max_bytes", _MAX_BYTES_LIMIT)
    size = path.stat().st_size
    if size > limit:
        raise ConfigurationError(f"{path.name} has {size} bytes, above max_bytes={limit}.")
    try:
        ts = tskit.load(str(path))
    except (tskit.TskitException, tskit.FileFormatError, OSError, ValueError) as exc:
        raise ConfigurationError(f"tskit could not read {path.name}: {exc}") from exc
    return path, ts


def _precision(value: object) -> int:
    return _integer(value, "precision", minimum=0, maximum=16)


def _count_rows(path: Path) -> int:
    with path.open(encoding="utf8") as handle:
        lines = sum(1 for _ in handle)
    return max(lines - 1, 0)


def tree_file_info(
    tree_sequence_path: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarise a ``.trees`` file, as tskit's ``info`` command does.

    The result holds ``summary_text``, which is the text of ``str(TreeSequence)``, and the table
    counts: trees, samples, nodes, edges, sites, mutations, individuals, populations, migrations and
    provenances, with the sequence length.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return {
        "status": "completed",
        "operation": "tree_file_info",
        "source": str(source),
        "summary_text": str(ts),
        "sequence_length": float(ts.sequence_length),
        "num_trees": int(ts.num_trees),
        "num_samples": int(ts.num_samples),
        "num_nodes": int(ts.num_nodes),
        "num_edges": int(ts.num_edges),
        "num_sites": int(ts.num_sites),
        "num_mutations": int(ts.num_mutations),
        "num_individuals": int(ts.num_individuals),
        "num_populations": int(ts.num_populations),
        "num_migrations": int(ts.num_migrations),
        "num_provenances": int(ts.num_provenances),
        "provenance": _provenance("tskit", _INFO_CONVENTION),
    }


def tree_file_trees(
    tree_sequence_path: str | Path,
    *,
    precision: int = 6,
    draw: bool = False,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the trees of a ``.trees`` file with their index, interval and number of sites, as
    tskit's ``trees`` does.

    ``precision`` is the number of decimal places of the interval (6 in tskit's command).
    ``draw=True`` adds tskit's unicode drawing of each tree (``Tree.draw(format="unicode")``).
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    places = _precision(precision)
    if not isinstance(draw, bool):
        raise ConfigurationError("draw must be true or false.")
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    trees: list[dict[str, Any]] = []
    for tree in ts.trees():
        record: dict[str, Any] = {
            "index": int(tree.index),
            "num_sites": int(tree.num_sites),
            "interval": [
                float(f"{tree.interval.left:.{places}f}"),
                float(f"{tree.interval.right:.{places}f}"),
            ],
        }
        if draw:
            record["drawing"] = tree.draw(format="unicode")
        trees.append(record)
    return {
        "status": "completed",
        "operation": "tree_file_trees",
        "source": str(source),
        "precision": places,
        "num_trees": len(trees),
        "trees": trees,
        "provenance": _provenance("tskit", _TREES_CONVENTION),
    }


def _dump_table(
    tskit: Any,
    source: Path,
    ts: Any,
    destination: str | Path,
    table: str,
    **options: Any,
) -> dict[str, Any]:
    target = _new_file(destination, "destination")
    try:
        with target.open("w", encoding="utf8") as handle:
            ts.dump_text(**{table: handle}, **options)
    except (tskit.TskitException, ValueError, OSError) as exc:
        target.unlink(missing_ok=True)
        raise BackendExecutionError(f"tskit could not write the {table} table: {exc}") from exc
    return {
        "status": "completed",
        "operation": f"export_tree_file_{table}_table",
        "source": str(source),
        "destination": str(target),
        "table": table,
        "row_count": _count_rows(target),
        "bytes": target.stat().st_size,
        "provenance": _provenance("tskit", _TABLE_CONVENTION),
    }


def export_tree_file_individuals_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    precision: int = 6,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the individuals table of a ``.trees`` file as text (``individuals`` command)."""

    tskit = require("tskit", _TSKIT_EXTRA)
    places = _precision(precision)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return _dump_table(tskit, source, ts, destination, "individuals", precision=places)


def export_tree_file_nodes_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    precision: int = 6,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the nodes table of a ``.trees`` file as text (tskit's ``nodes`` command)."""

    tskit = require("tskit", _TSKIT_EXTRA)
    places = _precision(precision)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return _dump_table(tskit, source, ts, destination, "nodes", precision=places)


def export_tree_file_edges_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    precision: int = 6,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the edges table of a ``.trees`` file as text (tskit's ``edges`` command)."""

    tskit = require("tskit", _TSKIT_EXTRA)
    places = _precision(precision)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return _dump_table(tskit, source, ts, destination, "edges", precision=places)


def export_tree_file_sites_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    precision: int = 6,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the sites table of a ``.trees`` file as text (tskit's ``sites`` command)."""

    tskit = require("tskit", _TSKIT_EXTRA)
    places = _precision(precision)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return _dump_table(tskit, source, ts, destination, "sites", precision=places)


def export_tree_file_mutations_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    precision: int = 6,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the mutations table of a ``.trees`` file as text (tskit's ``mutations`` command)."""

    tskit = require("tskit", _TSKIT_EXTRA)
    places = _precision(precision)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return _dump_table(tskit, source, ts, destination, "mutations", precision=places)


def export_tree_file_populations_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the populations table of a ``.trees`` file as text (``populations`` command)."""

    tskit = require("tskit", _TSKIT_EXTRA)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return _dump_table(tskit, source, ts, destination, "populations")


def export_tree_file_migrations_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    precision: int = 6,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the migrations table of a ``.trees`` file as text (tskit's ``migrations`` command)."""

    tskit = require("tskit", _TSKIT_EXTRA)
    places = _precision(precision)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    return _dump_table(tskit, source, ts, destination, "migrations", precision=places)


def export_tree_file_provenances_table(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    human_readable: bool = False,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the provenance records of a ``.trees`` file (tskit's ``provenances`` command).

    The default is the tabular text of the provenance table. ``human_readable=True`` writes one
    block per record, ``id=..., timestamp=..., record=...``, with the record pretty-printed as JSON,
    as tskit's ``-H`` option does. A record that is not JSON is an error.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    if not isinstance(human_readable, bool):
        raise ConfigurationError("human_readable must be true or false.")
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    if not human_readable:
        return _dump_table(tskit, source, ts, destination, "provenances")
    target = _new_file(destination, "destination")
    blocks: list[str] = []
    try:
        for provenance in ts.provenances():
            record = json.loads(provenance.record)
            blocks.append(
                f"id={provenance.id}, timestamp={provenance.timestamp}, "
                f"record={json.dumps(record, indent=4)}"
            )
    except ValueError as exc:
        raise ConfigurationError(f"a provenance record is not JSON: {exc}") from exc
    try:
        target.write_text("".join(f"{block}\n" for block in blocks), encoding="utf8")
    except OSError as exc:
        raise BackendExecutionError(f"could not write {target.name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "export_tree_file_provenances_table",
        "source": str(source),
        "destination": str(target),
        "table": "provenances",
        "human_readable": True,
        "row_count": len(blocks),
        "bytes": target.stat().st_size,
        "provenance": _provenance("tskit", _PROVENANCE_CONVENTION),
    }


def export_tree_file_fasta(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    wrap_width: int = 60,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the sample haplotype sequences of a ``.trees`` file as FASTA, using tskit.

    This is the operation of tskit's ``run_fasta()`` handler, which calls
    ``TreeSequence.write_fasta()`` with the wrap width (60 by default, as in tskit).
    ``wrap_width=0`` writes each sequence on one line. The destination must not exist yet.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    width = _integer(wrap_width, "wrap_width", minimum=0)
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    target = _new_file(destination, "destination")
    ts.write_fasta(str(target), wrap_width=width)
    return {
        "status": "completed",
        "operation": "export_tree_file_fasta",
        "source": str(source),
        "destination": str(target),
        "file_format": "fasta",
        "wrap_width": width,
        "bytes": target.stat().st_size,
        "num_samples": int(ts.num_samples),
        "sequence_length": float(ts.sequence_length),
        "num_sites": int(ts.num_sites),
        "provenance": _provenance("tskit", _FASTA_CONVENTION),
    }


def export_tree_file_vcf(
    tree_sequence_path: str | Path,
    destination: str | Path,
    *,
    ploidy: int | None = None,
    contig_id: str = "1",
    allow_position_zero: bool = False,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the variants of a ``.trees`` file as a VCF file, using tskit's ``write_vcf()``.

    This is the operation of tskit's ``run_vcf()`` handler. ``ploidy`` null uses tskit's default;
    otherwise the number of samples must be a multiple of ``ploidy`` (samples are grouped into
    individuals). ``contig_id`` names the CHROM column. The destination must not exist yet.
    """

    tskit = require("tskit", _TSKIT_EXTRA)
    count = None if ploidy is None else _integer(ploidy, "ploidy", minimum=1)
    contig = _text(contig_id, "contig_id")
    if not isinstance(allow_position_zero, bool):
        raise ConfigurationError("allow_position_zero must be true or false.")
    source, ts = _load_tree_file(tskit, tree_sequence_path, max_bytes)
    target = _new_file(destination, "destination")
    try:
        with target.open("w", encoding="utf8") as handle:
            ts.write_vcf(
                handle,
                ploidy=count,
                contig_id=contig,
                allow_position_zero=allow_position_zero,
            )
    except (tskit.TskitException, ValueError) as exc:
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"tskit could not write VCF: {exc}") from exc
    with target.open(encoding="utf8") as handle:
        records = sum(1 for line in handle if not line.startswith("#"))
    return {
        "status": "completed",
        "operation": "export_tree_file_vcf",
        "source": str(source),
        "destination": str(target),
        "file_format": "vcf",
        "ploidy": count,
        "contig_id": contig,
        "allow_position_zero": allow_position_zero,
        "bytes": target.stat().st_size,
        "num_samples": int(ts.num_samples),
        "num_records": records,
        "provenance": _provenance("tskit", _VCF_CONVENTION),
    }


__all__ = [
    "export_tree_file_edges_table",
    "export_tree_file_fasta",
    "export_tree_file_individuals_table",
    "export_tree_file_migrations_table",
    "export_tree_file_mutations_table",
    "export_tree_file_nodes_table",
    "export_tree_file_populations_table",
    "export_tree_file_provenances_table",
    "export_tree_file_sites_table",
    "export_tree_file_vcf",
    "tree_file_info",
    "tree_file_trees",
]

_register_public_verbose(globals())
