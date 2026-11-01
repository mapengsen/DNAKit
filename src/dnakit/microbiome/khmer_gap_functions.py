"""khmer k-mer table operations as typed DNAKit functions: table sizing from a memory budget, and building
counting (countgraph) and presence (nodegraph) tables from reads, saved to disk for later queries.

The user installs khmer (``pip install "dnakit[bio-khmer]"``). Defaults follow khmer's own scripts: k = 32,
4 tables, table size 1e6 buckets. Reads are upper-cased; reads with symbols other than A, C, G, T or shorter
than k are not consumed and are counted in the result. Saved tables can be queried with
``dnakit.microbiome.query_saved_countgraph``.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import Any, Literal

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_KHMER_EXTRA = "bio-khmer"
_DEFAULT_K = 32
_DEFAULT_TABLE_SIZE = 1_000_000
_DEFAULT_N_TABLES = 4
_DEFAULT_MAX_READS = 1_000_000
_DEFAULT_MAX_BASES = 100_000_000
_COUNT_LIMIT = 10**12
_ACGT = frozenset("ACGT")
_KINDS = ("countgraph", "smallcountgraph", "nodegraph")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be between {minimum} and {maximum}.")
    return value


def _ksize(value: object) -> int:
    return _integer(value, "ksize", 1, 32)


def _n_tables(value: object) -> int:
    return _integer(value, "n_tables", 1, 20)


def _table_size(value: object) -> int:
    return _integer(value, "table_size", 1, _COUNT_LIMIT)


def _limits(max_reads: object, max_bases: object) -> tuple[int, int]:
    return _integer(max_reads, "max_reads", 1, _COUNT_LIMIT), _integer(max_bases, "max_bases", 1, _COUNT_LIMIT)


def _new_output(value: str | Path) -> Path:
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"output_path must not exist yet: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for output_path does not exist: {target.parent}")
    return target


def _reads(khmer: Any, path: Path, max_reads: int, max_bases: int) -> Iterator[tuple[str, str]]:
    """Yield (name, upper-case sequence) through khmer's parser, enforcing the explicit read and base limits."""
    reads = 0
    bases = 0
    try:
        for read in khmer.ReadParser(str(path)):
            reads += 1
            if reads > max_reads:
                raise ConfigurationError(f"read input exceeds max_reads={max_reads}; increase the explicit limit.")
            bases += len(read.sequence)
            if bases > max_bases:
                raise ConfigurationError(f"read input exceeds max_bases={max_bases}; increase the explicit limit.")
            yield read.name, read.sequence.upper()
    except ConfigurationError:
        raise
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"khmer could not parse {path.name}: {exc}") from exc


def khmer_table_size_from_memory(
    graph_type: Literal["countgraph", "smallcountgraph", "nodegraph"],
    *,
    max_memory_bytes: int | None = None,
    max_tablesize: int | None = None,
    n_tables: int = _DEFAULT_N_TABLES,
    multiplier: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Buckets per table and estimated memory of a khmer table, from a memory budget or a table size.

    Give exactly one of ``max_memory_bytes`` (the table size follows from the memory and the number of tables,
    as khmer's calculate_graphsize does) and ``max_tablesize`` (buckets per table). ``multiplier`` scales the
    memory-derived size, as in khmer. The memory estimate is ``n_tables * table_size / buckets_per_byte``.
    """
    if graph_type not in _KINDS:
        raise ConfigurationError(f"graph_type must be one of {', '.join(_KINDS)}.")
    if (max_memory_bytes is None) == (max_tablesize is None):
        raise ConfigurationError("give exactly one of max_memory_bytes and max_tablesize.")
    tables = _n_tables(n_tables)
    khmer = require("khmer", _KHMER_EXTRA)
    buckets_per_byte = int(khmer._buckets_per_byte[graph_type])
    if max_memory_bytes is not None:
        if isinstance(max_memory_bytes, bool) or not isinstance(max_memory_bytes, (int, float)) or max_memory_bytes <= 0:
            raise ConfigurationError("max_memory_bytes must be a positive number.")
        if multiplier <= 0:
            raise ConfigurationError("multiplier must be positive.")
        table_size = int(float(multiplier) * buckets_per_byte * float(max_memory_bytes) / tables)
    else:
        table_size = _table_size(max_tablesize)
    if table_size < 1:
        raise ConfigurationError("The memory budget is too small for one bucket per table.")
    return {
        "status": "completed",
        "operation": "khmer_table_size_from_memory",
        "graph_type": graph_type,
        "buckets_per_byte": buckets_per_byte,
        "n_tables": tables,
        "table_size": table_size,
        "estimated_memory_bytes": tables * table_size / buckets_per_byte,
    }


def _build(
    source: str,
    output_path: str,
    *,
    kind: Literal["countgraph", "smallcountgraph", "nodegraph"],
    ksize: int,
    table_size: int,
    n_tables: int,
    max_reads: int,
    max_bases: int,
) -> dict[str, Any]:
    path = input_path(source)
    k = _ksize(ksize)
    size = _table_size(table_size)
    tables = _n_tables(n_tables)
    read_limit, base_limit = _limits(max_reads, max_bases)
    target = _new_output(output_path)
    khmer = require("khmer", _KHMER_EXTRA)
    if kind == "countgraph":
        table = khmer.Countgraph(k, size, tables)
    elif kind == "smallcountgraph":
        table = khmer.SmallCountgraph(k, size, tables)
    else:
        table = khmer.Nodegraph(k, size, tables)
    totals = {"reads": 0, "used_reads": 0, "skipped_non_acgt_reads": 0, "skipped_short_reads": 0}
    kmers = 0
    for _name, sequence in _reads(khmer, path, read_limit, base_limit):
        totals["reads"] += 1
        if len(sequence) < k:
            totals["skipped_short_reads"] += 1
        elif not set(sequence) <= _ACGT:
            totals["skipped_non_acgt_reads"] += 1
        else:
            kmers += int(table.consume(sequence))
            totals["used_reads"] += 1
    table.save(str(target))
    return {
        "status": "completed",
        "source": str(path),
        "output": str(target),
        "kind": kind,
        "ksize": k,
        "table_size": size,
        "n_tables": tables,
        "hashsizes": [int(value) for value in table.hashsizes()],
        "occupied_bins": int(table.n_occupied()),
        "kmers_consumed": kmers,
        **totals,
    }


def khmer_build_countgraph(
    source: str,
    output_path: str,
    *,
    ksize: int = _DEFAULT_K,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    small_count: bool = False,
    max_reads: int = _DEFAULT_MAX_READS,
    max_bases: int = _DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count canonical k-mers of FASTA or FASTQ reads into a khmer countgraph and save it (khmer.Countgraph).

    ``small_count=True`` uses khmer's SmallCountgraph (smaller counters, saturating counts). The saved file can be
    loaded with khmer.load_countgraph or queried with ``query_saved_countgraph``.
    """
    kind: Literal["countgraph", "smallcountgraph"] = "smallcountgraph" if small_count else "countgraph"
    return _build(
        source, output_path, kind=kind, ksize=ksize, table_size=table_size, n_tables=n_tables,
        max_reads=max_reads, max_bases=max_bases,
    )


def khmer_build_nodegraph(
    source: str,
    output_path: str,
    *,
    ksize: int = _DEFAULT_K,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    max_reads: int = _DEFAULT_MAX_READS,
    max_bases: int = _DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Record which canonical k-mers occur in FASTA or FASTQ reads in a khmer nodegraph (presence table).

    A nodegraph stores one bit per k-mer position, so it answers presence questions with a small, approximate
    false-positive rate (see ``estimate_kmer_table_size`` for sizing). The saved file can be loaded with
    khmer.load_nodegraph.
    """
    return _build(
        source, output_path, kind="nodegraph", ksize=ksize, table_size=table_size, n_tables=n_tables,
        max_reads=max_reads, max_bases=max_bases,
    )


__all__ = [
    "khmer_build_countgraph",
    "khmer_build_nodegraph",
    "khmer_table_size_from_memory",
]
_register_public_verbose(globals())
