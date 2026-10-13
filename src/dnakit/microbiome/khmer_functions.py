"""khmer k-mer helpers for microbial DNA: hashing, trimming, graph inspection and table sizing.

The functions here complement ``dnakit.microbiome.query_kmer_counts()``,
``estimate_unique_kmers()``, ``kmer_abundance_histogram()`` and ``normalize_reads_by_median()``,
which count reads. They cover k-mer hashes and 2-bit codes, low-abundance trimming of local reads,
inspection and querying of saved khmer tables, table sizing and a collision estimate.

khmer silently gives wrong answers for lowercase DNA, so every sequence is uppercased before it
reaches khmer, and k is limited to 1 to 32. Reads with symbols other than A, C, G and T are
skipped, as in the existing backends. khmer is imported only when a function runs; a missing package
raises ``BackendUnavailableError`` with the ``pip install "dnakit[bio-khmer]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import math
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_KHMER_EXTRA = "bio-khmer"
_KSIZE_LIMIT = 32
_TABLE_SIZE_LIMITS = (1_000, 100_000_000)
_N_TABLES_LIMITS = (1, 8)
_DEFAULT_TABLE_SIZE = 1_000_000
_DEFAULT_N_TABLES = 4
_DEFAULT_MAX_READS = 1_000_000
_MAX_READS_LIMIT = 100_000_000
_DEFAULT_MAX_BASES = 100_000_000
_MAX_BASES_LIMIT = 10_000_000_000
_MAX_CUTOFF = 1_000_000
_MAX_QUERY_KMERS = 100_000
_COUNT_LIMIT = 10**12
_ACGT = frozenset("ACGT")
_FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
_FASTQ_SUFFIXES = (".fq", ".fastq")
_GRAPH_KINDS = ("countgraph", "nodegraph")


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _kmer_size(value: object) -> int:
    return _integer(value, "ksize", 1, _KSIZE_LIMIT)


def _table_size(value: object) -> int:
    return _integer(value, "table_size", *_TABLE_SIZE_LIMITS)


def _n_tables(value: object) -> int:
    return _integer(value, "n_tables", *_N_TABLES_LIMITS)


def _read_limits(max_reads: object, max_bases: object) -> tuple[int, int]:
    return (
        _integer(max_reads, "max_reads", 1, _MAX_READS_LIMIT),
        _integer(max_bases, "max_bases", 1, _MAX_BASES_LIMIT),
    )


def _kmer(value: object, name: str) -> str:
    """Uppercase a DNA k-mer and check its alphabet and length (1 to 32)."""

    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = value.upper()
    if set(text) - _ACGT:
        raise ConfigurationError(f"{name} must contain only A, C, G and T.")
    if len(text) > _KSIZE_LIMIT:
        raise ConfigurationError(f"{name} length must be at most {_KSIZE_LIMIT}.")
    return text


def _kmer_queries(value: object) -> list[str]:
    if isinstance(value, (str, bytes, Mapping)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("kmers must be a non-empty list of DNA strings.")
    if len(value) > _MAX_QUERY_KMERS:
        raise ConfigurationError(f"kmers can contain at most {_MAX_QUERY_KMERS} items.")
    queries = [_kmer(item, "each k-mer") for item in value]
    if len({len(item) for item in queries}) != 1:
        raise ConfigurationError("all k-mers must have the same length.")
    return queries


def _new_output(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if not any(path.name.endswith(suffix) for suffix in _FASTA_SUFFIXES + _FASTQ_SUFFIXES):
        raise ConfigurationError("output must end with .fa, .fasta, .fna, .fq or .fastq.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {path.parent}")
    if path.exists():
        raise ConfigurationError(f"output already exists and will not be overwritten: {path}")
    return path


def _records(
    khmer: Any, path: Path, max_reads: int, max_bases: int
) -> Iterator[tuple[str, str, str | None]]:
    """Yield ``(name, uppercase sequence, quality or None)`` through khmer's parser with limits."""

    reads = 0
    bases = 0
    try:
        for read in khmer.ReadParser(str(path)):
            reads += 1
            if reads > max_reads:
                raise ConfigurationError(
                    f"read input exceeds max_reads={max_reads}; increase the explicit limit."
                )
            bases += len(read.sequence)
            if bases > max_bases:
                raise ConfigurationError(
                    f"read input exceeds max_bases={max_bases}; increase the explicit limit."
                )
            yield read.name, read.sequence.upper(), getattr(read, "quality", None)
    except ConfigurationError:
        raise
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"khmer could not parse {path.name}: {exc}") from exc


def _classify(sequence: str, k: int) -> str:
    """Return ``usable``, ``skipped_non_acgt`` or ``skipped_short`` for one uppercase read."""

    if set(sequence) - _ACGT:
        return "skipped_non_acgt"
    if len(sequence) < k:
        return "skipped_short"
    return "usable"


def _count_usable(
    khmer: Any, path: Path, k: int, table: Any, max_reads: int, max_bases: int
) -> dict[str, int]:
    """Consume usable reads into ``table``; return read totals for the pass."""

    totals = {"reads": 0, "used_reads": 0, "skipped_non_acgt_reads": 0, "skipped_short_reads": 0}
    for _name, sequence, _quality in _records(khmer, path, max_reads, max_bases):
        totals["reads"] += 1
        status = _classify(sequence, k)
        if status == "usable":
            table.consume(sequence)
            totals["used_reads"] += 1
        else:
            totals[f"{status}_reads"] += 1
    return totals


def hash_kmer_murmur3(
    kmer: str,
    *,
    canonical: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """MurmurHash3 of a DNA k-mer, computed by khmer.

    ``canonical=True`` calls ``khmer.hash_murmur3``, which gives the same value for a k-mer and its
    reverse complement. ``canonical=False`` calls ``khmer.hash_no_rc_murmur3``, which hashes the
    k-mer as written. Lowercase input is uppercased first, because khmer hashes lowercase bases
    incorrectly without warning.
    """

    text = _kmer(kmer, "kmer")
    use_canonical = _flag(canonical, "canonical")
    khmer = require("khmer", _KHMER_EXTRA)
    hashed = khmer.hash_murmur3(text) if use_canonical else khmer.hash_no_rc_murmur3(text)
    return {
        "status": "completed",
        "operation": "hash_kmer_murmur3",
        "kmer": text,
        "ksize": len(text),
        "canonical": use_canonical,
        "hash": int(hashed),
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "hash": "MurmurHash3 of the k-mer (khmer.hash_murmur3 or khmer.hash_no_rc_murmur3)",
            "coordinate_system": "k-mer hash; no genomic coordinates",
            "convention": "input is uppercased; canonical hashes are equal for a k-mer and its "
            "reverse complement; values are khmer's own and differ from sourmash hashes",
        },
    }


def encode_kmer_2bit(
    kmer: str,
    *,
    canonical: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Exact 2-bit integer code of a DNA k-mer (k up to 32), computed by khmer.

    ``canonical=True`` calls ``khmer.forward_hash`` and returns the smaller code of the k-mer and
    its reverse complement. ``canonical=False`` calls ``khmer.forward_hash_no_rc``. khmer's base
    order is A=0, T=1, C=2, G=3. ``decode_kmer_2bit()`` reverses the code.
    """

    text = _kmer(kmer, "kmer")
    use_canonical = _flag(canonical, "canonical")
    khmer = require("khmer", _KHMER_EXTRA)
    k = len(text)
    code = khmer.forward_hash(text, k) if use_canonical else khmer.forward_hash_no_rc(text, k)
    return {
        "status": "completed",
        "operation": "encode_kmer_2bit",
        "kmer": text,
        "ksize": k,
        "canonical": use_canonical,
        "code": int(code),
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "encoding": "2 bits per base in khmer order A=0, T=1, C=2, G=3; first base in the "
            "highest bits",
            "coordinate_system": "k-mer code; no genomic coordinates",
            "convention": "canonical codes choose the smaller of the k-mer and its reverse "
            "complement; input is uppercased",
        },
    }


def decode_kmer_2bit(
    value: int,
    *,
    ksize: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Recover the DNA k-mer from a 2-bit code produced by ``encode_kmer_2bit()``.

    Calls ``khmer.reverse_hash``. A canonical code decodes to the canonical k-mer.
    """

    k = _kmer_size(ksize)
    code = _integer(value, "value", 0, 4**k - 1)
    khmer = require("khmer", _KHMER_EXTRA)
    text = str(khmer.reverse_hash(code, k))
    return {
        "status": "completed",
        "operation": "decode_kmer_2bit",
        "value": code,
        "ksize": k,
        "kmer": text,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "encoding": "inverse of encode_kmer_2bit (khmer order A=0, T=1, C=2, G=3)",
            "coordinate_system": "k-mer code; no genomic coordinates",
            "convention": "the decoded k-mer is uppercase",
        },
    }


def trim_low_abundance_reads(
    source: str | Path,
    output: str | Path,
    *,
    ksize: int = 21,
    cutoff: int = 2,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    max_reads: int = _DEFAULT_MAX_READS,
    max_bases: int = _DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Trim each read at its first k-mer whose count is below ``cutoff``, using khmer.

    Pass one counts canonical k-mers from all usable reads. Pass two calls
    ``khmer.trimming.trim_record``: the read keeps the bases before the first k-mer (from the 5'
    end) whose count is below ``cutoff``. A read that becomes shorter than ``ksize`` is dropped.
    The output is FASTA (``.fa``, ``.fasta``, ``.fna``) or FASTQ (``.fq``, ``.fastq``). FASTQ output
    needs FASTQ input and keeps quality values for the kept bases. Bases are written uppercase.
    """

    path = input_path(source)
    k = _kmer_size(ksize)
    limit = _integer(cutoff, "cutoff", 1, _MAX_CUTOFF)
    size = _table_size(table_size)
    tables = _n_tables(n_tables)
    read_limit, base_limit = _read_limits(max_reads, max_bases)
    target = _new_output(output)
    fastq_output = any(target.name.endswith(suffix) for suffix in _FASTQ_SUFFIXES)
    khmer = require("khmer", _KHMER_EXTRA)
    trimming = require("khmer.trimming", _KHMER_EXTRA)
    table = khmer.Countgraph(k, size, tables)
    totals = _count_usable(khmer, path, k, table, read_limit, base_limit)
    totals.update({"kept_untrimmed_reads": 0, "trimmed_reads": 0, "dropped_short_reads": 0})
    try:
        with target.open("w", encoding="utf-8") as handle:
            for name, sequence, quality in _records(khmer, path, read_limit, base_limit):
                if _classify(sequence, k) != "usable":
                    continue
                if fastq_output and quality is None:
                    raise ConfigurationError("FASTQ output needs FASTQ input with quality values.")
                record = SimpleNamespace(name=name, sequence=sequence, cleaned_seq=sequence)
                if quality is not None:
                    record.quality = quality
                try:
                    trimmed, changed = trimming.trim_record(table, record, limit)
                except (ValueError, OSError) as exc:
                    raise ConfigurationError(f"khmer could not trim read {name!r}: {exc}") from exc
                if trimmed is None:
                    totals["dropped_short_reads"] += 1
                    continue
                if changed:
                    totals["trimmed_reads"] += 1
                else:
                    totals["kept_untrimmed_reads"] += 1
                if fastq_output:
                    handle.write(f"@{trimmed.name}\n{trimmed.sequence}\n+\n{trimmed.quality}\n")
                else:
                    handle.write(f">{trimmed.name}\n{trimmed.sequence}\n")
    except Exception:
        target.unlink(missing_ok=True)
        raise
    kept = totals["kept_untrimmed_reads"] + totals["trimmed_reads"]
    return {
        "status": "completed",
        "operation": "trim_low_abundance_reads",
        "source": str(path),
        "output": str(target),
        "ksize": k,
        "cutoff": limit,
        "table_size": size,
        "n_tables": tables,
        **totals,
        "kept_reads": kept,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "algorithm": "count canonical k-mers, then trim at the first k-mer with count below "
            "cutoff (khmer.trimming.trim_record)",
            "coordinate_system": "read order in the input file; no genomic coordinates",
            "convention": "reads with symbols other than A, C, G and T and reads shorter than "
            "ksize are skipped; reads shorter than ksize after trimming are dropped",
        },
    }


def inspect_khmer_graph(
    source: str | Path,
    *,
    kind: str = "countgraph",
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the header of a saved khmer countgraph or nodegraph without loading its tables.

    ``kind="countgraph"`` calls ``khmer.extract_countgraph_info``; ``kind="nodegraph"`` calls
    ``khmer.extract_nodegraph_info`` and maps its tuple to named fields. The table is not loaded, so
    the file can be large. A file of the other kind is rejected.
    """

    if kind not in _GRAPH_KINDS:
        raise ConfigurationError(f"kind must be one of: {', '.join(_GRAPH_KINDS)}.")
    path = input_path(source)
    khmer = require("khmer", _KHMER_EXTRA)
    try:
        if kind == "countgraph":
            info = khmer.extract_countgraph_info(str(path))
            fields: dict[str, Any] = {
                "ksize": int(info.ksize),
                "table_size": int(info.table_size),
                "n_tables": int(info.n_tables),
                "use_bigcount": bool(info.use_bigcount),
                "format_version": int(info.version),
                "table_type_flag": int(info.ht_type),
                "n_occupied": int(info.n_occupied),
            }
        else:
            values = khmer.extract_nodegraph_info(str(path))
            if len(values) != 6:
                raise ConfigurationError(f"unexpected nodegraph header in {path.name}.")
            ksize, size, tables, version, flag, occupied = (int(value) for value in values)
            fields = {
                "ksize": ksize,
                "table_size": size,
                "n_tables": tables,
                "format_version": version,
                "table_type_flag": flag,
                "n_occupied": occupied,
            }
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"{path.name} is not a khmer {kind}: {exc}") from exc
    expected_flag = 1 if kind == "countgraph" else 2
    if fields["table_type_flag"] != expected_flag or fields["n_tables"] < 1:
        raise ConfigurationError(f"{path.name} is not a khmer {kind}.")
    return {
        "status": "completed",
        "operation": "inspect_khmer_graph",
        "source": str(path),
        "kind": kind,
        **fields,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "source_function": (
                "khmer.extract_countgraph_info"
                if kind == "countgraph"
                else "khmer.extract_nodegraph_info"
            ),
            "coordinate_system": "k-mer table; no genomic coordinates",
            "convention": "table_type_flag is khmer's own code (1 countgraph, 2 nodegraph); "
            "n_occupied counts occupied bins over all tables",
        },
    }


def query_saved_countgraph(
    source: str | Path,
    kmers: Sequence[str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up canonical k-mer counts in a saved khmer countgraph, using ``khmer.load_countgraph``.

    All query k-mers must have the table's ``ksize``. A k-mer and its reverse complement share one
    count, and a hash collision can only raise a count. The result is sorted by k-mer.
    """

    path = input_path(source)
    queries = _kmer_queries(kmers)
    khmer = require("khmer", _KHMER_EXTRA)
    try:
        table = khmer.load_countgraph(str(path))
    except (OSError, ValueError) as exc:
        raise ConfigurationError(f"{path.name} is not a khmer countgraph: {exc}") from exc
    k = int(table.ksize())
    if len(queries[0]) != k:
        raise ConfigurationError(f"query k-mers have length {len(queries[0])}; table ksize is {k}.")
    counts = {kmer: int(table.get(kmer)) for kmer in sorted(set(queries))}
    return {
        "status": "completed",
        "operation": "query_saved_countgraph",
        "source": str(path),
        "ksize": k,
        "hashsizes": [int(value) for value in table.hashsizes()],
        "counts": counts,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "counting": "canonical k-mers; a k-mer and its reverse complement share one count",
            "coordinate_system": "k-mer counts; no genomic coordinates",
            "convention": "queries are uppercased; hash collisions can only raise counts",
        },
    }


def estimate_kmer_table_size(
    num_kmers: int,
    *,
    mem_cap: int | None = None,
    fp_rate: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Table layout for ``num_kmers`` distinct k-mers, from a memory cap or a false-positive rate.

    Give exactly one of ``mem_cap`` and ``fp_rate``. Calls
    ``khmer.khmer_args.optimal_size``. ``mem_use`` and ``mem_cap`` are in khmer's table units (the
    total number of hash bins over all tables); khmer's own helper does not convert them to bytes.
    """

    count = _integer(num_kmers, "num_kmers", 1, _COUNT_LIMIT)
    if (mem_cap is None) == (fp_rate is None):
        raise ConfigurationError("give exactly one of mem_cap or fp_rate.")
    cap = None if mem_cap is None else _integer(mem_cap, "mem_cap", 1, _COUNT_LIMIT)
    rate = None
    if fp_rate is not None:
        if isinstance(fp_rate, bool) or not isinstance(fp_rate, (int, float)):
            raise ConfigurationError("fp_rate must be a number in (0, 1).")
        rate = float(fp_rate)
        if not 0.0 < rate < 1.0 or not math.isfinite(rate):
            raise ConfigurationError("fp_rate must be a number in (0, 1).")
    args = require("khmer.khmer_args", _KHMER_EXTRA)
    try:
        layout = args.optimal_size(count, mem_cap=cap, fp_rate=rate)
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"khmer could not size the tables: {exc}") from exc
    return {
        "status": "completed",
        "operation": "estimate_kmer_table_size",
        "num_kmers": count,
        "mode": "mem_cap" if cap is not None else "fp_rate",
        "mem_cap": cap,
        "requested_fp_rate": rate,
        "num_htables": int(layout.num_htables),
        "htable_size": int(layout.htable_size),
        "mem_use": int(layout.mem_use),
        "expected_fp_rate": float(layout.fp_rate),
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "estimator": "khmer_args.optimal_size: number of tables and bins per table",
            "coordinate_system": "k-mer table sizing; no genomic coordinates",
            "convention": "expected_fp_rate is the Bloom-style estimate khmer uses for these "
            "tables; mem units are khmer table bins",
        },
    }


def expected_kmer_table_collisions(
    source: str | Path,
    *,
    ksize: int = 21,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    max_reads: int = _DEFAULT_MAX_READS,
    max_bases: int = _DEFAULT_MAX_BASES,
    max_false_positive: float = 0.2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected false-positive rate of a khmer presence table filled with the reads of a file.

    Builds a ``khmer.Nodegraph`` from the usable reads, then calls
    ``khmer.calc_expected_collisions`` on it. That function uses the occupied bins divided by the
    smallest table, raised to the number of tables. The estimate is quick and approximate. khmer's
    version exits or prints when the rate is too high; this wrapper does not, and reports
    ``within_threshold`` against ``max_false_positive`` instead.
    """

    path = input_path(source)
    k = _kmer_size(ksize)
    size = _table_size(table_size)
    tables = _n_tables(n_tables)
    read_limit, base_limit = _read_limits(max_reads, max_bases)
    if isinstance(max_false_positive, bool) or not isinstance(max_false_positive, (int, float)):
        raise ConfigurationError("max_false_positive must be a number in (0, 1].")
    threshold = float(max_false_positive)
    if not 0.0 < threshold <= 1.0:
        raise ConfigurationError("max_false_positive must be a number in (0, 1].")
    khmer = require("khmer", _KHMER_EXTRA)
    graph = khmer.Nodegraph(k, size, tables)
    totals = _count_usable(khmer, path, k, graph, read_limit, base_limit)
    rate = float(khmer.calc_expected_collisions(graph, force=True, max_false_pos=math.inf))
    return {
        "status": "completed",
        "operation": "expected_kmer_table_collisions",
        "source": str(path),
        "ksize": k,
        "table_size": size,
        "n_tables": tables,
        "hashsizes": [int(value) for value in graph.hashsizes()],
        "occupied_bins": int(graph.n_occupied()),
        **totals,
        "expected_false_positive_rate": rate,
        "max_false_positive": threshold,
        "within_threshold": rate <= threshold,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "estimator": "khmer.calc_expected_collisions on a Nodegraph of canonical k-mers",
            "coordinate_system": "k-mer presence table; no genomic coordinates",
            "convention": "approximate; reads with symbols other than A, C, G and T are skipped",
        },
    }


__all__ = [
    "decode_kmer_2bit",
    "encode_kmer_2bit",
    "estimate_kmer_table_size",
    "expected_kmer_table_collisions",
    "hash_kmer_murmur3",
    "inspect_khmer_graph",
    "query_saved_countgraph",
    "trim_low_abundance_reads",
]

_register_public_verbose(globals(), tuple(__all__))
