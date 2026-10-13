"""Optional MinHash sketch and k-mer backends for microbial DNA sequences.

The functions in this module are explicit adapters over optional packages:

- ``sketch_fasta_minhash()`` sketches FASTA records with ``sourmash.MinHash``
  (scaled sketches over canonical DNA k-mers). ``compare_minhash_sketches()``
  reports Jaccard similarity and containment with sourmash's own methods.
  ``write_sourmash_signatures()`` and ``read_sourmash_signatures()`` exchange
  the same sketches as sourmash signature files.
- ``query_kmer_counts()``, ``estimate_unique_kmers()``,
  ``kmer_abundance_histogram()`` and ``normalize_reads_by_median()`` use khmer
  counting and presence tables (canonical k-mers, ksize 1 to 32). Reads are
  uppercased and streamed with explicit read and base limits.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import gzip
import importlib.metadata
import re
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_SOURMASH_EXTRA = "bio-sourmash"
_HASH_LIMIT = 2**64 - 1
_KSIZE_LIMIT = 100
_SCALED_LIMIT = 1_000_000_000
_DEFAULT_MAX_RECORDS = 1_000
_MAX_RECORDS_LIMIT = 100_000
_DEFAULT_MAX_BASES = 10_000_000
_MAX_BASES_LIMIT = 1_000_000_000
_DEFAULT_MAX_SIGNATURES = 1_000
_MAX_SIGNATURES_LIMIT = 100_000
_SIGNATURE_SUFFIXES = (".sig", ".sig.gz")
_SEQUENCE_PATTERN = re.compile(r"[A-Za-z]+")


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _new_output(value: str | Path, suffixes: tuple[str, ...]) -> Path:
    path = Path(value).expanduser().resolve()
    if not any(path.name.endswith(suffix) for suffix in suffixes):
        raise ConfigurationError(f"output must end with one of: {', '.join(suffixes)}.")
    if not path.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {path.parent}")
    if path.exists():
        raise ConfigurationError(f"output already exists and will not be overwritten: {path}")
    return path


def _fasta_records(path: Path, max_records: int, max_bases: int) -> list[tuple[str, str]]:
    """Read FASTA records with explicit limits. This parser is pure Python."""

    records: list[tuple[str, list[str]]] = []
    bases = 0
    try:
        with path.open("r", encoding="utf-8") as handle:
            for number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line:
                    continue
                if "\x00" in line:
                    raise ConfigurationError(f"FASTA line {number} contains a NUL character.")
                if line.startswith(">"):
                    identifier = line[1:].split()
                    if not identifier:
                        raise ConfigurationError(f"FASTA header on line {number} has no name.")
                    if len(records) >= max_records:
                        raise ConfigurationError(
                            f"FASTA input exceeds max_records={max_records}; "
                            "increase the explicit limit."
                        )
                    records.append((identifier[0], []))
                    continue
                if not records:
                    raise ConfigurationError("FASTA sequence data must follow a '>' header.")
                if not _SEQUENCE_PATTERN.fullmatch(line):
                    raise ConfigurationError(f"FASTA line {number} contains non-letter characters.")
                bases += len(line)
                if bases > max_bases:
                    raise ConfigurationError(
                        f"FASTA input exceeds max_bases={max_bases}; increase the explicit limit."
                    )
                records[-1][1].append(line.upper())
    except UnicodeDecodeError as exc:
        raise ConfigurationError(f"FASTA input is not UTF-8 text: {exc}") from exc
    if not records:
        raise ConfigurationError("FASTA input contains no records.")
    return [(identifier, "".join(chunks)) for identifier, chunks in records]


def _sketch_hashes(sourmash: Any, sequence: str, ksize: int, scaled: int) -> list[int]:
    try:
        minhash = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
        minhash.add_sequence(sequence, force=True)
        return sorted(int(value) for value in minhash.hashes)
    except (AttributeError, OverflowError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash sketching failed: {exc}") from exc


def _minhash(sourmash: Any, ksize: int, scaled: int, hashes: Sequence[int]) -> Any:
    """Rebuild a sourmash MinHash from explicit hashes so sourmash compares them."""

    try:
        minhash = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
        minhash.add_many(hashes)
    except (AttributeError, OverflowError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash could not build a MinHash: {exc}") from exc
    if len(minhash) != len(set(hashes)):
        raise ConfigurationError("some hashes lie outside the scaled range of their sketch.")
    return minhash


def _sketch_fields(sketch: object, name: str) -> tuple[int, int, list[int]]:
    if not isinstance(sketch, Mapping):
        raise ConfigurationError(
            f"{name} must be a sketch mapping, such as one item of "
            "sketch_fasta_minhash()['sketches']."
        )
    ksize = _integer(sketch.get("ksize"), f"{name} ksize", 1, _KSIZE_LIMIT)
    scaled = _integer(sketch.get("scaled"), f"{name} scaled", 1, _SCALED_LIMIT)
    hashes = sketch.get("hashes")
    if isinstance(hashes, (str, bytes)) or not isinstance(hashes, Sequence) or not hashes:
        raise ConfigurationError(f"{name} hashes must be a non-empty list of integers.")
    values: list[int] = []
    for value in hashes:
        if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= _HASH_LIMIT:
            raise ConfigurationError(f"{name} hashes must be unsigned 64-bit integers.")
        values.append(value)
    return ksize, scaled, sorted(set(values))


def _identifier(value: object) -> str | None:
    return None if value is None else str(value)


def sketch_fasta_minhash(
    source: str | Path,
    *,
    ksize: int = 31,
    scaled: int = 1_000,
    max_records: int = _DEFAULT_MAX_RECORDS,
    max_bases: int = _DEFAULT_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sketch each FASTA record with a scaled sourmash MinHash over canonical DNA k-mers."""

    path = input_path(source)
    k = _integer(ksize, "ksize", 1, _KSIZE_LIMIT)
    s = _integer(scaled, "scaled", 1, _SCALED_LIMIT)
    record_limit = _integer(max_records, "max_records", 1, _MAX_RECORDS_LIMIT)
    base_limit = _integer(max_bases, "max_bases", 1, _MAX_BASES_LIMIT)
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    sketches = []
    for identifier, sequence in _fasta_records(path, record_limit, base_limit):
        hashes = _sketch_hashes(sourmash, sequence, k, s)
        sketches.append(
            {
                "id": identifier,
                "length": len(sequence),
                "ksize": k,
                "scaled": s,
                "hash_count": len(hashes),
                "hashes": hashes,
            }
        )
    return {
        "status": "completed",
        "operation": "sketch_fasta_minhash",
        "source": str(path),
        "ksize": k,
        "scaled": s,
        "records": len(sketches),
        "sketches": sketches,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "sketch": "scaled MinHash, murmur64 hashes, default seed 42, canonical DNA k-mers",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": (
                "k-mers containing symbols other than A, C, G and T are skipped by sourmash; "
                "records shorter than ksize have hash_count 0; FASTA parsing is pure Python"
            ),
        },
    }


def compare_minhash_sketches(
    left: Mapping[str, Any], right: Mapping[str, Any], *, verbose: bool = True
) -> dict[str, Any]:
    """Jaccard similarity and containment between two scaled MinHash sketches, using sourmash.

    Both sketches must share ``ksize`` and ``scaled``. Containment is directional:
    ``left_contained_by_right`` is the fraction of the left sketch's hashes that the
    right sketch also contains.
    """

    sourmash = require("sourmash", _SOURMASH_EXTRA)
    left_k, left_s, left_hashes = _sketch_fields(left, "left")
    right_k, right_s, right_hashes = _sketch_fields(right, "right")
    if (left_k, left_s) != (right_k, right_s):
        raise ConfigurationError("sketches must share ksize and scaled to be compared.")
    left_minhash = _minhash(sourmash, left_k, left_s, left_hashes)
    right_minhash = _minhash(sourmash, right_k, right_s, right_hashes)
    try:
        jaccard = float(
            left_minhash.similarity(right_minhash, ignore_abundance=True, downsample=False)
        )
        left_contained = float(left_minhash.contained_by(right_minhash, downsample=False))
        right_contained = float(right_minhash.contained_by(left_minhash, downsample=False))
        shared = int(left_minhash.count_common(right_minhash, downsample=False))
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"sourmash comparison failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "compare_minhash_sketches",
        "left_id": _identifier(left.get("id")),
        "right_id": _identifier(right.get("id")),
        "ksize": left_k,
        "scaled": left_s,
        "left_hash_count": len(left_hashes),
        "right_hash_count": len(right_hashes),
        "shared_hash_count": shared,
        "jaccard": jaccard,
        "left_contained_by_right": left_contained,
        "right_contained_by_left": right_contained,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "estimator": "set Jaccard and containment from sourmash MinHash methods",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "containment is directional: |L and R| / |L| for left_contained_by_right",
        },
    }


def write_sourmash_signatures(
    sketches: Sequence[Mapping[str, Any]],
    output: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write sketches as one sourmash JSON signature file (``.sig`` or ``.sig.gz``)."""

    if (
        isinstance(sketches, (str, bytes, Mapping))
        or not isinstance(sketches, Sequence)
        or not sketches
    ):
        raise ConfigurationError("sketches must be a non-empty list of sketch mappings.")
    if len(sketches) > _MAX_SIGNATURES_LIMIT:
        raise ConfigurationError(f"at most {_MAX_SIGNATURES_LIMIT} sketches can be written.")
    path = _new_output(output, _SIGNATURE_SUFFIXES)
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    signatures = []
    for index, sketch in enumerate(sketches):
        ksize, scaled, hashes = _sketch_fields(sketch, f"sketches[{index}]")
        name = _identifier(sketch.get("id")) or ""
        minhash = _minhash(sourmash, ksize, scaled, hashes)
        signatures.append(sourmash.signature.SourmashSignature(minhash, name=name))
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8") as handle:
        sourmash.save_signatures_to_json(signatures, handle)
    return {
        "status": "completed",
        "operation": "write_sourmash_signatures",
        "output": str(path),
        "signatures": len(signatures),
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "format": "sourmash JSON signatures; compressed when the name ends with .gz",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
        },
    }


def read_sourmash_signatures(
    source: str | Path,
    *,
    max_signatures: int = _DEFAULT_MAX_SIGNATURES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read DNA scaled MinHash signatures from a sourmash ``.sig`` or ``.sig.gz`` file."""

    path = input_path(source)
    limit = _integer(max_signatures, "max_signatures", 1, _MAX_SIGNATURES_LIMIT)
    sourmash = require("sourmash", _SOURMASH_EXTRA)
    loaded: list[Any] = []
    try:
        for item in sourmash.load_file_as_signatures(str(path)):
            if len(loaded) >= limit:
                raise ConfigurationError(
                    f"signature file exceeds max_signatures={limit}; increase the explicit limit."
                )
            loaded.append(item)
    except ConfigurationError:
        raise
    except Exception as exc:  # sourmash raises generic exceptions for malformed files
        raise ConfigurationError(f"sourmash could not read {path.name}: {exc}") from exc
    sketches = []
    for item in loaded:
        minhash = item.minhash
        if minhash.moltype != "DNA":
            raise ConfigurationError(
                f"signature {item.name!r} uses {minhash.moltype}; "
                "only DNA signatures are supported."
            )
        if int(minhash.scaled) == 0:
            raise ConfigurationError(
                f"signature {item.name!r} is num-based; only scaled signatures are supported."
            )
        hashes = sorted(int(value) for value in minhash.hashes)
        sketches.append(
            {
                "id": item.name,
                "length": None,
                "ksize": int(minhash.ksize),
                "scaled": int(minhash.scaled),
                "hash_count": len(hashes),
                "hashes": hashes,
                "md5sum": item.md5sum(),
            }
        )
    return {
        "status": "completed",
        "operation": "read_sourmash_signatures",
        "source": str(path),
        "signatures": len(sketches),
        "sketches": sketches,
        "provenance": {
            "package": f"sourmash {_version('sourmash')}",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "abundance values are dropped; signatures store no sequence length",
        },
    }


_KHMER_EXTRA = "bio-khmer"
_KHMER_KSIZE_LIMIT = 32
_TABLE_SIZE_LIMITS = (1_000, 100_000_000)
_N_TABLES_LIMITS = (1, 8)
_DEFAULT_TABLE_SIZE = 1_000_000
_DEFAULT_N_TABLES = 4
_DEFAULT_KHMER_MAX_READS = 1_000_000
_MAX_READS_LIMIT = 100_000_000
_DEFAULT_KHMER_MAX_BASES = 100_000_000
_MAX_BASES_LIMIT = 10_000_000_000
_MAX_QUERY_KMERS = 100_000
_MAX_CUTOFF = 1_000_000
_ACGT = frozenset("ACGT")
_FASTA_SUFFIXES = (".fa", ".fasta", ".fna")


def _kmer_size(value: object) -> int:
    return _integer(value, "ksize", 1, _KHMER_KSIZE_LIMIT)


def _table_size(value: object) -> int:
    return _integer(value, "table_size", *_TABLE_SIZE_LIMITS)


def _n_tables(value: object) -> int:
    return _integer(value, "n_tables", *_N_TABLES_LIMITS)


def _read_limits(max_reads: object, max_bases: object) -> tuple[int, int]:
    return (
        _integer(max_reads, "max_reads", 1, _MAX_READS_LIMIT),
        _integer(max_bases, "max_bases", 1, _MAX_BASES_LIMIT),
    )


def _kmer_queries(value: object) -> list[str]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("kmers must be a non-empty list of DNA strings.")
    if len(value) > _MAX_QUERY_KMERS:
        raise ConfigurationError(f"kmers can contain at most {_MAX_QUERY_KMERS} items.")
    queries: list[str] = []
    length: int | None = None
    for item in value:
        if not isinstance(item, str) or not item:
            raise ConfigurationError("each k-mer must be a non-empty string.")
        if not 1 <= len(item) <= _KHMER_KSIZE_LIMIT:
            raise ConfigurationError(f"k-mer length must be in [1, {_KHMER_KSIZE_LIMIT}].")
        if length is None:
            length = len(item)
        if len(item) != length:
            raise ConfigurationError("all k-mers must have the same length.")
        if set(item) - _ACGT:
            raise ConfigurationError("k-mers must contain only uppercase A, C, G and T.")
        queries.append(item)
    return queries


def _reads(
    khmer: Any, path: Path, max_reads: int, max_bases: int
) -> Iterator[tuple[Any, str | None]]:
    """Stream reads through khmer's parser under explicit limits.

    Yields ``(read, sequence)``. ``sequence`` is the uppercased read, or ``None``
    when it is empty or contains symbols other than A, C, G and T.
    """

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
            upper = read.sequence.upper()
            yield read, upper if upper and set(upper) <= _ACGT else None
    except OSError as exc:
        raise ConfigurationError(f"khmer could not parse {path.name}: {exc}") from exc


def _count_reads(
    khmer: Any, path: Path, ksize: int, table: Any, max_reads: int, max_bases: int
) -> dict[str, int]:
    """Consume usable reads into ``table`` and count the reads that were skipped."""

    totals = {"reads": 0, "used_reads": 0, "skipped_non_acgt_reads": 0, "skipped_short_reads": 0}
    positions = 0
    for _read, sequence in _reads(khmer, path, max_reads, max_bases):
        totals["reads"] += 1
        if sequence is None:
            totals["skipped_non_acgt_reads"] += 1
        elif len(sequence) < ksize:
            totals["skipped_short_reads"] += 1
        else:
            positions += int(table.consume(sequence))
            totals["used_reads"] += 1
    totals["kmer_positions"] = positions
    return totals


def query_kmer_counts(
    source: str | Path,
    kmers: Sequence[str],
    *,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    max_reads: int = _DEFAULT_KHMER_MAX_READS,
    max_bases: int = _DEFAULT_KHMER_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count how often each query k-mer occurs in a FASTA or FASTQ file, using khmer.

    All query k-mers must share one length from 1 to 32 and use uppercase A, C, G and T.
    A k-mer and its reverse complement share one count. Hash collisions can only raise
    counts, so use a larger ``table_size`` when counts look inflated.
    """

    path = input_path(source)
    queries = _kmer_queries(kmers)
    size = _table_size(table_size)
    tables = _n_tables(n_tables)
    read_limit, base_limit = _read_limits(max_reads, max_bases)
    k = len(queries[0])
    khmer = require("khmer", _KHMER_EXTRA)
    table = khmer.Countgraph(k, size, tables)
    totals = _count_reads(khmer, path, k, table, read_limit, base_limit)
    counts = {kmer: int(table.get(kmer)) for kmer in queries}
    return {
        "status": "completed",
        "operation": "query_kmer_counts",
        "source": str(path),
        "ksize": k,
        "table_size": size,
        "n_tables": tables,
        **totals,
        "counts": counts,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "counting": "canonical k-mers; a k-mer and its reverse complement share one count",
            "coordinate_system": "k-mer counts; no genomic coordinates",
            "convention": (
                "reads with symbols other than A, C, G and T are skipped; "
                "reads are uppercased before counting"
            ),
        },
    }


def estimate_unique_kmers(
    source: str | Path,
    *,
    ksize: int = 21,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    max_reads: int = _DEFAULT_KHMER_MAX_READS,
    max_bases: int = _DEFAULT_KHMER_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Estimate the number of distinct canonical k-mers in a FASTA or FASTQ file, using khmer.

    The estimate comes from a presence table (``Nodegraph``) and khmer's occupancy
    formula, so it is approximate. Reads with symbols other than A, C, G and T are skipped.
    """

    path = input_path(source)
    k = _kmer_size(ksize)
    size = _table_size(table_size)
    tables = _n_tables(n_tables)
    read_limit, base_limit = _read_limits(max_reads, max_bases)
    khmer = require("khmer", _KHMER_EXTRA)
    graph = khmer.Nodegraph(k, size, tables)
    totals = _count_reads(khmer, path, k, graph, read_limit, base_limit)
    return {
        "status": "completed",
        "operation": "estimate_unique_kmers",
        "source": str(path),
        "ksize": k,
        "table_size": size,
        "n_tables": tables,
        "hashsizes": [int(value) for value in graph.hashsizes()],
        **totals,
        "occupied_bins": int(graph.n_occupied()),
        "unique_kmers_estimate": int(graph.n_unique_kmers()),
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "estimator": "khmer Nodegraph occupancy formula; approximate",
            "counting": "canonical k-mers; a k-mer and its reverse complement count once",
            "coordinate_system": "k-mer hashes; no genomic coordinates",
            "convention": "reads with symbols other than A, C, G and T are skipped",
        },
    }


def kmer_abundance_histogram(
    source: str | Path,
    *,
    ksize: int = 21,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    max_reads: int = _DEFAULT_KHMER_MAX_READS,
    max_bases: int = _DEFAULT_KHMER_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Histogram of how many distinct canonical k-mers occur each number of times, using khmer.

    khmer rescans the file for the histogram, so every read must be uppercase and made
    only of A, C, G and T. Other inputs are rejected instead of being skipped.
    """

    path = input_path(source)
    k = _kmer_size(ksize)
    size = _table_size(table_size)
    tables = _n_tables(n_tables)
    read_limit, base_limit = _read_limits(max_reads, max_bases)
    khmer = require("khmer", _KHMER_EXTRA)
    counts = khmer.Countgraph(k, size, tables)
    for read, sequence in _reads(khmer, path, read_limit, base_limit):
        if sequence is None or read.sequence != sequence:
            raise ConfigurationError(
                "abundance histograms need reads made only of uppercase A, C, G and T; "
                f"found an unsupported read named {read.name!r}."
            )
        if len(sequence) < k:
            raise ConfigurationError(f"read {read.name!r} is shorter than ksize={k}.")
        counts.consume(sequence)
    histogram = counts.abundance_distribution(str(path), khmer.Nodegraph(k, size, tables))
    bins = [{"count": index, "kmers": int(value)} for index, value in enumerate(histogram) if value]
    return {
        "status": "completed",
        "operation": "kmer_abundance_histogram",
        "source": str(path),
        "ksize": k,
        "table_size": size,
        "n_tables": tables,
        "distinct_kmers": sum(item["kmers"] for item in bins),
        "histogram": bins,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "counting": "canonical k-mers; a k-mer and its reverse complement share one count",
            "coordinate_system": "k-mer counts; no genomic coordinates",
            "convention": "bins list only counts that occur; approximate under collisions",
        },
    }


def normalize_reads_by_median(
    source: str | Path,
    output: str | Path,
    *,
    ksize: int = 20,
    cutoff: int = 20,
    table_size: int = _DEFAULT_TABLE_SIZE,
    n_tables: int = _DEFAULT_N_TABLES,
    max_reads: int = _DEFAULT_KHMER_MAX_READS,
    max_bases: int = _DEFAULT_KHMER_MAX_BASES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Digital normalization: keep a read when its median k-mer count is below ``cutoff``.

    Reads are visited in file order. The median is taken before a read is counted, and
    kept reads are then counted, so coverage above ``cutoff`` is discarded. Kept reads are
    written as FASTA with their original sequence. Reads with symbols other than A, C, G
    and T, and reads shorter than ``ksize``, are skipped.
    """

    path = input_path(source)
    k = _kmer_size(ksize)
    limit = _integer(cutoff, "cutoff", 1, _MAX_CUTOFF)
    size = _table_size(table_size)
    tables = _n_tables(n_tables)
    read_limit, base_limit = _read_limits(max_reads, max_bases)
    target = _new_output(output, _FASTA_SUFFIXES)
    khmer = require("khmer", _KHMER_EXTRA)
    table = khmer.Countgraph(k, size, tables)
    totals = {"reads": 0, "kept_reads": 0, "discarded_high_median_reads": 0}
    totals.update({"skipped_non_acgt_reads": 0, "skipped_short_reads": 0})
    try:
        with target.open("w", encoding="utf-8") as handle:
            for read, sequence in _reads(khmer, path, read_limit, base_limit):
                totals["reads"] += 1
                if sequence is None:
                    totals["skipped_non_acgt_reads"] += 1
                    continue
                if len(sequence) < k:
                    totals["skipped_short_reads"] += 1
                    continue
                median = int(table.get_median_count(sequence)[0])
                if median < limit:
                    table.consume(sequence)
                    handle.write(f">{read.name}\n{read.sequence}\n")
                    totals["kept_reads"] += 1
                else:
                    totals["discarded_high_median_reads"] += 1
    except Exception:
        target.unlink(missing_ok=True)
        raise
    return {
        "status": "completed",
        "operation": "normalize_reads_by_median",
        "source": str(path),
        "output": str(target),
        "ksize": k,
        "cutoff": limit,
        "table_size": size,
        "n_tables": tables,
        **totals,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "algorithm": "digital normalization by median canonical k-mer count",
            "coordinate_system": "read order in the input file; no genomic coordinates",
            "convention": "output is FASTA; qualities and FASTQ headers are not written",
        },
    }


__all__ = [
    "compare_minhash_sketches",
    "estimate_unique_kmers",
    "kmer_abundance_histogram",
    "normalize_reads_by_median",
    "query_kmer_counts",
    "read_sourmash_signatures",
    "sketch_fasta_minhash",
    "write_sourmash_signatures",
]

_register_public_verbose(globals(), tuple(__all__))
