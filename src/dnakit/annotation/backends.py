"""Optional-backend interval, annotation, and region-set adapters for DNA and gene coordinates.

Each public function wraps one upstream operation on local inputs:

- ``bioframe``, ``pybedtools``, ``pyranges1``, and ``polars-bio``: merge, closest or nearest,
  coverage, intersection, complement, and FASTA extraction for BED intervals.
- ``gtfparse`` and ``gffutils``: GTF tables with attribute columns, gene and transcript
  selection, GFF/GTF database creation, and child queries.
- ``pyBigWig``: per-base values, intervals, and interval statistics from bigWig tracks.
- ``gtars``: base-pair Jaccard similarity, gaps between regions, and region tokenization.
- ``ncls`` and ``genomicranges``: interval index queries and range arithmetic.

Coordinates are BED-style 0-based half-open intervals unless a docstring says otherwise.
gffutils reports GFF feature coordinates as stored (1-based, inclusive). Packages are imported
only when a function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import gzip
import itertools
import math
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_DEFAULT_ROWS = 100_000
_MAX_ROWS = 1_000_000
_COMMENT_PREFIXES = ("#", "track", "browser")
_BIGWIG_STATISTICS = frozenset({"coverage", "max", "mean", "min", "std", "sum"})


# Validation and conversion helpers.


def _integer(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _text(value: object, name: str, *, maximum: int = 512) -> str:
    if not isinstance(value, str) or not value or len(value) > maximum:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {maximum} characters."
        )
    if "\x00" in value:
        raise ConfigurationError(f"{name} cannot contain NUL characters.")
    return value


def _text_list(value: object, name: str, *, maximum: int = 10_000) -> tuple[str, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise ConfigurationError(f"{name} must be a list of strings.")
    items = tuple(_text(item, name) for item in value)
    if len(items) > maximum:
        raise ConfigurationError(f"{name} accepts at most {maximum} values.")
    return items


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _rows_limit(value: object) -> int:
    return _integer(value, "max_rows", 1, _MAX_ROWS)


def _lengths(value: object, name: str) -> dict[str, int]:
    if not isinstance(value, Mapping) or not value or len(value) > 100_000:
        raise ConfigurationError(
            f"{name} must be a non-empty mapping from chromosome name to length."
        )
    return {
        _text(chrom, f"{name} chromosome"): _integer(length, f"{name} length", 1, 10**12)
        for chrom, length in value.items()
    }


def _plain(value: Any) -> Any:
    """Convert numpy, pandas, and missing values to JSON-compatible Python objects."""

    if value is None or type(value).__name__ in {"NAType", "NaTType"}:
        return None
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def _package(distribution: str) -> str:
    try:
        return f"{distribution} {version(distribution)}"
    except PackageNotFoundError:
        return distribution


def _checked(rows: list[Any], limit: int, what: str) -> list[Any]:
    if len(rows) > limit:
        raise ConfigurationError(f"{what} has {len(rows)} rows; max_rows={limit} is the limit.")
    return rows


def _record_order(record: Mapping[str, Any]) -> tuple[tuple[int, Any], ...]:
    """Type-safe sort key so rows come back in the same order regardless of hash seed."""

    keys: list[tuple[int, Any]] = []
    for value in record.values():
        if value is None:
            keys.append((3, ""))
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            keys.append((0, float(value)))
        else:
            keys.append((2, str(value)))
    return tuple(keys)


def _ordered(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    # Grouped upstream results can follow hash-dependent order (for example, per chromosome).
    return sorted(records, key=_record_order)


def _pandas_records(
    frame: Any, limit: int, what: str, *, ordered: bool = False
) -> list[dict[str, Any]]:
    _checked(frame, limit, what)
    records = [
        {str(key): _plain(value) for key, value in row.items()} for row in frame.to_dict("records")
    ]
    return _ordered(records) if ordered else records


def _gtf_rows(frame: Any, limit: int, what: str) -> list[dict[str, Any]]:
    """gtfparse fills absent attribute keys with empty strings; report them as None."""

    return [
        {key: (None if value == "" else value) for key, value in row.items()}
        for row in _pandas_records(frame, limit, what)
    ]


def _polars_records(
    frame: Any, limit: int, what: str, *, ordered: bool = False
) -> list[dict[str, Any]]:
    _checked(frame, limit, what)
    records = [{str(key): _plain(value) for key, value in row.items()} for row in frame.to_dicts()]
    return _ordered(records) if ordered else records


def _output_lines(text: str) -> list[list[str]]:
    return [line.split("\t") for line in text.splitlines() if line.strip()]


def _bedtools_call(operation: str, call: Any) -> Any:
    try:
        return call()
    except Exception as exc:  # pybedtools raises BEDToolsError with the BEDTools message
        raise ConfigurationError(f"pybedtools {operation} failed: {exc}") from exc


# Bounded file readers.


def _data_lines(path: Path) -> Iterator[tuple[int, list[str]]]:
    """Yield (line number, tab-separated fields) for every data line, skipping comments."""

    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                text = line.rstrip("\r\n")
                if not text.strip() or text.startswith(_COMMENT_PREFIXES):
                    continue
                yield number, text.split("\t")
    except UnicodeDecodeError as exc:
        raise ConfigurationError(f"{path.name} is not UTF-8 text.") from exc


def _scan_bed(path: Path, limit: int) -> int:
    """Check the row limit and BED column count, returning the number of columns (3 to 6)."""

    width = 0
    count = 0
    for _, fields in _data_lines(path):
        if not width:
            width = len(fields)
            if not 3 <= width <= 6:
                raise ConfigurationError(
                    f"{path.name} must be BED3 to BED6; found {width} columns."
                )
        count += 1
        if count > limit:
            raise ConfigurationError(f"{path.name} has more than max_rows={limit} data lines.")
    if not count:
        raise ConfigurationError(f"{path.name} contains no BED intervals.")
    return width


def _bed_rows(path: Path, limit: int) -> list[dict[str, Any]]:
    """Parse a bounded BED3 to BED6 file into chrom, start, end, name, and strand."""

    _scan_bed(path, limit)
    rows: list[dict[str, Any]] = []
    for number, fields in _data_lines(path):
        if not fields[0]:
            raise ConfigurationError(f"BED line {number} has an empty chromosome name.")
        try:
            start, end = int(fields[1]), int(fields[2])
        except (IndexError, ValueError) as exc:
            raise ConfigurationError(f"BED line {number} needs integer start and end.") from exc
        if not 0 <= start < end:
            raise ConfigurationError(f"BED line {number} must satisfy 0 <= start < end.")
        rows.append(
            {
                "chrom": fields[0],
                "start": start,
                "end": end,
                "name": fields[3] if len(fields) > 3 else None,
                "strand": fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else None,
            }
        )
    return rows


def _gtf_records_frame(gtfparse: Any, path: Path, features: tuple[str, ...], normalize: bool):
    kwargs: dict[str, Any] = {
        "features": set(features) if features else None,
        "attribute_aliases": gtfparse.GENCODE_BIOTYPE_ALIASES if normalize else None,
    }
    try:
        return gtfparse.read_gtf(str(path), **kwargs)
    except Exception as exc:  # gtfparse raises ParsingError and pandas errors for bad files
        raise ConfigurationError(f"gtfparse could not read {path.name}: {exc}") from exc


def _gtf_limit(path: Path, limit: int) -> None:
    """Reject GTF files with more data lines than the caller's max_rows."""

    for count, _ in enumerate(_data_lines(path), start=1):
        if count > limit:
            raise ConfigurationError(f"{path.name} has more than max_rows={limit} data lines.")


# Public functions: bioframe.


def merge_bioframe_intervals(
    source: str | Path,
    *,
    min_distance: int = 0,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge overlapping or nearby BED intervals with bioframe."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    distance = _integer(min_distance, "min_distance", 0, 1_000_000_000)
    bioframe = require("bioframe", "bio-bioframe")
    width = _scan_bed(path, limit)
    frame = bioframe.read_table(str(path), schema=f"bed{width}")
    merged = bioframe.merge(frame, min_dist=distance)
    intervals = _pandas_records(merged, limit, "Merged intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "merge_bioframe_intervals",
        "min_distance": distance,
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "0-based half-open; min_distance=0 also joins bookended intervals",
        },
    }


def closest_bioframe_intervals(
    source: str | Path,
    other: str | Path,
    *,
    k: int = 1,
    ignore_overlaps: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the k closest intervals in another BED file with bioframe."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    neighbours = _integer(k, "k", 1, 100)
    skip_overlaps = _flag(ignore_overlaps, "ignore_overlaps")
    bioframe = require("bioframe", "bio-bioframe")
    left = bioframe.read_table(str(path), schema=f"bed{_scan_bed(path, limit)}")
    right = bioframe.read_table(str(other_path), schema=f"bed{_scan_bed(other_path, limit)}")
    closest = bioframe.closest(left, right, k=neighbours, ignore_overlaps=skip_overlaps)
    records = _pandas_records(closest, limit, "Closest intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "closest_bioframe_intervals",
        "k": neighbours,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "distance is start_b - end_a for separated intervals, 0 when overlapping",
        },
    }


def coverage_bioframe_intervals(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count bases of another BED file covered by each interval, with bioframe."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    bioframe = require("bioframe", "bio-bioframe")
    left = bioframe.read_table(str(path), schema=f"bed{_scan_bed(path, limit)}")
    right = bioframe.read_table(str(other_path), schema=f"bed{_scan_bed(other_path, limit)}")
    coverage = bioframe.coverage(left, right)
    records = _pandas_records(coverage, limit, "Coverage table", ordered=True)
    return {
        "status": "completed",
        "operation": "coverage_bioframe_intervals",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("bioframe"),
            "convention": "coverage is the number of bases of other covered by each interval",
        },
    }


# Public functions: pybedtools.


def merge_pybedtools_intervals(
    source: str | Path,
    *,
    distance: int = 0,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge BED intervals with BEDTools merge through pybedtools."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    gap = _integer(distance, "distance", 0, 1_000_000_000)
    _scan_bed(path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    text = _bedtools_call(
        "merge",
        lambda: str(pybedtools.BedTool(str(path)).sort().merge(d=gap)),
    )
    intervals = [
        {"chrom": fields[0], "start": int(fields[1]), "end": int(fields[2])}
        for fields in _output_lines(text)
    ]
    _checked(intervals, limit, "Merged intervals")
    return {
        "status": "completed",
        "operation": "merge_pybedtools_intervals",
        "distance": gap,
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("pybedtools"),
            "backend_scope": "pybedtools wrapper plus the separately installed bedtools executable",
            "convention": "0-based half-open; bookended intervals are joined at distance 0",
        },
    }


def closest_pybedtools_intervals(
    source: str | Path,
    other: str | Path,
    *,
    k: int = 1,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the k closest intervals in another BED file with BEDTools closest."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    neighbours = _integer(k, "k", 1, 100)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    left = _bedtools_call("closest", lambda: pybedtools.BedTool(str(path)).sort())
    right = _bedtools_call("closest", lambda: pybedtools.BedTool(str(other_path)).sort())
    width_a = left.field_count()
    width_b = right.field_count()
    text = _bedtools_call(
        "closest",
        lambda: str(left.closest(right, k=neighbours, d=True, t="first")),
    )
    records: list[dict[str, Any]] = []
    for fields in _output_lines(text):
        if len(fields) != width_a + width_b + 1:
            raise ConfigurationError("Unexpected BEDTools closest output width.")
        missing = fields[width_a] == "."
        records.append(
            {
                "a": fields[:width_a],
                "b": None if missing else fields[width_a : width_a + width_b],
                "distance": None if missing else int(fields[-1]),
            }
        )
    _checked(records, limit, "Closest intervals")
    return {
        "status": "completed",
        "operation": "closest_pybedtools_intervals",
        "k": neighbours,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "backend_scope": "pybedtools wrapper plus the separately installed bedtools executable",
            "convention": "BEDTools distance counts one extra base for separated intervals",
        },
    }


def coverage_pybedtools_intervals(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count overlapping features, covered bases, and covered fraction with BEDTools coverage."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    pybedtools = require("pybedtools", "bio-pybedtools")
    left = _bedtools_call("coverage", lambda: pybedtools.BedTool(str(path)).sort())
    right = _bedtools_call("coverage", lambda: pybedtools.BedTool(str(other_path)).sort())
    width = left.field_count()
    text = _bedtools_call("coverage", lambda: str(left.coverage(right)))
    records: list[dict[str, Any]] = []
    for fields in _output_lines(text):
        if len(fields) != width + 4:
            raise ConfigurationError("Unexpected BEDTools coverage output width.")
        records.append(
            {
                "a": fields[:width],
                "overlaps": int(fields[width]),
                "bases_covered": int(fields[width + 1]),
                "length": int(fields[width + 2]),
                "fraction": float(fields[width + 3]),
            }
        )
    _checked(records, limit, "Coverage table")
    return {
        "status": "completed",
        "operation": "coverage_pybedtools_intervals",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pybedtools"),
            "backend_scope": "pybedtools wrapper plus the separately installed bedtools executable",
            "convention": "fraction is bases_covered divided by length",
        },
    }


def extract_pybedtools_sequences(
    source: str | Path,
    fasta: str | Path,
    *,
    use_strand: bool = False,
    max_intervals: int = 1_000,
    max_sequence_length: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Extract FASTA sequences for BED intervals with BEDTools getfasta through pybedtools."""

    path = input_path(source)
    fasta_path = input_path(fasta)
    strand = _flag(use_strand, "use_strand")
    limit = _integer(max_intervals, "max_intervals", 1, 100_000)
    span = _integer(max_sequence_length, "max_sequence_length", 1, 10_000_000)
    rows = _bed_rows(path, limit)
    for row in rows:
        if row["end"] - row["start"] > span:
            raise ConfigurationError(
                f"Interval {row['chrom']}:{row['start']}-{row['end']} exceeds "
                f"max_sequence_length={span}."
            )
    pybedtools = require("pybedtools", "bio-pybedtools")
    result = _bedtools_call(
        "getfasta",
        lambda: pybedtools.BedTool(str(path)).sequence(fi=str(fasta_path), tab=True, s=strand),
    )
    sequences: list[dict[str, Any]] = []
    for name, sequence in (
        line.split("\t", 1) for line in Path(result.seqfn).read_text(encoding="utf-8").splitlines()
    ):
        sequences.append({"name": name, "sequence": sequence, "length": len(sequence)})
    if len(sequences) != len(rows):
        raise ConfigurationError(
            "BEDTools returned fewer sequences than intervals; every interval must lie inside "
            "its FASTA sequence."
        )
    return {
        "status": "completed",
        "operation": "extract_pybedtools_sequences",
        "use_strand": strand,
        "sequences": sequences,
        "sequence_count": len(sequences),
        "provenance": {
            "package": _package("pybedtools"),
            "backend_scope": "pybedtools wrapper plus the separately installed bedtools executable",
            "convention": "names are chrom:start-end; bedtools creates <fasta>.fai if it is absent",
        },
    }


# Public functions: gtfparse.


def read_gtfparse_table(
    source: str | Path,
    *,
    features: Iterable[str] = (),
    normalize_gencode_biotypes: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse a GTF file into records with one column per attribute key, using gtfparse."""

    path = input_path(source)
    feature_names = _text_list(features, "features", maximum=100)
    normalize = _flag(normalize_gencode_biotypes, "normalize_gencode_biotypes")
    limit = _rows_limit(max_rows)
    _gtf_limit(path, limit)
    gtfparse = require("gtfparse", "bio-gtfparse")
    frame = _gtf_records_frame(gtfparse, path, feature_names, normalize)
    records = _gtf_rows(frame, limit, "GTF table")
    return {
        "status": "completed",
        "operation": "read_gtfparse_table",
        "columns": [str(column) for column in frame.columns],
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "GTF coordinates are 1-based and inclusive, as stored in the file",
            "gencode_biotype_aliases": normalize,
        },
    }


def select_gtfparse_records(
    source: str | Path,
    *,
    gene_ids: Iterable[str] = (),
    gene_names: Iterable[str] = (),
    transcript_ids: Iterable[str] = (),
    features: Iterable[str] = (),
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Select GTF records whose gene ID, gene name, or transcript ID is in the given lists."""

    path = input_path(source)
    ids = _text_list(gene_ids, "gene_ids")
    names = _text_list(gene_names, "gene_names")
    transcripts = _text_list(transcript_ids, "transcript_ids")
    if not (ids or names or transcripts):
        raise ConfigurationError("Give at least one of gene_ids, gene_names, or transcript_ids.")
    feature_names = _text_list(features, "features", maximum=100)
    limit = _rows_limit(max_rows)
    _gtf_limit(path, limit)
    gtfparse = require("gtfparse", "bio-gtfparse")
    frame = _gtf_records_frame(gtfparse, path, feature_names, normalize=False)
    rows = _gtf_rows(frame, _MAX_ROWS, "GTF table")
    wanted = {"gene_id": set(ids), "gene_name": set(names), "transcript_id": set(transcripts)}
    selected = [
        row
        for row in rows
        if any(
            row.get(column) is not None and row[column] in values
            for column, values in wanted.items()
        )
    ]
    _checked(selected, limit, "Selected GTF records")
    return {
        "status": "completed",
        "operation": "select_gtfparse_records",
        "records": selected,
        "record_count": len(selected),
        "gene_ids": sorted({str(row["gene_id"]) for row in selected if row.get("gene_id")}),
        "transcript_ids": sorted(
            {str(row["transcript_id"]) for row in selected if row.get("transcript_id")}
        ),
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "a record is kept when a gene ID, gene name, or transcript ID is listed",
        },
    }


def summarize_gtfparse_annotation(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count GTF features, genes, transcripts, gene biotypes, and chromosomes using gtfparse."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    _gtf_limit(path, limit)
    gtfparse = require("gtfparse", "bio-gtfparse")
    frame = _gtf_records_frame(gtfparse, path, (), normalize=False)
    rows = _gtf_rows(frame, _MAX_ROWS, "GTF table")
    biotypes: dict[str, str] = {}
    for row in rows:
        gene = row.get("gene_id")
        biotype = row.get("gene_biotype")
        if gene and biotype and gene not in biotypes:
            biotypes[str(gene)] = str(biotype)
    return {
        "status": "completed",
        "operation": "summarize_gtfparse_annotation",
        "record_count": len(rows),
        "feature_counts": dict(Counter(str(row["feature"]) for row in rows if row.get("feature"))),
        "gene_count": len({row["gene_id"] for row in rows if row.get("gene_id")}),
        "transcript_count": len({row["transcript_id"] for row in rows if row.get("transcript_id")}),
        "gene_biotype_counts": dict(Counter(biotypes.values())),
        "chromosome_counts": dict(
            Counter(str(row["seqname"]) for row in rows if row.get("seqname"))
        ),
        "provenance": {
            "package": _package("gtfparse"),
            "convention": "gene and transcript counts are distinct IDs; biotypes count genes",
        },
    }


# Public functions: gffutils.


def _feature_db(gffutils: Any, path: Path) -> Any:
    try:
        return gffutils.FeatureDB(str(path))
    except Exception as exc:  # sqlite errors for files that are not gffutils databases
        raise ConfigurationError(f"gffutils could not open {path.name}: {exc}") from exc


def create_gffutils_database(
    source: str | Path,
    database: str | Path,
    *,
    max_input_bytes: int = 100_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Create a gffutils SQLite database from a GFF3 or GTF file; the database path must be new."""

    path = input_path(source)
    limit = _integer(max_input_bytes, "max_input_bytes", 1, 10_000_000_000)
    if path.stat().st_size > limit:
        raise ConfigurationError(f"{path.name} exceeds max_input_bytes={limit}.")
    target = Path(_text(str(database), "database", maximum=4096)).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"Database already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"Database directory does not exist: {target.parent}")
    gffutils = require("gffutils", "bio-gffutils")
    try:
        db = gffutils.create_db(str(path), dbfn=str(target), force=False, verbose=False)
    except Exception as exc:  # gffutils raises DuplicateIDError and parser errors
        raise ConfigurationError(f"gffutils could not create a database: {exc}") from exc
    counts = {str(kind): int(db.count_features_of_type(kind)) for kind in sorted(db.featuretypes())}
    return {
        "status": "completed",
        "operation": "create_gffutils_database",
        "database": str(target),
        "input_format": str(db.dialect.get("fmt")),
        "featuretype_counts": counts,
        "feature_count": sum(counts.values()),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "feature coordinates are 1-based and inclusive, as stored",
        },
    }


def list_gffutils_children(
    database: str | Path,
    parent_id: str,
    *,
    featuretype: str | None = None,
    level: int = 1,
    max_features: int = 1_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """List child features of one feature in an existing gffutils database."""

    path = input_path(database)
    parent = _text(parent_id, "parent_id")
    kind = None if featuretype is None else _text(featuretype, "featuretype", maximum=128)
    depth = _integer(level, "level", 1, 10)
    limit = _integer(max_features, "max_features", 1, 5_000)
    gffutils = require("gffutils", "bio-gffutils")
    db = _feature_db(gffutils, path)
    try:
        db[parent]
    except gffutils.FeatureNotFoundError as exc:
        raise ConfigurationError(f"Feature {parent!r} is not in {path.name}.") from exc
    # gffutils reports a GTF-inferred transcript as its own child; that self-link is dropped.
    children = db.children(parent, level=depth, featuretype=kind, order_by="start")
    selected = [
        feature for feature in itertools.islice(children, limit + 2) if feature.id != parent
    ]
    _checked(selected, limit, "Child features")
    records = [
        {
            "id": feature.id,
            "featuretype": feature.featuretype,
            "seqid": feature.seqid,
            "start": int(feature.start),
            "end": int(feature.end),
            "strand": feature.strand,
        }
        for feature in selected
    ]
    return {
        "status": "completed",
        "operation": "list_gffutils_children",
        "parent_id": parent,
        "level": depth,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "feature coordinates are 1-based and inclusive, as stored",
        },
    }


def count_gffutils_featuretypes(
    database: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count features of each feature type in an existing gffutils database."""

    path = input_path(database)
    gffutils = require("gffutils", "bio-gffutils")
    db = _feature_db(gffutils, path)
    counts = {str(kind): int(db.count_features_of_type(kind)) for kind in sorted(db.featuretypes())}
    return {
        "status": "completed",
        "operation": "count_gffutils_featuretypes",
        "featuretype_counts": counts,
        "feature_count": sum(counts.values()),
        "provenance": {"package": _package("gffutils")},
    }


# Public functions: pyBigWig.


def _open_bigwig(pybigwig: Any, path: Path) -> Any:
    handle = pybigwig.open(str(path))
    if handle is None:
        raise ConfigurationError(f"pyBigWig could not open {path.name} as a bigWig file.")
    return handle


def _bigwig_region(
    handle: Any, chromosome: str, start: int, end: int, max_bases: int | None
) -> None:
    lengths = handle.chroms()
    if chromosome not in lengths:
        raise ConfigurationError(f"The bigWig file has no chromosome {chromosome!r}.")
    if not 0 <= start < end <= lengths[chromosome]:
        raise ConfigurationError(
            f"Region {chromosome}:{start}-{end} must satisfy 0 <= start < end <= "
            f"{lengths[chromosome]}."
        )
    if max_bases is not None and end - start > max_bases:
        raise ConfigurationError(f"Region spans {end - start} bases; max_bases={max_bases}.")


def read_pybigwig_values(
    source: str | Path,
    chromosome: str,
    start: int,
    end: int,
    *,
    max_bases: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the bigWig value at every base of a bounded region, with None where no data exists."""

    path = input_path(source)
    chrom = _text(chromosome, "chromosome")
    begin = _integer(start, "start", 0, 10**12)
    stop = _integer(end, "end", 1, 10**12)
    span = _integer(max_bases, "max_bases", 1, 1_000_000)
    pybigwig = require("pyBigWig", "bio-pybigwig")
    handle = _open_bigwig(pybigwig, path)
    try:
        _bigwig_region(handle, chrom, begin, stop, span)
        raw = handle.values(chrom, begin, stop)
    finally:
        handle.close()
    values = [None if value is None or math.isnan(value) else float(value) for value in raw]
    return {
        "status": "completed",
        "operation": "read_pybigwig_values",
        "chromosome": chrom,
        "start": begin,
        "end": stop,
        "values": values,
        "length": len(values),
        "covered_bases": sum(value is not None for value in values),
        "provenance": {
            "package": _package("pyBigWig"),
            "convention": "0-based half-open; one value per base, None where the track has no data",
        },
    }


def read_pybigwig_intervals(
    source: str | Path,
    chromosome: str,
    start: int,
    end: int,
    *,
    max_bases: int = 1_000_000,
    max_intervals: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the stored bigWig intervals that overlap a bounded region."""

    path = input_path(source)
    chrom = _text(chromosome, "chromosome")
    begin = _integer(start, "start", 0, 10**12)
    stop = _integer(end, "end", 1, 10**12)
    span = _integer(max_bases, "max_bases", 1, 10_000_000)
    limit = _integer(max_intervals, "max_intervals", 1, 1_000_000)
    pybigwig = require("pyBigWig", "bio-pybigwig")
    handle = _open_bigwig(pybigwig, path)
    try:
        _bigwig_region(handle, chrom, begin, stop, span)
        raw = handle.intervals(chrom, begin, stop) or ()
    finally:
        handle.close()
    records = [
        {"start": int(interval_start), "end": int(interval_end), "value": float(value)}
        for interval_start, interval_end, value in _checked(list(raw), limit, "bigWig intervals")
    ]
    return {
        "status": "completed",
        "operation": "read_pybigwig_intervals",
        "chromosome": chrom,
        "start": begin,
        "end": stop,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pyBigWig"),
            "convention": "intervals overlap the region and are returned as stored, not clipped",
        },
    }


def summarize_pybigwig_bed_intervals(
    source: str | Path,
    bed: str | Path,
    *,
    statistic: str = "mean",
    max_intervals: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute one bigWig statistic for every interval of a BED file."""

    path = input_path(source)
    bed_path = input_path(bed)
    if not isinstance(statistic, str) or statistic not in _BIGWIG_STATISTICS:
        raise ConfigurationError(f"statistic must be one of {sorted(_BIGWIG_STATISTICS)}.")
    limit = _integer(max_intervals, "max_intervals", 1, 100_000)
    rows = _bed_rows(bed_path, limit)
    pybigwig = require("pyBigWig", "bio-pybigwig")
    handle = _open_bigwig(pybigwig, path)
    records: list[dict[str, Any]] = []
    try:
        for row in rows:
            _bigwig_region(handle, row["chrom"], row["start"], row["end"], None)
            value = handle.stats(row["chrom"], row["start"], row["end"], type=statistic)[0]
            records.append(
                {
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "name": row["name"],
                    "value": _plain(value),
                }
            )
    finally:
        handle.close()
    return {
        "status": "completed",
        "operation": "summarize_pybigwig_bed_intervals",
        "statistic": statistic,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pyBigWig"),
            "convention": "value is None when the interval has no bigWig data",
        },
    }


# Public functions: pyranges1.


def intersect_pyranges1_intervals(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the parts of source intervals that overlap other intervals, using pyranges1."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    pyranges1 = require("pyranges1", "bio-pyranges1")
    intersection = pyranges1.read_bed(path).intersect_overlaps(pyranges1.read_bed(other_path))
    records = _pandas_records(intersection, limit, "Intersected intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "intersect_pyranges1_intervals",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pyranges1"),
            "convention": "pieces are clipped to the overlap; columns use pyranges names",
        },
    }


def nearest_pyranges1_intervals(
    source: str | Path,
    other: str | Path,
    *,
    k: int = 1,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the k nearest other-intervals for each source interval, using pyranges1."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    neighbours = _integer(k, "k", 1, 100)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    pyranges1 = require("pyranges1", "bio-pyranges1")
    nearest = pyranges1.read_bed(path).nearest_ranges(
        pyranges1.read_bed(other_path), k=neighbours, dist_col="Distance"
    )
    records = _pandas_records(nearest, limit, "Nearest intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "nearest_pyranges1_intervals",
        "k": neighbours,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pyranges1"),
            "convention": "Distance follows the pyranges1 convention, the same as BEDTools closest",
        },
    }


def complement_pyranges1_intervals(
    source: str | Path,
    *,
    chrom_sizes: Mapping[str, int],
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the unannotated gaps of each chromosome, from 0 to its length, using pyranges1."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    sizes = _lengths(chrom_sizes, "chrom_sizes")
    _scan_bed(path, limit)
    pyranges1 = require("pyranges1", "bio-pyranges1")
    gaps = pyranges1.read_bed(path).complement_ranges(
        chromsizes=sizes, use_strand=False, include_first_interval=True
    )
    records = _pandas_records(gaps, limit, "Complement intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "complement_pyranges1_intervals",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("pyranges1"),
            "convention": "strands are ignored; 0-based half-open coordinates",
        },
    }


# Public functions: polars-bio.


def merge_polars_bio_intervals(
    source: str | Path,
    *,
    min_distance: int = 0,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge BED intervals with polars-bio, keeping BED coordinates."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    distance = _integer(min_distance, "min_distance", 0, 1_000_000_000)
    _scan_bed(path, limit)
    polars_bio = require("polars_bio", "bio-polarsbio")
    frame = polars_bio.read_bed(str(path), use_zero_based=True)
    # polars-bio joins gaps strictly below min_dist, so the wrapper passes distance + 1.
    merged = polars_bio.merge(frame, min_dist=distance + 1, output_type="polars.DataFrame")
    intervals = _polars_records(merged, limit, "Merged intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "merge_polars_bio_intervals",
        "min_distance": distance,
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("polars-bio"),
            "convention": "0-based half-open; gaps up to min_distance are joined",
        },
    }


def nearest_polars_bio_intervals(
    source: str | Path,
    other: str | Path,
    *,
    k: int = 1,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find the k nearest other-intervals for each source interval, using polars-bio."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    neighbours = _integer(k, "k", 1, 100)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    polars_bio = require("polars_bio", "bio-polarsbio")
    nearest = polars_bio.nearest(
        polars_bio.read_bed(str(path), use_zero_based=True),
        polars_bio.read_bed(str(other_path), use_zero_based=True),
        k=neighbours,
        overlap=True,
        distance=True,
        output_type="polars.DataFrame",
    )
    records = _polars_records(nearest, limit, "Nearest intervals", ordered=True)
    return {
        "status": "completed",
        "operation": "nearest_polars_bio_intervals",
        "k": neighbours,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("polars-bio"),
            "convention": "distance is start_b - end_a for separated intervals, 0 when overlapping",
        },
    }


def coverage_polars_bio_intervals(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count bases of other intervals covered by each source interval, using polars-bio."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    polars_bio = require("polars_bio", "bio-polarsbio")
    coverage = polars_bio.coverage(
        polars_bio.read_bed(str(path), use_zero_based=True),
        polars_bio.read_bed(str(other_path), use_zero_based=True),
        output_type="polars.DataFrame",
    )
    records = _polars_records(coverage, limit, "Coverage table", ordered=True)
    return {
        "status": "completed",
        "operation": "coverage_polars_bio_intervals",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("polars-bio"),
            "convention": "coverage is the number of bases of other covered by each interval",
        },
    }


# Public functions: gtars.


def jaccard_gtars_regions(
    source: str | Path,
    other: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compute the base-pair Jaccard index of two BED region sets, using gtars."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    models = require("gtars.models", "bio-gtars")
    left = models.RegionSet(str(path))
    right = models.RegionSet(str(other_path))
    value = float(left.jaccard(right))
    if not math.isfinite(value):
        raise ConfigurationError("gtars returned a non-finite Jaccard index.")
    return {
        "status": "completed",
        "operation": "jaccard_gtars_regions",
        "jaccard": value,
        "source_regions": len(left),
        "other_regions": len(right),
        "provenance": {
            "package": _package("gtars"),
            "statistic": "base-pair Jaccard: shared covered bases over covered bases in either set",
        },
    }


def gaps_gtars_regions(
    source: str | Path,
    chrom_sizes: Mapping[str, int],
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the gaps between BED regions on each chromosome, bounded by chromosome sizes."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    sizes = _lengths(chrom_sizes, "chrom_sizes")
    _scan_bed(path, limit)
    models = require("gtars.models", "bio-gtars")
    regions = models.RegionSet(str(path))
    gaps = [
        {"chrom": str(region.chr), "start": int(region.start), "end": int(region.end)}
        for region in regions.gaps(sizes)
    ]
    _checked(gaps, limit, "Gaps")
    return {
        "status": "completed",
        "operation": "gaps_gtars_regions",
        "gaps": gaps,
        "gap_count": len(gaps),
        "provenance": {
            "package": _package("gtars"),
            "convention": "0-based half-open; chromosomes without regions give one full-length gap",
        },
    }


def tokenize_gtars_regions(
    source: str | Path,
    universe: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Map BED regions onto a BED3 universe of regions, returning one token per overlap."""

    path = input_path(source)
    universe_path = input_path(universe)
    limit = _rows_limit(max_rows)
    _scan_bed(path, limit)
    if _scan_bed(universe_path, limit) != 3:
        raise ConfigurationError("The tokenization universe must be a BED3 file.")
    tokenizers = require("gtars.tokenizers", "bio-gtars")
    try:
        tokenizer = tokenizers.Tokenizer.from_bed(str(universe_path))
        tokens = [str(token) for token in tokenizer.tokenize(str(path))]
    except Exception as exc:  # gtars raises RuntimeError for unreadable inputs
        raise ConfigurationError(f"gtars could not tokenize the regions: {exc}") from exc
    _checked(tokens, limit, "Tokens")
    token_ids = [int(token_id) for token_id in tokenizer.convert_tokens_to_ids(tokens)]
    return {
        "status": "completed",
        "operation": "tokenize_gtars_regions",
        "tokens": tokens,
        "token_ids": token_ids,
        "token_count": len(tokens),
        "provenance": {
            "package": _package("gtars"),
            "convention": "one token per overlapping universe region; duplicates are kept",
        },
    }


# Public functions: ncls.


def _ncls_groups(numpy: Any, rows: list[dict[str, Any]]) -> dict[str, tuple[Any, Any, Any]]:
    """Group intervals by chromosome as int64 arrays; the ids are the row positions."""

    groups: dict[str, list[int]] = {}
    for position, row in enumerate(rows):
        groups.setdefault(row["chrom"], []).append(position)
    return {
        chrom: (
            numpy.array([rows[i]["start"] for i in members], dtype=numpy.int64),
            numpy.array([rows[i]["end"] for i in members], dtype=numpy.int64),
            numpy.array(members, dtype=numpy.int64),
        )
        for chrom, members in groups.items()
    }


def _ncls_indexes(numpy: Any, ncls: Any, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Build one NCLS64 index per chromosome; each interval id is its row position."""

    return {
        chrom: ncls.NCLS64(starts, ends, ids)
        for chrom, (starts, ends, ids) in _ncls_groups(numpy, rows).items()
    }


def _ncls_pairs(
    numpy: Any,
    ncls: Any,
    target_rows: list[dict[str, Any]],
    query_rows: list[dict[str, Any]],
    operation: str,
    max_pairs: int,
) -> list[tuple[int, int]]:
    """Return (query row, target row) pairs for one NCLS query operation."""

    indexes = _ncls_indexes(numpy, ncls, target_rows)
    pairs: list[tuple[int, int]] = []
    for chrom, (starts, ends, ids) in _ncls_groups(numpy, query_rows).items():
        index = indexes.get(chrom)
        if index is None:
            continue
        query_ids, target_ids = getattr(index, operation)(starts, ends, ids)
        pairs.extend(zip(query_ids.tolist(), target_ids.tolist(), strict=True))
        if len(pairs) > max_pairs:
            raise ConfigurationError(f"Interval pairs exceed max_pairs={max_pairs}.")
    return sorted(pairs)


def _pair_records(
    pairs: list[tuple[int, int]],
    target_rows: list[dict[str, Any]],
    query_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    return [
        {
            "chrom": query_rows[query]["chrom"],
            "query_index": query,
            "query_start": query_rows[query]["start"],
            "query_end": query_rows[query]["end"],
            "target_index": target,
            "target_start": target_rows[target]["start"],
            "target_end": target_rows[target]["end"],
        }
        for query, target in pairs
    ]


def query_ncls_overlaps(
    source: str | Path,
    queries: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    max_pairs: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return every (query, source) pair of overlapping intervals, using an NCLS index."""

    path = input_path(source)
    query_path = input_path(queries)
    limit = _rows_limit(max_rows)
    pair_limit = _integer(max_pairs, "max_pairs", 1, _MAX_ROWS)
    target_rows = _bed_rows(path, limit)
    query_rows = _bed_rows(query_path, limit)
    numpy = require("numpy", "bio-ncls")
    ncls = require("ncls", "bio-ncls")
    pairs = _ncls_pairs(numpy, ncls, target_rows, query_rows, "all_overlaps_both", pair_limit)
    records = _pair_records(pairs, target_rows, query_rows)
    return {
        "status": "completed",
        "operation": "query_ncls_overlaps",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("ncls"),
            "convention": "0-based half-open; intervals on different chromosomes never pair",
        },
    }


def count_ncls_overlaps(
    source: str | Path,
    queries: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count the source intervals that overlap each query interval, using an NCLS index."""

    path = input_path(source)
    query_path = input_path(queries)
    limit = _rows_limit(max_rows)
    target_rows = _bed_rows(path, limit)
    query_rows = _bed_rows(query_path, limit)
    numpy = require("numpy", "bio-ncls")
    ncls = require("ncls", "bio-ncls")
    indexes = _ncls_indexes(numpy, ncls, target_rows)
    counts = numpy.zeros(len(query_rows), dtype=numpy.int64)
    for chrom, (starts, ends, ids) in _ncls_groups(numpy, query_rows).items():
        index = indexes.get(chrom)
        if index is not None:
            query_ids, _ = index.all_overlaps_both(starts, ends, ids)
            counts += numpy.bincount(query_ids, minlength=len(query_rows))
    records = [
        {
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "name": row["name"],
            "overlap_count": int(count),
        }
        for row, count in zip(query_rows, counts.tolist(), strict=True)
    ]
    return {
        "status": "completed",
        "operation": "count_ncls_overlaps",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("ncls"),
            "convention": "one count per query row, in input order",
        },
    }


def find_ncls_containing_intervals(
    source: str | Path,
    queries: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    max_pairs: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the source intervals that fully contain each query interval, using NCLS."""

    path = input_path(source)
    query_path = input_path(queries)
    limit = _rows_limit(max_rows)
    pair_limit = _integer(max_pairs, "max_pairs", 1, _MAX_ROWS)
    target_rows = _bed_rows(path, limit)
    query_rows = _bed_rows(query_path, limit)
    numpy = require("numpy", "bio-ncls")
    ncls = require("ncls", "bio-ncls")
    pairs = _ncls_pairs(numpy, ncls, target_rows, query_rows, "all_containments_both", pair_limit)
    records = _pair_records(pairs, target_rows, query_rows)
    return {
        "status": "completed",
        "operation": "find_ncls_containing_intervals",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("ncls"),
            "convention": "a source interval contains a query when it spans the query endpoints",
        },
    }


# Public functions: genomicranges.


def _granges(genomicranges: Any, pandas: Any, rows: list[dict[str, Any]]) -> Any:
    # genomicranges stores 1-based inclusive ranges; BED starts are shifted by one.
    frame = pandas.DataFrame(
        {
            "seqnames": [row["chrom"] for row in rows],
            "starts": [row["start"] + 1 for row in rows],
            "ends": [row["end"] for row in rows],
            "strand": [row["strand"] or "*" for row in rows],
        }
    )
    return genomicranges.GenomicRanges.from_pandas(frame)


def _granges_records(granges: Any, limit: int) -> list[dict[str, Any]]:
    frame = granges.to_pandas()
    _checked(frame, limit, "GenomicRanges result")
    records = []
    for row in frame.to_dict("records"):
        start = int(row["starts"]) - 1
        if start < 0:
            raise ConfigurationError("A genomic range extends before position 0.")
        strand = _plain(row["strand"])
        records.append(
            {
                "chrom": _plain(row["seqnames"]),
                "start": start,
                "end": int(row["ends"]),
                "strand": None if strand in {None, "*"} else strand,
            }
        )
    return records


def reduce_genomicranges_intervals(
    source: str | Path,
    *,
    min_gap_width: int = 1,
    use_strand: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge overlapping or adjacent ranges with genomicranges reduce."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    gap = _integer(min_gap_width, "min_gap_width", 0, 1_000_000)
    strand = _flag(use_strand, "use_strand")
    _scan_bed(path, limit)
    genomicranges = require("genomicranges", "bio-genomicranges")
    pandas = require("pandas", "bio-genomicranges")
    reduced = _granges(genomicranges, pandas, _bed_rows(path, limit)).reduce(
        min_gap_width=gap, ignore_strand=not strand
    )
    intervals = _granges_records(reduced, limit)
    return {
        "status": "completed",
        "operation": "reduce_genomicranges_intervals",
        "min_gap_width": gap,
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("genomicranges"),
            "convention": "min_gap_width=1 joins overlapping and bookended ranges; 0-based output",
        },
    }


def intersect_genomicranges_intervals(
    source: str | Path,
    other: str | Path,
    *,
    use_strand: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the positions covered by both BED interval sets, using genomicranges."""

    path = input_path(source)
    other_path = input_path(other)
    limit = _rows_limit(max_rows)
    strand = _flag(use_strand, "use_strand")
    _scan_bed(path, limit)
    _scan_bed(other_path, limit)
    genomicranges = require("genomicranges", "bio-genomicranges")
    pandas = require("pandas", "bio-genomicranges")
    left = _granges(genomicranges, pandas, _bed_rows(path, limit))
    right = _granges(genomicranges, pandas, _bed_rows(other_path, limit))
    intersection = left.intersect(right, ignore_strand=not strand)
    intervals = _granges_records(intersection, limit)
    return {
        "status": "completed",
        "operation": "intersect_genomicranges_intervals",
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("genomicranges"),
            "convention": "set intersection of covered positions; 0-based half-open output",
        },
    }


def flank_genomicranges_intervals(
    source: str | Path,
    width: int,
    *,
    side: str = "start",
    use_strand: bool = False,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return flanking ranges of a given width on the start or end side, using genomicranges."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    flank_width = _integer(width, "width", 1, 100_000_000)
    if side not in {"start", "end"}:
        raise ConfigurationError("side must be 'start' or 'end'.")
    strand = _flag(use_strand, "use_strand")
    _scan_bed(path, limit)
    genomicranges = require("genomicranges", "bio-genomicranges")
    pandas = require("pandas", "bio-genomicranges")
    flanks = _granges(genomicranges, pandas, _bed_rows(path, limit)).flank(
        width=flank_width, start=side == "start", ignore_strand=not strand
    )
    intervals = _granges_records(flanks, limit)
    return {
        "status": "completed",
        "operation": "flank_genomicranges_intervals",
        "width": flank_width,
        "side": side,
        "intervals": intervals,
        "interval_count": len(intervals),
        "provenance": {
            "package": _package("genomicranges"),
            "convention": "flanks are not clipped to chromosome length; 0-based half-open output",
        },
    }


__all__ = [
    "closest_bioframe_intervals",
    "closest_pybedtools_intervals",
    "complement_pyranges1_intervals",
    "count_gffutils_featuretypes",
    "count_ncls_overlaps",
    "coverage_bioframe_intervals",
    "coverage_polars_bio_intervals",
    "coverage_pybedtools_intervals",
    "create_gffutils_database",
    "extract_pybedtools_sequences",
    "find_ncls_containing_intervals",
    "flank_genomicranges_intervals",
    "gaps_gtars_regions",
    "intersect_genomicranges_intervals",
    "intersect_pyranges1_intervals",
    "jaccard_gtars_regions",
    "list_gffutils_children",
    "merge_bioframe_intervals",
    "merge_polars_bio_intervals",
    "merge_pybedtools_intervals",
    "nearest_polars_bio_intervals",
    "nearest_pyranges1_intervals",
    "query_ncls_overlaps",
    "read_gtfparse_table",
    "read_pybigwig_intervals",
    "read_pybigwig_values",
    "reduce_genomicranges_intervals",
    "select_gtfparse_records",
    "summarize_gtfparse_annotation",
    "summarize_pybigwig_bed_intervals",
    "tokenize_gtars_regions",
]

_register_public_verbose(globals(), tuple(__all__))
