"""Offline genomepy helpers for annotation tables, FASTA filters, gap and size files, and sampling.

The functions call genomepy's annotation, FASTA and genome code on local files. Input files are
never modified. Output paths must not exist yet, so existing files are never overwritten.
pyfaidx and genomepy write index sidecars (``.fai``, ``.sizes``, ``.gaps.bed``) beside the input
FASTA; pass a copy when the source directory is read-only. Downloads, provider catalogues and
aligner index builders are not wrapped here; ``gene_backends`` covers the genome installer.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-genomepy]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import random
import re
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-genomepy"
_DISTRIBUTION = "genomepy"
_FORMATS = ("gtf", "bed12")
_DEFAULT_MAX_BYTES = 64_000_000
_MAX_BYTES = 1_000_000_000
_MAX_RECORDS = 1_000_000
_MAX_REGIONS = 100_000
_MAX_LENGTH = 1_000_000
_MAX_CHROMS = 10_000
_MAX_PATTERN = 256
_MAX_SEED = 2**32 - 1


def _text(value: object, name: str, *, limit: int = 256) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string.")
    text = value.strip()
    if not text or len(text) > limit or "\x00" in text:
        raise ConfigurationError(
            f"{name} must be a non-empty string of at most {limit} characters."
        )
    return text


def _int(value: object, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer from {minimum} to {maximum}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _sequence_of(values: object, name: str, *, limit: int) -> list[object]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ConfigurationError(f"{name} must be a list or tuple.")
    if len(values) > limit:
        raise ConfigurationError(f"{name} has {len(values)} items; the limit is {limit}.")
    return list(values)


def _plain(value: object) -> Any:
    """Convert pandas and numpy scalars to JSON-compatible Python objects."""
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        value = value.item()
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return value if value == value and abs(value) != float("inf") else None
    return str(value)


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version(_DISTRIBUTION)
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    return {"package": f"{_DISTRIBUTION} {version}", **notes}


def _document(operation: str, provenance: dict[str, Any], **fields: object) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": provenance}


def _upstream(operation: str, call: Callable[[], Any]) -> Any:
    try:
        return call()
    except ConfigurationError:
        raise
    except Exception as exc:
        raise BackendExecutionError(f"genomepy {operation} failed: {exc}") from exc


def _new_output(value: object, name: str, source: Path | None = None) -> Path:
    """Return a fresh output path in an existing directory; never overwrite an input or file."""
    path = Path(_text(str(value), name, limit=4096)).expanduser().resolve()
    if source is not None and path == source:
        raise ConfigurationError(f"{name} must differ from the input file.")
    if path.exists():
        raise ConfigurationError(f"{name} already exists; choose a new file: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"Directory of {name} does not exist: {path.parent}")
    return path


def _columns(file_format: str) -> list[str]:
    if file_format not in _FORMATS:
        raise ConfigurationError(f"file_format must be one of {', '.join(_FORMATS)}.")
    annotation = require("genomepy.annotation.utils", _EXTRA)
    return list(annotation.GTF_FORMAT if file_format == "gtf" else annotation.BED12_FORMAT)


def _annotation_frame(
    records: object, file_format: str, *, limit: int, contig_text: bool = False
) -> Any:
    """Validate records against the column layout of one annotation format and build a frame."""
    columns = _columns(file_format)
    items = _sequence_of(records, "records", limit=limit)
    rows = []
    for position, item in enumerate(items):
        if not isinstance(item, Mapping) or set(item) != set(columns):
            raise ConfigurationError(f"records[{position}] must have exactly the keys {columns}.")
        row = []
        for column in columns:
            value = item[column]
            if isinstance(value, bool) or not isinstance(value, (str, int)):
                raise ConfigurationError(f"records[{position}].{column} must be a str or int.")
            if isinstance(value, str) and ("\t" in value or "\n" in value or "\x00" in value):
                raise ConfigurationError(f"records[{position}].{column} contains a tab or newline.")
            if contig_text and column == columns[0] and not isinstance(value, str):
                raise ConfigurationError(
                    f"records[{position}].{column} must be a string to filter."
                )
            row.append(value)
        rows.append(row)
    pandas = require("pandas", _EXTRA)
    return pandas.DataFrame(rows, columns=columns)


def read_annotation_records(
    annotation_path: str | Path,
    *,
    max_bytes: int = _DEFAULT_MAX_BYTES,
    max_records: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read a GTF (9 columns) or BED12 (12 columns) file into one record per line."""
    path = input_path(annotation_path)
    size_limit = _int(max_bytes, "max_bytes", 1, _MAX_BYTES)
    record_limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    size = path.stat().st_size
    if size > size_limit:
        raise ConfigurationError(f"Annotation file has {size} bytes; max_bytes is {size_limit}.")
    annotation = require("genomepy.annotation.utils", _EXTRA)
    try:
        frame = annotation.read_annot(str(path))
    except ValueError as exc:
        raise ConfigurationError(f"Annotation layout is not GTF or BED12: {exc}") from exc
    except Exception as exc:
        raise BackendExecutionError(f"genomepy read_annot failed: {exc}") from exc
    if len(frame) > record_limit:
        raise ConfigurationError(
            f"Annotation has {len(frame)} records; max_records is {record_limit}."
        )
    file_format = "gtf" if len(frame.columns) == 9 else "bed12"
    records = _plain(frame.to_dict("records"))
    return _document(
        "read_annotation_records",
        _provenance(input_format=file_format, coordinate_system="as stored in the file"),
        file_format=file_format,
        columns=[str(column) for column in frame.columns],
        record_count=len(records),
        records=records,
    )


def write_annotation_records(
    records: Sequence[Mapping[str, object]],
    output_path: str | Path,
    *,
    file_format: str = "gtf",
    max_records: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write GTF or BED12 records to a new tab-separated file without a header line."""
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    frame = _annotation_frame(records, file_format, limit=limit)
    target = _new_output(output_path, "output_path")
    annotation = require("genomepy.annotation.utils", _EXTRA)
    _upstream("write_annot", lambda: annotation.write_annot(frame, str(target)))
    return _document(
        "write_annotation_records",
        _provenance(output_format="tab-separated, no header"),
        output_path=str(target),
        file_format=file_format,
        record_count=len(frame),
    )


def filter_annotation_contigs(
    records: Sequence[Mapping[str, object]],
    pattern: str,
    *,
    file_format: str = "gtf",
    invert_match: bool = False,
    max_records: int = _MAX_RECORDS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep annotation records whose contig name matches a regular expression at its start."""
    text = _text(pattern, "pattern", limit=_MAX_PATTERN)
    try:
        re.compile(text)
    except re.error as exc:
        raise ConfigurationError(f"pattern is not a valid regular expression: {exc}") from exc
    flag = _flag(invert_match, "invert_match")
    limit = _int(max_records, "max_records", 1, _MAX_RECORDS)
    frame = _annotation_frame(records, file_format, limit=limit, contig_text=True)
    annotation = require("genomepy.annotation", _EXTRA)
    kept = _upstream(
        "filter_regex", lambda: annotation.filter_regex(frame, text, invert_match=flag)
    )
    kept_records = _plain(kept.to_dict("records"))
    return _document(
        "filter_annotation_contigs",
        _provenance(match="re.match on the contig column, as genomepy does"),
        pattern=text,
        invert_match=flag,
        kept_count=len(kept_records),
        removed_count=len(frame) - len(kept_records),
        records=kept_records,
    )


def write_gap_bed(
    fasta_path: str | Path,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a BED file of N-runs (gaps) in a FASTA file, with 0-based half-open coordinates."""
    source = input_path(fasta_path)
    target = _new_output(output_path, "output_path", source)
    genome = require("genomepy.genome", _EXTRA)
    _upstream("generate_gap_bed", lambda: genome.generate_gap_bed(str(source), str(target)))
    gap_count = sum(1 for line in target.read_text(encoding="utf-8").splitlines() if line)
    return _document(
        "write_gap_bed",
        _provenance(coordinate_system="BED 0-based half-open", side_effect="writes FASTA index"),
        output_path=str(target),
        gap_count=gap_count,
    )


def write_fasta_sizes(
    fasta_path: str | Path,
    output_path: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a two-column sizes file (sequence name and length) for a FASTA file."""
    source = input_path(fasta_path)
    target = _new_output(output_path, "output_path", source)
    genome = require("genomepy.genome", _EXTRA)
    _upstream("generate_fa_sizes", lambda: genome.generate_fa_sizes(str(source), str(target)))
    sequence_count = sum(1 for line in target.read_text(encoding="utf-8").splitlines() if line)
    return _document(
        "write_fasta_sizes",
        _provenance(columns="name, length", side_effect="writes FASTA index"),
        output_path=str(target),
        sequence_count=sequence_count,
    )


def filter_fasta_records(
    fasta_path: str | Path,
    output_path: str | Path,
    *,
    pattern: str = ".*",
    invert_match: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the FASTA records whose header matches a regular expression to a new file."""
    source = input_path(fasta_path)
    target = _new_output(output_path, "output_path", source)
    text = _text(pattern, "pattern", limit=_MAX_PATTERN)
    try:
        re.compile(text)
    except re.error as exc:
        raise ConfigurationError(f"pattern is not a valid regular expression: {exc}") from exc
    flag = _flag(invert_match, "invert_match")
    files = require("genomepy.files", _EXTRA)
    removed = _upstream(
        "filter_fasta",
        lambda: files.filter_fasta(str(source), str(target), regex=text, invert_match=flag),
    )
    kept_count = sum(
        1 for line in target.read_text(encoding="utf-8").splitlines() if line.startswith(">")
    )
    return _document(
        "filter_fasta_records",
        _provenance(header_match="regular expression on the header text"),
        output_path=str(target),
        pattern=text,
        invert_match=flag,
        kept_count=kept_count,
        removed=sorted(str(name) for name in removed),
    )


def sample_random_regions(
    fasta_path: str | Path,
    *,
    n: int = 10,
    length: int = 200,
    chroms: Sequence[str] | None = None,
    max_n: float = 0.1,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample fixed-length regions from a FASTA with genomepy, seeded and free of heavy N runs."""
    source = input_path(fasta_path)
    count = _int(n, "n", 1, _MAX_REGIONS)
    span = _int(length, "length", 1, _MAX_LENGTH)
    seed_value = _int(seed, "seed", 0, _MAX_SEED)
    if isinstance(max_n, bool) or not isinstance(max_n, (int, float)) or not 0 <= max_n <= 1:
        raise ConfigurationError("max_n must be a number from 0 to 1.")
    names = None
    if chroms is not None:
        names = [
            _text(item, "chroms item", limit=256)
            for item in _sequence_of(chroms, "chroms", limit=_MAX_CHROMS)
        ]
    genomepy = require("genomepy", _EXTRA)
    genome = _upstream("Genome", lambda: genomepy.Genome(str(source)))
    if names is not None:
        known = set(genome.keys())
        missing = sorted(name for name in names if name not in known)
        if missing:
            raise ConfigurationError(f"chroms not in the FASTA: {', '.join(missing)}")
    state = random.getstate()
    random.seed(seed_value)
    try:
        sampled = genome.get_random_sequences(
            n=count, length=span, chroms=names, max_n=float(max_n), outtype="list"
        )
    except Exception as exc:
        raise ConfigurationError(f"cannot sample the requested regions: {exc}") from exc
    finally:
        random.setstate(state)
    regions = [{"chrom": chrom, "start": start, "end": end} for chrom, start, end in sampled]
    return _document(
        "sample_random_regions",
        _provenance(
            sampler="genomepy Genome.get_random_sequences",
            seed="global random state is restored afterwards",
            side_effect="writes FASTA index and sizes",
        ),
        seed=seed_value,
        length=span,
        max_n=float(max_n),
        region_count=len(regions),
        regions=regions,
    )


__all__ = [
    "filter_annotation_contigs",
    "filter_fasta_records",
    "read_annotation_records",
    "sample_random_regions",
    "write_annotation_records",
    "write_fasta_sizes",
    "write_gap_bed",
]

_register_public_verbose(globals(), tuple(__all__))
