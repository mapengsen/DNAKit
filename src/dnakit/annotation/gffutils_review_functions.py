"""UCSC bins, GFF iteration, Biopython and pybedtools conversion, canonical transcripts, GFF
sanitising, feature merging, and dialect inference through gffutils (review group G18g).

Functions read local GFF3 or GTF files, gffutils databases (.db), and local FASTA files. A GFF
or GTF path is loaded into an in-memory gffutils database, so no database file is created. The
sanitising function writes one new GFF file and never modifies its input. Coordinates are 1-based
inclusive as stored in GFF, except where a function returns BED or Biopython output (0-based
half-open).
gffutils is imported only when a function runs; a missing package raises BackendUnavailableError
with the install hint.
"""

from __future__ import annotations

import contextlib
import io
import shutil
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import (
    _DEFAULT_ROWS,
    _checked,
    _package,
    _plain,
    _rows_limit,
    _text,
    _text_list,
)
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_MAX_COORDINATE = 2_000_000_000
_DATABASE_SUFFIXES = frozenset({".db", ".sqlite", ".sqlite3"})
_MERGE_CRITERIA = frozenset(
    {
        "exact_coordinates_only",
        "feature_type",
        "overlap_any_inclusive",
        "overlap_any_threshold",
        "overlap_end_inclusive",
        "overlap_end_threshold",
        "overlap_start_inclusive",
        "overlap_start_threshold",
        "seqid",
        "strand",
    }
)
_THRESHOLD_CRITERIA = frozenset(
    {"overlap_any_threshold", "overlap_end_threshold", "overlap_start_threshold"}
)
_BIO_STRAND = {"+": 1, "-": -1, ".": None, "?": 0}
_RESERVED_QUALIFIERS = frozenset({"seqid", "source", "score", "frame"})


def _open_database(gffutils: Any, path: Path) -> Any:
    """Open a .db file, or load a GFF or GTF file into an in-memory database."""

    try:
        if path.suffix in _DATABASE_SUFFIXES:
            return gffutils.FeatureDB(str(path))
        return gffutils.create_db(str(path), ":memory:", verbose=False)
    except Exception as exc:  # sqlite and gffutils errors for files that are neither
        raise ConfigurationError(f"gffutils could not open {path.name}: {exc}") from exc


def _feature_record(feature: Any) -> dict[str, Any]:
    attributes = {
        str(key): sorted(str(item) for item in values)
        for key, values in sorted(dict(feature.attributes).items())
    }
    return {
        "seqid": feature.seqid,
        "source": feature.source,
        "featuretype": feature.featuretype,
        "start": int(feature.start),
        "end": int(feature.end),
        "score": _plain(feature.score),
        "strand": feature.strand,
        "frame": _plain(feature.frame),
        "attributes": attributes,
    }


def _interval_record(interval: Any) -> dict[str, Any]:
    return {
        "chrom": interval.chrom,
        "start": int(interval.start),
        "end": int(interval.stop),
        "name": interval.name,
        "score": _plain(interval.score),
        "strand": interval.strand,
        "featuretype": interval.fields[2] if len(interval.fields) > 2 else None,
    }


def _dialect_value(value: Any) -> Any:
    if isinstance(value, (set, frozenset)):
        return sorted(_plain(item) for item in value)
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return _plain(value)


def read_gffutils_records(
    source: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Iterate the features of a GFF or GTF file with gffutils.DataIterator, in file order."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    gffutils = require("gffutils", "bio-gffutils")
    records: list[dict[str, Any]] = []
    try:
        for feature in gffutils.DataIterator(str(path)):
            records.append(_feature_record(feature))
            if len(records) > limit:
                break
    except Exception as exc:  # gffutils raises ValueError for malformed feature lines
        raise ConfigurationError(f"gffutils could not read {path.name}: {exc}") from exc
    _checked(records, limit, "GFF features")
    return {
        "status": "completed",
        "operation": "read_gffutils_records",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "start and end are 1-based inclusive; attribute values are sorted lists; "
            "records keep file order",
        },
    }


def bin_gffutils_interval(
    start: int,
    stop: int,
    *,
    fmt: str = "gff",
    one: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return UCSC genomic bins for an interval with gffutils.bins.bins."""

    if fmt not in {"gff", "bed"}:
        raise ConfigurationError("fmt must be 'gff' (1-based start) or 'bed' (0-based start).")
    minimum = 1 if fmt == "gff" else 0
    first = _integer_at_least(start, "start", minimum)
    last = _integer_at_least(stop, "stop", first)
    single = _flag(one, "one")
    bins_module = require("gffutils.bins", "bio-gffutils")
    found = bins_module.bins(first, last, fmt=fmt, one=single)
    values = [int(found)] if single else sorted(int(value) for value in found)
    return {
        "status": "completed",
        "operation": "bin_gffutils_interval",
        "start": first,
        "stop": last,
        "fmt": fmt,
        "one": single,
        "bins": values,
        "provenance": {
            "package": _package("gffutils"),
            "convention": "UCSC binning; one=True gives the smallest bin holding the interval; "
            "stop is used as given, as upstream does",
        },
    }


def _integer_at_least(value: object, name: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value:
        raise ConfigurationError(f"{name} must be an integer of at least {minimum}.")
    if value > _MAX_COORDINATE:
        raise ConfigurationError(f"{name} must be at most {_MAX_COORDINATE}.")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def asinterval_gffutils_feature(
    database: str | Path,
    feature_id: str,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert one feature to a pybedtools Interval with gffutils.helpers.asinterval."""

    path = input_path(database)
    wanted = _text(feature_id, "feature_id")
    gffutils = require("gffutils", "bio-gffutils")
    helpers = require("gffutils.helpers", "bio-gffutils")
    require("pybedtools", "bio-pybedtools")
    db = _open_database(gffutils, path)
    try:
        feature = db[wanted]
    except Exception as exc:  # gffutils raises FeatureNotFoundError for unknown IDs
        raise ConfigurationError(f"feature_id not found in {path.name}: {wanted}") from exc
    interval = helpers.asinterval(feature)
    record = _interval_record(interval)
    return {
        "status": "completed",
        "operation": "asinterval_gffutils_feature",
        "feature_id": wanted,
        **record,
        "fields": [str(field) for field in interval.fields],
        "provenance": {
            "package": _package("gffutils"),
            "convention": "pybedtools Interval for the GFF line: start is 0-based, end is the "
            "GFF end, and name is the feature ID",
        },
    }


def features_gffutils_to_bedtool(
    database: str | Path,
    *,
    featuretype: str | None = None,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert database features (optionally one feature type) to intervals with to_bedtool."""

    path = input_path(database)
    limit = _rows_limit(max_rows)
    kind = None if featuretype is None else _text(featuretype, "featuretype", maximum=64)
    gffutils = require("gffutils", "bio-gffutils")
    integration = require("gffutils.pybedtools_integration", "bio-gffutils")
    db = _open_database(gffutils, path)
    if kind is None:
        features = db.all_features(order_by=("seqid", "start"))
    else:
        features = db.features_of_type(kind, order_by=("seqid", "start"))
    records: list[dict[str, Any]] = []
    try:
        for interval in integration.to_bedtool(features):
            records.append(_interval_record(interval))
            if len(records) > limit:
                break
    except Exception as exc:  # pybedtools and gffutils raise for malformed features
        raise ConfigurationError(f"gffutils could not convert {path.name}: {exc}") from exc
    _checked(records, limit, "features")
    records.sort(key=lambda r: (r["chrom"], r["start"], r["end"], str(r["name"])))
    return {
        "status": "completed",
        "operation": "features_gffutils_to_bedtool",
        "featuretype": kind,
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "BED-style start is 0-based and end is exclusive; name is the feature ID",
        },
    }


def infer_gffutils_attribute_dialect(attributes: str, *, verbose: bool = True) -> dict[str, Any]:
    """Infer the GFF3 or GTF dialect from the attribute column text with gffutils.helpers."""

    text = _text(attributes, "attributes", maximum=65_536)
    helpers = require("gffutils.helpers", "bio-gffutils")
    try:
        dialect = helpers.infer_dialect(text)
    except Exception as exc:  # gffutils raises for attribute text it cannot split
        raise ConfigurationError(f"gffutils could not infer a dialect: {exc}") from exc
    plain = {str(key): _dialect_value(value) for key, value in sorted(dialect.items())}
    return {
        "status": "completed",
        "operation": "infer_gffutils_attribute_dialect",
        "fmt": plain.get("fmt"),
        "dialect": plain,
        "provenance": {
            "package": _package("gffutils"),
            "convention": "fmt is 'gff3' or 'gtf'; separators are the literal text gffutils "
            "uses to write the attributes",
        },
    }


def canonical_gffutils_transcript_sequences(
    database: str | Path,
    fasta: str | Path,
    *,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pick one canonical transcript per gene and return its exon sequence with gffutils."""

    db_path = input_path(database)
    fasta_path = input_path(fasta)
    limit = _rows_limit(max_rows)
    gffutils = require("gffutils", "bio-gffutils")
    helpers = require("gffutils.helpers", "bio-gffutils")
    require("pyfaidx", "bio-pyfaidx")
    db = _open_database(gffutils, db_path)
    genes = list(db.features_of_type("gene", order_by=("seqid", "start")))
    if not genes:
        raise ConfigurationError(f"{db_path.name} has no gene features.")
    for gene in genes:
        if not any(True for _ in db.children(gene, level=1)):
            raise ConfigurationError(
                f"gene {gene.id} has no transcript children; canonical selection needs them."
            )
    # Upstream prints each chosen transcript; the output is discarded and stdout is restored.
    chatter = io.StringIO()
    try:
        with contextlib.redirect_stdout(chatter):
            pairs = list(helpers.canonical_transcripts(db, str(fasta_path)))
    except Exception as exc:  # pyfaidx raises for FASTA files it cannot index
        raise ConfigurationError(f"canonical transcript selection failed: {exc}") from exc
    _checked(pairs, limit, "transcripts")
    records = []
    for transcript, sequence in pairs:
        records.append(
            {
                "transcript_id": transcript.id,
                "gene_ids": sorted(parent.id for parent in db.parents(transcript, level=1)),
                "seqid": transcript.seqid,
                "start": int(transcript.start),
                "end": int(transcript.end),
                "strand": transcript.strand,
                "length": len(sequence),
                "sequence": sequence,
            }
        )
    records.sort(key=lambda r: r["transcript_id"])
    return {
        "status": "completed",
        "operation": "canonical_gffutils_transcript_sequences",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "the transcript with the longest CDS wins, else the longest transcript; "
            "sequence joins CDS and UTR exons in transcript order (reverse complement on -)",
        },
    }


def gffutils_line_to_seqfeature(line: str, *, verbose: bool = True) -> dict[str, Any]:
    """Convert one nine-column GFF line to a Biopython SeqFeature summary."""

    text = _text(line, "line", maximum=65_536)
    if text.count("\t") != 8:
        raise ConfigurationError("line must be a tab-separated GFF record with nine columns.")
    require("Bio.SeqFeature", "bio-formats")
    integration = require("gffutils.biopython_integration", "bio-gffutils")
    try:
        feature = integration.to_seqfeature(text)
    except Exception as exc:  # gffutils raises ValueError for malformed lines
        raise ConfigurationError(f"gffutils could not convert the line: {exc}") from exc
    location = feature.location
    return {
        "status": "completed",
        "operation": "gffutils_line_to_seqfeature",
        "type": feature.type,
        "id": feature.id,
        "location": {
            "start": int(location.start),
            "end": int(location.end),
            "strand": _plain(location.strand),
        },
        "qualifiers": {
            str(key): [str(item) for item in values]
            for key, values in sorted(feature.qualifiers.items())
        },
        "provenance": {
            "package": _package("gffutils"),
            "convention": "location start is 0-based and end is exclusive, as in Biopython; "
            "strand is 1, -1, 0 or None",
        },
    }


def seqfeature_to_gffutils_line(
    seqid: str,
    feature_type: str,
    start: int,
    end: int,
    *,
    strand: str = "+",
    feature_id: str | None = None,
    source: str = ".",
    score: str = ".",
    frame: str = ".",
    attributes: Mapping[str, str | list[str]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a Biopython SeqFeature and convert it back to one GFF line with gffutils."""

    chrom = _text(seqid, "seqid", maximum=256)
    kind = _text(feature_type, "feature_type", maximum=128)
    begin = _integer_at_least(start, "start", 0)
    stop = _integer_at_least(end, "end", 1)
    if begin >= stop:
        raise ConfigurationError("start must be less than end (start is 0-based).")
    if strand not in _BIO_STRAND:
        raise ConfigurationError("strand must be one of '+', '-', '.' or '?'.")
    qualifiers: dict[str, list[str]] = {
        "seqid": [chrom],
        "source": [_text(source, "source", maximum=256)],
        "score": [_text(score, "score", maximum=64)],
        "frame": [_text(frame, "frame", maximum=8)],
    }
    if feature_id is not None:
        qualifiers["ID"] = [_text(feature_id, "feature_id")]
    for key, value in (attributes or {}).items():
        name = _text(key, "attribute name", maximum=128)
        if name in _RESERVED_QUALIFIERS:
            raise ConfigurationError(f"attribute {name!r} is a GFF column; use its own argument.")
        values = [value] if isinstance(value, str) else list(_text_list(value, name))
        qualifiers[name] = [_text(item, name) for item in values]
    bio = require("Bio.SeqFeature", "bio-formats")
    integration = require("gffutils.biopython_integration", "bio-gffutils")
    kwargs = {} if feature_id is None else {"id": feature_id}
    feature = bio.SeqFeature(
        bio.FeatureLocation(begin, stop, strand=_BIO_STRAND[strand]),
        type=kind,
        qualifiers=qualifiers,
        **kwargs,
    )
    try:
        converted = integration.from_seqfeature(feature)
    except Exception as exc:  # gffutils rejects SeqFeature values it cannot write as GFF
        raise ConfigurationError(f"gffutils could not write the feature: {exc}") from exc
    return {
        "status": "completed",
        "operation": "seqfeature_to_gffutils_line",
        "line": str(converted),
        "seqid": chrom,
        "featuretype": kind,
        "start": begin + 1,
        "end": stop,
        "strand": strand,
        "provenance": {
            "package": _package("gffutils"),
            "convention": "input start is 0-based and the GFF line start is 1-based; end is "
            "inclusive in both",
        },
    }


def sanitize_gffutils_file(
    source: str | Path,
    output: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a sanitised copy of a GFF file with gffutils.helpers.sanitize_gff_file."""

    path = input_path(source)
    target = Path(output).expanduser().resolve()
    if target == path:
        raise ConfigurationError("output must differ from the source file.")
    if target.exists():
        raise ConfigurationError(f"output already exists: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"output directory does not exist: {target.parent}")
    helpers = require("gffutils.helpers", "bio-gffutils")
    shutil.copyfile(path, target)
    try:
        # in_place rewrites the copy only; the source file is never opened for writing.
        helpers.sanitize_gff_file(str(target), in_memory=True, in_place=True)
    except Exception as exc:  # gffutils raises for files it cannot load
        target.unlink(missing_ok=True)
        raise ConfigurationError(f"gffutils could not sanitise {path.name}: {exc}") from exc
    rows = [
        fields
        for fields in (line.split("\t") for line in target.read_text(encoding="utf-8").splitlines())
        if len(fields) == 9 and not fields[0].startswith("#")
    ]
    return {
        "status": "completed",
        "operation": "sanitize_gffutils_file",
        "path": str(target),
        "feature_count": len(rows),
        "gene_count": sum(1 for fields in rows if fields[2] == "gene"),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "features outside a gene are dropped; reversed coordinates are swapped; "
            "each record gets a gid attribute with its gene ID; the header has a timestamp",
        },
    }


def merge_gffutils_features(
    source: str | Path,
    *,
    featuretype: str,
    criteria: tuple[str, ...] | list[str] = (
        "seqid",
        "overlap_end_inclusive",
        "strand",
        "feature_type",
    ),
    threshold: int = 0,
    max_rows: int = _DEFAULT_ROWS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge features of one type whose merge criteria all hold, with gffutils merge_criteria."""

    path = input_path(source)
    limit = _rows_limit(max_rows)
    kind = _text(featuretype, "featuretype", maximum=64)
    names = _text_list(criteria, "criteria", maximum=10)
    if not names:
        raise ConfigurationError("criteria must name at least one merge criterion.")
    if len(set(names)) != len(names):
        raise ConfigurationError("criteria must not repeat a merge criterion.")
    for name in names:
        if name not in _MERGE_CRITERIA:
            raise ConfigurationError(f"unknown merge criterion: {name}")
    gap = _integer_at_least(threshold, "threshold", 0)
    gffutils = require("gffutils", "bio-gffutils")
    merge_criteria = require("gffutils.merge_criteria", "bio-gffutils")
    callbacks = [
        getattr(merge_criteria, name)(gap)
        if name in _THRESHOLD_CRITERIA
        else getattr(merge_criteria, name)
        for name in names
    ]
    db = _open_database(gffutils, path)
    features = db.features_of_type(kind, order_by=("seqid", "start"))
    merged: list[dict[str, Any]] = []
    try:
        for feature in db.merge(features, merge_criteria=callbacks, multiline=False):
            merged.append(
                {
                    "id": feature.id,
                    "seqid": feature.seqid,
                    "featuretype": feature.featuretype,
                    "start": int(feature.start),
                    "end": int(feature.end),
                    "strand": feature.strand,
                    "member_ids": sorted(child.id for child in feature.children),
                }
            )
            if len(merged) > limit:
                break
    except Exception as exc:  # gffutils raises for merge criteria it cannot evaluate
        raise ConfigurationError(f"gffutils could not merge {kind} features: {exc}") from exc
    _checked(merged, limit, "merged features")
    merged.sort(key=lambda r: (r["seqid"], r["start"], r["end"], r["id"]))
    return {
        "status": "completed",
        "operation": "merge_gffutils_features",
        "featuretype": kind,
        "criteria": list(names),
        "threshold": gap,
        "merged": merged,
        "merged_count": len(merged),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "features are visited in seqid and start order; a feature joins the "
            "current group when every criterion holds; member_ids lists the grouped features",
        },
    }


__all__ = [
    "asinterval_gffutils_feature",
    "bin_gffutils_interval",
    "canonical_gffutils_transcript_sequences",
    "features_gffutils_to_bedtool",
    "gffutils_line_to_seqfeature",
    "infer_gffutils_attribute_dialect",
    "merge_gffutils_features",
    "read_gffutils_records",
    "sanitize_gffutils_file",
    "seqfeature_to_gffutils_line",
]

_register_public_verbose(globals(), tuple(__all__))
