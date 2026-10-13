"""Transcript BED12 export, TSS intervals, GFF line parsing, and file inspection through gffutils.

Functions read local GFF3 or GTF files and gffutils SQLite databases only. Feature coordinates are
kept as stored (1-based, inclusive) except where a function returns BED output, which is 0-based
half-open. gffutils is imported only when a function runs; a missing package raises
BackendUnavailableError with the install hint.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.backends import _bedtools_call, _package, _plain, _text
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError


def _open_database(gffutils: Any, path: Path) -> Any:
    try:
        return gffutils.FeatureDB(str(path))
    except Exception as exc:  # sqlite and gffutils errors for files that are not databases
        raise ConfigurationError(f"gffutils could not open {path.name}: {exc}") from exc


def to_bed12_gffutils_transcript(
    database: str | Path,
    transcript_id: str,
    *,
    child_type: str = "exon",
    name_field: str = "ID",
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert one transcript and its child features to a BED12 line with gffutils."""

    path = input_path(database)
    feature_id = _text(transcript_id, "transcript_id")
    kind = _text(child_type, "child_type", maximum=64)
    name = _text(name_field, "name_field", maximum=64)
    gffutils = require("gffutils", "bio-gffutils")
    convert = require("gffutils.convert", "bio-gffutils")
    db = _open_database(gffutils, path)
    try:
        feature = db[feature_id]
    except Exception as exc:  # gffutils raises FeatureNotFoundError for unknown IDs
        raise ConfigurationError(f"transcript_id not found in {path.name}: {feature_id}") from exc
    line = _bedtools_call(
        "to_bed12",
        lambda: convert.to_bed12(feature, db, child_type=kind, name_field=name),
    )
    return {
        "status": "completed",
        "operation": "to_bed12_gffutils_transcript",
        "transcript_id": feature_id,
        "bed12": str(line).rstrip("\n"),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "BED12 start is 0-based; blocks come from child_type features",
        },
    }


def tss_gffutils_intervals(
    database: str | Path,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """List one-base transcription start sites from a gffutils database as BED6 rows."""

    path = input_path(database)
    gffutils = require("gffutils", "bio-gffutils")
    integration = require("gffutils.pybedtools_integration", "bio-gffutils")
    db = _open_database(gffutils, path)
    sites = _bedtools_call("tsses", lambda: integration.tsses(db, as_bed6=True))
    records = []
    for line in sites:
        fields = str(line).rstrip("\n").split("\t")
        records.append(
            {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "name": fields[3],
                "strand": fields[5],
            }
        )
    records.sort(key=lambda row: (row["chrom"], row["start"], row["name"]))
    return {
        "status": "completed",
        "operation": "tss_gffutils_intervals",
        "records": records,
        "record_count": len(records),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "BED6 one-base sites; start is 0-based, one site per transcript",
        },
    }


def parse_gffutils_feature_line(line: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse one nine-column GFF3 line into a feature record with gffutils.feature_from_line."""

    text = _text(line, "line", maximum=65_536)
    if text.count("\t") != 8:
        raise ConfigurationError("line must be a tab-separated GFF record with nine columns.")
    feature_module = require("gffutils.feature", "bio-gffutils")
    try:
        feature = feature_module.feature_from_line(text, strict=True)
    except Exception as exc:  # gffutils raises ValueError for malformed attributes or columns
        raise ConfigurationError(f"gffutils could not parse the line: {exc}") from exc
    attributes = {
        str(key): sorted(str(item) for item in values)
        for key, values in sorted(dict(feature.attributes).items())
    }
    return {
        "status": "completed",
        "operation": "parse_gffutils_feature_line",
        "seqid": feature.seqid,
        "source": feature.source,
        "featuretype": feature.featuretype,
        "start": int(feature.start),
        "end": int(feature.end),
        "score": _plain(feature.score),
        "strand": feature.strand,
        "frame": _plain(feature.frame),
        "attributes": attributes,
        "provenance": {
            "package": _package("gffutils"),
            "convention": "start and end are 1-based inclusive as stored in GFF",
        },
    }


def inspect_gffutils_source(path: str | Path, *, verbose: bool = True) -> dict[str, Any]:
    """Count feature types, chromosomes, attribute keys, and features in a GFF or GTF file."""

    source = input_path(path)
    inspector = require("gffutils.inspect", "bio-gffutils")
    counts = _bedtools_call("inspect", lambda: inspector.inspect(str(source), verbose=False))

    def ordered(mapping: Any) -> dict[str, int]:
        return {str(key): int(value) for key, value in sorted(dict(mapping).items())}

    return {
        "status": "completed",
        "operation": "inspect_gffutils_source",
        "featuretype": ordered(counts["featuretype"]),
        "chrom": ordered(counts["chrom"]),
        "attribute_keys": ordered(counts["attribute_keys"]),
        "feature_count": int(counts["feature_count"]),
        "provenance": {
            "package": _package("gffutils"),
            "convention": "counts are per feature line; keys are sorted",
        },
    }


__all__ = [
    "inspect_gffutils_source",
    "parse_gffutils_feature_line",
    "to_bed12_gffutils_transcript",
    "tss_gffutils_intervals",
]

_register_public_verbose(globals(), tuple(__all__))
