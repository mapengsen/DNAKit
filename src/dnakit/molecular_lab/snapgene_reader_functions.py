"""SnapGene ``.dna`` record summaries through snapgene_reader (optional, read only).

snapgene_reader is imported only inside the function. The file is read and never
written. ``snapgene_file_to_dict`` and ``snapgene_file_to_gbk`` are already exposed
by :mod:`dnakit.molecular_lab.backends` as ``read_snapgene_sequence``,
``read_snapgene_features`` and ``convert_snapgene_to_genbank``.
"""

from __future__ import annotations

import importlib.metadata
import struct
from pathlib import Path
from typing import Any
from xml.parsers.expat import ExpatError

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_MAX_FILE_BYTES = 50_000_000
_MAX_BASES = 2_000_000
_MAX_FEATURES = 100_000


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        version = importlib.metadata.version("snapgene-reader")
    except importlib.metadata.PackageNotFoundError:
        version = None
    return {"package": "snapgene_reader", "version": version, **notes}


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _label(value: object) -> str | None:
    if isinstance(value, (list, tuple)):
        value = value[0] if value else None
    return None if value is None else str(value)


def _record_summary(
    record: Any, source: Path, *, max_bases: int, max_features: int
) -> dict[str, Any]:
    sequence = str(record.seq)
    if len(sequence) > max_bases:
        raise ConfigurationError(f"SnapGene sequence exceeds max_bases={max_bases}.")
    if len(record.features) > max_features:
        raise ConfigurationError(f"SnapGene record exceeds max_features={max_features}.")
    rows = []
    for feature in record.features:
        location = feature.location
        rows.append(
            {
                "type": str(feature.type),
                "start": int(location.start),
                "end": int(location.end),
                "strand": None if location.strand is None else int(location.strand),
                "label": _label(feature.qualifiers.get("label")),
            }
        )
    rows.sort(key=lambda row: (row["start"], row["end"], row["type"], row["label"] or ""))
    return {
        "status": "completed",
        "operation": "snapgene_seqrecord_summary",
        "source": str(source),
        "record_id": str(record.id),
        "name": str(record.name),
        "description": str(record.description),
        "length": len(sequence),
        "sequence": sequence,
        "feature_count": len(rows),
        "features": rows,
        "provenance": _provenance(format="SnapGene .dna, read only", view="Biopython SeqRecord"),
    }


def snapgene_seqrecord_summary(
    source: str | Path,
    *,
    max_bases: int = _MAX_BASES,
    max_features: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Summarize a SnapGene ``.dna`` file through snapgene_reader's SeqRecord view.

    The result holds the record identifiers, the DNA sequence and the features
    sorted by position. A missing file, or a file that snapgene_reader cannot
    parse, raises ``ConfigurationError``.
    """

    snapgene = require("snapgene_reader", "bio-snapgene")
    if not isinstance(source, (str, Path)):
        raise ConfigurationError("source must be a local SnapGene .dna file path.")
    path = input_path(source)
    if path.stat().st_size > _MAX_FILE_BYTES:
        raise ConfigurationError(f"Input file exceeds {_MAX_FILE_BYTES} bytes: {path.name}")
    limit = _integer(max_bases, "max_bases", minimum=1, maximum=_MAX_BASES)
    feature_limit = _integer(max_features, "max_features", minimum=0, maximum=_MAX_FEATURES)
    try:
        record = snapgene.snapgene_file_to_seqrecord(filepath=str(path))
    except (ExpatError, IndexError, KeyError, OSError, TypeError, ValueError, struct.error) as exc:
        raise ConfigurationError(f"snapgene_reader could not parse {path.name}: {exc}") from exc
    return _record_summary(record, path, max_bases=limit, max_features=feature_limit)


__all__ = ["snapgene_seqrecord_summary"]

_register_public_verbose(globals(), tuple(__all__))
