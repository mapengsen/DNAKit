"""Sequence records with the detected alphabet through gtars (family-variant check, group H4).

``gtars.refget.digest_sequence`` and ``gtars.refget.load_fasta`` return sequence records that
carry the alphabet gtars detects (``dna2bit``, ``dna3bit``, ``dnaio``, ``protein`` or ``ascii``)
and, for ``load_fasta``, the loaded residues. The existing ``sequence_digest`` and
``digest_fasta_sequences`` wrappers return the digests and lengths but not the alphabet.
gtars is imported only when a function runs; a missing package raises
``BackendUnavailableError`` with the install hint.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.annotation.gtars_functions import (
    _DEFAULT_ROWS,
    _checked,
    _integer,
    _provenance,
    _rows_limit,
    _text,
)
from dnakit.bioinformatics._common import input_path, require, sequence_text
from dnakit.exceptions import ConfigurationError

_MAX_BASES = 1_000_000_000


def _record(metadata: Any) -> dict[str, Any]:
    description = metadata.description
    return {
        "name": str(metadata.name),
        "description": None if description in (None, "") else str(description),
        "length": int(metadata.length),
        "alphabet": str(metadata.alphabet),
        "sha512t24u": str(metadata.sha512t24u),
        "md5": str(metadata.md5),
    }


def sequence_digest_record(
    sequence: str,
    *,
    name: str | None = None,
    description: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the GA4GH sha512t24u and MD5 digests, length and detected alphabet of a sequence.

    The sequence is uppercased first, as gtars does. ``name`` defaults to an empty string.
    """

    text = sequence_text(sequence)
    label = "" if name is None else _text(name, "name", maximum=256)
    note = None if description is None else _text(description, "description", maximum=1024)
    refget = require("gtars.refget", "bio-gtars")
    try:
        record = refget.digest_sequence(text.encode("ascii"), name=label, description=note)
    except Exception as exc:  # gtars reports invalid input with runtime errors
        raise ConfigurationError(f"gtars could not digest the sequence: {exc}") from exc
    return {
        "status": "completed",
        "operation": "sequence_digest_record",
        **_record(record.metadata),
        "provenance": _provenance(
            "sequence is uppercased before digesting; alphabet is detected by gtars"
        ),
    }


def load_fasta_sequence_collection(
    fasta: str | Path,
    *,
    include_sequences: bool = True,
    max_rows: int = _DEFAULT_ROWS,
    max_bases: int = 100_000_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a FASTA file with gtars and return each sequence with its digests and alphabet.

    With ``include_sequences=True`` each record also holds the loaded residues (uppercase,
    line breaks removed). The collection digest equals the one from ``digest_fasta_sequences``.
    """

    source = input_path(fasta)
    limit = _rows_limit(max_rows)
    base_limit = _integer(max_bases, "max_bases", 1, _MAX_BASES)
    if not isinstance(include_sequences, bool):
        raise ConfigurationError("include_sequences must be a boolean.")
    refget = require("gtars.refget", "bio-gtars")
    try:
        collection = refget.load_fasta(str(source))
    except Exception as exc:  # gtars reports unreadable FASTA with runtime errors
        raise ConfigurationError(f"gtars could not load {source.name}: {exc}") from exc
    rows = []
    total = 0
    for item in collection.sequences:
        row = _record(item.metadata)
        total += row["length"]
        if total > base_limit:
            raise ConfigurationError(f"FASTA has more than max_bases={base_limit} bases.")
        if include_sequences:
            row["sequence"] = bytes(item.sequence).decode("ascii")
        rows.append(row)
    _checked(rows, limit, "FASTA file")
    rows.sort(key=lambda row: row["name"])
    return {
        "status": "completed",
        "operation": "load_fasta_sequence_collection",
        "sequence_count": len(rows),
        "total_bases": total,
        "collection_digest": str(collection.digest),
        "sequences": rows,
        "provenance": _provenance(
            "digests follow the GA4GH refget sequence specification; sequences sorted by name"
        ),
    }


__all__ = ["load_fasta_sequence_collection", "sequence_digest_record"]

_register_public_verbose(globals(), tuple(__all__))
