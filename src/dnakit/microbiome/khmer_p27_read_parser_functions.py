"""Bounded, JSON-compatible FASTA/FASTQ reads and pairs using khmer.ReadParser."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

from . import _khmer_p27_common as c


def _record(read: Any) -> dict[str, str]:
    item = {"name": read.name, "sequence": read.sequence}
    quality = getattr(read, "quality", None)
    if quality is not None:
        item["quality"] = quality
    return item


def _parse(path: str | Path, paired: bool, max_records: int, max_bases: int) -> dict[str, Any]:
    source = input_path(path)
    limit = c.integer(max_records, "max_records", 1, 100_000)
    base_limit = c.integer(max_bases, "max_bases", 1, 10_000_000)
    khmer = require("khmer", c.KHMER_EXTRA)
    records: list[Any] = []
    bases = 0
    truncated = False
    try:
        parser = khmer.ReadParser(str(source))
        iterator = parser.iter_read_pairs() if paired else iter(parser)
        for entry in iterator:
            reads = list(entry) if paired else [entry]
            size = sum(len(read.sequence) for read in reads)
            if len(records) >= limit or bases + size > base_limit:
                truncated = True
                break
            converted = [_record(read) for read in reads]
            records.append(converted if paired else converted[0])
            bases += size
    except (OSError, ValueError, RuntimeError) as exc:
        raise ConfigurationError(f"khmer could not parse {'paired ' if paired else ''}reads: {exc}") from exc
    return {"status": "completed", "source": str(source), "paired": paired,
            "records": records, "record_count": len(records), "base_count": bases,
            "truncated": truncated, "provenance": c.provenance(source_function=
                "khmer.ReadParser.iter_read_pairs" if paired else "khmer.ReadParser")}


def khmer_parse_reads(path: str | Path, *, max_records: int = 1000, max_bases: int = 1_000_000,
                       verbose: bool = True) -> dict[str, Any]:
    """Read FASTA/FASTQ records, preserving names, sequences and FASTQ qualities.

    ``truncated`` explicitly signals that more records exist beyond the record or
    base limit. Returned sequences preserve the upstream parser's letter case.
    """
    return _parse(path, False, max_records, max_bases)


def khmer_parse_read_pairs(path: str | Path, *, max_records: int = 1000,
                            max_bases: int = 1_000_000, verbose: bool = True) -> dict[str, Any]:
    """Read interleaved FASTA/FASTQ pairs with khmer's name-pairing validation.

    Each returned record is a two-read list; ``max_records`` counts pairs. Invalid
    pairs raise ConfigurationError rather than silently discarding orphan reads.
    """
    return _parse(path, True, max_records, max_bases)


__all__ = ["khmer_parse_reads", "khmer_parse_read_pairs"]
_register_public_verbose(globals(), tuple(__all__))
