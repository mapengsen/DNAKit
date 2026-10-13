"""khmer read-cleaning and paired-read helpers for DNA, using local sequences or a local file.

``clean_khmer_read_sequences()`` wraps khmer's own cleaning of reads, which uppercases bases and
replaces N with A. The replaced positions are counted and reported, because that replacement
changes the sequence. ``read_paired_fastx()`` pairs the records of a local FASTA or FASTQ file
with khmer's broken-paired reader and reports singletons. Sequences are uppercased before they
reach khmer, since khmer silently gives wrong results for lowercase bases.

khmer is imported only when a function runs. A missing package raises ``BackendUnavailableError``
with the ``pip install "dnakit[bio-khmer]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
from collections.abc import Iterator, Mapping
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_KHMER_EXTRA = "bio-khmer"
_ACGTN = frozenset("ACGTN")
_READ_LENGTH_LIMIT = 100_000
_CLEAN_READ_LIMIT = 10_000
_DEFAULT_MAX_READS = 100_000
_MAX_READS_LIMIT = 10_000_000


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


def _read_id(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError("read IDs must be non-empty strings.")
    return value


def _read_sequence(value: object, name: str) -> str:
    """Uppercase a read and check that it has only A, C, G, T and N."""

    if not isinstance(value, str) or not value:
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = value.upper()
    if len(text) > _READ_LENGTH_LIMIT:
        raise ConfigurationError(f"{name} length must be at most {_READ_LENGTH_LIMIT}.")
    if set(text) - _ACGTN:
        raise ConfigurationError(f"{name} must contain only A, C, G, T and N.")
    return text


def _limited_records(khmer: Any, path: Path, limit: int) -> Iterator[Any]:
    """Stream khmer's parser output and stop with a clear error after ``limit`` records."""

    for count, record in enumerate(khmer.ReadParser(str(path)), start=1):
        if count > limit:
            raise ConfigurationError(
                f"read input exceeds max_reads={limit}; increase the explicit limit."
            )
        yield record


def _read_entry(record: Any) -> dict[str, str]:
    return {"name": str(record.name), "sequence": str(record.sequence).upper()}


def clean_khmer_read_sequences(
    sequences: Mapping[str, str],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Uppercase reads and replace N with A, as khmer's ``clean_input_reads`` does.

    ``sequences`` maps a read ID to a DNA string. Calls ``khmer.utils.clean_input_reads``, which
    sets ``cleaned_seq`` to the uppercase read with N replaced by A. The replacement is lossy: an
    N becomes an A, so the cleaned read is not the original sequence. Each read reports its
    ``replaced_n_count``. Only A, C, G, T and N are accepted; other symbols raise an error. Reads
    are returned sorted by ID.
    """

    if not isinstance(sequences, Mapping) or not sequences:
        raise ConfigurationError("sequences must be a non-empty mapping of read ID to DNA string.")
    if len(sequences) > _CLEAN_READ_LIMIT:
        raise ConfigurationError(f"at most {_CLEAN_READ_LIMIT} reads can be cleaned at once.")
    texts = {
        _read_id(key): _read_sequence(value, f"sequence {key!r}")
        for key, value in sequences.items()
    }
    utils = require("khmer.utils", _KHMER_EXTRA)
    records = [SimpleNamespace(name=read_id, sequence=texts[read_id]) for read_id in sorted(texts)]
    reads = []
    for record in utils.clean_input_reads(records):
        reads.append(
            {
                "id": record.name,
                "sequence": record.sequence,
                "cleaned_sequence": record.cleaned_seq,
                "replaced_n_count": record.sequence.count("N"),
            }
        )
    replaced_total = sum(item["replaced_n_count"] for item in reads)
    return {
        "status": "completed",
        "operation": "clean_khmer_read_sequences",
        "read_count": len(reads),
        "reads_with_n": sum(1 for item in reads if item["replaced_n_count"]),
        "replaced_n_total": replaced_total,
        "reads": reads,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "function": "khmer.utils.clean_input_reads",
            "coordinate_system": "read order by ID; no genomic coordinates",
            "convention": "input is uppercased, then N is replaced by A (khmer's cleaned_seq); "
            "the replaced N count is reported, and the cleaned read is not the original sequence",
        },
    }


def read_paired_fastx(
    source: str | Path,
    *,
    min_length: int | None = None,
    require_paired: bool = False,
    force_single: bool = False,
    max_reads: int = _DEFAULT_MAX_READS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pair the records of a local FASTA or FASTQ file, and keep unpaired reads as singletons.

    Calls ``khmer.utils.broken_paired_reader`` over ``khmer.ReadParser``. Two consecutive records
    form a pair when their names are a mate pair, such as ``name/1`` and ``name/2``, or the
    ``1:``/``2:`` and ``seq/1``/``seq/2`` forms. ``min_length`` drops reads shorter than it, and
    ``require_paired`` turns any orphan into an error. ``force_single`` returns every record as a
    singleton. ``units`` lists each pair or singleton in file order, with ``index`` counted as
    khmer counts records. Sequences are uppercased. ``max_reads`` bounds the number of records.
    """

    path = input_path(source)
    minimum = (
        None if min_length is None else _integer(min_length, "min_length", 1, _READ_LENGTH_LIMIT)
    )
    paired_required = _flag(require_paired, "require_paired")
    single_only = _flag(force_single, "force_single")
    if paired_required and single_only:
        raise ConfigurationError("require_paired and force_single cannot both be true.")
    limit = _integer(max_reads, "max_reads", 1, _MAX_READS_LIMIT)
    khmer = require("khmer", _KHMER_EXTRA)
    utils = require("khmer.utils", _KHMER_EXTRA)
    units: list[dict[str, Any]] = []
    try:
        for index, is_pair, first, second in utils.broken_paired_reader(
            _limited_records(khmer, path, limit),
            min_length=minimum,
            force_single=single_only,
            require_paired=paired_required,
        ):
            units.append(
                {
                    "index": int(index),
                    "is_pair": bool(is_pair),
                    "read1": _read_entry(first),
                    "read2": None if second is None else _read_entry(second),
                }
            )
    except ConfigurationError:
        raise
    except (OSError, ValueError) as exc:
        raise ConfigurationError(
            f"khmer could not read paired records from {path.name}: {exc}"
        ) from exc
    pair_count = sum(1 for unit in units if unit["is_pair"])
    return {
        "status": "completed",
        "operation": "read_paired_fastx",
        "source": str(path),
        "min_length": minimum,
        "require_paired": paired_required,
        "force_single": single_only,
        "unit_count": len(units),
        "pair_count": pair_count,
        "singleton_count": len(units) - pair_count,
        "units": units,
        "provenance": {
            "package": f"khmer {_version('khmer')}",
            "reader": "khmer.ReadParser with khmer.utils.broken_paired_reader",
            "coordinate_system": "record order in the input file; no genomic coordinates",
            "convention": "sequences are uppercase; index counts records like enumerate and "
            "advances by 2 for a pair; a pair holds two records of the same fragment",
        },
    }


__all__ = [
    "clean_khmer_read_sequences",
    "read_paired_fastx",
]

_register_public_verbose(globals(), tuple(__all__))
