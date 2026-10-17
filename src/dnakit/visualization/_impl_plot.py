"""Sequence parsing and statistics shared by the ``dnakit.visualization`` candidate functions."""

from __future__ import annotations

import os
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit.exceptions import InputFormatError

MAX_SEQUENCES = 100_000
MAX_TOTAL_BASES = 500_000_000


def _is_file(path: Path) -> bool:
    try:
        return path.is_file()
    except (OSError, ValueError):
        return False


def read_sequences(value: object) -> list[tuple[str, str]]:
    """``[(id, sequence)]`` (upper-case, whitespace removed, U read as T) from FASTA text/path, a plain sequence,
    a list of sequences / ``(id, sequence)`` pairs / mappings with ``id`` and ``sequence``, or ``{id: sequence}``."""

    if isinstance(value, Mapping) and "records" in value:
        value = value["records"]
    entries: list[tuple[str, object]] = []
    if isinstance(value, Mapping):
        if "sequence" in value:
            entries = [(str(value.get("id", "seq1")), value["sequence"])]
        else:
            entries = [(str(key), item) for key, item in value.items()]
    elif isinstance(value, (str, os.PathLike)):
        text: str
        if isinstance(value, os.PathLike) or ("\n" not in value and _is_file(Path(value))):
            path = Path(value)
            if not _is_file(path):
                raise InputFormatError(f"File not found: {path}", code="FILE_NOT_FOUND")
            text = path.read_text(encoding="utf-8", errors="replace")
        else:
            text = value
        if text.lstrip().startswith(">"):
            identifier: str | None = None
            chunks: list[str] = []
            for line in text.splitlines():
                line = line.strip()
                if line.startswith(">"):
                    if identifier is not None:
                        entries.append((identifier, "".join(chunks)))
                    header = line[1:].split(None, 1)
                    identifier, chunks = (header[0] if header else f"seq{len(entries) + 1}"), []
                elif line and identifier is not None:
                    chunks.append(line)
            if identifier is not None:
                entries.append((identifier, "".join(chunks)))
        else:
            entries = [("seq1", text)]
    elif isinstance(value, Sequence):
        for index, item in enumerate(value, start=1):
            if isinstance(item, Mapping):
                entries.append((str(item.get("id", f"seq{index}")), item.get("sequence", "")))
            elif isinstance(item, Sequence) and not isinstance(item, str) and len(item) == 2:
                entries.append((str(item[0]), item[1]))
            else:
                entries.append((f"seq{index}", item))
    else:
        raise InputFormatError("Provide FASTA text, a sequence, or a list/mapping of sequences.", code="INVALID_SEQUENCE_INPUT")
    if not entries or len(entries) > MAX_SEQUENCES:
        raise InputFormatError(f"Provide between 1 and {MAX_SEQUENCES} sequences.", code="INVALID_SEQUENCE_INPUT")
    result: list[tuple[str, str]] = []
    total = 0
    for identifier, sequence in entries:
        if not isinstance(sequence, str):
            raise InputFormatError(f"Sequence {identifier!r} must be text.", code="INVALID_SEQUENCE")
        cleaned = "".join(sequence.split()).upper().replace("U", "T")
        invalid = sorted({c for c in cleaned if not ("A" <= c <= "Z" or c in "-*.")})
        if invalid:
            raise InputFormatError(f"Sequence {identifier!r} has invalid symbols {invalid[:5]}.", code="INVALID_SEQUENCE")
        total += len(cleaned)
        if total > MAX_TOTAL_BASES:
            raise InputFormatError("The sequences are too large.", code="INPUT_TOO_LARGE")
        result.append((identifier, cleaned))
    return result


def gc_percent(sequence: str) -> float | None:
    """100 * (G + C) / (A + C + G + T); ``None`` when the sequence has no unambiguous base."""

    counts = Counter(sequence)
    acgt = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    return 100.0 * (counts["G"] + counts["C"]) / acgt if acgt else None


def n_statistic(lengths: Sequence[int], fraction: float = 0.5) -> int:
    """Smallest length ``L`` such that the sequences of length >= ``L`` hold at least ``fraction`` of all bases."""

    total, running = sum(lengths), 0
    for value in sorted(lengths, reverse=True):
        running += value
        if running >= total * fraction:
            return value
    return 0


def median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    return float(ordered[middle]) if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2


def summarize(records: Sequence[tuple[str, str]]) -> dict[str, object]:
    lengths = [len(sequence) for _, sequence in records]
    counts: Counter[str] = Counter()
    for _, sequence in records:
        counts.update(sequence)
    acgt = sum(counts[base] for base in "ACGT")
    return {
        "sequence_count": len(records), "total_bases": sum(lengths),
        "length": {"min": min(lengths), "max": max(lengths), "mean": sum(lengths) / len(lengths), "median": median(lengths), "n50": n_statistic(lengths)},
        "gc_percent": 100.0 * (counts["G"] + counts["C"]) / acgt if acgt else None,
        "base_counts": dict(sorted(counts.items())),
        "non_acgt_bases": sum(counts.values()) - acgt,
        "n_fraction": counts["N"] / sum(lengths) if sum(lengths) else 0.0,
    }
