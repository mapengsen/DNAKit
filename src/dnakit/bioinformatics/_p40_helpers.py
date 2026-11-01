"""Private helpers shared by the pack P40 modules (explicit local files, bounded results)."""

from __future__ import annotations

import gzip
import math
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError

MAX_ROWS = 100_000
DEFAULT_ROWS = 1_000


def plain(value: Any) -> Any:
    """Numpy/pandas scalars to plain Python; NaN and infinity become None."""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item") and not isinstance(value, (str, dict, list, tuple)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def existing(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ConfigurationError(f"{name} must be an existing file: {path}")
    return path


def new_file(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    if not path.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {path.parent}")
    return path


def new_dir(value: str | Path, name: str) -> Path:
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} must not exist yet: {path}")
    path.mkdir(parents=True)
    return path


def integer(value: object, name: str, *, minimum: int | None = None, maximum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(f"{name} must be an integer.")
    if minimum is not None and value < minimum:
        raise ConfigurationError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return value


def choice(value: object, options: Sequence[str], name: str) -> str:
    if value not in options:
        raise ConfigurationError(f"{name} must be one of {', '.join(options)}.")
    return str(value)


def text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")
    if "\x00" in value:
        raise ConfigurationError(f"{name} must not contain NUL characters.")
    return value


def open_text(path: Path) -> Any:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def read_fasta(path: Path, *, max_bases: int = 500_000_000) -> list[tuple[str, str, str]]:
    """Plain FASTA reader returning (id, description, sequence) with the sequence as in the file."""
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    total = 0

    def flush() -> None:
        if header is not None:
            parts = header.split(None, 1)
            records.append((parts[0] if parts else "", parts[1] if len(parts) > 1 else "", "".join(chunks)))

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                flush()
                header, chunks = line[1:], []
            elif line.strip():
                if header is None:
                    raise ConfigurationError(f"{path} is not FASTA: sequence before the first '>' line.")
                chunks.append(line.strip())
                total += len(line)
                if total > max_bases:
                    raise ConfigurationError(f"{path} is larger than {max_bases} bases.")
    flush()
    if not records:
        raise ConfigurationError(f"{path} contains no FASTA records.")
    return records


def iter_fastx(path: Path) -> Iterator[tuple[str, str, str | None]]:
    """(name line without marker, sequence, qualities or None) for FASTA or FASTQ files."""
    with open_text(path) as handle:
        first = handle.readline()
        if not first:
            return
        if first.startswith("@"):
            line = first
            while line:
                name = line.rstrip("\r\n")[1:]
                seq = handle.readline().rstrip("\r\n")
                handle.readline()
                qual = handle.readline().rstrip("\r\n")
                yield name, seq, qual
                line = handle.readline()
            return
        if not first.startswith(">"):
            raise ConfigurationError(f"{path} is neither FASTA nor FASTQ.")
    for rid, desc, seq in read_fasta(path):
        yield (rid + (" " + desc if desc else "")), seq, None
