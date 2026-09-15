"""Private staging writes: the enclosing operation commits and hashes final outputs."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from dnakit.core import DNARecord
from dnakit.io import WriteConfig, WriteResult
from dnakit.io.api import _infer_format, _writer


def write_fragment(records: Iterable[DNARecord], path: Path, config: WriteConfig) -> WriteResult:
    fmt = _infer_format(str(path), None)
    with path.open("w", encoding=config.encoding, newline="", buffering=1024**2) as handle:
        count = _writer(handle, records, fmt, config)
    return WriteResult(fmt, count)
