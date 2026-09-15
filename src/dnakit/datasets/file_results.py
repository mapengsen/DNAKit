"""Lazy file-backed subsets for bounded-memory random and K-fold file splitting."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.core import DNARecord

from .results import SplitResult


@dataclass(frozen=True)
class FileRecords:
    """Re-open a generated subset when iterated; no sequence payloads are retained."""

    path: Path
    ids: tuple[str, ...]
    read_config: Any = None

    def __len__(self) -> int:
        return len(self.ids)

    def __iter__(self) -> Iterator[DNARecord]:
        from dnakit.io import read_inputs

        with read_inputs(self.path, read_config=self.read_config) as records:
            yield from records


@dataclass(frozen=True)
class FileSplitSubset:
    name: str
    records: FileRecords


@dataclass(frozen=True)
class FileSplitResult:
    """Split audit plus lazy subsets; ``to_dict`` matches the in-memory split audit."""

    subsets: tuple[FileSplitSubset, ...]
    output_path: Path
    _plan: SplitResult

    def get(self, name: str) -> FileRecords:
        for subset in self.subsets:
            if subset.name == name:
                return subset.records
        raise KeyError(name)

    def __getattr__(self, name: str) -> Any:
        return getattr(object.__getattribute__(self, "_plan"), name)

    def to_dict(self) -> dict[str, Any]:
        return self._plan.to_dict()
