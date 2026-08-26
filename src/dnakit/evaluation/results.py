"""Immutable, JSON-compatible results for DNA evaluation workflows."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core import DNASet, Issue, Provenance
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError

_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _text(value: object, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")


@dataclass(frozen=True, slots=True)
class EvaluationEntry:
    """One subject-level evaluation row in stable input order."""

    subject_id: str
    input_index: int
    metrics: FrozenDict
    issues: tuple[Issue, ...] = ()

    def __post_init__(self) -> None:
        _text(self.subject_id, "subject_id")
        if (
            isinstance(self.input_index, bool)
            or not isinstance(self.input_index, int)
            or self.input_index < 0
        ):
            raise ConfigurationError("input_index must be a non-negative integer.")
        if not isinstance(self.metrics, FrozenDict):
            raise ConfigurationError("EvaluationEntry metrics must be FrozenDict.")
        if any(not isinstance(issue, Issue) for issue in self.issues):
            raise ConfigurationError("EvaluationEntry issues must contain Issue objects.")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class EvaluationReport:
    """A deterministic evaluation report with full method and provenance context."""

    name: str
    method: str
    algorithm_version: str
    parameters: FrozenDict
    metrics: FrozenDict
    entries: tuple[EvaluationEntry, ...]
    provenance: Provenance
    issues: tuple[Issue, ...] = ()

    def __post_init__(self) -> None:
        for name in ("name", "method", "algorithm_version"):
            _text(getattr(self, name), name)
        if not isinstance(self.parameters, FrozenDict) or not isinstance(self.metrics, FrozenDict):
            raise ConfigurationError("Evaluation report mappings must be FrozenDict objects.")
        if any(not isinstance(entry, EvaluationEntry) for entry in self.entries):
            raise ConfigurationError("Evaluation report entries must be EvaluationEntry objects.")
        if any(entry.input_index != index for index, entry in enumerate(self.entries)):
            raise ConfigurationError("Evaluation entry indices must be contiguous input order.")
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError("Evaluation report provenance must be Provenance.")
        if any(not isinstance(issue, Issue) for issue in self.issues):
            raise ConfigurationError("Evaluation report issues must contain Issue objects.")

    @property
    def score(self) -> float | None:
        value = self.metrics.get("score")
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return None
        return float(value)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class ReferenceLibrary:
    """A local, versioned reference collection with a deterministic content digest."""

    records: DNASet
    name: str
    version: str
    source: str
    digest: str
    digest_scope: str
    date: str | None
    filters: FrozenDict
    index_parameters: FrozenDict

    def __post_init__(self) -> None:
        if not isinstance(self.records, DNASet):
            raise ConfigurationError("ReferenceLibrary records must be DNASet.")
        for name in ("name", "version", "source", "digest_scope"):
            _text(getattr(self, name), name)
        if self.date is not None:
            _text(self.date, "date")
        if not isinstance(self.digest, str) or not _SHA256.fullmatch(self.digest):
            raise ConfigurationError("ReferenceLibrary digest must be a lowercase SHA-256.")
        if len(set(self.records.ids)) != len(self.records):
            raise ConfigurationError("ReferenceLibrary record IDs must be unique.")
        if not isinstance(self.filters, FrozenDict) or not isinstance(
            self.index_parameters, FrozenDict
        ):
            raise ConfigurationError("ReferenceLibrary metadata mappings must be FrozenDict.")

    def to_dict(self) -> dict[str, Any]:
        """Serialize provenance without duplicating reference sequence content."""

        return cast(
            dict[str, Any],
            to_json_compatible(
                {
                    "name": self.name,
                    "version": self.version,
                    "source": self.source,
                    "digest": self.digest,
                    "digest_scope": self.digest_scope,
                    "date": self.date,
                    "filters": self.filters,
                    "index_parameters": self.index_parameters,
                    "record_count": len(self.records),
                    "record_ids": self.records.ids,
                }
            ),
        )


__all__ = ["EvaluationEntry", "EvaluationReport", "ReferenceLibrary"]
