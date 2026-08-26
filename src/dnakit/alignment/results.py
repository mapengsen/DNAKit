"""Immutable pairwise-alignment result objects."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Literal, cast

from dnakit.core import (
    ExecutionMode,
    ImplementationInfo,
    ImplementationLabel,
    Issue,
    OriginClass,
    Provenance,
)
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError


def _alignment_provenance() -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.REIMPLEMENTATION,
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=OriginClass.PUBLISHED_ALGORITHM,
        )
    )


@dataclass(frozen=True, slots=True)
class AlignmentColumn:
    """One alignment column and its optional source coordinates."""

    query_symbol: str
    target_symbol: str
    query_position: int | None
    target_position: int | None
    operation: Literal["match", "mismatch", "insertion", "deletion"]

    def __post_init__(self) -> None:
        for name in ("query_symbol", "target_symbol"):
            value = getattr(self, name)
            if not isinstance(value, str) or len(value) != 1:
                raise ConfigurationError(f"Alignment column {name} must be one character.")
        for name in ("query_position", "target_position"):
            value = getattr(self, name)
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, int) or value < 0
            ):
                raise ConfigurationError(f"Alignment column {name} must be non-negative or None.")
        if self.operation not in ("match", "mismatch", "insertion", "deletion"):
            raise ConfigurationError("Alignment column operation is invalid.")
        if self.operation == "insertion":
            valid = (
                self.query_symbol == "-"
                and self.query_position is None
                and self.target_symbol != "-"
                and self.target_position is not None
            )
        elif self.operation == "deletion":
            valid = (
                self.target_symbol == "-"
                and self.target_position is None
                and self.query_symbol != "-"
                and self.query_position is not None
            )
        else:
            valid = (
                self.query_symbol != "-"
                and self.target_symbol != "-"
                and self.query_position is not None
                and self.target_position is not None
                and ((self.query_symbol == self.target_symbol) == (self.operation == "match"))
            )
        if not valid:
            raise ConfigurationError("Alignment column fields conflict with its operation.")


@dataclass(frozen=True, slots=True)
class AlignmentResult:
    """Traceable global or local alignment using alignment-only '-' symbols."""

    name: str
    method: Literal["global", "local", "semi_global"]
    algorithm_version: str
    score: float
    aligned_query: str
    aligned_target: str
    query_id: str | None
    target_id: str | None
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    matches: int
    mismatches: int
    insertions: int
    deletions: int
    identity: float | None
    query_coverage: float
    target_coverage: float
    columns: tuple[AlignmentColumn, ...]
    parameters: FrozenDict
    provenance: Provenance = field(default_factory=_alignment_provenance)
    issues: tuple[Issue, ...] = ()

    def __post_init__(self) -> None:
        for name in ("name", "algorithm_version"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Alignment {name} must be non-empty.")
        if self.method not in ("global", "local", "semi_global"):
            raise ConfigurationError("Alignment method must be global, local, or semi_global.")
        for name in ("query_id", "target_id"):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError(f"Alignment {name} must be non-empty or None.")
        if (
            isinstance(self.score, bool)
            or not isinstance(self.score, (int, float))
            or not math.isfinite(self.score)
        ):
            raise ConfigurationError("Alignment score must be finite.")
        if not isinstance(self.aligned_query, str) or not isinstance(self.aligned_target, str):
            raise ConfigurationError("Aligned values must be strings.")
        if len(self.aligned_query) != len(self.aligned_target):
            raise ConfigurationError("Aligned strings must have equal lengths.")
        if not isinstance(self.columns, tuple) or any(
            not isinstance(column, AlignmentColumn) for column in self.columns
        ):
            raise ConfigurationError("Alignment columns must be AlignmentColumn objects.")
        if len(self.columns) != len(self.aligned_query) or any(
            column.query_symbol != query_symbol or column.target_symbol != target_symbol
            for column, query_symbol, target_symbol in zip(
                self.columns, self.aligned_query, self.aligned_target, strict=True
            )
        ):
            raise ConfigurationError("Alignment columns must match aligned-string length.")
        for name in ("query_start", "query_end", "target_start", "target_end"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"Alignment {name} must be non-negative.")
        if self.query_end < self.query_start or self.target_end < self.target_start:
            raise ConfigurationError("Alignment source coordinates must be ordered.")
        for name in ("matches", "mismatches", "insertions", "deletions"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"Alignment {name} must be non-negative.")
        counts = {
            "match": self.matches,
            "mismatch": self.mismatches,
            "insertion": self.insertions,
            "deletion": self.deletions,
        }
        if any(
            sum(column.operation == operation for column in self.columns) != count
            for operation, count in counts.items()
        ):
            raise ConfigurationError("Alignment operation counts do not match its columns.")
        if not 0.0 <= self.query_coverage <= 1.0 or not 0.0 <= self.target_coverage <= 1.0:
            raise ConfigurationError("Alignment coverage must be between zero and one.")
        if self.identity is not None and not 0.0 <= self.identity <= 1.0:
            raise ConfigurationError("Alignment identity must be between zero and one.")
        expected_identity = self.matches / len(self.columns) if self.columns else None
        if self.identity is None:
            if expected_identity is not None:
                raise ConfigurationError("Non-empty alignment must report identity.")
        elif expected_identity is None or not math.isclose(self.identity, expected_identity):
            raise ConfigurationError("Alignment identity does not match its columns.")
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("Alignment parameters must be FrozenDict.")
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError("Alignment provenance must be Provenance.")
        if not isinstance(self.issues, tuple) or any(
            not isinstance(issue, Issue) for issue in self.issues
        ):
            raise ConfigurationError("Alignment issues must contain Issue objects.")

    @property
    def alignment_length(self) -> int:
        return len(self.columns)

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = ["AlignmentColumn", "AlignmentResult"]
