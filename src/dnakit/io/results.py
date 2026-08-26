"""Structured results produced by DNAKit writers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any, cast

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.core.enums import ExecutionMode, ImplementationLabel, OriginClass
from dnakit.core.issues import Issue
from dnakit.core.provenance import ArtifactRef, ImplementationInfo, Provenance
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True, slots=True)
class GeneratedID:
    """Stable ID assigned to an anonymous sequence at one input position."""

    input_index: int
    generated_id: str

    def __post_init__(self) -> None:
        if (
            isinstance(self.input_index, bool)
            or not isinstance(self.input_index, int)
            or self.input_index < 0
        ):
            raise ConfigurationError(
                "GeneratedID input_index must be a non-negative integer.",
                code="INVALID_GENERATED_ID",
            )
        if not isinstance(self.generated_id, str) or not self.generated_id.strip():
            raise ConfigurationError(
                "GeneratedID generated_id must be a non-empty string.",
                code="INVALID_GENERATED_ID",
            )


@dataclass(frozen=True, init=False)
class WriteResult:
    """Auditable summary of a completed serialization operation."""

    format: str
    record_count: int
    byte_count: int | None
    generated_ids: tuple[GeneratedID, ...]
    target_artifact: ArtifactRef | None
    parameters: FrozenDict
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __init__(
        self,
        format: str,
        record_count: int,
        *,
        byte_count: int | None = None,
        generated_ids: Iterable[GeneratedID] = (),
        target_artifact: ArtifactRef | None = None,
        parameters: Mapping[str, object] | None = None,
        provenance: Provenance | None = None,
        issues: Iterable[Issue] = (),
    ) -> None:
        if not isinstance(format, str) or not format.strip():
            raise ConfigurationError(
                "WriteResult format must be a non-empty string.",
                code="INVALID_WRITE_RESULT",
            )
        if isinstance(record_count, bool) or not isinstance(record_count, int) or record_count < 0:
            raise ConfigurationError(
                "WriteResult record_count must be a non-negative integer.",
                code="INVALID_WRITE_RESULT",
            )
        if byte_count is not None and (
            isinstance(byte_count, bool) or not isinstance(byte_count, int) or byte_count < 0
        ):
            raise ConfigurationError(
                "WriteResult byte_count must be non-negative or None.",
                code="INVALID_WRITE_RESULT",
            )
        generated_tuple = tuple(generated_ids)
        issue_tuple = tuple(issues)
        if any(not isinstance(item, GeneratedID) for item in generated_tuple):
            raise ConfigurationError(
                "WriteResult generated_ids must contain GeneratedID objects.",
                code="INVALID_WRITE_RESULT",
            )
        if any(not isinstance(item, Issue) for item in issue_tuple):
            raise ConfigurationError(
                "WriteResult issues must contain Issue objects.",
                code="INVALID_WRITE_RESULT",
            )
        if target_artifact is not None and not isinstance(target_artifact, ArtifactRef):
            raise ConfigurationError(
                "WriteResult target_artifact must be ArtifactRef or None.",
                code="INVALID_WRITE_RESULT",
            )
        if provenance is not None and not isinstance(provenance, Provenance):
            raise ConfigurationError(
                "WriteResult provenance must be Provenance or None.",
                code="INVALID_WRITE_RESULT",
            )
        object.__setattr__(self, "format", format)
        object.__setattr__(self, "record_count", record_count)
        object.__setattr__(self, "byte_count", byte_count)
        object.__setattr__(self, "generated_ids", generated_tuple)
        object.__setattr__(self, "target_artifact", target_artifact)
        object.__setattr__(self, "parameters", freeze_mapping(parameters))
        object.__setattr__(
            self,
            "provenance",
            provenance
            if provenance is not None
            else Provenance(
                implementation=ImplementationInfo(
                    ImplementationLabel.REIMPLEMENTATION,
                    ExecutionMode.INTERNAL,
                    OriginClass.STANDARD,
                )
            ),
        )
        object.__setattr__(self, "issues", issue_tuple)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible representation."""

        return cast(dict[str, Any], to_json_compatible(self))


__all__ = ["GeneratedID", "WriteResult"]
