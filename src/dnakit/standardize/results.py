"""Structured result objects returned by standardization and validation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal, cast

from dnakit.core._json import to_json_compatible
from dnakit.core.gap import Gap
from dnakit.core.issues import Issue
from dnakit.core.provenance import Provenance
from dnakit.core.sequence import DNASequence

if TYPE_CHECKING:
    from .config import DatasetValidationConfig, NormalizationConfig, ValidationConfig


@dataclass(frozen=True, slots=True)
class InputPosition:
    """A zero-based character position in the original, possibly multipart, input."""

    part_index: int
    offset: int
    absolute_offset: int


@dataclass(frozen=True, slots=True)
class NormalizationChange:
    """One auditable edit tied to an original input coordinate."""

    operation: str
    position: InputPosition
    before: str
    after: str
    normalized_offset: int
    reason: str


@dataclass(frozen=True, slots=True)
class NormalizationStep:
    """Summary of one deterministic step in the normalization pipeline."""

    name: str
    status: Literal["applied", "no-op", "failed"]
    change_count: int
    parameters: tuple[tuple[str, str], ...] = ()


@dataclass(frozen=True, slots=True)
class RawInputSnapshot:
    """Full in-memory input snapshot plus a persistence-safe digest."""

    input_type: Literal["str", "bytes", "parts", "DNASequence"]
    content: str | bytes | tuple[str | Gap, ...] | DNASequence
    sha256: str
    character_count: int

    def to_dict(self, *, include_content: bool = False) -> dict[str, Any]:
        """Return metadata and include sensitive raw content only when requested."""
        result: dict[str, Any] = {
            "input_type": self.input_type,
            "sha256": self.sha256,
            "character_count": self.character_count,
        }
        if include_content:
            result["content"] = _raw_content_to_json(self.content)
        return result


@dataclass(frozen=True, slots=True)
class SymbolOccurrence:
    """A symbol occurrence in normalized sequence coordinates."""

    symbol: str
    symbol_offset: int
    part_index: int
    part_offset: int


@dataclass(frozen=True, slots=True)
class SymbolCount:
    """Count and positions for one symbol."""

    symbol: str
    count: int
    positions: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ProbabilityResolution:
    """Per-position base probabilities for an IUPAC ambiguity symbol."""

    symbol: str
    symbol_offset: int
    probabilities: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class AmbiguityReport:
    """Counts, fraction, positions and optional probability resolutions."""

    total_count: int
    denominator: int | None
    fraction: float | None
    by_symbol: tuple[SymbolCount, ...]
    occurrences: tuple[SymbolOccurrence, ...]
    probability_resolutions: tuple[ProbabilityResolution, ...] = ()

    def count(self, symbol: str) -> int:
        """Return the count for ``symbol`` without exposing mutable mappings."""
        upper = symbol.upper()
        return next((item.count for item in self.by_symbol if item.symbol == upper), 0)


@dataclass(frozen=True, slots=True)
class InvalidSymbol:
    """An invalid Unicode code point and every original occurrence."""

    symbol: str
    codepoint: str
    positions: tuple[InputPosition, ...]


@dataclass(frozen=True, slots=True)
class NormalizationResult:
    """Complete, non-mutating normalization output and audit trail."""

    raw_input: RawInputSnapshot
    sequence: DNASequence | None
    config: NormalizationConfig
    algorithm_version: str
    provenance: Provenance
    normalized_parts: tuple[str | Gap, ...]
    steps: tuple[NormalizationStep, ...]
    changes: tuple[NormalizationChange, ...]
    issues: tuple[Issue, ...]
    ambiguity: AmbiguityReport
    invalid_symbols: tuple[InvalidSymbol, ...]
    u_positions: tuple[InputPosition, ...]

    @property
    def is_valid(self) -> bool:
        """Whether a valid :class:`DNASequence` was produced."""
        return self.sequence is not None and not _has_error(self.issues)

    @property
    def was_modified(self) -> bool:
        """Whether at least one input character was edited or removed."""
        return bool(self.changes)

    def to_dict(self, *, include_raw_content: bool = False) -> dict[str, Any]:
        """Serialize safely; raw sequence content is opt-in."""
        result = cast(dict[str, Any], to_json_compatible(self))
        result["raw_input"] = self.raw_input.to_dict(include_content=include_raw_content)
        return result


@dataclass(frozen=True, slots=True)
class QualitySummary:
    """Basic numeric summary for a PHRED quality annotation."""

    count: int
    minimum: float | None
    maximum: float | None
    mean: float | None


@dataclass(frozen=True, slots=True)
class ValidationReport:
    """Validation and QC findings for one sequence or record."""

    record_id: str | None
    sequence: DNASequence
    config: ValidationConfig
    algorithm_version: str
    provenance: Provenance
    is_valid: bool
    symbol_length: int
    coordinate_span: int | None
    ambiguity: AmbiguityReport
    quality: QualitySummary | None
    issues: tuple[Issue, ...]

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible representation."""
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class DuplicateID:
    """A repeated record identifier and its stable input indices."""

    id: str
    indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class DatasetValidationReport:
    """Collection-level validation outcome."""

    record_count: int
    config: DatasetValidationConfig
    algorithm_version: str
    provenance: Provenance
    ids_unique: bool
    duplicate_ids: tuple[DuplicateID, ...]
    is_valid: bool
    issues: tuple[Issue, ...]
    record_reports: tuple[ValidationReport, ...] | None

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible representation."""
        return cast(dict[str, Any], to_json_compatible(self))


def _has_error(issues: tuple[Issue, ...]) -> bool:
    return any(getattr(issue.severity, "value", issue.severity) == "error" for issue in issues)


def _raw_content_to_json(content: object) -> object:
    if isinstance(content, bytes):
        return {"encoding": "hex", "data": content.hex()}
    if isinstance(content, DNASequence):
        return _parts_to_json(content.parts)
    if isinstance(content, tuple):
        return _parts_to_json(content)
    return content


def _parts_to_json(parts: tuple[str | Gap, ...]) -> list[object]:
    result: list[object] = []
    for part in parts:
        if isinstance(part, str):
            result.append(part)
        else:
            result.append(
                {
                    "gap": {
                        "length": part.length,
                        "kind": getattr(part.kind, "value", part.kind),
                        "crossable": part.crossable,
                        "evidence": list(part.evidence),
                        "metadata": to_json_compatible(part.metadata),
                    }
                }
            )
    return result
