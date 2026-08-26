"""Deterministic origin operations for circular DNA sequences."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias

from dnakit.core import DNASequence, Topology
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError, UnsupportedGapOperationError
from dnakit.ops._common import require_sequence, resolved_span
from dnakit.ops.edit import _split_at

CircularOperation: TypeAlias = Literal["rotate", "canonical_origin"]


@dataclass(frozen=True, slots=True)
class CircularOperationResult:
    """A circular sequence with an auditable origin-selection decision."""

    sequence: DNASequence
    operation: CircularOperation
    requested_offset: int | None
    effective_offset: int
    sequence_span: int
    rule: str
    parameters: FrozenDict

    def __post_init__(self) -> None:
        if self.operation not in {"rotate", "canonical_origin"}:
            raise ConfigurationError("Unknown circular operation result kind.")
        if self.sequence.topology is not Topology.CIRCULAR:
            raise ConfigurationError("A circular operation result must remain circular.")
        if self.requested_offset is not None and (
            isinstance(self.requested_offset, bool) or not isinstance(self.requested_offset, int)
        ):
            raise ConfigurationError("requested_offset must be an integer or None.")
        if self.sequence_span <= 0:
            raise ConfigurationError("A circular operation requires a positive sequence span.")
        if not 0 <= self.effective_offset < self.sequence_span:
            raise ConfigurationError("effective_offset must lie within the circular span.")
        if self.sequence.coordinate_span != self.sequence_span:
            raise ConfigurationError("sequence_span must match the transformed sequence.")
        if not isinstance(self.rule, str) or not self.rule.strip():
            raise ConfigurationError("A circular operation result rule must be non-empty.")
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("Circular operation parameters must be FrozenDict.")

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-compatible audit payload without flattening explicit gaps."""

        payload = {
            "operation": self.operation,
            "requested_offset": self.requested_offset,
            "effective_offset": self.effective_offset,
            "sequence_span": self.sequence_span,
            "rule": self.rule,
            "parameters": self.parameters,
        }
        converted = to_json_compatible(payload)
        assert isinstance(converted, dict)
        return converted


def _require_circular(sequence: DNASequence, *, operation: str) -> DNASequence:
    source = require_sequence(sequence)
    if source.topology is not Topology.CIRCULAR:
        raise ConfigurationError(
            f"{operation} requires topology='circular'.",
            code="CIRCULAR_TOPOLOGY_REQUIRED",
            context={"operation": operation, "topology": source.topology.value},
        )
    return source


def rotate(sequence: DNASequence, offset: int) -> CircularOperationResult:
    """Left-rotate a circular sequence so ``offset`` becomes coordinate zero.

    Offsets are normalized modulo the resolved coordinate span. A rotation may
    split a nucleotide fragment, but never an explicit Gap. Unknown-length Gaps
    are rejected because the effective coordinate cannot be determined.
    """

    source = _require_circular(sequence, operation="rotate")
    if isinstance(offset, bool) or not isinstance(offset, int):
        raise ConfigurationError(
            "Circular rotation offset must be an integer.",
            code="INVALID_ROTATION_OFFSET",
            context={"offset": offset},
        )
    span = resolved_span(source, operation="rotate")
    effective = offset % span
    try:
        left, right = _split_at(source, effective)
    except UnsupportedGapOperationError as exc:
        raise UnsupportedGapOperationError(
            "A circular origin cannot be placed inside a Gap.",
            code="ROTATION_ORIGIN_INSIDE_GAP",
            context={"requested_offset": offset, "effective_offset": effective},
            hint="Choose a nucleotide or Gap boundary as the circular origin.",
        ) from exc
    rotated = DNASequence(
        (*right, *left),
        alphabet=source.alphabet,
        topology=Topology.CIRCULAR,
        strandedness=source.strandedness,
    )
    return CircularOperationResult(
        rotated,
        "rotate",
        offset,
        effective,
        span,
        "left_rotation_offset_modulo_coordinate_span",
        FrozenDict(
            {
                "direction": "left",
                "offset_normalization": "modulo_coordinate_span",
                "gap_interior_policy": "error",
            }
        ),
    )


def _minimal_rotation_offset(text: str) -> int:
    """Return the smallest offset having the lexicographically minimal rotation.

    This is Booth's linear-time algorithm with deterministic smallest-index
    tie-breaking for periodic strings.
    """

    length = len(text)
    if length < 2:
        return 0
    doubled = text + text
    first, second, matched = 0, 1, 0
    while first < length and second < length and matched < length:
        left = doubled[first + matched]
        right = doubled[second + matched]
        if left == right:
            matched += 1
            continue
        if left > right:
            first = first + matched + 1
            if first <= second:
                first = second + 1
        else:
            second = second + matched + 1
            if second <= first:
                second = first + 1
        matched = 0
    return min(first, second)


def canonical_origin(sequence: DNASequence) -> CircularOperationResult:
    """Choose the lexicographically minimal forward rotation as coordinate zero.

    The rule compares normalized DNA symbols using their ordinary code-point
    order and preserves strand orientation. Equal rotations select the smallest
    offset. Explicit gaps are rejected because omitting or serializing them into
    a comparison token would change the sequence model.
    """

    source = _require_circular(sequence, operation="canonical_origin")
    if source.is_gapped:
        raise UnsupportedGapOperationError(
            "Canonical lexicographic origin selection does not accept explicit Gaps.",
            code="CANONICAL_ORIGIN_GAPPED_SEQUENCE",
            hint="Resolve the origin explicitly with rotate() at a Gap boundary.",
        )
    span = resolved_span(source, operation="canonical_origin")
    offset = _minimal_rotation_offset(source.symbols)
    rotated = rotate(source, offset).sequence
    return CircularOperationResult(
        rotated,
        "canonical_origin",
        None,
        offset,
        span,
        "lexicographically_minimal_forward_rotation_then_smallest_offset",
        FrozenDict(
            {
                "comparison": "normalized_symbol_code_point_order",
                "orientation": "forward_only",
                "tie_break": "smallest_offset",
                "algorithm": "booth",
                "algorithm_complexity": "O(n)",
            }
        ),
    )


__all__ = ["CircularOperation", "CircularOperationResult", "canonical_origin", "rotate"]
