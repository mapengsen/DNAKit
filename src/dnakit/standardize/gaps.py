"""Auditable conversion of ambiguous runs and strict AGP object assembly."""

from __future__ import annotations

import hashlib
import itertools
import json
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any, Literal, TypeAlias, cast

from dnakit.core._json import FrozenDict, freeze_mapping, to_json_compatible
from dnakit.core.coordinates import Interval
from dnakit.core.enums import (
    DNAAlphabet,
    ExecutionMode,
    GapKind,
    ImplementationLabel,
    IssueSeverity,
    OriginClass,
    Strandedness,
    Topology,
)
from dnakit.core.facade import DNA
from dnakit.core.gap import Gap
from dnakit.core.issues import Issue
from dnakit.core.provenance import ImplementationInfo, Provenance
from dnakit.core.record import DNARecord
from dnakit.core.sequence import STRICT_SYMBOLS, DNASequence
from dnakit.exceptions import ConfigurationError, SequenceError, UnsupportedGapOperationError
from dnakit.io.annotations import AGPComponent, AGPDocument, AGPEntry, AGPGap

GapNormalizable: TypeAlias = DNA | DNASequence | DNARecord
AGPComponentValue: TypeAlias = DNA | DNASequence | DNARecord
CircularBoundaryPolicy = Literal["separate", "error"]

_N_RUN = re.compile(r"N+")
_GAP_AUDIT_KEY = "dnakit_gap_normalization"
_AGP_AUDIT_KEY = "dnakit_agp_assembly"
_AGP_KIND_BY_TYPE = {
    "scaffold": GapKind.SCAFFOLD,
    "contig": GapKind.CONTIG,
    "centromere": GapKind.CENTROMERE,
    "short_arm": GapKind.SHORT_ARM,
    "heterochromatin": GapKind.HETEROCHROMATIN,
    "telomere": GapKind.TELOMERE,
    "repeat": GapKind.REPEAT,
    "contamination": GapKind.CONTAMINATION,
}


def _native_provenance() -> Provenance:
    return Provenance(
        implementation=ImplementationInfo(
            label=ImplementationLabel.NATIVE,
            execution_mode=ExecutionMode.INTERNAL,
            origin_class=OriginClass.DNAKIT,
        )
    )


def _positive_integer(value: object, name: str, code: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigurationError(
            f"{name} must be a positive integer.",
            code=code,
            context={"field": name, "value": value},
        )
    return value


def _non_empty_optional_text(value: object, name: str, code: str) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            f"{name} must be a non-empty string or None.",
            code=code,
            context={"field": name},
        )
    return value


def _coerce_evidence(value: object, code: str) -> tuple[str, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Iterable):
        raise ConfigurationError(
            "evidence must be an iterable of non-empty strings.",
            code=code,
        )
    evidence = tuple(value)
    if any(not isinstance(item, str) or not item.strip() for item in evidence):
        raise ConfigurationError(
            "evidence must contain only non-empty strings.",
            code=code,
        )
    if len(set(evidence)) != len(evidence):
        raise ConfigurationError("evidence must not contain duplicates.", code=code)
    return cast(tuple[str, ...], evidence)


@dataclass(frozen=True, slots=True)
class GapNormalizationConfig:
    """Strict policy and resource limits for :func:`normalize_gaps`."""

    min_run_length: int = 10
    kind: GapKind | str = GapKind.UNKNOWN
    crossable: bool | None = None
    evidence: tuple[str, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)
    circular_boundary_policy: CircularBoundaryPolicy = "separate"
    max_input_symbols: int = 10_000_000
    max_n_run_length: int = 1_000_000
    max_converted_runs: int = 100_000
    max_output_parts: int = 200_001

    def __post_init__(self) -> None:
        for name in (
            "min_run_length",
            "max_input_symbols",
            "max_n_run_length",
            "max_converted_runs",
            "max_output_parts",
        ):
            _positive_integer(getattr(self, name), name, "INVALID_GAP_NORMALIZATION_LIMIT")
        if self.min_run_length > self.max_n_run_length:
            raise ConfigurationError(
                "min_run_length must not exceed max_n_run_length.",
                code="INVALID_GAP_NORMALIZATION_LIMIT",
            )
        try:
            kind = self.kind if isinstance(self.kind, GapKind) else GapKind(self.kind)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "kind is not a supported GapKind.",
                code="INVALID_GAP_NORMALIZATION_KIND",
                context={"kind": self.kind},
            ) from exc
        if self.crossable is not None and not isinstance(self.crossable, bool):
            raise ConfigurationError(
                "crossable must be True, False, or None.",
                code="INVALID_GAP_NORMALIZATION_POLICY",
            )
        evidence = _coerce_evidence(self.evidence, "INVALID_GAP_NORMALIZATION_EVIDENCE")
        try:
            metadata = freeze_mapping(self.metadata)
        except (AttributeError, ConfigurationError, TypeError, ValueError) as exc:
            raise ConfigurationError(
                "metadata must be a JSON-compatible mapping.",
                code="INVALID_GAP_NORMALIZATION_METADATA",
            ) from exc
        if _GAP_AUDIT_KEY in metadata:
            raise ConfigurationError(
                f"metadata key {_GAP_AUDIT_KEY!r} is reserved for source audit data.",
                code="RESERVED_GAP_NORMALIZATION_METADATA",
            )
        if self.circular_boundary_policy not in {"separate", "error"}:
            raise ConfigurationError(
                "circular_boundary_policy must be 'separate' or 'error'.",
                code="INVALID_GAP_NORMALIZATION_POLICY",
            )
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "evidence", evidence)
        object.__setattr__(self, "metadata", metadata)


@dataclass(frozen=True, slots=True)
class GapNormalizationChange:
    """One N-run replacement with original zero-based source coordinates."""

    source_part_index: int
    part_interval: Interval
    symbol_interval: Interval
    coordinate_interval: Interval | None
    original_symbol: str
    original_length: int
    original_sha256: str
    replacement: Gap

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class GapNormalizationResult:
    """Immutable output and audit trail for an N-run normalization operation."""

    sequence: DNASequence
    record: DNARecord | None
    source_kind: Literal["sequence", "record"]
    source_sha256: str
    normalized_parts: tuple[str | Gap, ...]
    changes: tuple[GapNormalizationChange, ...]
    preserved_gap_count: int
    parameters: FrozenDict
    algorithm_version: str
    provenance: Provenance
    issues: tuple[Issue, ...]

    @property
    def was_modified(self) -> bool:
        return bool(self.changes)

    @property
    def output(self) -> DNASequence | DNARecord:
        return self.record if self.record is not None else self.sequence

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def normalize_gaps(
    value: GapNormalizable,
    *,
    config: GapNormalizationConfig | None = None,
) -> GapNormalizationResult:
    """Convert qualifying N-runs to known gaps without mutating the source object.

    Runs are evaluated independently inside each text part. Existing explicit
    :class:`Gap` objects are retained by identity and are never merged with a
    newly created gap. Circular origin-spanning runs are deliberately split at
    the stored origin unless ``circular_boundary_policy='error'`` is selected.
    """

    if not isinstance(value, (DNA, DNASequence, DNARecord)):
        raise TypeError("value must be DNA, DNASequence, or DNARecord.")
    resolved = GapNormalizationConfig() if config is None else config
    if not isinstance(resolved, GapNormalizationConfig):
        raise ConfigurationError(
            "config must be GapNormalizationConfig or None.",
            code="INVALID_GAP_NORMALIZATION_CONFIG",
        )
    source_record: DNARecord | None
    source_sequence: DNASequence
    if isinstance(value, DNA):
        source_record = value.record
        source_sequence = value.sequence
    elif isinstance(value, DNARecord):
        source_record = value
        source_sequence = value.sequence
    else:
        source_record = None
        source_sequence = value
    if source_sequence.symbol_length > resolved.max_input_symbols:
        raise SequenceError(
            "Sequence exceeds max_input_symbols.",
            code="GAP_NORMALIZATION_INPUT_LIMIT_EXCEEDED",
            context={
                "symbol_length": source_sequence.symbol_length,
                "max_input_symbols": resolved.max_input_symbols,
            },
        )

    candidate_runs = _candidate_n_runs(source_sequence, resolved)
    if candidate_runs and source_record is not None and source_record.letter_annotations:
        raise UnsupportedGapOperationError(
            "Converting N-runs removes symbol positions used by letter annotations.",
            code="GAP_NORMALIZATION_LETTER_ANNOTATIONS_UNSUPPORTED",
            context={
                "record_id": source_record.id,
                "annotation_names": tuple(source_record.letter_annotations),
            },
            hint="Remove or remap per-symbol annotations explicitly before converting runs.",
        )

    issues = _circular_boundary_issues(source_sequence, resolved)
    output_parts: list[str | Gap] = []
    changes: list[GapNormalizationChange] = []
    symbol_cursor = 0
    coordinate_cursor: int | None = 0
    preserved_gap_count = 0

    for part_index, part in enumerate(source_sequence.parts):
        if isinstance(part, Gap):
            output_parts.append(part)
            _check_output_part_limit(output_parts, resolved.max_output_parts)
            preserved_gap_count += 1
            if coordinate_cursor is not None:
                coordinate_cursor = None if part.length is None else coordinate_cursor + part.length
            continue

        part_cursor = 0
        for match in _N_RUN.finditer(part):
            start, end = match.span()
            run_length = end - start
            if run_length < resolved.min_run_length:
                continue
            if run_length > resolved.max_n_run_length:
                raise SequenceError(
                    "An N-run exceeds max_n_run_length.",
                    code="GAP_NORMALIZATION_RUN_LIMIT_EXCEEDED",
                    context={
                        "source_part_index": part_index,
                        "part_start": start,
                        "run_length": run_length,
                        "max_n_run_length": resolved.max_n_run_length,
                    },
                )
            if len(changes) >= resolved.max_converted_runs:
                raise SequenceError(
                    "N-run count exceeds max_converted_runs.",
                    code="GAP_NORMALIZATION_CHANGE_LIMIT_EXCEEDED",
                    context={"max_converted_runs": resolved.max_converted_runs},
                )
            if start > part_cursor:
                output_parts.append(part[part_cursor:start])
                _check_output_part_limit(output_parts, resolved.max_output_parts)
            coordinate_interval = (
                None
                if coordinate_cursor is None
                else Interval(coordinate_cursor + start, coordinate_cursor + end)
            )
            audit = {
                "source": "N-run",
                "original_symbol": "N",
                "original_length": run_length,
                "source_part_index": part_index,
                "part_interval": {"start": start, "end": end},
                "symbol_interval": {
                    "start": symbol_cursor + start,
                    "end": symbol_cursor + end,
                },
            }
            gap_metadata = dict(resolved.metadata)
            gap_metadata[_GAP_AUDIT_KEY] = audit
            replacement = Gap(
                run_length,
                kind=resolved.kind,
                crossable=resolved.crossable,
                evidence=resolved.evidence,
                metadata=gap_metadata,
            )
            output_parts.append(replacement)
            _check_output_part_limit(output_parts, resolved.max_output_parts)
            changes.append(
                GapNormalizationChange(
                    source_part_index=part_index,
                    part_interval=Interval(start, end),
                    symbol_interval=Interval(symbol_cursor + start, symbol_cursor + end),
                    coordinate_interval=coordinate_interval,
                    original_symbol="N",
                    original_length=run_length,
                    original_sha256=hashlib.sha256(match.group().encode("ascii")).hexdigest(),
                    replacement=replacement,
                )
            )
            part_cursor = end
        if part_cursor < len(part):
            output_parts.append(part[part_cursor:])
            _check_output_part_limit(output_parts, resolved.max_output_parts)
        symbol_cursor += len(part)
        if coordinate_cursor is not None:
            coordinate_cursor += len(part)

    if not changes:
        output_sequence = source_sequence
        output_record = source_record
        normalized_parts = source_sequence.parts
    else:
        remaining_symbols = sum(len(part) for part in output_parts if isinstance(part, str))
        if source_sequence.topology is Topology.CIRCULAR and remaining_symbols == 0:
            raise UnsupportedGapOperationError(
                "Gap normalization would leave a circular sequence with no nucleotide symbols.",
                code="GAP_NORMALIZATION_EMPTY_CIRCULAR_OUTPUT",
            )
        output_sequence = DNASequence(
            output_parts,
            alphabet=source_sequence.alphabet,
            topology=source_sequence.topology,
            strandedness=source_sequence.strandedness,
        )
        normalized_parts = output_sequence.parts
        output_record = (
            None
            if source_record is None
            else DNARecord(
                output_sequence,
                source_record.id,
                description=source_record.description,
                features=source_record.features,
                metadata=source_record.metadata,
            )
        )

    return GapNormalizationResult(
        sequence=output_sequence,
        record=output_record,
        source_kind="record" if source_record is not None else "sequence",
        source_sha256=_sequence_sha256(source_sequence),
        normalized_parts=normalized_parts,
        changes=tuple(changes),
        preserved_gap_count=preserved_gap_count,
        parameters=freeze_mapping(
            {
                "min_run_length": resolved.min_run_length,
                "kind": cast(GapKind, resolved.kind).value,
                "crossable": resolved.crossable,
                "evidence": resolved.evidence,
                "metadata": resolved.metadata,
                "circular_boundary_policy": resolved.circular_boundary_policy,
                "max_input_symbols": resolved.max_input_symbols,
                "max_n_run_length": resolved.max_n_run_length,
                "max_converted_runs": resolved.max_converted_runs,
                "max_output_parts": resolved.max_output_parts,
                "coordinate_system": "0-based-half-open",
            }
        ),
        algorithm_version="dnakit-gap-normalization-v1",
        provenance=_native_provenance(),
        issues=tuple(issues),
    )


def _check_output_part_limit(parts: list[str | Gap], max_output_parts: int) -> None:
    if len(parts) > max_output_parts:
        raise SequenceError(
            "Normalized sequence exceeds max_output_parts.",
            code="GAP_NORMALIZATION_OUTPUT_LIMIT_EXCEEDED",
            context={"max_output_parts": max_output_parts},
        )


def _sequence_sha256(sequence: DNASequence) -> str:
    payload = json.dumps(
        to_json_compatible(sequence),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _candidate_n_runs(sequence: DNASequence, config: GapNormalizationConfig) -> int:
    count = 0
    for part in sequence.parts:
        if not isinstance(part, str):
            continue
        for match in _N_RUN.finditer(part):
            run_length = len(match.group())
            if run_length < config.min_run_length:
                continue
            if run_length > config.max_n_run_length:
                raise SequenceError(
                    "An N-run exceeds max_n_run_length.",
                    code="GAP_NORMALIZATION_RUN_LIMIT_EXCEEDED",
                    context={
                        "run_length": run_length,
                        "max_n_run_length": config.max_n_run_length,
                    },
                )
            count += 1
            if count > config.max_converted_runs:
                raise SequenceError(
                    "N-run count exceeds max_converted_runs.",
                    code="GAP_NORMALIZATION_CHANGE_LIMIT_EXCEEDED",
                    context={"max_converted_runs": config.max_converted_runs},
                )
    return count


def _circular_boundary_issues(
    sequence: DNASequence,
    config: GapNormalizationConfig,
) -> list[Issue]:
    if sequence.topology is not Topology.CIRCULAR or len(sequence.parts) == 0:
        return []
    first = sequence.parts[0]
    last = sequence.parts[-1]
    if not isinstance(first, str) or not isinstance(last, str):
        return []
    leading = len(first) - len(first.lstrip("N"))
    trailing = len(last) - len(last.rstrip("N"))
    if leading == 0 or trailing == 0 or leading + trailing < config.min_run_length:
        return []
    if config.circular_boundary_policy == "error":
        raise UnsupportedGapOperationError(
            "A qualifying N-run spans the stored origin of a circular sequence.",
            code="GAP_NORMALIZATION_CIRCULAR_BOUNDARY_UNRESOLVED",
            context={"leading_n": leading, "trailing_n": trailing},
            hint="Rotate the sequence origin or use circular_boundary_policy='separate'.",
        )
    return [
        Issue(
            "STD_CIRCULAR_N_RUN_SPLIT_AT_ORIGIN",
            IssueSeverity.INFO,
            "A circular origin-spanning N-run was evaluated as two stored linear runs.",
            details={
                "leading_n": leading,
                "trailing_n": trailing,
                "policy": "separate",
            },
        )
    ]


@dataclass(frozen=True, slots=True)
class AGPAssemblyConfig:
    """Strict output choices and resource limits for :func:`sequence_from_agp`."""

    output_alphabet: DNAAlphabet | str | None = None
    output_strandedness: Strandedness | str = Strandedness.SINGLE
    record_description: str = ""
    record_metadata: Mapping[str, object] = field(default_factory=dict)
    max_entries: int = 1_000_000
    max_components: int = 1_000_000
    max_component_symbols: int = 100_000_000
    max_output_symbols: int = 1_000_000_000
    max_output_span: int = 1_000_000_000
    max_header_lines: int = 10_000
    max_header_length: int = 1_000_000

    def __post_init__(self) -> None:
        for name in (
            "max_entries",
            "max_components",
            "max_component_symbols",
            "max_output_symbols",
            "max_output_span",
            "max_header_lines",
            "max_header_length",
        ):
            _positive_integer(getattr(self, name), name, "INVALID_AGP_ASSEMBLY_LIMIT")
        try:
            alphabet = (
                None
                if self.output_alphabet is None
                else self.output_alphabet
                if isinstance(self.output_alphabet, DNAAlphabet)
                else DNAAlphabet(self.output_alphabet)
            )
            strandedness = (
                self.output_strandedness
                if isinstance(self.output_strandedness, Strandedness)
                else Strandedness(self.output_strandedness)
            )
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Invalid AGP output alphabet or strandedness.",
                code="INVALID_AGP_ASSEMBLY_CONFIG",
            ) from exc
        if not isinstance(self.record_description, str):
            raise ConfigurationError(
                "record_description must be a string.",
                code="INVALID_AGP_ASSEMBLY_CONFIG",
            )
        try:
            metadata = freeze_mapping(self.record_metadata)
        except (AttributeError, ConfigurationError, TypeError, ValueError) as exc:
            raise ConfigurationError(
                "record_metadata must be a JSON-compatible mapping.",
                code="INVALID_AGP_ASSEMBLY_METADATA",
            ) from exc
        if _AGP_AUDIT_KEY in metadata:
            raise ConfigurationError(
                f"record_metadata key {_AGP_AUDIT_KEY!r} is reserved.",
                code="RESERVED_AGP_ASSEMBLY_METADATA",
            )
        object.__setattr__(self, "output_alphabet", alphabet)
        object.__setattr__(self, "output_strandedness", strandedness)
        object.__setattr__(self, "record_metadata", metadata)


@dataclass(frozen=True, slots=True)
class AGPAssemblySegment:
    """Auditable source coordinates for one assembled AGP row."""

    object_id: str
    object_interval: Interval
    part_number: int
    segment_type: Literal["component", "known_gap", "unknown_gap"]
    agp_component_type: str
    component_id: str | None = None
    component_interval: Interval | None = None
    orientation: str | None = None
    gap_type: str | None = None
    linkage: bool | None = None
    linkage_evidence: tuple[str, ...] = ()
    gap: Gap | None = None

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True, slots=True)
class AGPAssemblyResult:
    """A strict AGP assembly record with retained source-row coordinates."""

    record: DNARecord
    object_id: str
    headers: tuple[str, ...]
    segments: tuple[AGPAssemblySegment, ...]
    used_components: tuple[str, ...]
    parameters: FrozenDict
    algorithm_version: str
    provenance: Provenance
    issues: tuple[Issue, ...]

    @property
    def sequence(self) -> DNASequence:
        return self.record.sequence

    @property
    def coordinate_span_unresolved(self) -> bool:
        return self.sequence.coordinate_span is None

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def sequence_from_agp(
    document: AGPDocument | Iterable[AGPEntry],
    components: Mapping[str, AGPComponentValue],
    *,
    object_id: str | None = None,
    record_id: str | None = None,
    config: AGPAssemblyConfig | None = None,
) -> AGPAssemblyResult:
    """Build one linear record from AGP rows without inferring missing semantics.

    Only the selected object's component values are accessed. Consequently a
    lazy mapping may resolve component sequences on demand without materializing
    an entire assembly database.
    """

    resolved = AGPAssemblyConfig() if config is None else config
    if not isinstance(resolved, AGPAssemblyConfig):
        raise ConfigurationError(
            "config must be AGPAssemblyConfig or None.",
            code="INVALID_AGP_ASSEMBLY_CONFIG",
        )
    selected_object_id = _non_empty_optional_text(object_id, "object_id", "INVALID_AGP_OBJECT_ID")
    selected_record_id = _non_empty_optional_text(record_id, "record_id", "INVALID_AGP_RECORD_ID")
    if not isinstance(components, Mapping):
        raise ConfigurationError(
            "components must map component IDs to DNASequence or DNARecord objects.",
            code="INVALID_AGP_COMPONENT_MAPPING",
        )
    entries, headers = _materialize_agp_entries(document, resolved)
    object_ids = _agp_object_blocks(entries)
    if selected_object_id is None:
        if len(object_ids) != 1:
            raise ConfigurationError(
                "object_id is required when AGP entries contain multiple objects.",
                code="AGP_OBJECT_ID_REQUIRED",
                context={"object_ids": object_ids},
            )
        selected_object_id = object_ids[0]
    selected = tuple(entry for entry in entries if entry.object_id == selected_object_id)
    if not selected:
        raise SequenceError(
            "The requested AGP object does not exist.",
            code="AGP_OBJECT_NOT_FOUND",
            context={"object_id": selected_object_id},
        )

    _validate_agp_continuity(selected, resolved.max_output_span)
    output_parts: list[str | Gap] = []
    segments: list[AGPAssemblySegment] = []
    used_components: list[str] = []
    seen_components: set[str] = set()
    inferred_alphabet = DNAAlphabet.STRICT
    output_symbol_count = 0
    component_count = 0
    unknown_gaps: list[AGPGap] = []

    for entry in selected:
        if isinstance(entry, AGPGap):
            _validate_agp_gap(entry)
            output_parts.append(entry.gap)
            segment_type: Literal["known_gap", "unknown_gap"] = (
                "known_gap" if entry.component_type == "N" else "unknown_gap"
            )
            if entry.component_type == "U":
                unknown_gaps.append(entry)
            segments.append(
                AGPAssemblySegment(
                    entry.object_id,
                    entry.object_interval,
                    entry.part_number,
                    segment_type,
                    entry.component_type,
                    gap_type=entry.gap_type,
                    linkage=entry.linkage,
                    linkage_evidence=entry.linkage_evidence,
                    gap=entry.gap,
                )
            )
            continue

        component_count += 1
        if component_count > resolved.max_components:
            raise SequenceError(
                "Selected AGP object exceeds max_components.",
                code="AGP_COMPONENT_LIMIT_EXCEEDED",
                context={"max_components": resolved.max_components},
            )
        remaining_output_symbols = resolved.max_output_symbols - output_symbol_count
        if len(entry.component_interval) > remaining_output_symbols:
            raise SequenceError(
                "Selected AGP object exceeds max_output_symbols.",
                code="AGP_OUTPUT_SYMBOL_LIMIT_EXCEEDED",
                context={
                    "part_number": entry.part_number,
                    "requested_fragment_symbols": len(entry.component_interval),
                    "remaining_output_symbols": remaining_output_symbols,
                    "max_output_symbols": resolved.max_output_symbols,
                },
            )
        fragment, component_alphabet = _resolve_agp_component(entry, components, resolved)
        output_symbol_count += len(fragment)
        if output_symbol_count > resolved.max_output_symbols:
            raise SequenceError(
                "Selected AGP object exceeds max_output_symbols.",
                code="AGP_OUTPUT_SYMBOL_LIMIT_EXCEEDED",
                context={"max_output_symbols": resolved.max_output_symbols},
            )
        output_parts.append(fragment)
        if component_alphabet is DNAAlphabet.IUPAC:
            inferred_alphabet = DNAAlphabet.IUPAC
        if entry.component_id not in seen_components:
            used_components.append(entry.component_id)
            seen_components.add(entry.component_id)
        segments.append(
            AGPAssemblySegment(
                entry.object_id,
                entry.object_interval,
                entry.part_number,
                "component",
                entry.component_type,
                component_id=entry.component_id,
                component_interval=entry.component_interval,
                orientation=entry.orientation,
            )
        )

    alphabet = cast(DNAAlphabet | None, resolved.output_alphabet) or inferred_alphabet
    if alphabet is DNAAlphabet.STRICT and any(
        isinstance(part, str) and set(part) - STRICT_SYMBOLS for part in output_parts
    ):
        raise SequenceError(
            "An assembled component fragment is incompatible with strict output alphabet.",
            code="AGP_OUTPUT_ALPHABET_MISMATCH",
        )
    output_sequence = DNASequence(
        output_parts,
        alphabet=alphabet,
        topology=Topology.LINEAR,
        strandedness=cast(Strandedness, resolved.output_strandedness),
    )
    metadata = dict(resolved.record_metadata)
    metadata[_AGP_AUDIT_KEY] = {
        "object_id": selected_object_id,
        "coordinate_system": "0-based-half-open",
        "entry_count": len(selected),
        "component_count": component_count,
        "known_gap_count": sum(
            isinstance(entry, AGPGap) and entry.component_type == "N" for entry in selected
        ),
        "unknown_gap_count": len(unknown_gaps),
    }
    output_record = DNARecord(
        output_sequence,
        selected_record_id or selected_object_id,
        description=resolved.record_description,
        metadata=metadata,
    )
    issues: list[Issue] = []
    if unknown_gaps:
        issues.append(
            Issue(
                "STD_AGP_COORDINATE_SPAN_UNRESOLVED",
                IssueSeverity.INFO,
                "AGP U rows retain Gap(length=None); core coordinate_span is therefore unresolved.",
                details={
                    "unknown_gap_count": len(unknown_gaps),
                    "agp_placeholder_intervals": [
                        {
                            "start": entry.object_interval.start,
                            "end": entry.object_interval.end,
                        }
                        for entry in unknown_gaps
                    ],
                },
            )
        )
    return AGPAssemblyResult(
        record=output_record,
        object_id=selected_object_id,
        headers=headers,
        segments=tuple(segments),
        used_components=tuple(used_components),
        parameters=freeze_mapping(
            {
                "object_id": selected_object_id,
                "record_id": selected_record_id or selected_object_id,
                "output_alphabet": alphabet.value,
                "output_strandedness": cast(Strandedness, resolved.output_strandedness).value,
                "output_topology": Topology.LINEAR.value,
                "coordinate_system": "0-based-half-open",
                "record_description": resolved.record_description,
                "record_metadata": resolved.record_metadata,
                "max_entries": resolved.max_entries,
                "max_components": resolved.max_components,
                "max_component_symbols": resolved.max_component_symbols,
                "max_output_symbols": resolved.max_output_symbols,
                "max_output_span": resolved.max_output_span,
                "max_header_lines": resolved.max_header_lines,
                "max_header_length": resolved.max_header_length,
            }
        ),
        algorithm_version="dnakit-agp-sequence-assembly-v1",
        provenance=Provenance(
            implementation=ImplementationInfo(
                label=ImplementationLabel.REIMPLEMENTATION,
                execution_mode=ExecutionMode.INTERNAL,
                origin_class=OriginClass.STANDARD,
            )
        ),
        issues=tuple(issues),
    )


def _materialize_agp_entries(
    document: AGPDocument | Iterable[AGPEntry], config: AGPAssemblyConfig
) -> tuple[tuple[AGPEntry, ...], tuple[str, ...]]:
    if isinstance(document, AGPDocument):
        if len(document.entries) > config.max_entries:
            raise SequenceError(
                "AGP document exceeds max_entries.",
                code="AGP_ENTRY_LIMIT_EXCEEDED",
                context={"max_entries": config.max_entries},
            )
        if len(document.headers) > config.max_header_lines:
            raise SequenceError(
                "AGP document exceeds max_header_lines.",
                code="AGP_HEADER_LIMIT_EXCEEDED",
                context={"max_header_lines": config.max_header_lines},
            )
        oversized_header = next(
            (
                index
                for index, header in enumerate(document.headers)
                if len(header) > config.max_header_length
            ),
            None,
        )
        if oversized_header is not None:
            raise SequenceError(
                "An AGP header exceeds max_header_length.",
                code="AGP_HEADER_LENGTH_LIMIT_EXCEEDED",
                context={
                    "header_index": oversized_header,
                    "max_header_length": config.max_header_length,
                },
            )
        entries = document.entries
        headers = document.headers
    else:
        raw_document: object = document
        if isinstance(raw_document, (str, bytes)) or not isinstance(raw_document, Iterable):
            raise TypeError("document must be AGPDocument or an iterable of AGP entries.")
        materialized = tuple(
            itertools.islice(iter(cast(Iterable[object], raw_document)), config.max_entries + 1)
        )
        if len(materialized) > config.max_entries:
            raise SequenceError(
                "AGP entry iterable exceeds max_entries.",
                code="AGP_ENTRY_LIMIT_EXCEEDED",
                context={"max_entries": config.max_entries},
            )
        validated_entries: list[AGPEntry] = []
        for index, entry in enumerate(materialized):
            if not isinstance(entry, (AGPComponent, AGPGap)):
                raise SequenceError(
                    "AGP iterable contains an unsupported entry.",
                    code="INVALID_AGP_ASSEMBLY_ENTRY",
                    context={"entry_index": index, "entry_type": type(entry).__name__},
                )
            validated_entries.append(entry)
        entries = tuple(validated_entries)
        headers = ()
    if not entries:
        raise SequenceError("AGP input is empty.", code="EMPTY_AGP_ASSEMBLY")
    for index, entry in enumerate(entries):
        if not isinstance(entry, (AGPComponent, AGPGap)):
            raise SequenceError(
                "AGP iterable contains an unsupported entry.",
                code="INVALID_AGP_ASSEMBLY_ENTRY",
                context={"entry_index": index, "entry_type": type(entry).__name__},
            )
    return entries, headers


def _agp_object_blocks(entries: tuple[AGPEntry, ...]) -> tuple[str, ...]:
    object_ids: list[str] = []
    completed: set[str] = set()
    active: str | None = None
    for entry in entries:
        if entry.object_id == active:
            continue
        if active is not None:
            completed.add(active)
        if entry.object_id in completed:
            raise SequenceError(
                "AGP object rows must form one contiguous block.",
                code="AGP_OBJECT_BLOCK_DISCONTIGUOUS",
                context={"object_id": entry.object_id},
            )
        active = entry.object_id
        object_ids.append(active)
    return tuple(object_ids)


def _validate_agp_continuity(entries: tuple[AGPEntry, ...], max_output_span: int) -> None:
    expected_start = 0
    expected_part = 1
    for entry in entries:
        if entry.object_interval.start != expected_start or entry.part_number != expected_part:
            raise SequenceError(
                "AGP object coordinates and part numbers must be contiguous from zero/one.",
                code="INVALID_AGP_ASSEMBLY_CONTINUITY",
                context={
                    "part_number": entry.part_number,
                    "object_start": entry.object_interval.start,
                    "expected_part_number": expected_part,
                    "expected_object_start": expected_start,
                },
            )
        if entry.object_interval.end > max_output_span:
            raise SequenceError(
                "AGP object exceeds max_output_span.",
                code="AGP_OUTPUT_SPAN_LIMIT_EXCEEDED",
                context={"max_output_span": max_output_span},
            )
        expected_start = entry.object_interval.end
        expected_part += 1


def _validate_agp_gap(entry: AGPGap) -> None:
    span = len(entry.object_interval)
    expected_kind = _AGP_KIND_BY_TYPE.get(entry.gap_type)
    if expected_kind is None or entry.gap.kind is not expected_kind:
        raise SequenceError(
            "AGP gap_type and Gap.kind do not agree.",
            code="INVALID_AGP_ASSEMBLY_GAP_KIND",
            context={"gap_type": entry.gap_type, "kind": entry.gap.kind.value},
        )
    if entry.gap.crossable is not entry.linkage or entry.gap.evidence != entry.linkage_evidence:
        raise SequenceError(
            "AGP linkage fields and embedded Gap do not agree.",
            code="INVALID_AGP_ASSEMBLY_GAP_LINKAGE",
            context={"part_number": entry.part_number},
        )
    metadata_type = entry.gap.metadata.get("agp_gap_type")
    if metadata_type is not None and metadata_type != entry.gap_type:
        raise SequenceError(
            "AGP gap metadata conflicts with gap_type.",
            code="INVALID_AGP_ASSEMBLY_GAP_METADATA",
            context={"part_number": entry.part_number},
        )
    valid_length = (entry.component_type == "N" and entry.gap.length == span) or (
        entry.component_type == "U" and entry.gap.length is None and span == 100
    )
    if not valid_length:
        raise SequenceError(
            "AGP known/unknown gap representation does not match its object span.",
            code="INVALID_AGP_ASSEMBLY_GAP_LENGTH",
            context={"part_number": entry.part_number, "component_type": entry.component_type},
        )


def _resolve_agp_component(
    entry: AGPComponent,
    components: Mapping[str, AGPComponentValue],
    config: AGPAssemblyConfig,
) -> tuple[str, DNAAlphabet]:
    if len(entry.component_interval) != len(entry.object_interval):
        raise SequenceError(
            "AGP component and object interval lengths do not match.",
            code="INVALID_AGP_ASSEMBLY_COMPONENT_SPAN",
            context={"part_number": entry.part_number},
        )
    try:
        value = components[entry.component_id]
    except KeyError as exc:
        raise SequenceError(
            "AGP component sequence is missing.",
            code="AGP_COMPONENT_NOT_FOUND",
            context={"component_id": entry.component_id},
        ) from exc
    if not isinstance(value, (DNA, DNASequence, DNARecord)):
        raise SequenceError(
            "AGP component mapping value must be DNA, DNASequence, or DNARecord.",
            code="INVALID_AGP_COMPONENT_VALUE",
            context={"component_id": entry.component_id, "value_type": type(value).__name__},
        )
    sequence = value.sequence if isinstance(value, (DNA, DNARecord)) else value
    if sequence.topology is Topology.CIRCULAR:
        raise UnsupportedGapOperationError(
            "AGP component extraction from a circular sequence requires an explicit origin policy.",
            code="AGP_CIRCULAR_COMPONENT_UNSUPPORTED",
            context={"component_id": entry.component_id},
        )
    if sequence.is_gapped:
        raise UnsupportedGapOperationError(
            "AGP component coordinates cannot be applied to an already-gapped sequence.",
            code="AGP_GAPPED_COMPONENT_UNSUPPORTED",
            context={"component_id": entry.component_id},
        )
    if sequence.symbol_length > config.max_component_symbols:
        raise SequenceError(
            "AGP component exceeds max_component_symbols.",
            code="AGP_COMPONENT_SYMBOL_LIMIT_EXCEEDED",
            context={
                "component_id": entry.component_id,
                "symbol_length": sequence.symbol_length,
                "max_component_symbols": config.max_component_symbols,
            },
        )
    if entry.component_interval.end > sequence.symbol_length:
        raise SequenceError(
            "AGP component interval exceeds its sequence.",
            code="AGP_COMPONENT_INTERVAL_OUT_OF_BOUNDS",
            context={
                "component_id": entry.component_id,
                "component_end": entry.component_interval.end,
                "symbol_length": sequence.symbol_length,
            },
        )
    if entry.orientation not in {"+", "-"}:
        raise UnsupportedGapOperationError(
            "AGP component orientation is unresolved; DNAKit will not guess it.",
            code="AGP_COMPONENT_ORIENTATION_UNRESOLVED",
            context={"component_id": entry.component_id, "orientation": entry.orientation},
        )
    fragment = sequence.symbols[entry.component_interval.start : entry.component_interval.end]
    if len(fragment) != len(entry.object_interval):
        raise SequenceError(
            "Extracted AGP component length does not equal its object span.",
            code="INVALID_AGP_ASSEMBLY_COMPONENT_SPAN",
            context={"component_id": entry.component_id},
        )
    if entry.orientation == "-":
        fragment = DNASequence(fragment, alphabet=sequence.alphabet).reverse_complement().symbols
    return fragment, sequence.alphabet


__all__ = [
    "AGPAssemblyConfig",
    "AGPAssemblyResult",
    "AGPAssemblySegment",
    "AGPComponentValue",
    "CircularBoundaryPolicy",
    "GapNormalizable",
    "GapNormalizationChange",
    "GapNormalizationConfig",
    "GapNormalizationResult",
    "normalize_gaps",
    "sequence_from_agp",
]
