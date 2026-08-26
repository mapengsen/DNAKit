"""Sequence-, record- and collection-level DNA validation."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable
from typing import overload

from dnakit.core.collection import DNASet
from dnakit.core.enums import DNAAlphabet, IssueSeverity
from dnakit.core.facade import DNA
from dnakit.core.gap import Gap
from dnakit.core.issues import Issue
from dnakit.core.provenance import Provenance
from dnakit.core.record import DNARecord
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import ConfigurationError

from ._shared import CANONICAL_BASES, IUPAC_BASES, ambiguity_report
from .config import DatasetValidationConfig, ValidationConfig
from .results import (
    DatasetValidationReport,
    DuplicateID,
    QualitySummary,
    ValidationReport,
)


@overload
def validate(  # type: ignore[overload-overlap]
    value: DNA,
    *,
    config: ValidationConfig | DatasetValidationConfig | None = None,
) -> ValidationReport | DatasetValidationReport: ...


@overload
def validate(
    value: DNASequence | DNARecord,
    *,
    config: ValidationConfig | DatasetValidationConfig | None = None,
) -> ValidationReport: ...


@overload
def validate(
    value: DNASet | Iterable[DNARecord],
    *,
    config: ValidationConfig | DatasetValidationConfig | None = None,
) -> DatasetValidationReport: ...


def validate(
    value: DNA | DNASequence | DNARecord | DNASet | Iterable[DNARecord],
    *,
    config: ValidationConfig | DatasetValidationConfig | None = None,
) -> ValidationReport | DatasetValidationReport:
    """Validate one sequence/record or a record collection with one entry point.

    The returned report type follows the input type. A collection is invalid
    when any record violates its configured rules, when a required unique-ID
    rule is violated, or when the collection is empty.
    """
    if isinstance(value, DNA):
        if value.is_single:
            return _validate_one(value.record, config=_resolve_record_config(config))
        return _validate_collection(value.records, config=_resolve_dataset_config(config))
    if isinstance(value, (DNASequence, DNARecord)):
        return _validate_one(value, config=_resolve_record_config(config))
    if isinstance(value, Iterable):
        return _validate_collection(value, config=_resolve_dataset_config(config))
    raise TypeError(
        "validate expects a DNASequence, DNARecord, DNASet, or an iterable of DNARecord."
    )


def _resolve_record_config(
    config: ValidationConfig | DatasetValidationConfig | None,
) -> ValidationConfig:
    if config is None:
        return ValidationConfig()
    if isinstance(config, ValidationConfig):
        return config
    if isinstance(config, DatasetValidationConfig):
        return config.record
    raise ConfigurationError(
        "config must be ValidationConfig, DatasetValidationConfig, or None.",
        code="INVALID_VALIDATION_CONFIG",
    )


def _resolve_dataset_config(
    config: ValidationConfig | DatasetValidationConfig | None,
) -> DatasetValidationConfig:
    if config is None:
        return DatasetValidationConfig()
    if isinstance(config, DatasetValidationConfig):
        return config
    if isinstance(config, ValidationConfig):
        return DatasetValidationConfig(record=config)
    raise ConfigurationError(
        "config must be ValidationConfig, DatasetValidationConfig, or None.",
        code="INVALID_DATASET_VALIDATION_CONFIG",
    )


def _validate_one(
    sequence: DNASequence | DNARecord,
    *,
    config: ValidationConfig,
) -> ValidationReport:
    resolved = config
    if isinstance(sequence, DNARecord):
        record = sequence
        dna = record.sequence
        record_id: str | None = record.id
    elif isinstance(sequence, DNASequence):
        record = None
        dna = sequence
        record_id = None
    else:
        raise TypeError("validate expects a DNASequence or DNARecord.")

    issues: list[Issue] = []
    alphabet = resolved.alphabet or dna.alphabet
    allowed = CANONICAL_BASES if alphabet is DNAAlphabet.STRICT else IUPAC_BASES
    invalid_positions: dict[str, list[int]] = defaultdict(list)
    symbol_offset = 0
    unknown_gap_count = 0
    gap_count = 0

    for part_index, part in enumerate(dna.parts):
        if isinstance(part, Gap):
            gap_count += 1
            if not resolved.allow_gaps:
                issues.append(
                    _issue(
                        "STD_GAP_NOT_ALLOWED",
                        IssueSeverity.ERROR,
                        "The sequence contains a Gap but gaps are disabled.",
                        {"part_index": part_index},
                    )
                )
            if part.length is None:
                unknown_gap_count += 1
            continue
        for symbol in part:
            if symbol not in allowed:
                invalid_positions[symbol].append(symbol_offset)
            symbol_offset += 1

    for symbol, positions in sorted(invalid_positions.items()):
        issues.append(
            _issue(
                "STD_INVALID_SYMBOL",
                IssueSeverity.ERROR,
                f"Character {symbol!r} is not valid for the {alphabet.value} DNA alphabet.",
                {
                    "symbol": symbol,
                    "codepoint": f"U+{ord(symbol):04X}",
                    "positions": positions,
                },
            )
        )

    if unknown_gap_count and not resolved.allow_unknown_gap_length:
        issues.append(
            _issue(
                "STD_UNKNOWN_GAP_LENGTH",
                IssueSeverity.ERROR,
                f"Found {unknown_gap_count} Gap object(s) with unknown length.",
                {"count": unknown_gap_count},
            )
        )

    length = dna.symbol_length
    if length == 0 and not resolved.allow_empty:
        issues.append(
            _issue(
                "STD_EMPTY_SEQUENCE",
                IssueSeverity.ERROR,
                "The sequence contains no nucleotide symbols.",
                {},
            )
        )
    if resolved.sequence_length is not None and length != resolved.sequence_length:
        issues.append(
            _issue(
                "STD_SEQUENCE_LENGTH_MISMATCH",
                IssueSeverity.ERROR,
                f"Sequence length {length} does not match the configured exact length "
                f"{resolved.sequence_length}.",
                {"length": length, "expected_length": resolved.sequence_length},
            )
        )
    if (
        resolved.min_length is not None
        and length < resolved.min_length
        and not (length == 0 and (resolved.allow_empty or resolved.min_length > 0))
    ):
        issues.append(
            _issue(
                "STD_SEQUENCE_TOO_SHORT",
                IssueSeverity.ERROR,
                f"Sequence length {length} is below the configured minimum {resolved.min_length}.",
                {"length": length, "minimum": resolved.min_length},
            )
        )
    if resolved.max_length is not None and length > resolved.max_length:
        issues.append(
            _issue(
                "STD_SEQUENCE_TOO_LONG",
                IssueSeverity.ERROR,
                f"Sequence length {length} exceeds the configured maximum {resolved.max_length}.",
                {"length": length, "maximum": resolved.max_length},
            )
        )

    ambiguity = ambiguity_report(
        dna.parts,
        denominator_includes_gap=resolved.ambiguity_denominator_includes_gap,
    )
    if resolved.max_ambiguity_fraction is not None:
        if ambiguity.fraction is None and ambiguity.denominator is None:
            issues.append(
                _issue(
                    "STD_AMBIGUITY_DENOMINATOR_UNKNOWN",
                    IssueSeverity.WARNING,
                    "Ambiguity fraction cannot include an unknown-length Gap.",
                    {"gap_count": gap_count, "unknown_gap_count": unknown_gap_count},
                )
            )
        elif (
            ambiguity.fraction is not None and ambiguity.fraction > resolved.max_ambiguity_fraction
        ):
            issues.append(
                _issue(
                    "STD_AMBIGUITY_FRACTION_HIGH",
                    IssueSeverity.ERROR,
                    f"Ambiguity fraction {ambiguity.fraction:.6g} exceeds the "
                    f"configured maximum {resolved.max_ambiguity_fraction:.6g}.",
                    {
                        "fraction": ambiguity.fraction,
                        "maximum": resolved.max_ambiguity_fraction,
                    },
                )
            )

    quality: QualitySummary | None = None
    if record is not None:
        _validate_metadata(record, resolved, issues)
        quality = _validate_letter_annotations(record, resolved, issues)

    issue_tuple = tuple(issues)
    return ValidationReport(
        record_id=record_id,
        sequence=dna,
        config=resolved,
        algorithm_version="std-validate-v2",
        provenance=Provenance(),
        is_valid=not _has_error(issue_tuple),
        symbol_length=length,
        coordinate_span=dna.coordinate_span,
        ambiguity=ambiguity,
        quality=quality,
        issues=issue_tuple,
    )


def _validate_collection(
    records: DNASet | Iterable[DNARecord],
    *,
    config: DatasetValidationConfig,
) -> DatasetValidationReport:
    """Validate records once, including stable duplicate-ID detection."""
    resolved = config
    seen_indices: dict[str, list[int]] = defaultdict(list)
    reports: list[ValidationReport] = []
    invalid_indices: list[int] = []
    count = 0

    for index, record in enumerate(records):
        if not isinstance(record, DNARecord):
            raise TypeError(
                f"records item {index} has type {type(record).__name__}; expected DNARecord."
            )
        seen_indices[record.id].append(index)
        report = _validate_one(record, config=resolved.record)
        if not report.is_valid:
            invalid_indices.append(index)
        if resolved.collect_record_reports:
            reports.append(report)
        count += 1

    duplicates = tuple(
        DuplicateID(id=record_id, indices=tuple(indices))
        for record_id, indices in seen_indices.items()
        if len(indices) > 1
    )
    issues: list[Issue] = []
    if count == 0:
        issues.append(
            _issue(
                "STD_EMPTY_DATASET",
                IssueSeverity.ERROR,
                "The dataset contains no records.",
                {},
            )
        )
    if duplicates:
        severity = IssueSeverity.ERROR if resolved.require_unique_ids else IssueSeverity.WARNING
        issues.append(
            _issue(
                "STD_DUPLICATE_RECORD_ID",
                severity,
                f"Found {len(duplicates)} duplicated record identifier(s).",
                {
                    "duplicates": [
                        {"id": item.id, "indices": list(item.indices)} for item in duplicates
                    ]
                },
            )
        )
    if invalid_indices:
        issues.append(
            _issue(
                "STD_INVALID_RECORDS",
                IssueSeverity.ERROR,
                f"Found {len(invalid_indices)} record(s) that failed validation.",
                {"indices": invalid_indices},
            )
        )

    issue_tuple = tuple(issues)
    return DatasetValidationReport(
        record_count=count,
        config=resolved,
        algorithm_version="std-validate-set-v2",
        provenance=Provenance(),
        ids_unique=not duplicates,
        duplicate_ids=duplicates,
        is_valid=not invalid_indices and not _has_error(issue_tuple),
        issues=issue_tuple,
        record_reports=tuple(reports) if resolved.collect_record_reports else None,
    )


def validate_set(
    records: DNA | DNASet | Iterable[DNARecord],
    *,
    config: DatasetValidationConfig | None = None,
) -> DatasetValidationReport:
    """Backward-compatible alias for ``validate(records, config=config)``."""
    return _validate_collection(records, config=_resolve_dataset_config(config))


def _validate_metadata(
    record: DNARecord,
    config: ValidationConfig,
    issues: list[Issue],
) -> None:
    missing = tuple(
        field
        for field in config.required_metadata_fields
        if field not in record.metadata or record.metadata[field] is None
    )
    if missing:
        issues.append(
            _issue(
                "STD_REQUIRED_METADATA_MISSING",
                IssueSeverity.ERROR,
                f"Record {record.id!r} is missing required metadata.",
                {"fields": list(missing)},
            )
        )


def _validate_letter_annotations(
    record: DNARecord,
    config: ValidationConfig,
    issues: list[Issue],
) -> QualitySummary | None:
    missing = tuple(
        field
        for field in config.required_letter_annotations
        if field not in record.letter_annotations
    )
    if missing:
        issues.append(
            _issue(
                "STD_REQUIRED_LETTER_ANNOTATION_MISSING",
                IssueSeverity.ERROR,
                f"Record {record.id!r} is missing required letter annotations.",
                {"fields": list(missing)},
            )
        )

    if not config.check_phred_quality or "phred_quality" not in record.letter_annotations:
        return None

    raw_values = record.letter_annotations["phred_quality"]
    if len(raw_values) != record.sequence.symbol_length:
        issues.append(
            _issue(
                "STD_QUALITY_LENGTH_MISMATCH",
                IssueSeverity.ERROR,
                "PHRED quality count does not match sequence symbol length.",
                {
                    "quality_count": len(raw_values),
                    "symbol_length": record.sequence.symbol_length,
                },
            )
        )

    values: list[float] = []
    invalid_indices: list[int] = []
    out_of_range: list[int] = []
    for index, value in enumerate(raw_values):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            invalid_indices.append(index)
            continue
        numeric = float(value)
        if not math.isfinite(numeric):
            invalid_indices.append(index)
            continue
        values.append(numeric)
        if not config.minimum_phred <= numeric <= config.maximum_phred:
            out_of_range.append(index)

    if invalid_indices:
        issues.append(
            _issue(
                "STD_QUALITY_NON_NUMERIC",
                IssueSeverity.ERROR,
                "PHRED quality contains non-finite or non-numeric values.",
                {"indices": invalid_indices},
            )
        )
    if out_of_range:
        issues.append(
            _issue(
                "STD_QUALITY_OUT_OF_RANGE",
                IssueSeverity.ERROR,
                "PHRED quality contains values outside the configured range.",
                {
                    "indices": out_of_range,
                    "minimum": config.minimum_phred,
                    "maximum": config.maximum_phred,
                },
            )
        )

    mean = sum(values) / len(values) if values else None
    if (
        config.minimum_mean_phred is not None
        and mean is not None
        and mean < config.minimum_mean_phred
    ):
        issues.append(
            _issue(
                "STD_QUALITY_MEAN_LOW",
                IssueSeverity.ERROR,
                f"Mean PHRED quality {mean:.6g} is below the configured minimum "
                f"{config.minimum_mean_phred:.6g}.",
                {"mean": mean, "minimum": config.minimum_mean_phred},
            )
        )
    return QualitySummary(
        count=len(raw_values),
        minimum=min(values) if values else None,
        maximum=max(values) if values else None,
        mean=mean,
    )


def _has_error(issues: tuple[Issue, ...]) -> bool:
    return any(issue.severity is IssueSeverity.ERROR for issue in issues)


def _issue(
    code: str,
    severity: IssueSeverity,
    message: str,
    details: dict[str, object],
) -> Issue:
    return Issue(code=code, severity=severity, message=message, details=details)
