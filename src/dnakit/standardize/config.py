"""Configuration objects for deterministic sequence standardization and validation."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from itertools import islice
from types import MappingProxyType

from dnakit.core.enums import DNAAlphabet
from dnakit.exceptions import ConfigurationError

_MAX_REQUIRED_FIELDS = 10_000


class UPolicy(str, Enum):
    """How :func:`normalize` handles uracil in a DNA input.

    ``WARN`` and ``KEEP`` retain U for audit purposes, but U is not a valid
    DNA alphabet symbol.  Consequently those policies produce no
    :class:`~dnakit.core.sequence.DNASequence`; ``WARN`` additionally emits a
    warning.  ``DELETE`` removes U and ``REPLACE`` converts U to valid DNA.
    """

    DELETE = "delete"
    ERROR = "error"
    WARN = "warn"
    REPLACE = "replace"
    KEEP = "keep"


class AmbiguityPolicy(str, Enum):
    """How IUPAC ambiguity symbols are represented during normalization.

    ``IGNORE`` means that the symbols are retained and excluded from any
    canonical-base-only calculation.  It does not delete sequence positions.
    ``PROBABILITY`` also retains the symbols and records a deterministic base
    probability distribution in the normalization result.  ``DELETE`` removes
    ambiguity symbols while retaining an audit change for every source position.
    """

    DELETE = "delete"
    ERROR = "error"
    IGNORE = "ignore"
    MASK = "mask"
    PROBABILITY = "probability"


@dataclass(frozen=True, slots=True)
class NormalizationConfig:
    """Control a normalization run and retain an audit entry for every edit.

    The three simple keep flags define the default character policy.  Explicit
    ``u_policy`` or ``ambiguity_policy`` values remain available for advanced
    behavior and take precedence over their corresponding keep flag.
    """

    alphabet: DNAAlphabet | None = None
    uppercase: bool = True
    remove_whitespace: bool = True
    remove_invisible: bool = True
    removable_separators: tuple[str, ...] = ()
    keep_ambiguous: bool = True
    keep_u: bool = False
    keep_other: bool = False
    u_policy: UPolicy | None = None
    ambiguity_policy: AmbiguityPolicy | None = None
    ambiguity_mask: str = "N"
    base_priors: Mapping[str, float] = field(
        default_factory=lambda: MappingProxyType({"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
    )
    allow_gaps: bool = True
    raise_on_error: bool = False
    operator: str | None = None
    run_id: str | None = None

    def __post_init__(self) -> None:
        for field_name in (
            "uppercase",
            "remove_whitespace",
            "remove_invisible",
            "keep_ambiguous",
            "keep_u",
            "keep_other",
            "allow_gaps",
            "raise_on_error",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ConfigurationError(
                    f"{field_name} must be a boolean.",
                    code="INVALID_NORMALIZATION_FLAG",
                    context={"field": field_name},
                )
        for field_name in ("operator", "run_id"):
            value = getattr(self, field_name)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ConfigurationError(
                    f"{field_name} must be a non-empty string or None.",
                    code="INVALID_NORMALIZATION_AUDIT_FIELD",
                    context={"field": field_name},
                )
        raw_separators: object = self.removable_separators
        if isinstance(raw_separators, (str, bytes)):
            raise ConfigurationError(
                "removable_separators must be an iterable of one-character strings.",
                code="INVALID_NORMALIZATION_SEPARATOR",
            )
        separators = tuple(self.removable_separators)
        if any(not isinstance(separator, str) or len(separator) != 1 for separator in separators):
            raise ConfigurationError(
                "Each removable separator must contain exactly one character.",
                code="INVALID_NORMALIZATION_SEPARATOR",
            )
        if len(set(separators)) != len(separators):
            raise ConfigurationError(
                "removable_separators must not contain duplicates.",
                code="DUPLICATE_NORMALIZATION_SEPARATOR",
            )
        if any(separator.isspace() for separator in separators):
            raise ConfigurationError(
                "Whitespace must be controlled with remove_whitespace.",
                code="WHITESPACE_SEPARATOR_CONFLICT",
            )

        if not isinstance(self.ambiguity_mask, str):
            raise ConfigurationError(
                "ambiguity_mask must be a string.", code="INVALID_AMBIGUITY_MASK"
            )
        mask = self.ambiguity_mask.upper()
        if len(mask) != 1 or mask not in "RYSWKMBDHVN":
            raise ConfigurationError(
                "ambiguity_mask must be one IUPAC ambiguity symbol.",
                code="INVALID_AMBIGUITY_MASK",
            )

        try:
            raw_priors = dict(self.base_priors.items())
        except (AttributeError, TypeError, ValueError) as exc:
            raise ConfigurationError(
                "base_priors must map A, C, G and T to finite numbers.",
                code="INVALID_BASE_PRIORS",
            ) from exc
        if any(
            not isinstance(base, str) or base not in {"A", "C", "G", "T"} for base in raw_priors
        ):
            raise ConfigurationError(
                "base_priors keys must be exactly uppercase A, C, G and T.",
                code="INVALID_BASE_PRIORS",
            )
        if any(
            isinstance(value, bool) or not isinstance(value, (int, float))
            for value in raw_priors.values()
        ):
            raise ConfigurationError(
                "base_priors must map A, C, G and T to finite numbers.",
                code="INVALID_BASE_PRIORS",
            )
        priors = {base: float(value) for base, value in raw_priors.items()}
        if set(priors) != {"A", "C", "G", "T"}:
            raise ConfigurationError(
                "base_priors must define exactly A, C, G and T.",
                code="INVALID_BASE_PRIORS",
            )
        if (
            any(not math.isfinite(value) or value < 0.0 for value in priors.values())
            or sum(priors.values()) <= 0.0
        ):
            raise ConfigurationError(
                "base_priors must be non-negative and have a positive sum.",
                code="INVALID_BASE_PRIORS",
            )

        try:
            alphabet = (
                None
                if self.alphabet is None
                else self.alphabet
                if isinstance(self.alphabet, DNAAlphabet)
                else DNAAlphabet(self.alphabet)
            )
            u_policy = (
                None
                if self.u_policy is None
                else self.u_policy
                if isinstance(self.u_policy, UPolicy)
                else UPolicy(self.u_policy)
            )
            ambiguity_policy = (
                None
                if self.ambiguity_policy is None
                else self.ambiguity_policy
                if isinstance(self.ambiguity_policy, AmbiguityPolicy)
                else AmbiguityPolicy(self.ambiguity_policy)
            )
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Unknown normalization policy or DNA alphabet.",
                code="INVALID_NORMALIZATION_POLICY",
            ) from exc

        object.__setattr__(self, "alphabet", alphabet)
        object.__setattr__(self, "u_policy", u_policy)
        object.__setattr__(self, "ambiguity_policy", ambiguity_policy)
        object.__setattr__(self, "removable_separators", separators)
        object.__setattr__(self, "ambiguity_mask", mask)
        object.__setattr__(self, "base_priors", MappingProxyType(priors))

    @property
    def effective_u_policy(self) -> UPolicy:
        """Return the advanced U policy or the simple keep/delete policy."""

        if self.u_policy is not None:
            return self.u_policy
        return UPolicy.KEEP if self.keep_u else UPolicy.DELETE

    @property
    def effective_ambiguity_policy(self) -> AmbiguityPolicy:
        """Return the advanced ambiguity policy or the simple keep/delete policy."""

        if self.ambiguity_policy is not None:
            return self.ambiguity_policy
        return AmbiguityPolicy.IGNORE if self.keep_ambiguous else AmbiguityPolicy.DELETE


@dataclass(frozen=True, slots=True)
class ValidationConfig:
    """Thresholds and policies used by :func:`validate`.

    ``sequence_length`` optionally requires an exact ``symbol_length``.  It is
    independent of the inclusive ``min_length`` and ``max_length`` bounds.
    """

    alphabet: DNAAlphabet | None = None
    allow_empty: bool = False
    allow_gaps: bool = True
    allow_unknown_gap_length: bool = True
    min_length: int | None = 1
    max_length: int | None = None
    max_ambiguity_fraction: float | None = 1.0
    ambiguity_denominator_includes_gap: bool = False
    required_metadata_fields: tuple[str, ...] = ()
    required_letter_annotations: tuple[str, ...] = ()
    check_phred_quality: bool = True
    minimum_phred: float = 0.0
    maximum_phred: float = 93.0
    minimum_mean_phred: float | None = None
    sequence_length: int | None = None

    def __post_init__(self) -> None:
        try:
            alphabet = (
                None
                if self.alphabet is None
                else self.alphabet
                if isinstance(self.alphabet, DNAAlphabet)
                else DNAAlphabet(self.alphabet)
            )
        except (TypeError, ValueError) as exc:
            raise ConfigurationError(
                "Unknown validation DNA alphabet.", code="INVALID_VALIDATION_ALPHABET"
            ) from exc
        for field_name in (
            "allow_empty",
            "allow_gaps",
            "allow_unknown_gap_length",
            "ambiguity_denominator_includes_gap",
            "check_phred_quality",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ConfigurationError(
                    f"{field_name} must be a boolean.",
                    code="INVALID_VALIDATION_FLAG",
                    context={"field": field_name},
                )
        for field_name, threshold in (
            ("min_length", self.min_length),
            ("max_length", self.max_length),
            ("sequence_length", self.sequence_length),
        ):
            if threshold is not None and (
                isinstance(threshold, bool) or not isinstance(threshold, int) or threshold < 0
            ):
                raise ConfigurationError(
                    f"{field_name} must be a non-negative integer or None.",
                    code="INVALID_LENGTH_THRESHOLD",
                    context={"field": field_name},
                )
        if (
            self.min_length is not None
            and self.max_length is not None
            and self.min_length > self.max_length
        ):
            raise ConfigurationError(
                "min_length must not exceed max_length.", code="INVALID_LENGTH_THRESHOLD"
            )
        ambiguity_fraction = self.max_ambiguity_fraction
        if ambiguity_fraction is not None and (
            isinstance(ambiguity_fraction, bool)
            or not isinstance(ambiguity_fraction, (int, float))
            or not math.isfinite(ambiguity_fraction)
            or not 0.0 <= ambiguity_fraction <= 1.0
        ):
            raise ConfigurationError(
                "max_ambiguity_fraction must be a finite number between 0 and 1.",
                code="INVALID_AMBIGUITY_THRESHOLD",
            )
        quality_values: dict[str, float | None] = {}
        quality_thresholds: tuple[tuple[str, float | None, bool], ...] = (
            ("minimum_phred", self.minimum_phred, False),
            ("maximum_phred", self.maximum_phred, False),
            ("minimum_mean_phred", self.minimum_mean_phred, True),
        )
        for field_name, quality_threshold, optional in quality_thresholds:
            if quality_threshold is None and optional:
                quality_values[field_name] = None
                continue
            if (
                isinstance(quality_threshold, bool)
                or not isinstance(quality_threshold, (int, float))
                or not math.isfinite(quality_threshold)
            ):
                raise ConfigurationError(
                    f"{field_name} must be a finite number" + (" or None." if optional else "."),
                    code="INVALID_QUALITY_THRESHOLD",
                    context={"field": field_name},
                )
            quality_values[field_name] = float(quality_threshold)
        minimum_phred = quality_values["minimum_phred"]
        maximum_phred = quality_values["maximum_phred"]
        assert minimum_phred is not None and maximum_phred is not None
        if minimum_phred > maximum_phred:
            raise ConfigurationError(
                "minimum_phred must not exceed maximum_phred.",
                code="INVALID_QUALITY_THRESHOLD",
            )

        resolved_required_fields: dict[str, tuple[str, ...]] = {}
        for field_name, values in (
            ("required_metadata_fields", self.required_metadata_fields),
            ("required_letter_annotations", self.required_letter_annotations),
        ):
            raw_values: object = values
            if isinstance(raw_values, (str, bytes)) or not isinstance(raw_values, Iterable):
                raise ConfigurationError(
                    f"{field_name} must be an iterable of non-empty strings.",
                    code="INVALID_REQUIRED_FIELD",
                    context={"field": field_name},
                )
            resolved = tuple(islice(iter(values), _MAX_REQUIRED_FIELDS + 1))
            if len(resolved) > _MAX_REQUIRED_FIELDS:
                raise ConfigurationError(
                    f"{field_name} exceeds the supported field limit.",
                    code="REQUIRED_FIELD_LIMIT",
                    context={"field": field_name, "max_fields": _MAX_REQUIRED_FIELDS},
                )
            if any(not isinstance(value, str) or not value.strip() for value in resolved):
                raise ConfigurationError(
                    f"{field_name} must contain only non-empty strings.",
                    code="INVALID_REQUIRED_FIELD",
                    context={"field": field_name},
                )
            if len(set(resolved)) != len(resolved):
                raise ConfigurationError(
                    f"{field_name} must not contain duplicates.",
                    code="DUPLICATE_REQUIRED_FIELD",
                    context={"field": field_name},
                )
            resolved_required_fields[field_name] = resolved

        object.__setattr__(self, "alphabet", alphabet)
        object.__setattr__(
            self,
            "max_ambiguity_fraction",
            None if ambiguity_fraction is None else float(ambiguity_fraction),
        )
        object.__setattr__(self, "minimum_phred", minimum_phred)
        object.__setattr__(self, "maximum_phred", maximum_phred)
        object.__setattr__(self, "minimum_mean_phred", quality_values["minimum_mean_phred"])
        object.__setattr__(
            self, "required_metadata_fields", resolved_required_fields["required_metadata_fields"]
        )
        object.__setattr__(
            self,
            "required_letter_annotations",
            resolved_required_fields["required_letter_annotations"],
        )


@dataclass(frozen=True, slots=True)
class DatasetValidationConfig:
    """Collection-level validation rules."""

    record: ValidationConfig = field(default_factory=ValidationConfig)
    require_unique_ids: bool = True
    collect_record_reports: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.record, ValidationConfig):
            raise ConfigurationError(
                "record must be a ValidationConfig object.",
                code="INVALID_DATASET_VALIDATION_CONFIG",
            )
        for field_name in ("require_unique_ids", "collect_record_reports"):
            if not isinstance(getattr(self, field_name), bool):
                raise ConfigurationError(
                    f"{field_name} must be a boolean.",
                    code="INVALID_VALIDATION_FLAG",
                    context={"field": field_name},
                )
