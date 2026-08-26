"""Immutable structured results returned by DNAKit fingerprint functions."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, TypeAlias, cast

from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.exceptions import ConfigurationError
from dnakit.fingerprints._shared import (
    FingerprintAmbiguityPolicy,
    FingerprintRepresentation,
    GapEncodingPolicy,
    IntegerAmbiguityPolicy,
    KmerFingerprintMode,
    OneHotAmbiguityPolicy,
)

Numeric: TypeAlias = int | float
FingerprintValues: TypeAlias = tuple[Numeric, ...] | FrozenDict
BitFingerprintValues: TypeAlias = tuple[int, ...] | FrozenDict


def _validate_common(result: RepresentationResult) -> None:
    for name in ("name", "method", "schema_version"):
        value = getattr(result, name)
        if not isinstance(value, str) or not value.strip():
            raise ConfigurationError(f"Fingerprint result {name} must be non-empty.")
    if result.sequence_id is not None and (
        not isinstance(result.sequence_id, str) or not result.sequence_id.strip()
    ):
        raise ConfigurationError("Fingerprint result sequence_id must be non-empty or None.")
    for name in ("symbol_length", "gap_count", "unknown_gap_count"):
        value = getattr(result, name)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ConfigurationError(f"Fingerprint result {name} must be non-negative.")
    if result.unknown_gap_count > result.gap_count:
        raise ConfigurationError("unknown_gap_count cannot exceed gap_count.")


def _is_non_negative_number(value: object) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(value)
        and value >= 0
    )


@dataclass(frozen=True)
class RepresentationResult:
    """Common auditable context for encodings and fingerprints."""

    name: str
    method: str
    schema_version: str
    sequence_id: str | None
    symbol_length: int
    gap_count: int
    unknown_gap_count: int

    def __post_init__(self) -> None:
        _validate_common(self)

    def to_dict(self) -> dict[str, Any]:
        """Return a plain JSON-compatible representation."""

        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class IntegerEncodingResult(RepresentationResult):
    """Position-preserving integer encoding with an explicit codebook."""

    values: tuple[int, ...]
    codebook: FrozenDict
    ambiguity_policy: IntegerAmbiguityPolicy
    gap_policy: GapEncodingPolicy
    output_length: int
    encoded_ambiguity_count: int
    expanded_gap_length: int
    omitted_gap_count: int
    max_output_length: int

    def __post_init__(self) -> None:
        super().__post_init__()
        if not isinstance(self.codebook, FrozenDict):
            raise ConfigurationError("Integer encoding codebook must be a FrozenDict.")
        if not isinstance(self.ambiguity_policy, IntegerAmbiguityPolicy) or not isinstance(
            self.gap_policy, GapEncodingPolicy
        ):
            raise ConfigurationError("Integer encoding policies must use their public enums.")
        code_values = tuple(self.codebook.values())
        if any(isinstance(value, bool) or not isinstance(value, int) for value in code_values):
            raise ConfigurationError("Integer encoding codebook values must all be integers.")
        if len(set(cast(tuple[int, ...], code_values))) != len(code_values):
            raise ConfigurationError("Integer encoding codebook values must be unique.")
        if any(isinstance(value, bool) or not isinstance(value, int) for value in self.values):
            raise ConfigurationError("Integer encoding values must all be integers.")
        if isinstance(self.output_length, bool) or not isinstance(self.output_length, int):
            raise ConfigurationError("Integer encoding output_length must be an integer.")
        if self.output_length != len(self.values):
            raise ConfigurationError("Integer encoding output_length does not match its values.")
        for name in ("encoded_ambiguity_count", "expanded_gap_length", "omitted_gap_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"Integer encoding {name} must be non-negative.")
        if (
            isinstance(self.max_output_length, bool)
            or not isinstance(self.max_output_length, int)
            or self.max_output_length <= 0
            or self.output_length > self.max_output_length
        ):
            raise ConfigurationError("Integer encoding max_output_length is inconsistent.")
        if self.output_length != self.symbol_length + self.expanded_gap_length:
            raise ConfigurationError(
                "Integer encoding output length must equal symbols plus expanded Gap positions."
            )
        if self.encoded_ambiguity_count > self.symbol_length:
            raise ConfigurationError("encoded_ambiguity_count cannot exceed symbol_length.")
        if self.omitted_gap_count > self.gap_count:
            raise ConfigurationError("omitted_gap_count cannot exceed gap_count.")

    @property
    def dimension(self) -> int:
        """Number of distinct codes in the active immutable codebook."""

        return len(self.codebook)


@dataclass(frozen=True)
class OneHotEncodingResult(RepresentationResult):
    """A/C/G/T one-hot rows with explicit ambiguity and Gap policies."""

    values: tuple[tuple[float, ...], ...]
    feature_names: tuple[str, ...]
    ambiguity_policy: OneHotAmbiguityPolicy
    gap_policy: GapEncodingPolicy
    output_length: int
    encoded_ambiguity_count: int
    expanded_gap_length: int
    omitted_gap_count: int
    max_output_length: int

    def __post_init__(self) -> None:
        super().__post_init__()
        if not isinstance(self.ambiguity_policy, OneHotAmbiguityPolicy) or not isinstance(
            self.gap_policy, GapEncodingPolicy
        ):
            raise ConfigurationError("One-hot encoding policies must use their public enums.")
        if isinstance(self.output_length, bool) or not isinstance(self.output_length, int):
            raise ConfigurationError("One-hot output_length must be an integer.")
        if self.output_length != len(self.values):
            raise ConfigurationError("One-hot output_length does not match its rows.")
        if len(self.feature_names) != 4 or set(self.feature_names) != set("ACGT"):
            raise ConfigurationError("One-hot feature_names must be a permutation of A/C/G/T.")
        for row in self.values:
            if len(row) != self.dimension or any(
                not _is_non_negative_number(value) or value > 1 for value in row
            ):
                raise ConfigurationError("One-hot rows must contain four finite values in [0, 1].")
            if not math.isclose(sum(row), 0.0) and not math.isclose(sum(row), 1.0):
                raise ConfigurationError("Each one-hot row must sum to zero or one.")
        for name in ("encoded_ambiguity_count", "expanded_gap_length", "omitted_gap_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"One-hot encoding {name} must be non-negative.")
        if (
            isinstance(self.max_output_length, bool)
            or not isinstance(self.max_output_length, int)
            or self.max_output_length <= 0
            or self.output_length > self.max_output_length
        ):
            raise ConfigurationError("One-hot encoding max_output_length is inconsistent.")
        if self.output_length != self.symbol_length + self.expanded_gap_length:
            raise ConfigurationError(
                "One-hot output length must equal symbols plus expanded Gap positions."
            )
        if self.encoded_ambiguity_count > self.symbol_length:
            raise ConfigurationError("encoded_ambiguity_count cannot exceed symbol_length.")
        if self.omitted_gap_count > self.gap_count:
            raise ConfigurationError("omitted_gap_count cannot exceed gap_count.")

    @property
    def dimension(self) -> int:
        """Number of columns in each encoded row."""

        return len(self.feature_names)


@dataclass(frozen=True)
class FingerprintResult(RepresentationResult):
    """Fixed-schema exact k-mer fingerprint in dense or sparse form."""

    k: int
    canonical: bool
    mode: KmerFingerprintMode
    representation: FingerprintRepresentation
    ambiguity_policy: FingerprintAmbiguityPolicy
    overlapping: bool
    cross_gaps: bool
    feature_names: tuple[str, ...]
    values: FingerprintValues
    observation_count: int
    ignored_ambiguity_count: int
    max_dimension: int

    def __post_init__(self) -> None:
        super().__post_init__()
        if not isinstance(self.mode, KmerFingerprintMode):
            raise ConfigurationError("Fingerprint mode must be KmerFingerprintMode.")
        if not isinstance(self.representation, FingerprintRepresentation):
            raise ConfigurationError(
                "Fingerprint representation must be FingerprintRepresentation."
            )
        if not isinstance(self.ambiguity_policy, FingerprintAmbiguityPolicy):
            raise ConfigurationError(
                "Fingerprint ambiguity_policy must be FingerprintAmbiguityPolicy."
            )
        if isinstance(self.k, bool) or not isinstance(self.k, int) or self.k <= 0:
            raise ConfigurationError("Fingerprint k must be a positive integer.")
        if (
            isinstance(self.max_dimension, bool)
            or not isinstance(self.max_dimension, int)
            or self.max_dimension <= 0
        ):
            raise ConfigurationError("Fingerprint max_dimension must be a positive integer.")
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("Fingerprint canonical must be a boolean.")
        if not isinstance(self.overlapping, bool) or not isinstance(self.cross_gaps, bool):
            raise ConfigurationError("Fingerprint traversal flags must be booleans.")
        if len(set(self.feature_names)) != len(self.feature_names) or any(
            not isinstance(name, str) or not name for name in self.feature_names
        ):
            raise ConfigurationError("Fingerprint feature_names must be unique non-empty strings.")
        if self.dimension > self.max_dimension:
            raise ConfigurationError("Fingerprint dimension exceeds max_dimension.")
        if isinstance(self.values, tuple):
            if self.representation is not FingerprintRepresentation.DENSE:
                raise ConfigurationError("Tuple fingerprint values require dense representation.")
            if len(self.values) != self.dimension:
                raise ConfigurationError("Dense fingerprint length does not match its schema.")
            numeric_values: tuple[object, ...] = self.values
        elif isinstance(self.values, FrozenDict):
            if self.representation is not FingerprintRepresentation.SPARSE:
                raise ConfigurationError(
                    "Mapping fingerprint values require sparse representation."
                )
            unknown_features = set(self.values) - set(self.feature_names)
            if unknown_features:
                raise ConfigurationError(
                    "Sparse fingerprint contains features outside its schema.",
                    context={"features": sorted(unknown_features)},
                )
            numeric_values = tuple(self.values.values())
        else:
            raise ConfigurationError("Fingerprint values must be a tuple or FrozenDict.")
        if any(not _is_non_negative_number(value) for value in numeric_values):
            raise ConfigurationError("Fingerprint values must be finite non-negative numbers.")
        if isinstance(self.values, FrozenDict) and any(value == 0 for value in numeric_values):
            raise ConfigurationError("Sparse fingerprint values must omit zero entries.")
        for name in ("observation_count", "ignored_ambiguity_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"Fingerprint {name} must be non-negative.")
        if self.mode is KmerFingerprintMode.COUNT and any(
            not isinstance(value, int) or isinstance(value, bool) for value in numeric_values
        ):
            raise ConfigurationError("Count fingerprint values must be integers.")
        if (
            self.mode is KmerFingerprintMode.COUNT
            and sum(cast(int, value) for value in numeric_values) != self.observation_count
        ):
            raise ConfigurationError("Count fingerprint values must sum to observation_count.")
        if self.mode is KmerFingerprintMode.BINARY and any(
            not isinstance(value, int) or isinstance(value, bool) or value not in (0, 1)
            for value in numeric_values
        ):
            raise ConfigurationError("Binary fingerprint values must be zero or one.")
        if self.mode is KmerFingerprintMode.FREQUENCY and any(
            not isinstance(value, (int, float)) or value > 1 for value in numeric_values
        ):
            raise ConfigurationError("Frequency fingerprint values must be in [0, 1].")
        if self.mode is KmerFingerprintMode.FREQUENCY:
            expected_total = 1.0 if self.observation_count else 0.0
            actual_total = sum(float(cast(Numeric, value)) for value in numeric_values)
            if not math.isclose(actual_total, expected_total):
                raise ConfigurationError(
                    "Frequency fingerprint values must sum to one, or zero with no observations."
                )

    @property
    def dimension(self) -> int:
        """Number of features in the immutable schema."""

        return len(self.feature_names)

    def dense_values(self) -> tuple[Numeric, ...]:
        """Return values ordered exactly like :attr:`feature_names`."""

        if isinstance(self.values, tuple):
            return self.values
        return tuple(cast(Numeric, self.values.get(name, 0)) for name in self.feature_names)

    def sparse_values(self) -> FrozenDict:
        """Return only non-zero values keyed by schema feature name."""

        if isinstance(self.values, FrozenDict):
            return self.values
        return FrozenDict(
            {
                name: value
                for name, value in zip(self.feature_names, self.values, strict=True)
                if value != 0
            }
        )


KmerFingerprintResult = FingerprintResult


@dataclass(frozen=True)
class BitFingerprintResult(RepresentationResult):
    """Fixed-schema binary fingerprint in dense or sparse form."""

    representation: FingerprintRepresentation
    feature_names: tuple[str, ...]
    values: BitFingerprintValues
    parameters: FrozenDict
    observation_count: int
    ignored_ambiguity_count: int
    max_dimension: int

    def __post_init__(self) -> None:
        super().__post_init__()
        if not isinstance(self.representation, FingerprintRepresentation):
            raise ConfigurationError(
                "Bit fingerprint representation must be FingerprintRepresentation."
            )
        if not self.feature_names or len(set(self.feature_names)) != len(self.feature_names):
            raise ConfigurationError("Bit fingerprint feature_names must be unique and non-empty.")
        if any(not isinstance(name, str) or not name for name in self.feature_names):
            raise ConfigurationError(
                "Bit fingerprint feature_names must be unique non-empty strings."
            )
        if (
            isinstance(self.max_dimension, bool)
            or not isinstance(self.max_dimension, int)
            or self.max_dimension <= 0
            or self.dimension > self.max_dimension
        ):
            raise ConfigurationError("Bit fingerprint max_dimension is inconsistent.")
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("Bit fingerprint parameters must be a FrozenDict.")
        for name in ("observation_count", "ignored_ambiguity_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"Bit fingerprint {name} must be non-negative.")
        if isinstance(self.values, tuple):
            if self.representation is not FingerprintRepresentation.DENSE:
                raise ConfigurationError(
                    "Tuple bit fingerprint values require dense representation."
                )
            if len(self.values) != self.dimension:
                raise ConfigurationError("Dense bit fingerprint length does not match its schema.")
            numeric_values: tuple[object, ...] = self.values
        elif isinstance(self.values, FrozenDict):
            if self.representation is not FingerprintRepresentation.SPARSE:
                raise ConfigurationError(
                    "Mapping bit fingerprint values require sparse representation."
                )
            unknown_features = set(self.values) - set(self.feature_names)
            if unknown_features:
                raise ConfigurationError(
                    "Sparse bit fingerprint contains features outside its schema.",
                    context={"features": sorted(unknown_features)},
                )
            numeric_values = tuple(self.values.values())
        else:
            raise ConfigurationError("Bit fingerprint values must be a tuple or FrozenDict.")
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value not in (0, 1)
            for value in numeric_values
        ):
            raise ConfigurationError("Bit fingerprint values must be zero or one.")
        if isinstance(self.values, FrozenDict) and any(value == 0 for value in numeric_values):
            raise ConfigurationError("Sparse bit fingerprint values must omit zero entries.")

    @property
    def dimension(self) -> int:
        """Number of bits in the immutable schema."""

        return len(self.feature_names)

    @property
    def set_bit_count(self) -> int:
        """Number of bits set to one."""

        if isinstance(self.values, FrozenDict):
            return len(self.values)
        return sum(self.values)

    def dense_values(self) -> tuple[int, ...]:
        """Return bits ordered exactly like :attr:`feature_names`."""

        if isinstance(self.values, tuple):
            return self.values
        return tuple(cast(int, self.values.get(name, 0)) for name in self.feature_names)

    def sparse_values(self) -> FrozenDict:
        """Return only set bits keyed by schema feature name."""

        if isinstance(self.values, FrozenDict):
            return self.values
        return FrozenDict(
            {
                name: value
                for name, value in zip(self.feature_names, self.values, strict=True)
                if value
            }
        )


@dataclass(frozen=True)
class SketchResult:
    """Versioned set of deterministic 64-bit k-mer hashes."""

    name: str
    method: str
    schema_version: str
    sequence_id: str | None
    symbol_length: int
    k: int
    canonical: bool
    seed: int
    selection: str
    hashes: tuple[int, ...]
    num_hashes: int | None
    scaled: int | None
    threshold: int | None
    observation_count: int
    unique_hash_count: int
    max_hashes: int
    parameters: FrozenDict

    def __post_init__(self) -> None:
        for name in ("name", "method", "schema_version"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"Sketch {name} must be non-empty.")
        if self.sequence_id is not None and (
            not isinstance(self.sequence_id, str) or not self.sequence_id.strip()
        ):
            raise ConfigurationError("Sketch sequence_id must be non-empty or None.")
        for name in ("symbol_length", "observation_count", "unique_hash_count"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ConfigurationError(f"Sketch {name} must be a non-negative integer.")
        for name in ("k", "max_hashes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ConfigurationError(f"Sketch {name} must be a positive integer.")
        if not isinstance(self.canonical, bool):
            raise ConfigurationError("Sketch canonical must be boolean.")
        if (
            isinstance(self.seed, bool)
            or not isinstance(self.seed, int)
            or not 0 <= self.seed < 2**64
        ):
            raise ConfigurationError("Sketch seed must be an unsigned 64-bit integer.")
        if self.selection not in {"bottom_k", "scaled"}:
            raise ConfigurationError("Sketch selection must be bottom_k or scaled.")
        if len(self.hashes) > self.max_hashes or len(set(self.hashes)) != len(self.hashes):
            raise ConfigurationError("Sketch hashes violate uniqueness or size limits.")
        if tuple(sorted(self.hashes)) != self.hashes or any(
            isinstance(item, bool) or not isinstance(item, int) or not 0 <= item < 2**64
            for item in self.hashes
        ):
            raise ConfigurationError("Sketch hashes must be sorted unsigned 64-bit integers.")
        if self.selection == "bottom_k":
            if self.num_hashes is None or self.scaled is not None or self.threshold is not None:
                raise ConfigurationError("Bottom-k sketch selection fields are inconsistent.")
            if (
                isinstance(self.num_hashes, bool)
                or not isinstance(self.num_hashes, int)
                or not 0 < self.num_hashes <= self.max_hashes
                or len(self.hashes) > self.num_hashes
            ):
                raise ConfigurationError("Bottom-k num_hashes is invalid or exceeded.")
        elif self.scaled is None or self.threshold is None or self.num_hashes is not None:
            raise ConfigurationError("Scaled sketch selection fields are inconsistent.")
        else:
            if (
                isinstance(self.scaled, bool)
                or not isinstance(self.scaled, int)
                or not 0 < self.scaled <= 2**64
                or isinstance(self.threshold, bool)
                or not isinstance(self.threshold, int)
                or self.threshold != 2**64 // self.scaled
                or any(item >= self.threshold for item in self.hashes)
            ):
                raise ConfigurationError("Scaled sketch threshold fields are invalid.")
        if self.unique_hash_count < len(self.hashes):
            raise ConfigurationError("unique_hash_count cannot be smaller than retained hashes.")
        if self.unique_hash_count > self.observation_count:
            raise ConfigurationError("unique_hash_count cannot exceed observation_count.")
        expected_observations = max(0, self.symbol_length - self.k + 1)
        if self.observation_count != expected_observations:
            raise ConfigurationError(
                "Sketch observation_count does not match sequence length and k."
            )
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("Sketch parameters must be FrozenDict.")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


__all__ = [
    "BitFingerprintResult",
    "BitFingerprintValues",
    "FingerprintResult",
    "FingerprintValues",
    "IntegerEncodingResult",
    "KmerFingerprintResult",
    "Numeric",
    "OneHotEncodingResult",
    "RepresentationResult",
    "SketchResult",
]
