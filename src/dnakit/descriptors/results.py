"""Immutable result objects returned by native descriptor functions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

from dnakit.core import Provenance
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.descriptors._shared import DescriptorAmbiguityPolicy
from dnakit.descriptors.schema import DESCRIPTOR_NAMES_V1, DESCRIPTOR_SCHEMA_VERSION
from dnakit.exceptions import ConfigurationError


@dataclass(frozen=True)
class DescriptorResult:
    """Common auditable context shared by all descriptor result types."""

    name: str
    method: str
    sequence_id: str | None
    ambiguity_policy: DescriptorAmbiguityPolicy | None
    cross_gaps: bool
    gap_count: int
    unknown_gap_count: int

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class AllDescriptorsResult:
    """One ordered 240-field vector plus explicit unavailable-value reasons."""

    schema_version: str
    sequence_id: str | None
    values: FrozenDict
    unavailable_reasons: FrozenDict
    conditions: FrozenDict
    provenance: Provenance

    def __post_init__(self) -> None:
        if self.schema_version != DESCRIPTOR_SCHEMA_VERSION:
            raise ConfigurationError("Unknown all-descriptor schema version.")
        if self.sequence_id is not None and (
            not isinstance(self.sequence_id, str) or not self.sequence_id.strip()
        ):
            raise ConfigurationError("AllDescriptorsResult sequence_id must be non-empty or None.")
        if tuple(self.values) != DESCRIPTOR_NAMES_V1:
            raise ConfigurationError(
                "AllDescriptorsResult values must exactly match descriptor_schema_v1 order."
            )
        unavailable = {name for name, value in self.values.items() if value is None}
        if set(self.unavailable_reasons) != unavailable or any(
            not isinstance(reason, str) or not reason.strip()
            for reason in self.unavailable_reasons.values()
        ):
            raise ConfigurationError(
                "Every unavailable descriptor must have exactly one non-empty reason."
            )
        if not isinstance(self.conditions, FrozenDict):
            raise ConfigurationError("AllDescriptorsResult conditions must be FrozenDict.")
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError("AllDescriptorsResult provenance must be Provenance.")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class LengthResult(DescriptorResult):
    """Explicit symbol, coordinate, ambiguity, and gap lengths."""

    symbol_length: int
    coordinate_span: int | None
    canonical_base_count: int
    ambiguity_length: int
    known_gap_length: int


@dataclass(frozen=True)
class CompositionResult(DescriptorResult):
    """Canonical base counts and fractions with an explicit denominator."""

    counts: FrozenDict
    fractions: FrozenDict
    denominator: int
    ignored_ambiguity_count: int


@dataclass(frozen=True)
class ContentResult(DescriptorResult):
    """GC and AT composition calculated over canonical bases."""

    gc_count: int
    at_count: int
    gc_fraction: float | None
    at_fraction: float | None
    denominator: int
    ignored_ambiguity_count: int


@dataclass(frozen=True)
class SkewResult(DescriptorResult):
    """GC=(G-C)/(G+C) and AT=(A-T)/(A+T), or None at zero denominator."""

    gc_skew: float | None
    at_skew: float | None
    gc_denominator: int
    at_denominator: int
    ignored_ambiguity_count: int


@dataclass(frozen=True)
class CpGResult(DescriptorResult):
    """CpG count, adjacency density, and observed/expected ratio."""

    cpg_count: int
    density: float | None
    observed_expected: float | None
    adjacent_pair_denominator: int
    expected_length_denominator: int
    c_count: int
    g_count: int
    ignored_ambiguity_count: int
    density_formula: str
    observed_expected_formula: str


@dataclass(frozen=True)
class KmerResult(DescriptorResult):
    """Sparse k-mer counts, frequencies, and presence values."""

    k: int
    overlapping: bool
    canonical: bool
    counts: FrozenDict
    frequencies: FrozenDict
    presence: tuple[str, ...]
    denominator: int
    ignored_ambiguity_count: int


@dataclass(frozen=True)
class EntropyResult(DescriptorResult):
    """Shannon entropy of canonical bases or valid canonical k-mers."""

    entropy: float
    unit: str
    k: int
    log_base: float
    observation_count: int
    category_count: int
    ignored_ambiguity_count: int


@dataclass(frozen=True)
class HomopolymerRun:
    """One canonical-base run in zero-based, half-open coordinates."""

    base: str
    length: int
    symbol_start: int
    symbol_end: int
    coordinate_start: int | None
    coordinate_end: int | None
    crossed_gap_count: int
    crossed_unknown_gap: bool


@dataclass(frozen=True)
class HomopolymerResult(DescriptorResult):
    """Longest runs and all runs meeting the requested minimum length."""

    min_run_length: int
    longest_length: int
    longest_by_base: FrozenDict
    runs: tuple[HomopolymerRun, ...]
    ignored_ambiguity_count: int


@dataclass(frozen=True)
class WindowResult:
    """Descriptors for one zero-based, half-open symbol window."""

    symbol_start: int
    symbol_end: int
    coordinate_start: int | None
    coordinate_end: int | None
    is_partial: bool
    crossed_gap_count: int
    crossed_unknown_gap: bool
    values: FrozenDict


@dataclass(frozen=True)
class WindowDescriptorResult(DescriptorResult):
    """Ordered position-by-feature rows for a sliding-window calculation."""

    window_size: int
    step: int
    include_partial: bool
    entropy_log_base: float
    descriptors: tuple[str, ...]
    windows: tuple[WindowResult, ...]


@dataclass(frozen=True)
class CodonResult(DescriptorResult):
    """In-frame codon counts, frequencies, and start/stop densities."""

    frame: int
    genetic_code: int
    counts: FrozenDict
    frequencies: FrozenDict
    codon_count: int
    start_count: int
    stop_count: int
    start_density: float | None
    stop_density: float | None
    incomplete_base_count: int
    ignored_ambiguity_codon_count: int
    gap_interrupted_codon_count: int
    phase_coordinate_system: str = "sequence-coordinate"
    phase_unresolved_after_gap: bool = False
    unresolved_downstream_base_count: int = 0


@dataclass(frozen=True)
class ComplexityResult(DescriptorResult):
    """Linguistic vocabulary complexity with per-word-size audit."""

    score: float
    max_word_size: int
    by_k: FrozenDict
    observed_by_k: FrozenDict
    possible_by_k: FrozenDict
    observation_count: int
    max_observations: int
    formula: str


@dataclass(frozen=True)
class RepeatRun:
    unit: str
    unit_length: int
    repeat_count: int
    symbol_start: int
    symbol_end: int


@dataclass(frozen=True)
class ExactRepeatResult(DescriptorResult):
    """Union coverage and runs from an exact tandem-repeat scan."""

    repeat_fraction: float
    repeated_base_count: int
    denominator: int
    min_unit_length: int
    max_unit_length: int
    min_repeats: int
    runs: tuple[RepeatRun, ...]
    repeat_count_by_unit: FrozenDict
    comparisons: int
    max_comparisons: int


__all__ = [
    "AllDescriptorsResult",
    "CodonResult",
    "ComplexityResult",
    "CompositionResult",
    "ContentResult",
    "CpGResult",
    "DescriptorResult",
    "EntropyResult",
    "ExactRepeatResult",
    "HomopolymerResult",
    "HomopolymerRun",
    "KmerResult",
    "LengthResult",
    "RepeatRun",
    "SkewResult",
    "WindowDescriptorResult",
    "WindowResult",
]
