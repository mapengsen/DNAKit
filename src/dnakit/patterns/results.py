"""Auditable immutable results for sequence-pattern analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Generic, TypeVar, cast

from dnakit.core import (
    CompoundLocation,
    DNASequence,
    Interval,
    Location,
    Provenance,
    Strand,
    UnresolvedLocation,
)
from dnakit.core._json import FrozenDict, to_json_compatible
from dnakit.core.issues import Issue
from dnakit.exceptions import ConfigurationError

HitT = TypeVar("HitT")
_IUPAC = frozenset("ACGTRYSWKMBDHVN")


def _require_non_negative(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConfigurationError(f"{name} must be a non-negative integer.")


def _require_positive(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ConfigurationError(f"{name} must be a positive integer.")


def _require_text(value: str, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty string.")


def _require_location(value: object, name: str) -> None:
    if not isinstance(value, (Interval, CompoundLocation, UnresolvedLocation)):
        raise ConfigurationError(f"{name} must be an internal Location.")


def _require_strand(value: object, *, allow_both: bool = False) -> None:
    accepted = (
        (Strand.FORWARD, Strand.REVERSE, Strand.BOTH)
        if allow_both
        else (
            Strand.FORWARD,
            Strand.REVERSE,
        )
    )
    if value not in accepted:
        raise ConfigurationError("Result strand is outside the supported strand values.")


def _require_iupac(value: str, name: str) -> None:
    _require_text(value, name)
    if set(value) - _IUPAC:
        raise ConfigurationError(f"{name} must contain uppercase DNA IUPAC symbols.")


@dataclass(frozen=True)
class PatternResult(Generic[HitT]):
    """A bounded pattern-analysis result with complete algorithm context."""

    name: str
    method: str
    algorithm_version: str
    sequence_id: str | None
    parameters: FrozenDict
    hits: tuple[HitT, ...]
    inspected_symbol_count: int
    gap_count: int
    unknown_gap_count: int
    max_matches: int
    truncated: bool
    coordinate_system: str
    gap_policy: str
    topology: str
    provenance: Provenance
    issues: tuple[Issue, ...]

    def __post_init__(self) -> None:
        for name in ("name", "method", "algorithm_version", "coordinate_system", "gap_policy"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ConfigurationError(f"PatternResult {name} must be non-empty.")
        if self.sequence_id is not None and (
            not isinstance(self.sequence_id, str) or not self.sequence_id.strip()
        ):
            raise ConfigurationError("PatternResult sequence_id must be non-empty or None.")
        for name in ("inspected_symbol_count", "gap_count", "unknown_gap_count"):
            _require_non_negative(getattr(self, name), name)
        if not isinstance(self.parameters, FrozenDict):
            raise ConfigurationError("PatternResult parameters must be FrozenDict.")
        if not isinstance(self.hits, tuple):
            raise ConfigurationError("PatternResult hits must be a tuple.")
        if isinstance(self.max_matches, bool) or not isinstance(self.max_matches, int):
            raise ConfigurationError("max_matches must be a positive integer.")
        if self.max_matches <= 0 or len(self.hits) > self.max_matches:
            raise ConfigurationError("PatternResult hits exceed max_matches.")
        if not isinstance(self.truncated, bool):
            raise ConfigurationError("PatternResult truncated must be boolean.")
        if self.coordinate_system != "0-based-half-open":
            raise ConfigurationError("PatternResult coordinates must be 0-based half-open.")
        if self.topology not in ("linear", "circular"):
            raise ConfigurationError("PatternResult topology must be linear or circular.")
        if not isinstance(self.provenance, Provenance):
            raise ConfigurationError("PatternResult provenance must be Provenance.")
        if any(not isinstance(issue, Issue) for issue in self.issues):
            raise ConfigurationError("PatternResult issues must contain Issue objects.")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


@dataclass(frozen=True)
class MotifHit:
    """One exact, IUPAC, regular-expression, or PWM match."""

    motif_name: str
    matched_sequence: str
    strand: Strand
    symbol_location: Location
    coordinate_location: Location
    score: float | None
    threshold: float | None
    wraps_origin: bool

    def __post_init__(self) -> None:
        _require_text(self.motif_name, "motif_name")
        _require_text(self.matched_sequence, "matched_sequence")
        _require_strand(self.strand, allow_both=True)
        _require_location(self.symbol_location, "symbol_location")
        _require_location(self.coordinate_location, "coordinate_location")
        validate_finite_score(self.score, "score")
        validate_finite_score(self.threshold, "threshold")
        if (self.score is None) != (self.threshold is None):
            raise ConfigurationError(
                "MotifHit score and threshold must both be set or both be None."
            )
        if not isinstance(self.wraps_origin, bool):
            raise ConfigurationError("wraps_origin must be boolean.")


@dataclass(frozen=True)
class CodonSite:
    """One in-frame start or stop codon."""

    kind: str
    codon: str
    strand: Strand
    frame: int
    symbol_location: Location
    coordinate_location: Location

    def __post_init__(self) -> None:
        if self.kind not in ("start", "stop"):
            raise ConfigurationError("CodonSite kind must be start or stop.")
        if len(self.codon) != 3 or set(self.codon) - set("ACGT"):
            raise ConfigurationError("CodonSite codon must be a canonical DNA triplet.")
        _require_strand(self.strand)
        if self.frame not in (-3, -2, -1, 1, 2, 3):
            raise ConfigurationError("CodonSite frame must be one of +/-1, +/-2, or +/-3.")
        _require_location(self.symbol_location, "symbol_location")
        _require_location(self.coordinate_location, "coordinate_location")


@dataclass(frozen=True)
class ORFHit:
    """One start-anchored open reading frame."""

    strand: Strand
    frame: int
    symbol_location: Location
    coordinate_location: Location
    nucleotide_length: int
    codon_count: int
    start_codon: str
    stop_codon: str | None
    complete: bool
    translation: str

    def __post_init__(self) -> None:
        _require_strand(self.strand)
        if self.frame not in (-3, -2, -1, 1, 2, 3):
            raise ConfigurationError("ORFHit frame must be one of +/-1, +/-2, or +/-3.")
        _require_location(self.symbol_location, "symbol_location")
        _require_location(self.coordinate_location, "coordinate_location")
        _require_positive(self.nucleotide_length, "nucleotide_length")
        _require_positive(self.codon_count, "codon_count")
        if self.nucleotide_length != 3 * self.codon_count:
            raise ConfigurationError("ORF nucleotide_length must equal three times codon_count.")
        if len(self.start_codon) != 3 or set(self.start_codon) - set("ACGT"):
            raise ConfigurationError("ORF start_codon must be a canonical DNA triplet.")
        if self.stop_codon is not None and (
            len(self.stop_codon) != 3 or set(self.stop_codon) - set("ACGT")
        ):
            raise ConfigurationError("ORF stop_codon must be a canonical triplet or None.")
        if not isinstance(self.complete, bool) or self.complete != (self.stop_codon is not None):
            raise ConfigurationError("ORF complete must agree with stop_codon presence.")
        if not isinstance(self.translation, str) or len(self.translation) != self.codon_count:
            raise ConfigurationError("ORF translation length must equal codon_count.")


@dataclass(frozen=True)
class RestrictionSiteHit:
    """One recognition site and its top/bottom strand cut coordinates."""

    enzyme: str
    recognition_sequence: str
    strand: Strand
    symbol_location: Location
    coordinate_location: Location
    top_cut: int | None
    bottom_cut: int | None
    wraps_origin: bool

    def __post_init__(self) -> None:
        _require_text(self.enzyme, "enzyme")
        _require_iupac(self.recognition_sequence, "recognition_sequence")
        _require_strand(self.strand, allow_both=True)
        _require_location(self.symbol_location, "symbol_location")
        _require_location(self.coordinate_location, "coordinate_location")
        for name in ("top_cut", "bottom_cut"):
            value = getattr(self, name)
            if value is not None:
                _require_non_negative(value, name)
        if not isinstance(self.wraps_origin, bool):
            raise ConfigurationError("wraps_origin must be boolean.")


@dataclass(frozen=True)
class GuideCandidate:
    """A rule-matched CRISPR guide and adjacent PAM; no efficiency prediction."""

    nuclease: str
    strand: Strand
    guide_sequence: str
    pam_sequence: str
    guide_symbol_location: Location
    pam_symbol_location: Location
    guide_coordinate_location: Location
    pam_coordinate_location: Location
    gc_fraction: float
    wraps_origin: bool

    def __post_init__(self) -> None:
        _require_text(self.nuclease, "nuclease")
        _require_strand(self.strand)
        _require_iupac(self.guide_sequence, "guide_sequence")
        _require_iupac(self.pam_sequence, "pam_sequence")
        for name in (
            "guide_symbol_location",
            "pam_symbol_location",
            "guide_coordinate_location",
            "pam_coordinate_location",
        ):
            _require_location(getattr(self, name), name)
        if (
            isinstance(self.gc_fraction, bool)
            or not isinstance(self.gc_fraction, (int, float))
            or not math.isfinite(self.gc_fraction)
            or not 0 <= self.gc_fraction <= 1
        ):
            raise ConfigurationError("GuideCandidate gc_fraction must be between 0 and 1.")
        if not isinstance(self.wraps_origin, bool):
            raise ConfigurationError("wraps_origin must be boolean.")


@dataclass(frozen=True)
class RegionHit:
    """One region identified by a deterministic window or composition rule."""

    kind: str
    symbol_location: Location
    coordinate_location: Location
    length: int
    score: float | None
    attributes: FrozenDict

    def __post_init__(self) -> None:
        _require_text(self.kind, "kind")
        _require_location(self.symbol_location, "symbol_location")
        _require_location(self.coordinate_location, "coordinate_location")
        _require_positive(self.length, "length")
        validate_finite_score(self.score, "score")
        if not isinstance(self.attributes, FrozenDict):
            raise ConfigurationError("RegionHit attributes must be FrozenDict.")


@dataclass(frozen=True)
class PalindromeHit:
    """One reverse-complement palindrome."""

    sequence: str
    symbol_location: Location
    coordinate_location: Location
    length: int

    def __post_init__(self) -> None:
        _require_iupac(self.sequence, "sequence")
        _require_location(self.symbol_location, "symbol_location")
        _require_location(self.coordinate_location, "coordinate_location")
        _require_positive(self.length, "length")
        if len(self.sequence) != self.length:
            raise ConfigurationError("Palindrome sequence length must equal length.")


@dataclass(frozen=True)
class InvertedRepeatHit:
    """Two reverse-complement arms separated by a loop."""

    left_arm: str
    right_arm: str
    arm_length: int
    loop_length: int
    symbol_location: Location
    coordinate_location: Location
    left_symbol_location: Location
    loop_symbol_location: Location
    right_symbol_location: Location

    def __post_init__(self) -> None:
        _require_iupac(self.left_arm, "left_arm")
        _require_iupac(self.right_arm, "right_arm")
        _require_positive(self.arm_length, "arm_length")
        _require_non_negative(self.loop_length, "loop_length")
        if len(self.left_arm) != self.arm_length or len(self.right_arm) != self.arm_length:
            raise ConfigurationError("Inverted-repeat arm lengths must equal arm_length.")
        for name in (
            "symbol_location",
            "coordinate_location",
            "left_symbol_location",
            "loop_symbol_location",
            "right_symbol_location",
        ):
            _require_location(getattr(self, name), name)


@dataclass(frozen=True)
class TandemRepeatHit:
    """One maximal exact tandem repeat using its smallest accepted period."""

    unit: str
    unit_length: int
    repeat_count: int
    sequence: str
    symbol_location: Location
    coordinate_location: Location

    def __post_init__(self) -> None:
        _require_iupac(self.unit, "unit")
        _require_positive(self.unit_length, "unit_length")
        _require_positive(self.repeat_count, "repeat_count")
        if self.repeat_count < 2:
            raise ConfigurationError("TandemRepeatHit repeat_count must be at least 2.")
        if len(self.unit) != self.unit_length or self.sequence != self.unit * self.repeat_count:
            raise ConfigurationError(
                "Tandem-repeat sequence must equal unit repeated repeat_count times."
            )
        _require_location(self.symbol_location, "symbol_location")
        _require_location(self.coordinate_location, "coordinate_location")


@dataclass(frozen=True)
class LowComplexityResult:
    """Low-complexity regions and a gap-preserving IUPAC masked sequence."""

    analysis: PatternResult[RegionHit]
    masked_sequence: DNASequence
    mask_symbol: str

    def __post_init__(self) -> None:
        if not isinstance(self.analysis, PatternResult):
            raise ConfigurationError("LowComplexityResult analysis must be PatternResult.")
        if not isinstance(self.masked_sequence, DNASequence):
            raise ConfigurationError("LowComplexityResult masked_sequence must be DNASequence.")
        if not isinstance(self.mask_symbol, str) or len(self.mask_symbol) != 1:
            raise ConfigurationError("LowComplexityResult mask_symbol must be one character.")
        _require_iupac(self.mask_symbol, "mask_symbol")

    def to_dict(self) -> dict[str, Any]:
        return cast(dict[str, Any], to_json_compatible(self))


def validate_finite_score(value: float | None, name: str) -> None:
    """Validate optional scores used by public pattern definitions."""

    if value is not None and (
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
    ):
        raise ConfigurationError(f"{name} must be a finite number or None.")


__all__ = [
    "CodonSite",
    "GuideCandidate",
    "InvertedRepeatHit",
    "LowComplexityResult",
    "MotifHit",
    "ORFHit",
    "PalindromeHit",
    "PatternResult",
    "RegionHit",
    "RestrictionSiteHit",
    "TandemRepeatHit",
]
