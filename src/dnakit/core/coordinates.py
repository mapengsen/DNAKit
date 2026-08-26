"""Internal and external DNA coordinate value objects and conversions."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TypeAlias

from dnakit.core.enums import CoordinateSystem, Strand
from dnakit.exceptions import CoordinateError


def _require_non_negative_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise CoordinateError(
            f"{name} must be a non-negative integer.",
            context={name: value},
        )


def _coerce_coordinate_system(value: CoordinateSystem | str) -> CoordinateSystem:
    try:
        return value if isinstance(value, CoordinateSystem) else CoordinateSystem(value)
    except (TypeError, ValueError) as exc:
        raise CoordinateError(
            "Unknown coordinate system.",
            context={"system": value},
            hint=f"Choose one of: {', '.join(item.value for item in CoordinateSystem)}.",
        ) from exc


def _coerce_strand(value: Strand | str) -> Strand:
    try:
        return value if isinstance(value, Strand) else Strand(value)
    except (TypeError, ValueError) as exc:
        raise CoordinateError(
            "Unknown strand value.",
            context={"strand": value},
            hint=f"Choose one of: {', '.join(item.value for item in Strand)}.",
        ) from exc


@dataclass(frozen=True, order=True)
class Interval:
    """An internal zero-based, half-open interval."""

    start: int
    end: int

    def __post_init__(self) -> None:
        _require_non_negative_integer(self.start, "start")
        _require_non_negative_integer(self.end, "end")
        if self.end < self.start:
            raise CoordinateError(
                "Internal interval end cannot be smaller than start.",
                context={"start": self.start, "end": self.end},
            )

    def __len__(self) -> int:
        return self.end - self.start


@dataclass(frozen=True, init=False)
class CompoundLocation:
    """An ordered set of internal intervals, including origin-spanning features."""

    parts: tuple[Interval, ...]

    def __init__(self, parts: Iterable[Interval]) -> None:
        resolved = tuple(parts)
        if not resolved:
            raise CoordinateError("A compound location requires at least one interval.")
        if any(not isinstance(part, Interval) for part in resolved):
            raise CoordinateError("Compound location parts must all be Interval objects.")
        object.__setattr__(self, "parts", resolved)

    def __len__(self) -> int:
        return sum(len(part) for part in self.parts)


@dataclass(frozen=True, init=False)
class UnresolvedLocation:
    """A location that cannot be resolved, usually because a gap length is unknown."""

    reason: str
    anchors: tuple[Interval, ...]

    def __init__(self, reason: str, anchors: Iterable[Interval] = ()) -> None:
        if not isinstance(reason, str) or not reason.strip():
            raise CoordinateError("An unresolved location requires a non-empty reason.")
        resolved_anchors = tuple(anchors)
        if any(not isinstance(anchor, Interval) for anchor in resolved_anchors):
            raise CoordinateError("Unresolved location anchors must be Interval objects.")
        object.__setattr__(self, "reason", reason)
        object.__setattr__(self, "anchors", resolved_anchors)


Location: TypeAlias = Interval | CompoundLocation | UnresolvedLocation


@dataclass(frozen=True, init=False)
class ExternalInterval:
    """An interval retaining its external coordinate convention and strand reference."""

    start: int
    end: int
    system: CoordinateSystem
    strand: Strand

    def __init__(
        self,
        start: int,
        end: int,
        system: CoordinateSystem | str,
        strand: Strand | str = Strand.UNKNOWN,
    ) -> None:
        resolved_system = _coerce_coordinate_system(system)
        minimum = 0 if resolved_system.value.startswith("0-") else 1
        if isinstance(start, bool) or not isinstance(start, int) or start < minimum:
            raise CoordinateError(
                "External interval start is outside its coordinate system.",
                context={"start": start, "system": resolved_system.value},
            )
        if isinstance(end, bool) or not isinstance(end, int) or end < minimum:
            raise CoordinateError(
                "External interval end is outside its coordinate system.",
                context={"end": end, "system": resolved_system.value},
            )
        object.__setattr__(self, "start", start)
        object.__setattr__(self, "end", end)
        object.__setattr__(self, "system", resolved_system)
        object.__setattr__(self, "strand", _coerce_strand(strand))


def _validate_sequence_length(sequence_length: int | None) -> None:
    if sequence_length is not None:
        _require_non_negative_integer(sequence_length, "sequence_length")


def _to_internal_bounds(interval: ExternalInterval) -> tuple[int, int]:
    if interval.system is CoordinateSystem.ZERO_BASED_HALF_OPEN:
        return interval.start, interval.end
    if interval.system is CoordinateSystem.ZERO_BASED_CLOSED:
        return interval.start, interval.end + 1
    if interval.system is CoordinateSystem.ONE_BASED_CLOSED:
        return interval.start - 1, interval.end
    if interval.system is CoordinateSystem.ONE_BASED_HALF_OPEN:
        return interval.start - 1, interval.end - 1
    raise AssertionError("Unhandled coordinate system.")


def _import_one(interval: ExternalInterval, sequence_length: int | None) -> tuple[Interval, ...]:
    start, end = _to_internal_bounds(interval)
    wraps_origin = interval.start > interval.end
    if not wraps_origin:
        if sequence_length is not None and (start > sequence_length or end > sequence_length):
            raise CoordinateError(
                "External interval exceeds the sequence length.",
                context={
                    "start": interval.start,
                    "end": interval.end,
                    "sequence_length": sequence_length,
                },
            )
        return (Interval(start, end),)

    if sequence_length is None:
        raise CoordinateError(
            "A wrapped interval requires sequence_length.",
            context={"start": interval.start, "end": interval.end},
            hint="Provide the circular sequence length or split the interval explicitly.",
        )
    if start > sequence_length or end > sequence_length:
        raise CoordinateError(
            "Wrapped interval exceeds the sequence length.",
            context={
                "start": interval.start,
                "end": interval.end,
                "sequence_length": sequence_length,
            },
        )
    parts = tuple(
        part for part in (Interval(start, sequence_length), Interval(0, end)) if len(part) > 0
    )
    if not parts:
        raise CoordinateError("Wrapped interval resolves to an empty location.")
    return parts


def import_location(
    external: ExternalInterval | Sequence[ExternalInterval],
    *,
    sequence_length: int | None = None,
) -> Location:
    """Normalize one or more external intervals to internal coordinates."""

    _validate_sequence_length(sequence_length)
    intervals = (external,) if isinstance(external, ExternalInterval) else tuple(external)
    if not intervals:
        raise CoordinateError("At least one external interval is required.")
    if any(not isinstance(interval, ExternalInterval) for interval in intervals):
        raise CoordinateError("All external locations must be ExternalInterval objects.")

    parts = tuple(part for interval in intervals for part in _import_one(interval, sequence_length))
    return parts[0] if len(parts) == 1 else CompoundLocation(parts)


def _from_internal_bounds(
    interval: Interval,
    target_system: CoordinateSystem,
) -> tuple[int, int]:
    if target_system is CoordinateSystem.ZERO_BASED_HALF_OPEN:
        return interval.start, interval.end
    if target_system is CoordinateSystem.ZERO_BASED_CLOSED:
        if len(interval) == 0:
            raise CoordinateError("A closed coordinate system cannot represent an empty interval.")
        return interval.start, interval.end - 1
    if target_system is CoordinateSystem.ONE_BASED_CLOSED:
        if len(interval) == 0:
            raise CoordinateError("A closed coordinate system cannot represent an empty interval.")
        return interval.start + 1, interval.end
    if target_system is CoordinateSystem.ONE_BASED_HALF_OPEN:
        return interval.start + 1, interval.end + 1
    raise AssertionError("Unhandled coordinate system.")


def export_location(
    location: Location,
    *,
    target_system: CoordinateSystem | str,
    sequence_length: int | None = None,
) -> tuple[ExternalInterval, ...]:
    """Export internal coordinates without merging compound-location parts."""

    _validate_sequence_length(sequence_length)
    system = _coerce_coordinate_system(target_system)
    if isinstance(location, UnresolvedLocation):
        raise CoordinateError(
            "An unresolved location cannot be exported.",
            context={"reason": location.reason},
        )
    if isinstance(location, Interval):
        parts: tuple[Interval, ...] = (location,)
    elif isinstance(location, CompoundLocation):
        parts = location.parts
    else:
        raise CoordinateError("location must be an internal Location object.")

    if sequence_length is not None and any(part.end > sequence_length for part in parts):
        raise CoordinateError(
            "Internal location exceeds the sequence length.",
            context={"sequence_length": sequence_length},
        )
    return tuple(
        ExternalInterval(*_from_internal_bounds(part, system), system=system) for part in parts
    )


def reverse_strand_location(location: Location, *, sequence_length: int) -> Location:
    """Map an internal location through reverse-complement coordinates."""

    _validate_sequence_length(sequence_length)
    if isinstance(location, Interval):
        return _reverse_interval(location, sequence_length)
    if isinstance(location, CompoundLocation):
        return CompoundLocation(
            _reverse_interval(part, sequence_length) for part in reversed(location.parts)
        )
    if isinstance(location, UnresolvedLocation):
        return UnresolvedLocation(
            location.reason,
            (_reverse_interval(anchor, sequence_length) for anchor in reversed(location.anchors)),
        )
    raise CoordinateError("location must be an internal Location object.")


def _reverse_interval(location: Interval, sequence_length: int) -> Interval:
    if location.end > sequence_length:
        raise CoordinateError(
            "Location exceeds the sequence length.",
            context={"end": location.end, "sequence_length": sequence_length},
        )
    return Interval(sequence_length - location.end, sequence_length - location.start)


__all__ = [
    "CompoundLocation",
    "ExternalInterval",
    "Interval",
    "Location",
    "UnresolvedLocation",
    "export_location",
    "import_location",
    "reverse_strand_location",
]
