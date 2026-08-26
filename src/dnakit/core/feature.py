"""Immutable DNA feature annotations."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

from dnakit.core._json import FrozenDict, freeze_mapping
from dnakit.core.coordinates import CompoundLocation, Interval, Location, UnresolvedLocation
from dnakit.core.enums import Strand
from dnakit.exceptions import FeatureError


@dataclass(frozen=True, init=False)
class DNAFeature:
    """A typed annotation referring to a location without copying sequence data."""

    type: str
    location: Location
    id: str | None
    strand: Strand
    label: str | None
    score: float | None
    phase: int | None
    qualifiers: FrozenDict
    source: str | None

    __hash__ = None  # type: ignore[assignment]

    def __init__(
        self,
        type: str,
        location: Location,
        *,
        id: str | None = None,
        strand: Strand | str = Strand.UNKNOWN,
        label: str | None = None,
        score: float | None = None,
        phase: int | None = None,
        qualifiers: Mapping[str, object] | None = None,
        source: str | None = None,
    ) -> None:
        if not isinstance(type, str) or not type.strip():
            raise FeatureError("Feature type must be a non-empty string.")
        if not isinstance(location, (Interval, CompoundLocation, UnresolvedLocation)):
            raise FeatureError("Feature location must be an internal Location object.")
        for field_name, value in (("id", id), ("label", label), ("source", source)):
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise FeatureError(
                    f"Feature {field_name} must be None or a non-empty string.",
                    context={field_name: value},
                )
        try:
            resolved_strand = strand if isinstance(strand, Strand) else Strand(strand)
        except (TypeError, ValueError) as exc:
            raise FeatureError(
                "Unknown feature strand.",
                context={"strand": strand},
                hint=f"Choose one of: {', '.join(item.value for item in Strand)}.",
            ) from exc
        if score is not None and (
            isinstance(score, bool)
            or not isinstance(score, (int, float))
            or not math.isfinite(score)
        ):
            raise FeatureError(
                "Feature score must be a finite number or None.",
                context={"score": score},
            )
        if phase is not None and (
            isinstance(phase, bool) or not isinstance(phase, int) or phase not in (0, 1, 2)
        ):
            raise FeatureError(
                "Feature phase must be 0, 1, 2, or None.",
                context={"phase": phase},
            )

        object.__setattr__(self, "type", type)
        object.__setattr__(self, "location", location)
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "strand", resolved_strand)
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "score", None if score is None else float(score))
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "qualifiers", freeze_mapping(qualifiers))
        object.__setattr__(self, "source", source)


__all__ = ["DNAFeature"]
