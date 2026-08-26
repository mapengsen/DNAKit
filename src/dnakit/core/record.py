"""Immutable DNA record value object."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import islice
from typing import cast

from dnakit.core._json import FrozenDict, freeze_mapping
from dnakit.core.coordinates import CompoundLocation, UnresolvedLocation
from dnakit.core.feature import DNAFeature
from dnakit.core.gap import Gap
from dnakit.core.sequence import DNASequence
from dnakit.exceptions import SequenceError


@dataclass(frozen=True, init=False)
class DNARecord:
    """A DNA sequence plus stable identity, annotations, and metadata."""

    sequence: DNASequence
    id: str
    description: str
    features: tuple[DNAFeature, ...]
    metadata: FrozenDict
    letter_annotations: Mapping[str, tuple[int | float, ...]]

    __hash__ = None  # type: ignore[assignment]

    def __init__(
        self,
        sequence: DNASequence,
        id: str,
        *,
        description: str = "",
        features: Iterable[DNAFeature] = (),
        metadata: Mapping[str, object] | None = None,
        letter_annotations: Mapping[str, Iterable[int | float]] | None = None,
    ) -> None:
        if not isinstance(sequence, DNASequence):
            raise SequenceError(
                "DNARecord sequence must be a DNASequence object.",
                code="INVALID_RECORD_SEQUENCE",
                context={"sequence_type": type(sequence).__name__},
            )
        if not isinstance(id, str) or not id.strip():
            raise SequenceError(
                "DNARecord id must be a non-empty string.",
                code="INVALID_RECORD_ID",
            )
        if not isinstance(description, str):
            raise SequenceError(
                "DNARecord description must be a string.",
                code="INVALID_RECORD_DESCRIPTION",
            )
        feature_tuple = tuple(features)
        if any(not isinstance(feature, DNAFeature) for feature in feature_tuple):
            raise SequenceError(
                "DNARecord features must all be DNAFeature objects.",
                code="INVALID_RECORD_FEATURE",
            )
        resolvable_span = 0
        for part in sequence.parts:
            if isinstance(part, str):
                resolvable_span += len(part)
            elif isinstance(part, Gap) and part.length is not None:
                resolvable_span += part.length
            else:
                break
        for index, feature in enumerate(feature_tuple):
            location = feature.location
            if isinstance(location, UnresolvedLocation):
                intervals = location.anchors
            else:
                intervals = (
                    location.parts if isinstance(location, CompoundLocation) else (location,)
                )
            if any(interval.end > resolvable_span for interval in intervals):
                unknown_span = sequence.coordinate_span is None
                raise SequenceError(
                    (
                        "A resolved feature extends beyond the resolvable sequence prefix."
                        if unknown_span
                        else "A feature location extends beyond the sequence coordinate span."
                    ),
                    code=(
                        "FEATURE_LOCATION_UNRESOLVED" if unknown_span else "FEATURE_OUT_OF_BOUNDS"
                    ),
                    context={
                        "feature_index": index,
                        "feature_type": feature.type,
                        "location_kind": type(location).__name__,
                        "resolvable_span": resolvable_span,
                        "sequence_span": sequence.coordinate_span,
                    },
                    hint=(
                        "Use UnresolvedLocation for coordinates at or after an unknown-length Gap."
                        if unknown_span
                        else None
                    ),
                )

        annotations_input: dict[str, object] = {}
        for name, values in (letter_annotations or {}).items():
            if not isinstance(name, str) or not name.strip():
                raise SequenceError(
                    "Letter annotation names must be non-empty strings.",
                    code="INVALID_LETTER_ANNOTATION_NAME",
                )
            if isinstance(values, (str, bytes)):
                raise SequenceError(
                    "Letter annotation values must be numeric iterables.",
                    code="INVALID_LETTER_ANNOTATION",
                    context={"name": name},
                )
            try:
                value_tuple = tuple(islice(iter(values), sequence.symbol_length + 1))
            except TypeError as exc:
                raise SequenceError(
                    "Letter annotation values must be numeric iterables.",
                    code="INVALID_LETTER_ANNOTATION",
                    context={"name": name},
                ) from exc
            if len(value_tuple) != sequence.symbol_length:
                raise SequenceError(
                    "Letter annotation length must equal sequence symbol_length.",
                    code="LETTER_ANNOTATION_LENGTH_MISMATCH",
                    context={
                        "name": name,
                        "annotation_length": len(value_tuple),
                        "symbol_length": sequence.symbol_length,
                    },
                )
            if any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                for value in value_tuple
            ):
                raise SequenceError(
                    "Letter annotations must contain only finite numeric values.",
                    code="INVALID_LETTER_ANNOTATION",
                    context={"name": name},
                )
            annotations_input[name] = value_tuple

        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "id", id)
        object.__setattr__(self, "description", description)
        object.__setattr__(self, "features", feature_tuple)
        object.__setattr__(self, "metadata", freeze_mapping(metadata))
        frozen_annotations = cast(
            Mapping[str, tuple[int | float, ...]],
            freeze_mapping(annotations_input),
        )
        object.__setattr__(self, "letter_annotations", frozen_annotations)


__all__ = ["DNARecord"]
