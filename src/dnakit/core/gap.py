"""Explicit known- and unknown-length sequence gaps."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from dnakit.core._json import FrozenDict, freeze_mapping
from dnakit.core.enums import GapKind
from dnakit.exceptions import SequenceError


@dataclass(frozen=True, init=False)
class Gap:
    """A gap between sequence fragments.

    ``length=None`` denotes an unknown coordinate span.  A gap never stores a
    second copy of its flanking sequence; the flanks are determined by its
    position in :class:`DNASequence.parts`.
    """

    length: int | None
    kind: GapKind
    crossable: bool | None
    evidence: tuple[str, ...]
    metadata: FrozenDict

    def __init__(
        self,
        length: int | None,
        kind: GapKind | str = GapKind.UNKNOWN,
        crossable: bool | None = None,
        evidence: Iterable[str] = (),
        metadata: Mapping[str, object] | None = None,
    ) -> None:
        if length is not None and (isinstance(length, bool) or not isinstance(length, int)):
            raise SequenceError(
                "Gap length must be a positive integer or None.",
                code="INVALID_GAP_LENGTH",
                context={"length": length},
            )
        if length is not None and length <= 0:
            raise SequenceError(
                "Gap length must be greater than zero.",
                code="INVALID_GAP_LENGTH",
                context={"length": length},
                hint="Use None when the gap length is unknown.",
            )
        try:
            resolved_kind = kind if isinstance(kind, GapKind) else GapKind(kind)
        except (TypeError, ValueError) as exc:
            raise SequenceError(
                "Unknown gap kind.",
                code="INVALID_GAP_KIND",
                context={"kind": kind},
                hint=f"Choose one of: {', '.join(item.value for item in GapKind)}.",
            ) from exc
        if crossable is not None and not isinstance(crossable, bool):
            raise SequenceError(
                "Gap crossable policy must be True, False, or None.",
                code="INVALID_GAP_POLICY",
                context={"crossable": crossable},
            )
        evidence_tuple = tuple(evidence)
        if any(not isinstance(item, str) or not item.strip() for item in evidence_tuple):
            raise SequenceError(
                "Gap evidence entries must be non-empty strings.",
                code="INVALID_GAP_EVIDENCE",
            )

        object.__setattr__(self, "length", length)
        object.__setattr__(self, "kind", resolved_kind)
        object.__setattr__(self, "crossable", crossable)
        object.__setattr__(self, "evidence", evidence_tuple)
        object.__setattr__(self, "metadata", freeze_mapping(metadata))

    @property
    def is_known_length(self) -> bool:
        return self.length is not None


__all__ = ["Gap"]
