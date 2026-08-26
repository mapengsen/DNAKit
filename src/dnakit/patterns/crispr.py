"""Rule-based CRISPR PAM and adjacent-guide candidate scanning."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import islice
from typing import Literal, cast

from dnakit.core import Location, Strand
from dnakit.exceptions import ConfigurationError
from dnakit.patterns._shared import (
    SequenceInput,
    build_result,
    circular_location,
    coerce_strands,
    iupac_compatible,
    pattern_provenance,
    require_ungapped_circular,
    resolve_sequence,
    reverse_complement_text,
    segment_location,
    segments,
    validate_bool,
    validate_iupac_text,
    validate_positive_int,
)
from dnakit.patterns.results import GuideCandidate, PatternResult

PAMSide = Literal["3prime", "5prime"]
PAM_CATALOG_VERSION = "common-nucleases-v1"


@dataclass(frozen=True, init=False)
class PAMRule:
    """A nuclease name, IUPAC PAM, guide orientation, and default guide length."""

    name: str
    pam: str
    pam_side: PAMSide
    guide_length: int
    source: str

    def __init__(
        self,
        name: str,
        pam: str,
        pam_side: PAMSide,
        guide_length: int,
        *,
        source: str = "user-provided",
    ) -> None:
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("PAM rule name must be non-empty.")
        resolved_pam = validate_iupac_text(pam, "pam")
        if pam_side not in ("3prime", "5prime"):
            raise ConfigurationError("pam_side must be '3prime' or '5prime'.")
        validate_positive_int(guide_length, "guide_length")
        if not isinstance(source, str) or not source.strip():
            raise ConfigurationError("PAM rule source must be non-empty.")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "pam", resolved_pam)
        object.__setattr__(self, "pam_side", pam_side)
        object.__setattr__(self, "guide_length", guide_length)
        object.__setattr__(self, "source", source)


_SOURCE = "DNAKit common PAM boundary catalog; verify experimental design against primary sources"
BUILTIN_PAM_RULES: Mapping[str, PAMRule] = {
    rule.name: rule
    for rule in (
        PAMRule("SpCas9", "NGG", "3prime", 20, source=_SOURCE),
        PAMRule("SaCas9", "NNGRRT", "3prime", 21, source=_SOURCE),
        PAMRule("AsCas12a", "TTTV", "5prime", 20, source=_SOURCE),
    )
}


def _resolve_rule(rule: str | PAMRule) -> PAMRule:
    if isinstance(rule, PAMRule):
        return rule
    if isinstance(rule, str):
        try:
            return BUILTIN_PAM_RULES[rule]
        except KeyError as exc:
            raise ConfigurationError(
                "Unknown built-in PAM rule.",
                context={"rule": rule, "available": sorted(BUILTIN_PAM_RULES)},
                hint="Pass a PAMRule for another nuclease.",
            ) from exc
    raise ConfigurationError("rule must be a built-in name or PAMRule.")


def _slice_circular(text: str, start: int, length: int) -> str:
    if length > len(text):
        return ""
    start %= len(text)
    doubled = text + text
    return doubled[start : start + length]


def _location(
    *,
    item: object,
    oriented_start: int,
    length: int,
    strand: Strand,
    circular: bool,
) -> tuple[Location, Location, bool]:
    from dnakit.patterns._shared import Segment

    if not isinstance(item, Segment):
        raise AssertionError("Internal segment type mismatch.")
    original_start = (
        oriented_start if strand is Strand.FORWARD else len(item.text) - (oriented_start + length)
    )
    if circular:
        return circular_location(original_start, length, len(item.text))
    symbol_location, coordinate_location = segment_location(
        item, original_start, original_start + length
    )
    return symbol_location, coordinate_location, False


def _contains_excluded(guide: str, motifs: tuple[str, ...]) -> bool:
    for motif in motifs:
        for start in range(len(guide) - len(motif) + 1):
            if all(
                iupac_compatible(guide[start + offset], symbol)
                for offset, symbol in enumerate(motif)
            ):
                return True
    return False


def scan_pam_candidates(
    value: SequenceInput,
    rule: str | PAMRule,
    *,
    guide_length: int | None = None,
    strand: Strand | str = Strand.BOTH,
    min_gc: float = 0.0,
    max_gc: float = 1.0,
    exclude_motifs: Iterable[str] = (),
    allow_ambiguous_guides: bool = False,
    max_matches: int = 100_000,
    max_scan_length: int = 1_000_000,
    max_pam_length: int = 1_000,
    max_scan_cells: int = 50_000_000,
    max_exclude_motifs: int = 1_000,
    max_filter_cells: int = 50_000_000,
) -> PatternResult[GuideCandidate]:
    """Find PAM-adjacent guide candidates without efficiency/off-target prediction."""

    sequence, sequence_id = resolve_sequence(value)
    definition = _resolve_rule(rule)
    resolved_guide_length = definition.guide_length if guide_length is None else guide_length
    validate_positive_int(resolved_guide_length, "guide_length")
    validate_positive_int(max_matches, "max_matches")
    validate_positive_int(max_scan_length, "max_scan_length")
    validate_positive_int(max_pam_length, "max_pam_length")
    validate_positive_int(max_scan_cells, "max_scan_cells")
    validate_positive_int(max_exclude_motifs, "max_exclude_motifs")
    validate_positive_int(max_filter_cells, "max_filter_cells")
    validate_bool(allow_ambiguous_guides, "allow_ambiguous_guides")
    if sequence.symbol_length > max_scan_length:
        raise ConfigurationError("Sequence exceeds max_scan_length.")
    if len(definition.pam) > max_pam_length:
        raise ConfigurationError("PAM definition exceeds max_pam_length.")
    for value_, name in ((min_gc, "min_gc"), (max_gc, "max_gc")):
        if (
            isinstance(value_, bool)
            or not isinstance(value_, (int, float))
            or not math.isfinite(value_)
            or value_ < 0
            or value_ > 1
        ):
            raise ConfigurationError(f"{name} must be finite and between 0 and 1.")
    if min_gc > max_gc:
        raise ConfigurationError("min_gc cannot exceed max_gc.")
    if isinstance(exclude_motifs, (str, bytes)):
        raise ConfigurationError("exclude_motifs must be an iterable of IUPAC motifs.")
    try:
        raw_excluded = tuple(islice(iter(exclude_motifs), max_exclude_motifs + 1))
    except TypeError as exc:
        raise ConfigurationError("exclude_motifs must be an iterable of IUPAC motifs.") from exc
    if len(raw_excluded) > max_exclude_motifs:
        raise ConfigurationError("exclude_motifs exceeds max_exclude_motifs.")
    excluded = tuple(validate_iupac_text(motif, "exclude_motif") for motif in raw_excluded)
    circular = sequence.topology.value == "circular"
    require_ungapped_circular(sequence, "CRISPR PAM scanning")
    total_footprint = resolved_guide_length + len(definition.pam)
    selected_strands = coerce_strands(strand)
    estimated_cells = sequence.symbol_length * len(definition.pam) * len(selected_strands)
    if estimated_cells > max_scan_cells:
        raise ConfigurationError(
            "CRISPR PAM scan exceeds max_scan_cells.",
            context={"estimated_cells": estimated_cells, "max_scan_cells": max_scan_cells},
        )
    estimated_filter_cells = (
        sequence.symbol_length
        * len(selected_strands)
        * sum(max(0, resolved_guide_length - len(motif) + 1) * len(motif) for motif in excluded)
    )
    if estimated_filter_cells > max_filter_cells:
        raise ConfigurationError(
            "CRISPR guide filters exceed max_filter_cells.",
            context={
                "estimated_filter_cells": estimated_filter_cells,
                "max_filter_cells": max_filter_cells,
            },
        )

    hits: list[GuideCandidate] = []
    truncated = False
    for item in segments(sequence):
        if circular and total_footprint > len(item.text):
            continue
        for resolved_strand in selected_strands:
            oriented = (
                item.text
                if resolved_strand is Strand.FORWARD
                else reverse_complement_text(item.text)
            )
            pam_length = len(definition.pam)
            searchable = oriented + oriented[: pam_length - 1] if circular else oriented
            limit = len(oriented) if circular else len(oriented) - pam_length + 1
            for pam_start in range(max(0, limit)):
                pam_sequence = searchable[pam_start : pam_start + pam_length]
                if not all(
                    iupac_compatible(target, pattern)
                    for target, pattern in zip(pam_sequence, definition.pam, strict=True)
                ):
                    continue
                guide_start = (
                    pam_start - resolved_guide_length
                    if definition.pam_side == "3prime"
                    else pam_start + pam_length
                )
                if not circular and (
                    guide_start < 0 or guide_start + resolved_guide_length > len(oriented)
                ):
                    continue
                guide_sequence = (
                    _slice_circular(oriented, guide_start, resolved_guide_length)
                    if circular
                    else oriented[guide_start : guide_start + resolved_guide_length]
                )
                if not allow_ambiguous_guides and set(guide_sequence) - set("ACGT"):
                    continue
                gc_fraction = sum(base in "GC" for base in guide_sequence) / len(guide_sequence)
                if not min_gc <= gc_fraction <= max_gc:
                    continue
                if _contains_excluded(guide_sequence, excluded):
                    continue
                if len(hits) >= max_matches:
                    truncated = True
                    break
                guide_symbol, guide_coordinate, guide_wrap = _location(
                    item=item,
                    oriented_start=guide_start % len(oriented),
                    length=resolved_guide_length,
                    strand=resolved_strand,
                    circular=circular,
                )
                pam_symbol, pam_coordinate, pam_wrap = _location(
                    item=item,
                    oriented_start=pam_start,
                    length=pam_length,
                    strand=resolved_strand,
                    circular=circular,
                )
                hits.append(
                    GuideCandidate(
                        nuclease=definition.name,
                        strand=resolved_strand,
                        guide_sequence=guide_sequence,
                        pam_sequence=pam_sequence,
                        guide_symbol_location=guide_symbol,
                        pam_symbol_location=pam_symbol,
                        guide_coordinate_location=guide_coordinate,
                        pam_coordinate_location=pam_coordinate,
                        gc_fraction=gc_fraction,
                        wraps_origin=guide_wrap or pam_wrap,
                    )
                )
            if truncated:
                break
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="crispr_pam_candidate_scan",
        method="iupac_pam_with_adjacent_guide",
        algorithm_version="1.0",
        parameters={
            "nuclease": definition.name,
            "pam": definition.pam,
            "pam_side": definition.pam_side,
            "guide_length": resolved_guide_length,
            "strand": str(strand),
            "min_gc": float(min_gc),
            "max_gc": float(max_gc),
            "exclude_motifs": excluded,
            "allow_ambiguous_guides": allow_ambiguous_guides,
            "source": definition.source,
            "catalog_version": PAM_CATALOG_VERSION,
            "efficiency_prediction": False,
            "off_target_prediction": False,
            "max_scan_length": max_scan_length,
            "max_pam_length": max_pam_length,
            "estimated_scan_cells": estimated_cells,
            "max_scan_cells": max_scan_cells,
            "max_exclude_motifs": max_exclude_motifs,
            "estimated_filter_cells": estimated_filter_cells,
            "max_filter_cells": max_filter_cells,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="DNAKit common PAM rule definitions",
            reference_version=PAM_CATALOG_VERSION,
        ),
    )
    return cast(PatternResult[GuideCandidate], result)


__all__ = [
    "BUILTIN_PAM_RULES",
    "PAM_CATALOG_VERSION",
    "PAMRule",
    "scan_pam_candidates",
]
