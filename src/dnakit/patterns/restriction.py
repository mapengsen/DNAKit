"""Small audited restriction-enzyme catalog and recognition-site scanning."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import cast

from dnakit.core import CompoundLocation, Interval, Strand, UnresolvedLocation
from dnakit.exceptions import ConfigurationError
from dnakit.patterns._shared import (
    SequenceInput,
    build_result,
    pattern_provenance,
    resolve_sequence,
    validate_iupac_text,
    validate_positive_int,
)
from dnakit.patterns.motif import scan_motif
from dnakit.patterns.results import PatternResult, RestrictionSiteHit

RESTRICTION_CATALOG_VERSION = "common-enzymes-v1"

_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)


@dataclass(frozen=True, init=False)
class RestrictionEnzyme:
    """Recognition sequence and 0-based cleavage offsets from its forward start."""

    name: str
    recognition_sequence: str
    top_cut_offset: int
    bottom_cut_offset: int
    source: str

    def __init__(
        self,
        name: str,
        recognition_sequence: str,
        top_cut_offset: int,
        bottom_cut_offset: int,
        *,
        source: str = "user-provided",
    ) -> None:
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("Restriction enzyme name must be non-empty.")
        site = validate_iupac_text(recognition_sequence, "recognition_sequence")
        for value, field_name in (
            (top_cut_offset, "top_cut_offset"),
            (bottom_cut_offset, "bottom_cut_offset"),
        ):
            if isinstance(value, bool) or not isinstance(value, int):
                raise ConfigurationError(f"{field_name} must be an integer.")
            if value < 0 or value > len(site):
                raise ConfigurationError(
                    f"{field_name} must lie inside the recognition-site boundary.",
                    hint=(
                        "Use an external enzyme definition adapter for Type IIS cuts "
                        "outside the site."
                    ),
                )
        if not isinstance(source, str) or not source.strip():
            raise ConfigurationError("Restriction enzyme source must be non-empty.")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "recognition_sequence", site)
        object.__setattr__(self, "top_cut_offset", top_cut_offset)
        object.__setattr__(self, "bottom_cut_offset", bottom_cut_offset)
        object.__setattr__(self, "source", source)


_SOURCE = "DNAKit common-enzyme boundary catalog; verify release work against REBASE/vendor data"
BUILTIN_RESTRICTION_ENZYMES: Mapping[str, RestrictionEnzyme] = {
    enzyme.name: enzyme
    for enzyme in (
        RestrictionEnzyme("BamHI", "GGATCC", 1, 5, source=_SOURCE),
        RestrictionEnzyme("EcoRI", "GAATTC", 1, 5, source=_SOURCE),
        RestrictionEnzyme("HaeIII", "GGCC", 2, 2, source=_SOURCE),
        RestrictionEnzyme("HindIII", "AAGCTT", 1, 5, source=_SOURCE),
        RestrictionEnzyme("NotI", "GCGGCCGC", 2, 6, source=_SOURCE),
        RestrictionEnzyme("SmaI", "CCCGGG", 3, 3, source=_SOURCE),
    )
}


def _resolve_enzymes(
    enzymes: Iterable[str | RestrictionEnzyme],
    *,
    max_enzymes: int,
) -> tuple[RestrictionEnzyme, ...]:
    if isinstance(enzymes, (str, bytes)):
        raise ConfigurationError("enzymes must be an iterable of names or definitions.")
    resolved: list[RestrictionEnzyme] = []
    for item in enzymes:
        if len(resolved) >= max_enzymes:
            raise ConfigurationError("Restriction panel exceeds max_enzymes.")
        if isinstance(item, RestrictionEnzyme):
            resolved.append(item)
        elif isinstance(item, str):
            try:
                resolved.append(BUILTIN_RESTRICTION_ENZYMES[item])
            except KeyError as exc:
                raise ConfigurationError(
                    "Restriction enzyme is absent from the small built-in catalog.",
                    context={"enzyme": item, "available": sorted(BUILTIN_RESTRICTION_ENZYMES)},
                    hint=(
                        "Pass a RestrictionEnzyme definition or use a future "
                        "Biopython/REBASE adapter."
                    ),
                ) from exc
        else:
            raise ConfigurationError("Each enzyme must be a name or RestrictionEnzyme.")
    if not resolved:
        raise ConfigurationError("At least one restriction enzyme is required.")
    names = [enzyme.name for enzyme in resolved]
    if len(names) != len(set(names)):
        raise ConfigurationError("Restriction enzyme names must be unique in one scan.")
    return tuple(resolved)


def _location_start(location: object) -> int | None:
    if isinstance(location, Interval):
        return location.start
    if isinstance(location, CompoundLocation):
        return location.parts[0].start
    if isinstance(location, UnresolvedLocation):
        return None
    raise AssertionError("Unexpected restriction-site location.")


def scan_restriction_sites(
    value: SequenceInput,
    enzymes: Iterable[str | RestrictionEnzyme],
    *,
    max_matches: int = 100_000,
    max_scan_length: int = 1_000_000,
    max_scan_cells: int = 50_000_000,
    max_enzymes: int = 10_000,
) -> PatternResult[RestrictionSiteHit]:
    """Scan a small built-in or user-provided restriction-enzyme panel.

    This is not a complete REBASE snapshot. The result records the catalog
    boundary and every user/built-in definition source.
    """

    sequence, sequence_id = resolve_sequence(value)
    validate_positive_int(max_matches, "max_matches")
    validate_positive_int(max_scan_length, "max_scan_length")
    validate_positive_int(max_scan_cells, "max_scan_cells")
    validate_positive_int(max_enzymes, "max_enzymes")
    definitions = _resolve_enzymes(enzymes, max_enzymes=max_enzymes)
    estimated_total_cells = (
        2 * sequence.symbol_length * sum(len(enzyme.recognition_sequence) for enzyme in definitions)
    )
    if estimated_total_cells > max_scan_cells:
        raise ConfigurationError(
            "Restriction panel exceeds max_scan_cells.",
            context={
                "estimated_total_cells": estimated_total_cells,
                "max_scan_cells": max_scan_cells,
            },
        )
    hits: list[RestrictionSiteHit] = []
    truncated = False
    for enzyme in definitions:
        remaining = max_matches - len(hits)
        if remaining <= 0:
            truncated = True
            break
        motif_result = scan_motif(
            sequence,
            enzyme.recognition_sequence,
            mode="iupac",
            name=enzyme.name,
            strand=Strand.BOTH,
            merge_strands=(
                enzyme.recognition_sequence
                == enzyme.recognition_sequence.translate(_COMPLEMENT)[::-1]
            ),
            max_matches=remaining,
            max_scan_length=max_scan_length,
            max_scan_cells=max_scan_cells,
        )
        for motif_hit in motif_result.hits:
            symbol_start = _location_start(motif_hit.symbol_location)
            coordinate_start = _location_start(motif_hit.coordinate_location)
            if symbol_start is None:
                raise AssertionError("Symbol locations are always resolved.")
            site_length = len(enzyme.recognition_sequence)
            if motif_hit.strand is Strand.REVERSE:
                top_offset = site_length - enzyme.bottom_cut_offset
                bottom_offset = site_length - enzyme.top_cut_offset
            else:
                top_offset = enzyme.top_cut_offset
                bottom_offset = enzyme.bottom_cut_offset
            top_cut = None if coordinate_start is None else coordinate_start + top_offset
            bottom_cut = None if coordinate_start is None else coordinate_start + bottom_offset
            if sequence.topology.value == "circular":
                span = sequence.coordinate_span
                if span is None:
                    raise AssertionError("A scanned circular sequence has a known span.")
                if top_cut is not None:
                    top_cut %= span
                if bottom_cut is not None:
                    bottom_cut %= span
            hits.append(
                RestrictionSiteHit(
                    enzyme=enzyme.name,
                    recognition_sequence=motif_hit.matched_sequence,
                    strand=motif_hit.strand,
                    symbol_location=motif_hit.symbol_location,
                    coordinate_location=motif_hit.coordinate_location,
                    top_cut=top_cut,
                    bottom_cut=bottom_cut,
                    wraps_origin=motif_hit.wraps_origin,
                )
            )
        truncated = truncated or motif_result.truncated
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="restriction_site_scan",
        method="iupac_recognition_and_offset_map",
        algorithm_version="1.0",
        parameters={
            "enzymes": [
                {
                    "name": enzyme.name,
                    "recognition_sequence": enzyme.recognition_sequence,
                    "top_cut_offset": enzyme.top_cut_offset,
                    "bottom_cut_offset": enzyme.bottom_cut_offset,
                    "source": enzyme.source,
                }
                for enzyme in definitions
            ],
            "builtin_catalog_version": RESTRICTION_CATALOG_VERSION,
            "complete_rebase_catalog": False,
            "type_iis_outside_site_supported": False,
            "max_scan_length": max_scan_length,
            "estimated_total_scan_cells": estimated_total_cells,
            "max_scan_cells": max_scan_cells,
            "max_enzymes": max_enzymes,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="DNAKit common restriction-enzyme definitions",
            reference_version=RESTRICTION_CATALOG_VERSION,
        ),
    )
    return cast(PatternResult[RestrictionSiteHit], result)


__all__ = [
    "BUILTIN_RESTRICTION_ENZYMES",
    "RESTRICTION_CATALOG_VERSION",
    "RestrictionEnzyme",
    "scan_restriction_sites",
]
