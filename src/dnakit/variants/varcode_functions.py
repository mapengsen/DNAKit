"""Offline varcode wrappers: splice-site scans, SV breakend parsing, genotype and allele helpers.

Every function here runs on local inputs and needs no reference genome, pyensembl cache or
network access. Effect prediction, transcript-aware germline phasing and reference-based
left-alignment are not wrapped because varcode needs Ensembl or pyensembl data for them.

Coordinates are 1-based and inclusive unless a function says otherwise. Packages are imported
only when a function runs; a missing ``varcode`` raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-varcode]"`` hint.
"""

from __future__ import annotations

import dataclasses
import re
from collections.abc import Mapping, Sequence
from importlib import metadata
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import ConfigurationError

_VARCODE_EXTRA = "bio-varcode"
_MAX_TEXT = 10_000
_MAX_SEQUENCE = 1_000_000
_MAX_POSITION = 2_000_000_000
_MAX_ROWS = 100_000
_BASES = re.compile(r"[ACGTN]*")
_GENOTYPE = re.compile(r"(\.|\d+)([/|](\.|\d+))*")
_SPLICE_WINDOWS = {"donor": 9, "acceptor": 4}
_SAMPLE_KINDS = ("germline", "somatic", "mixed")


def _text(value: object, name: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str) or (not value and not allow_empty):
        raise ConfigurationError(
            f"{name} must be a {'string' if allow_empty else 'non-empty string'}."
        )
    if len(value) > _MAX_TEXT or any(
        character.isspace() or character == "\x00" for character in value
    ):
        raise ConfigurationError(
            f"{name} must be at most {_MAX_TEXT} characters without whitespace."
        )
    return value


def _bases(value: object, name: str, *, allow_empty: bool = True) -> str:
    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a string of A, C, G, T or N.")
    text = value.upper()
    if (not text and not allow_empty) or len(text) > _MAX_SEQUENCE or not _BASES.fullmatch(text):
        raise ConfigurationError(f"{name} must contain only A, C, G, T or N.")
    return text


def _int(value: object, name: str, *, minimum: int = 0, maximum: int = _MAX_POSITION) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _unit(value: object, name: str, *, low: float = 0.0, high: float = 1.0) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be a number in [{low}, {high}].")
    return float(value)


def _mapping(value: object, name: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{name} must be a mapping.")
    return dict(value)


def _provenance() -> dict[str, str]:
    return {
        "package": f"varcode {metadata.version('varcode')}",
        "reference": "none: offline inputs only",
    }


def _variant_key(variant: Any) -> dict[str, Any]:
    return {
        "contig": str(variant.contig),
        "start": int(variant.start),
        "ref": str(variant.ref),
        "alt": str(variant.alt),
    }


def _variant(record: object, name: str) -> Any:
    """Build a varcode Variant from a mapping with contig, start, ref and alt."""

    item = _mapping(record, name)
    missing = sorted({"contig", "start", "ref", "alt"} - item.keys())
    if missing:
        raise ConfigurationError(f"{name} is missing fields: {', '.join(missing)}.")
    contig = _text(item["contig"], f"{name}.contig")
    start = _int(item["start"], f"{name}.start", minimum=1)
    ref = _bases(item["ref"], f"{name}.ref", allow_empty=False)
    alt = _bases(item["alt"], f"{name}.alt", allow_empty=True)
    module = require("varcode", _VARCODE_EXTRA)
    try:
        return module.Variant(contig, start, ref, alt)
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"{name} is not a valid variant: {exc}") from exc


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except ValueError as exc:
        raise ConfigurationError(str(exc)) from exc


def enumerate_cryptic_exon_candidates(
    contig: str,
    sequence: str,
    sequence_start: int = 1,
    *,
    donor_threshold: float = 0.7,
    acceptor_threshold: float = 0.7,
    min_intron_length: int = 30,
    max_intron_length: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan a local sequence for (acceptor, donor) pairs that could form a cryptic internal exon."""

    name = _text(contig, "contig")
    bases = _bases(sequence, "sequence", allow_empty=False)
    origin = _int(sequence_start, "sequence_start", minimum=1)
    donor = _unit(donor_threshold, "donor_threshold")
    acceptor = _unit(acceptor_threshold, "acceptor_threshold")
    minimum = _int(min_intron_length, "min_intron_length", maximum=_MAX_SEQUENCE)
    maximum = _int(max_intron_length, "max_intron_length", maximum=_MAX_SEQUENCE)
    if maximum < minimum:
        raise ConfigurationError("max_intron_length must be at least min_intron_length.")
    module = require("varcode.cryptic_exons", _VARCODE_EXTRA)
    found = module.enumerate_candidates(
        name,
        bases,
        origin,
        variant=None,
        donor_threshold=donor,
        acceptor_threshold=acceptor,
        min_intron_length=minimum,
        max_intron_length=maximum,
    )
    rows = [
        {
            "contig": str(item.contig),
            "interval_start": int(item.interval_start),
            "interval_end": int(item.interval_end),
            "donor_score": round(float(item.donor_score), 6),
            "acceptor_score": round(float(item.acceptor_score), 6),
        }
        for item in found
    ]
    rows.sort(
        key=lambda row: (
            -(row["donor_score"] + row["acceptor_score"]),
            row["interval_start"],
            row["interval_end"],
        )
    )
    return {
        "status": "completed",
        "operation": "enumerate_cryptic_exon_candidates",
        "sequence_length": len(bases),
        "count": len(rows),
        "candidates": rows,
        "provenance": _provenance(),
    }


def score_splice_site(window: str, *, site: str = "donor", verbose: bool = True) -> dict[str, Any]:
    """Score a donor (9 bp) or acceptor (4 bp) window with varcode's consensus motif model."""

    if site not in _SPLICE_WINDOWS:
        raise ConfigurationError("site must be 'donor' or 'acceptor'.")
    bases = _bases(window, "window", allow_empty=False)
    expected = _SPLICE_WINDOWS[site]
    if len(bases) != expected:
        raise ConfigurationError(f"A {site} window must be {expected} bp long.")
    module = require("varcode.cryptic_exons", _VARCODE_EXTRA)
    scorer = module.score_donor if site == "donor" else module.score_acceptor
    return {
        "status": "completed",
        "operation": "score_splice_site",
        "site": site,
        "window": bases,
        "window_length": expected,
        "score": round(float(scorer(bases)), 6),
        "provenance": _provenance(),
    }


def parse_symbolic_sv_alt(
    contig: str,
    start: int,
    ref: str,
    alt: str,
    *,
    info: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse one symbolic (<DEL>, <DUP>, <INV>, <CNV>, <INS>) or breakend ALT into SV fields."""

    name = _text(contig, "contig")
    position = _int(start, "start", minimum=1)
    reference = _bases(ref, "ref", allow_empty=False)
    allele = _text(alt, "alt")
    extra = None if info is None else _mapping(info, "info")
    module = require("varcode.sv_allele_parser", _VARCODE_EXTRA)
    parsed = _call(module.parse_symbolic_alt, name, position, reference, allele, info=extra)
    if parsed is None:
        return {
            "status": "completed",
            "operation": "parse_symbolic_sv_alt",
            "symbolic": False,
            "alt": allele,
            "provenance": _provenance(),
        }
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "parse_symbolic_sv_alt",
        "symbolic": True,
        "contig": str(parsed.contig),
        "start": int(parsed.start),
        "end": int(parsed.end),
        "affected_start": None if parsed.affected_start is None else int(parsed.affected_start),
        "affected_end": None if parsed.affected_end is None else int(parsed.affected_end),
        "sv_type": str(parsed.sv_type),
        "alt": str(parsed.alt),
        "mate_contig": getattr(parsed, "mate_contig", None),
        "mate_start": getattr(parsed, "mate_start", None),
        "provenance": _provenance(),
    }
    return result


def breakend_alt_sides(alt: str, *, verbose: bool = True) -> dict[str, Any]:
    """Report which side of each breakpoint a VCF breakend ALT keeps."""

    allele = _text(alt, "alt")
    module = require("varcode.sv_allele_parser", _VARCODE_EXTRA)
    sides = module.breakend_sides(allele)
    return {
        "status": "completed",
        "operation": "breakend_alt_sides",
        "alt": allele,
        "valid": sides is not None,
        "this_side": None if sides is None else sides[0],
        "mate_side": None if sides is None else sides[1],
        "provenance": _provenance(),
    }


def breakend_mate_position(alt: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the mate contig, mate position and bracket orientation named by a breakend ALT."""

    allele = _text(alt, "alt")
    module = require("varcode.sv_allele_parser", _VARCODE_EXTRA)
    mate = module.breakend_mate(allele)
    return {
        "status": "completed",
        "operation": "breakend_mate_position",
        "alt": allele,
        "valid": mate is not None,
        "mate_contig": None if mate is None else mate[0],
        "mate_position": None if mate is None else mate[1],
        "mate_orientation": None if mate is None else mate[2],
        "provenance": _provenance(),
    }


def breakend_insert_bases(ref: str, alt: str, *, verbose: bool = True) -> dict[str, Any]:
    """Return the bases a breakend ALT inserts next to a one-base REF anchor."""

    anchor = _bases(ref, "ref", allow_empty=False)
    allele = _text(alt, "alt")
    module = require("varcode.sv_allele_parser", _VARCODE_EXTRA)
    inserted = module.breakend_inserted_sequence(anchor, allele)
    return {
        "status": "completed",
        "operation": "breakend_insert_bases",
        "ref": anchor,
        "alt": allele,
        "known": inserted is not None,
        "inserted_bases": inserted,
        "provenance": _provenance(),
    }


def classify_breakend_pair(
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Type two paired breakend records as DEL, DUP or INV when they describe one event."""

    module = require("varcode.sv_allele_parser", _VARCODE_EXTRA)
    structural = require("varcode.structural_variant", _VARCODE_EXTRA)
    sides: list[Any] = []
    for label, record in (("first", first), ("second", second)):
        item = _mapping(record, label)
        missing = sorted({"contig", "start", "ref", "alt"} - item.keys())
        if missing:
            raise ConfigurationError(f"{label} is missing fields: {', '.join(missing)}.")
        info = item.get("info")
        parsed = _call(
            module.parse_symbolic_alt,
            _text(item["contig"], f"{label}.contig"),
            _int(item["start"], f"{label}.start", minimum=1),
            _bases(item["ref"], f"{label}.ref", allow_empty=False),
            _text(item["alt"], f"{label}.alt"),
            info=None if info is None else _mapping(info, f"{label}.info"),
        )
        if parsed is None or getattr(parsed, "mate_contig", None) is None:
            raise ConfigurationError(f"{label} must be a breakend ALT with a mate position.")
        sides.append(parsed)
    typed = structural.typed_event_from_breakends(sides[0], sides[1])
    if typed is None:
        return {
            "status": "completed",
            "operation": "classify_breakend_pair",
            "typed": False,
            "sv_type": None,
            "start": None,
            "end": None,
            "provenance": _provenance(),
        }
    sv_type, start, end = typed
    return {
        "status": "completed",
        "operation": "classify_breakend_pair",
        "typed": True,
        "sv_type": str(sv_type),
        "start": int(start),
        "end": int(end),
        "provenance": _provenance(),
    }


def _info_text(info: Mapping[str, Any]) -> str:
    """Write a mapping as VCF INFO text: flags as bare keys, lists comma-joined."""

    parts = []
    for key, value in info.items():
        if value is True:
            parts.append(str(key))
        elif isinstance(value, (list, tuple)):
            parts.append(f"{key}={','.join(map(str, value))}")
        else:
            parts.append(f"{key}={value}")
    return ";".join(parts)


def compare_sv_calls(
    rows: Sequence[Mapping[str, Any]],
    *,
    max_distance: int = 100,
    verbose: bool = True,
) -> dict[str, Any]:
    """Group structural-variant calls whose breakpoints lie within a distance of each other."""

    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)) or not rows:
        raise ConfigurationError("rows must be a non-empty list of call records.")
    if len(rows) > _MAX_ROWS:
        raise ConfigurationError(f"rows must contain at most {_MAX_ROWS} records.")
    distance = _int(max_distance, "max_distance", maximum=_MAX_POSITION)
    records = []
    for index, row in enumerate(rows):
        item = dict(_mapping(row, "rows item"))
        info = item.get("info")
        if isinstance(info, Mapping):
            # varcode parses VCF INFO text; a mapping is written back in that form.
            item["info"] = _info_text(info)
        elif info is not None and not isinstance(info, str):
            raise ConfigurationError(f"rows[{index}].info must be a mapping or VCF INFO text.")
        records.append(item)
    module = require("varcode.sv_comparison", _VARCODE_EXTRA)
    result = _call(module.compare_sv_calls, records, max_distance=distance)
    return {
        "status": "completed",
        "operation": "compare_sv_calls",
        "max_distance": distance,
        "call_count": len(records),
        "result": result,
        "provenance": _provenance(),
    }


def trim_shared_allele_flanks(ref: str, alt: str, *, verbose: bool = True) -> dict[str, Any]:
    """Remove the shared prefix and suffix of two alleles, as varcode does before classifying."""

    first = _bases(ref, "ref", allow_empty=True)
    second = _bases(alt, "alt", allow_empty=True)
    module = require("varcode.string_helpers", _VARCODE_EXTRA)
    trimmed_ref, trimmed_alt, prefix, suffix = module.trim_shared_flanking_strings(first, second)
    return {
        "status": "completed",
        "operation": "trim_shared_allele_flanks",
        "ref": first,
        "alt": second,
        "trimmed_ref": str(trimmed_ref),
        "trimmed_alt": str(trimmed_alt),
        "prefix": str(prefix),
        "suffix": str(suffix),
        "provenance": _provenance(),
    }


def parse_vcf_genotype(genotype: str, *, verbose: bool = True) -> dict[str, Any]:
    """Parse a VCF GT string such as ``0|1`` or ``./.`` into allele indices and phase."""

    text = _text(genotype, "genotype")
    if not _GENOTYPE.fullmatch(text):
        raise ConfigurationError(f"genotype {text!r} is not a VCF GT string.")
    module = require("varcode.genotype", _VARCODE_EXTRA)
    alleles, phased = _call(module.parse_gt_string, text)
    return {
        "status": "completed",
        "operation": "parse_vcf_genotype",
        "genotype": text,
        "alleles": [None if allele is None else int(allele) for allele in alleles],
        "phased": bool(phased),
        "ploidy": len(alleles),
        "provenance": _provenance(),
    }


def normalize_nucleotide_allele(
    allele: str,
    *,
    allow_extended_nucleotides: bool = False,
    empty_chars: str = ".-",
    treat_nan_as_empty: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Normalize an allele string: uppercase it and turn gap characters into an empty allele."""

    text = _text(allele, "allele", allow_empty=True)
    if not isinstance(allow_extended_nucleotides, bool) or not isinstance(treat_nan_as_empty, bool):
        raise ConfigurationError(
            "allow_extended_nucleotides and treat_nan_as_empty must be booleans."
        )
    if not isinstance(empty_chars, str) or len(empty_chars) > 8:
        raise ConfigurationError("empty_chars must be a short string of gap characters.")
    module = require("varcode.nucleotides", _VARCODE_EXTRA)
    normalized = _call(
        module.normalize_nucleotide_string,
        text,
        allow_extended_nucleotides=allow_extended_nucleotides,
        empty_chars=empty_chars,
        treat_nan_as_empty=treat_nan_as_empty,
    )
    return {
        "status": "completed",
        "operation": "normalize_nucleotide_allele",
        "input": text,
        "normalized": str(normalized),
        "is_empty": normalized == "",
        "provenance": _provenance(),
    }


def variant_overlaps_interval(
    variant_start: int,
    n_ref_bases: int,
    interval_start: int,
    interval_end: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Test whether a variant overlaps a 1-based inclusive interval on one contig."""

    start = _int(variant_start, "variant_start", minimum=1)
    width = _int(n_ref_bases, "n_ref_bases", maximum=_MAX_POSITION)
    left = _int(interval_start, "interval_start", minimum=1)
    right = _int(interval_end, "interval_end", minimum=1)
    if right < left:
        raise ConfigurationError("interval_end must be at least interval_start.")
    module = require("varcode.effects.effect_helpers", _VARCODE_EXTRA)
    # The upstream helper ends a variant at start + n_ref_bases, one base past its inclusive end.
    # Passing n_ref_bases - 1 gives the inclusive test; an SNV (width 1) maps to the insertion
    # branch, which tests start <= position <= end and is the same condition for one base.
    upstream_width = width - 1 if width else 0
    overlaps = bool(module.variant_overlaps_interval(start, upstream_width, left, right))
    return {
        "status": "completed",
        "operation": "variant_overlaps_interval",
        "variant_start": start,
        "n_ref_bases": width,
        "variant_end": start + width - 1 if width else None,
        "interval_start": left,
        "interval_end": right,
        "overlaps": overlaps,
        "provenance": _provenance(),
    }


def four_gamete_phase_test(
    both: int,
    first: int,
    second: int,
    neither: int,
    *,
    min_fragments: int = 2,
    error_rate: float = 0.01,
    max_p_value: float = 0.05,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decide whether two variants are in cis or trans from counts of supporting read fragments."""

    counts = {
        "both": _int(both, "both", maximum=_MAX_POSITION),
        "first": _int(first, "first", maximum=_MAX_POSITION),
        "second": _int(second, "second", maximum=_MAX_POSITION),
        "neither": _int(neither, "neither", maximum=_MAX_POSITION),
    }
    fragments = _int(min_fragments, "min_fragments", minimum=1)
    error = _unit(error_rate, "error_rate", low=0.0, high=0.5)
    if error == 0.0:
        raise ConfigurationError("error_rate must be above 0.")
    threshold = _unit(max_p_value, "max_p_value", low=0.0, high=1.0)
    if threshold == 0.0:
        raise ConfigurationError("max_p_value must be above 0.")
    module = require("varcode.rna_read_phasing", _VARCODE_EXTRA)
    answer = module.four_gamete_phase(
        counts["both"],
        counts["first"],
        counts["second"],
        counts["neither"],
        min_fragments=fragments,
        error_rate=error,
        max_p_value=threshold,
    )
    relation = "unknown" if answer is None else ("cis" if answer else "trans")
    return {
        "status": "completed",
        "operation": "four_gamete_phase_test",
        **counts,
        "relation": relation,
        "in_cis": None if answer is None else bool(answer),
        "provenance": _provenance(),
    }


class _EvidenceResolver:
    """Phase resolver backed by explicit in-cis evidence; unknown pairs return None."""

    def __init__(self, pairs: dict[frozenset, bool]) -> None:
        self._pairs = pairs

    def in_cis(self, first: Any, second: Any) -> bool | None:
        key = frozenset({_key_text(first), _key_text(second)})
        return self._pairs.get(key)


def _key_text(variant: Any) -> str:
    return f"{variant.contig}:{int(variant.start)}:{variant.ref}>{variant.alt}"


def _evidence(records: Sequence[Mapping[str, Any]] | None) -> _EvidenceResolver:
    pairs: dict[frozenset, bool] = {}
    for index, record in enumerate(records or ()):
        item = _mapping(record, f"in_cis_evidence[{index}]")
        if "in_cis" not in item or not isinstance(item["in_cis"], bool):
            raise ConfigurationError(f"in_cis_evidence[{index}].in_cis must be True or False.")
        first = _variant(item.get("first"), f"in_cis_evidence[{index}].first")
        second = _variant(item.get("second"), f"in_cis_evidence[{index}].second")
        pairs[frozenset({_key_text(first), _key_text(second)})] = item["in_cis"]
    return _EvidenceResolver(pairs)


def _variant_list(records: Sequence[Mapping[str, Any]] | None, name: str) -> list[Any]:
    if records is None:
        return []
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        raise ConfigurationError(f"{name} must be a list of variant mappings.")
    return [_variant(item, f"{name}[{index}]") for index, item in enumerate(records)]


def detect_germline_allele_overlap(
    somatic: Mapping[str, Any],
    germline: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report which germline variants share the contig, position, REF and ALT of a somatic one."""

    target = _variant(somatic, "somatic")
    candidates = _variant_list(germline, "germline")
    module = require("varcode.germline", _VARCODE_EXTRA)
    matched = [
        index
        for index, candidate in enumerate(candidates)
        if module.detect_germline_overlap(target, [candidate])
    ]
    return {
        "status": "completed",
        "operation": "detect_germline_allele_overlap",
        "somatic": _variant_key(target),
        "germline_count": len(candidates),
        "overlap": bool(matched),
        "matched_indices": matched,
        "provenance": _provenance(),
    }


def partition_germline_by_phase(
    somatic: Mapping[str, Any],
    germline: Sequence[Mapping[str, Any]],
    *,
    homozygous: Sequence[Mapping[str, Any]] | None = None,
    cis_variants: Sequence[Mapping[str, Any]] | None = None,
    in_cis_evidence: Sequence[Mapping[str, Any]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split germline variants into cis, trans and linked groups relative to a somatic variant."""

    target = _variant(somatic, "somatic")
    rest = _variant_list(germline, "germline")
    module = require("varcode.germline", _VARCODE_EXTRA)
    partition = module.partition_germline_by_phase(
        target,
        rest,
        phase_resolver=_evidence(in_cis_evidence),
        homozygous=tuple(_variant_list(homozygous, "homozygous")),
        cis_variants=tuple(_variant_list(cis_variants, "cis_variants")),
    )

    def keys(items: Sequence[Any]) -> list[dict[str, Any]]:
        return sorted(
            (_variant_key(item) for item in items),
            key=lambda row: (row["contig"], row["start"], row["ref"], row["alt"]),
        )

    placed = {_key_text(item) for group in (partition.cis, partition.trans) for item in group}
    for first, second in partition.linked:
        placed.update(_key_text(item) for item in (*first, *second))
    unphased = [item for item in rest if _key_text(item) not in placed]
    return {
        "status": "completed",
        "operation": "partition_germline_by_phase",
        "somatic": _variant_key(target),
        "germline_count": len(rest),
        "implicit": bool(partition.implicit),
        "cis": keys(partition.cis),
        "trans": keys(partition.trans),
        "linked": sorted(
            (
                {"group_a": keys(first), "group_b": keys(second)}
                for first, second in partition.linked
            ),
            key=lambda row: repr(row),
        ),
        "unphased": keys(unphased),
        "provenance": _provenance(),
    }


def enumerate_phase_hypotheses(
    somatic: Mapping[str, Any],
    germline: Sequence[Mapping[str, Any]],
    *,
    in_cis_evidence: Sequence[Mapping[str, Any]] | None = None,
    max_hypotheses: int = 8,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the phase configurations of a somatic variant against germline variants."""

    target = _variant(somatic, "somatic")
    rest = _variant_list(germline, "germline")
    limit = _int(max_hypotheses, "max_hypotheses", minimum=1, maximum=64)
    module = require("varcode.germline", _VARCODE_EXTRA)
    try:
        hypotheses = module.enumerate_phase_hypotheses(
            target, rest, _evidence(in_cis_evidence), max_hypotheses=limit
        )
    except Exception as exc:
        if type(exc).__name__ == "HypothesisLimitError":
            raise ConfigurationError(str(exc)) from exc
        raise
    rows = [
        {
            "phase_state": str(item.phase_state),
            "haplotype": str(item.haplotype),
            "cis": sorted((_variant_key(v) for v in item.cis), key=repr),
            "trans": sorted((_variant_key(v) for v in item.trans), key=repr),
        }
        for item in hypotheses
    ]
    rows.sort(key=lambda row: repr(row))
    return {
        "status": "completed",
        "operation": "enumerate_phase_hypotheses",
        "somatic": _variant_key(target),
        "count": len(rows),
        "hypotheses": rows,
        "provenance": _provenance(),
    }


def load_maf_variants(path: str, *, verbose: bool = True) -> dict[str, Any]:
    """Load a local TCGA-style MAF file as variant records with contig, start, ref and alt."""

    source = input_path(path)
    module = require("varcode.maf", _VARCODE_EXTRA)
    collection = _call(module.load_maf, str(source))
    variants = list(collection)
    references = sorted(
        {str(variant.reference_name) for variant in variants if variant.reference_name}
    )
    rows = []
    for variant in variants:
        row = _variant_key(variant)
        row["reference_name"] = (
            None if variant.reference_name is None else str(variant.reference_name)
        )
        rows.append(row)
    return {
        "status": "completed",
        "operation": "load_maf_variants",
        "reference": references[0] if len(references) == 1 else None,
        "count": len(rows),
        "variants": rows,
        "provenance": _provenance(),
    }


def check_vcf_sample_identity(
    path: str,
    *,
    kind: str = "germline",
    assembly: str | None = None,
    config: Mapping[str, Any] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a local VCF as sample fingerprints and report whether the samples match each other."""

    source = input_path(path)
    if kind not in _SAMPLE_KINDS:
        raise ConfigurationError("kind must be 'germline', 'somatic' or 'mixed'.")
    if assembly is not None and not isinstance(assembly, str):
        raise ConfigurationError("assembly must be a string such as 'GRCh38'.")
    module = require("varcode.sample_identity", _VARCODE_EXTRA)
    settings = None
    if config is not None:
        allowed = {field.name for field in dataclasses.fields(module.SampleCheckConfig)}
        overrides = _mapping(config, "config")
        unknown = sorted(set(overrides) - allowed)
        if unknown:
            raise ConfigurationError(f"Unknown config field(s): {', '.join(unknown)}.")
        settings = _call(module.SampleCheckConfig, **overrides)
    samples = _call(module.load_vcf_samples, source, kind=kind, assembly=assembly, config=settings)
    report = _call(module.check_sample_identity, samples, config=settings)
    return {
        "status": "completed",
        "operation": "check_vcf_sample_identity",
        "path": str(source),
        "kind": kind,
        "sample_count": len(samples),
        "report": report,
        "provenance": _provenance(),
    }


__all__ = [
    "breakend_alt_sides",
    "breakend_insert_bases",
    "breakend_mate_position",
    "check_vcf_sample_identity",
    "classify_breakend_pair",
    "compare_sv_calls",
    "detect_germline_allele_overlap",
    "enumerate_cryptic_exon_candidates",
    "enumerate_phase_hypotheses",
    "four_gamete_phase_test",
    "load_maf_variants",
    "normalize_nucleotide_allele",
    "parse_symbolic_sv_alt",
    "parse_vcf_genotype",
    "partition_germline_by_phase",
    "score_splice_site",
    "trim_shared_allele_flanks",
    "variant_overlaps_interval",
]

_register_public_verbose(globals(), tuple(__all__))
