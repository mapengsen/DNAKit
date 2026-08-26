"""Start/stop codon and six-frame ORF annotation."""

from __future__ import annotations

from collections.abc import Iterable
from itertools import islice
from typing import cast

from dnakit.core import Strand
from dnakit.exceptions import ConfigurationError
from dnakit.ops import translate
from dnakit.patterns._shared import (
    SequenceInput,
    build_result,
    coerce_strands,
    pattern_provenance,
    require_linear,
    resolve_sequence,
    reverse_complement_text,
    segment_location,
    segments,
    validate_bool,
    validate_positive_int,
)
from dnakit.patterns.results import CodonSite, ORFHit, PatternResult

_GENETIC_CODE_STARTS = {
    1: frozenset({"ATG"}),
    11: frozenset({"ATG", "GTG", "TTG"}),
}
_GENETIC_CODE_STOPS = {
    1: frozenset({"TAA", "TAG", "TGA"}),
    11: frozenset({"TAA", "TAG", "TGA"}),
}


def _codon_set(
    values: Iterable[str] | None,
    *,
    default: frozenset[str],
    name: str,
) -> frozenset[str]:
    if values is None:
        return default
    if isinstance(values, (str, bytes)):
        raise ConfigurationError(f"{name} must be an iterable of codon strings.")
    try:
        materialized = tuple(islice(iter(values), 65))
    except TypeError as exc:
        raise ConfigurationError(f"{name} must be an iterable of codon strings.") from exc
    if len(materialized) > 64:
        raise ConfigurationError(
            f"{name} cannot contain more than the 64 canonical triplets.",
            code="CODON_SET_SIZE_LIMIT",
        )
    resolved = frozenset(materialized)
    if not resolved:
        raise ConfigurationError(f"{name} must contain at least one codon.")
    if any(
        not isinstance(codon, str) or len(codon) != 3 or set(codon) - set("ACGT")
        for codon in resolved
    ):
        raise ConfigurationError(f"{name} must contain uppercase canonical DNA triplets.")
    return resolved


def _code_sets(
    genetic_code: int,
    start_codons: Iterable[str] | None,
    stop_codons: Iterable[str] | None,
) -> tuple[frozenset[str], frozenset[str]]:
    if (
        isinstance(genetic_code, bool)
        or not isinstance(genetic_code, int)
        or genetic_code not in _GENETIC_CODE_STARTS
    ):
        raise ConfigurationError(
            "Only NCBI genetic codes 1 and 11 are supported for ORF boundaries.",
            context={"genetic_code": genetic_code},
        )
    starts = _codon_set(
        start_codons,
        default=_GENETIC_CODE_STARTS[genetic_code],
        name="start_codons",
    )
    stops = _codon_set(
        stop_codons,
        default=_GENETIC_CODE_STOPS[genetic_code],
        name="stop_codons",
    )
    if starts & stops:
        raise ConfigurationError("start_codons and stop_codons must be disjoint.")
    return starts, stops


def _original_bounds(length: int, start: int, end: int, strand: Strand) -> tuple[int, int]:
    if strand is Strand.FORWARD:
        return start, end
    return length - end, length - start


def scan_codon_sites(
    value: SequenceInput,
    *,
    genetic_code: int = 1,
    start_codons: Iterable[str] | None = None,
    stop_codons: Iterable[str] | None = None,
    strand: Strand | str = Strand.BOTH,
    max_matches: int = 100_000,
    max_codon_checks: int = 3_000_000,
) -> PatternResult[CodonSite]:
    """Locate in-frame start and stop codons in the selected reading frames."""

    sequence, sequence_id = resolve_sequence(value)
    require_linear(sequence, "codon-site scanning")
    validate_positive_int(max_matches, "max_matches")
    validate_positive_int(max_codon_checks, "max_codon_checks")
    starts, stops = _code_sets(genetic_code, start_codons, stop_codons)
    selected_strands = coerce_strands(strand)
    estimated_checks = 2 * sequence.symbol_length * len(selected_strands)
    if estimated_checks > 3 * max_codon_checks:
        raise ConfigurationError(
            "Codon scan exceeds max_codon_checks.",
            context={"estimated_checks": estimated_checks // 3},
        )

    hits: list[CodonSite] = []
    checks = 0
    truncated = False
    for item in segments(sequence):
        for resolved_strand in selected_strands:
            oriented = (
                item.text
                if resolved_strand is Strand.FORWARD
                else reverse_complement_text(item.text)
            )
            for offset in range(3):
                signed_frame = offset + 1 if resolved_strand is Strand.FORWARD else -(offset + 1)
                for position in range(offset, len(oriented) - 2, 3):
                    checks += 1
                    if checks > max_codon_checks:
                        raise ConfigurationError("Codon scan exceeded max_codon_checks.")
                    codon = oriented[position : position + 3]
                    kind = "start" if codon in starts else "stop" if codon in stops else None
                    if kind is None:
                        continue
                    if len(hits) >= max_matches:
                        truncated = True
                        break
                    original_start, original_end = _original_bounds(
                        len(item.text), position, position + 3, resolved_strand
                    )
                    symbol_location, coordinate_location = segment_location(
                        item, original_start, original_end
                    )
                    hits.append(
                        CodonSite(
                            kind=kind,
                            codon=codon,
                            strand=resolved_strand,
                            frame=signed_frame,
                            symbol_location=symbol_location,
                            coordinate_location=coordinate_location,
                        )
                    )
                if truncated:
                    break
            if truncated:
                break
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="codon_site_scan",
        method="six_frame_triplet_scan",
        algorithm_version="1.0",
        parameters={
            "genetic_code": genetic_code,
            "start_codons": sorted(starts),
            "stop_codons": sorted(stops),
            "strand": str(strand),
            "codon_checks": checks,
            "max_codon_checks": max_codon_checks,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="NCBI genetic code",
            reference_version=str(genetic_code),
        ),
    )
    return cast(PatternResult[CodonSite], result)


def scan_orfs(
    value: SequenceInput,
    *,
    genetic_code: int = 1,
    start_codons: Iterable[str] | None = None,
    stop_codons: Iterable[str] | None = None,
    min_length: int = 0,
    require_complete: bool = True,
    strand: Strand | str = Strand.BOTH,
    max_matches: int = 100_000,
    max_codon_checks: int = 3_000_000,
) -> PatternResult[ORFHit]:
    """Find start-anchored ORFs on up to six frames.

    Every in-frame start is retained until the next stop, so nested starts yield
    distinct ORFs. In incomplete mode, an unclosed start extends to the final
    complete codon of its gap-delimited fragment. Circular inputs must first be
    linearized at a documented origin.
    """

    sequence, sequence_id = resolve_sequence(value)
    require_linear(sequence, "ORF scanning")
    if isinstance(min_length, bool) or not isinstance(min_length, int) or min_length < 0:
        raise ConfigurationError("min_length must be a non-negative integer.")
    validate_bool(require_complete, "require_complete")
    validate_positive_int(max_matches, "max_matches")
    validate_positive_int(max_codon_checks, "max_codon_checks")
    starts, stops = _code_sets(genetic_code, start_codons, stop_codons)
    selected_strands = coerce_strands(strand)

    hits: list[ORFHit] = []
    checks = 0
    truncated = False

    def add_orf(
        item: object,
        resolved_strand: Strand,
        signed_frame: int,
        oriented: str,
        start: int,
        end: int,
        stop_codon: str | None,
    ) -> bool:
        from dnakit.patterns._shared import Segment

        if not isinstance(item, Segment):
            raise AssertionError("Internal segment type mismatch.")
        nucleotide_length = end - start
        if nucleotide_length < min_length:
            return True
        if len(hits) >= max_matches:
            return False
        original_start, original_end = _original_bounds(len(item.text), start, end, resolved_strand)
        symbol_location, coordinate_location = segment_location(item, original_start, original_end)
        coding = oriented[start:end]
        translation_table = 1
        protein = translate(coding, table=translation_table)
        if genetic_code == 11 and protein and coding[:3] in starts:
            protein = "M" + protein[1:]
        hits.append(
            ORFHit(
                strand=resolved_strand,
                frame=signed_frame,
                symbol_location=symbol_location,
                coordinate_location=coordinate_location,
                nucleotide_length=nucleotide_length,
                codon_count=nucleotide_length // 3,
                start_codon=oriented[start : start + 3],
                stop_codon=stop_codon,
                complete=stop_codon is not None,
                translation=protein,
            )
        )
        return True

    for item in segments(sequence):
        for resolved_strand in selected_strands:
            oriented = (
                item.text
                if resolved_strand is Strand.FORWARD
                else reverse_complement_text(item.text)
            )
            for offset in range(3):
                signed_frame = offset + 1 if resolved_strand is Strand.FORWARD else -(offset + 1)
                open_starts: list[int] = []
                coding_end = offset + ((len(oriented) - offset) // 3) * 3
                for position in range(offset, coding_end, 3):
                    checks += 1
                    if checks > max_codon_checks:
                        raise ConfigurationError("ORF scan exceeded max_codon_checks.")
                    codon = oriented[position : position + 3]
                    if codon in starts:
                        open_starts.append(position)
                    if codon not in stops:
                        continue
                    for start in open_starts:
                        if not add_orf(
                            item,
                            resolved_strand,
                            signed_frame,
                            oriented,
                            start,
                            position + 3,
                            codon,
                        ):
                            truncated = True
                            break
                    open_starts.clear()
                    if truncated:
                        break
                if not truncated and not require_complete:
                    for start in open_starts:
                        if not add_orf(
                            item,
                            resolved_strand,
                            signed_frame,
                            oriented,
                            start,
                            coding_end,
                            None,
                        ):
                            truncated = True
                            break
                if truncated:
                    break
            if truncated:
                break
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="orf_scan",
        method="start_to_next_in_frame_stop",
        algorithm_version="1.0",
        parameters={
            "genetic_code": genetic_code,
            "translation_table": 1,
            "start_codons": sorted(starts),
            "stop_codons": sorted(stops),
            "min_length": min_length,
            "minimum_length_includes_terminal_stop": True,
            "require_complete": require_complete,
            "strand": str(strand),
            "nested_starts": "emit_each",
            "frame_origin": "each-gap-delimited-fragment",
            "terminal_stop_in_translation": True,
            "codon_checks": checks,
            "max_codon_checks": max_codon_checks,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="NCBI genetic code and ORF boundary rules",
            reference_version=str(genetic_code),
        ),
    )
    return cast(PatternResult[ORFHit], result)


__all__ = ["scan_codon_sites", "scan_orfs"]
