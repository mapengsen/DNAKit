"""Bounded exact, IUPAC, regular-expression, PWM, promoter, and TF motif scans."""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from itertools import islice
from typing import Literal, cast

from dnakit.core import Location, Strand
from dnakit.core._json import to_json_compatible
from dnakit.exceptions import ConfigurationError, SequenceError
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
from dnakit.patterns.results import MotifHit, PatternResult, validate_finite_score

MotifMode = Literal["exact", "iupac", "regex"]

PROMOTER_MOTIFS: Mapping[str, str] = {
    "eukaryotic_TATA_box_consensus": "TATAWAWR",
    "bacterial_minus_10_consensus": "TATAAT",
    "bacterial_minus_35_consensus": "TTGACA",
}
PROMOTER_CATALOG_VERSION = "builtin-boundary-v1"
MAX_PWM_LENGTH = 100_000


@dataclass(frozen=True, init=False)
class PWM:
    """A position-weight input matrix normalized independently by column."""

    name: str
    a: tuple[float, ...]
    c: tuple[float, ...]
    g: tuple[float, ...]
    t: tuple[float, ...]

    def __init__(self, name: str, matrix: Mapping[str, Iterable[float]]) -> None:
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("PWM name must be a non-empty string.")
        if not isinstance(matrix, Mapping) or set(matrix) != set("ACGT"):
            raise ConfigurationError("PWM matrix must contain exactly A, C, G, and T rows.")
        rows: dict[str, tuple[float, ...]] = {}
        for base in "ACGT":
            raw = matrix[base]
            if isinstance(raw, (str, bytes)):
                raise ConfigurationError("PWM rows must be numerical sequences.")
            try:
                raw_values = tuple(islice(iter(raw), MAX_PWM_LENGTH + 1))
            except TypeError as exc:
                raise ConfigurationError("PWM rows must be numerical sequences.") from exc
            if len(raw_values) > MAX_PWM_LENGTH:
                raise ConfigurationError(
                    "PWM row exceeds the supported length limit.",
                    code="PWM_LENGTH_LIMIT",
                    context={"max_pwm_length": MAX_PWM_LENGTH},
                )
            if any(
                isinstance(value, bool) or not isinstance(value, (int, float))
                for value in raw_values
            ):
                raise ConfigurationError("PWM rows must contain real numbers, not booleans.")
            try:
                values = tuple(float(value) for value in raw_values)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ConfigurationError("PWM rows must contain numbers.") from exc
            if not values or any(not math.isfinite(value) or value < 0 for value in values):
                raise ConfigurationError("PWM values must be finite non-negative numbers.")
            rows[base] = values
        lengths = {len(row) for row in rows.values()}
        if len(lengths) != 1:
            raise ConfigurationError("All PWM rows must have the same non-zero length.")
        for column in range(len(rows["A"])):
            if math.fsum(rows[base][column] for base in "ACGT") <= 0:
                raise ConfigurationError("Every PWM column must have positive total weight.")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "a", rows["A"])
        object.__setattr__(self, "c", rows["C"])
        object.__setattr__(self, "g", rows["G"])
        object.__setattr__(self, "t", rows["T"])

    @property
    def length(self) -> int:
        return len(self.a)

    def row(self, base: str) -> tuple[float, ...]:
        if base == "A":
            return self.a
        if base == "C":
            return self.c
        if base == "G":
            return self.g
        if base == "T":
            return self.t
        raise ConfigurationError("PWM scores are defined only for canonical DNA bases.")

    def to_dict(self) -> dict[str, object]:
        return cast(dict[str, object], to_json_compatible(self))


def _fixed_matches(
    text: str,
    motif: str,
    *,
    mode: Literal["exact", "iupac"],
    overlapping: bool,
    circular: bool,
) -> Iterator[tuple[int, int, str]]:
    length = len(motif)
    if circular and length > len(text):
        raise ConfigurationError("A circular motif cannot be longer than the sequence.")
    searchable = text + text[: length - 1] if circular and length > 1 else text
    step = 1 if overlapping else length
    start = 0
    limit = len(text) if circular else len(text) - length + 1
    while start < max(0, limit):
        candidate = searchable[start : start + length]
        matched = candidate == motif
        if mode == "iupac":
            matched = len(candidate) == length and all(
                iupac_compatible(left, right) for left, right in zip(candidate, motif, strict=True)
            )
        if matched:
            yield start, start + length, candidate
            start += step
        else:
            start += 1


def _regex_matches(
    text: str,
    pattern: str,
    *,
    overlapping: bool,
) -> Iterator[tuple[int, int, str]]:
    try:
        compiled = re.compile(pattern)
    except (re.error, OverflowError) as exc:
        raise ConfigurationError(
            "Invalid regular-expression motif.", context={"regex_error": str(exc)}
        ) from exc
    position = 0
    while position <= len(text):
        match = compiled.search(text, position)
        if match is None:
            return
        if match.end() == match.start():
            raise ConfigurationError("Regular-expression motifs must not produce empty matches.")
        yield match.start(), match.end(), match.group(0)
        position = match.start() + 1 if overlapping else match.end()


def _validate_safe_regex(pattern: str, scan_length: int) -> int:
    """Reject constructs with unbounded backtracking or non-DNA semantics."""

    allowed = frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn.[]^-?*+{},0123456789$")
    invalid = sorted(set(pattern) - allowed)
    if invalid or any(token in pattern for token in ("(", ")", "|", "\\")):
        raise ConfigurationError(
            "Regex motif is outside DNAKit's safe DNA-regex subset.",
            context={"invalid_characters": invalid},
            hint="Use literals, dot, character classes, anchors, and simple quantifiers only.",
        )
    if re.search(r"(?:[*+?]|\{\d+(?:,\d*)?\})\s*(?:[*+?]|\{)", pattern):
        raise ConfigurationError("Nested or adjacent regex quantifiers are not allowed.")
    quantifier_text = re.sub(r"\[[^\]]*\]", "X", pattern)
    if "{," in quantifier_text:
        raise ConfigurationError("Regex repetition bounds must include an explicit lower bound.")
    unbounded = re.findall(r"[*+]|\{\d+,\}", quantifier_text)
    if len(unbounded) > 1:
        raise ConfigurationError(
            "The safe DNA-regex subset allows at most one unbounded quantifier."
        )
    quantifiers = re.findall(r"[*+?]|\{\d+(?:,\d*)?\}", quantifier_text)
    quantifier_count = len(quantifiers)
    if quantifier_count > 32:
        raise ConfigurationError("Regex motif contains too many quantifiers for bounded scanning.")
    multiplier = 1
    for quantifier in quantifiers:
        if quantifier in ("*", "+") or quantifier.endswith(",}"):
            factor = max(1, scan_length)
        elif quantifier == "?":
            factor = 2
        else:
            bounds = quantifier[1:-1].split(",")
            lower = int(bounds[0])
            upper = lower if len(bounds) == 1 else int(bounds[1])
            if upper < lower:
                raise ConfigurationError(
                    "Regex repetition upper bound cannot be below its lower bound."
                )
            factor = max(1, upper, upper - lower + 1)
        multiplier *= factor
        if multiplier > 50_000_000:
            raise ConfigurationError("Regex quantifier expansion exceeds the safety limit.")
    return multiplier


def _mapped_location(
    *,
    local_start: int,
    local_end: int,
    segment_text_length: int,
    segment: object,
    strand: Strand,
    circular: bool,
) -> tuple[Location, Location, bool]:
    from dnakit.patterns._shared import Segment

    if not isinstance(segment, Segment):
        raise AssertionError("Internal segment type mismatch.")
    length = local_end - local_start
    original_start = local_start if strand is Strand.FORWARD else segment_text_length - local_end
    if circular:
        return circular_location(original_start, length, segment_text_length)
    symbol_location, coordinate_location = segment_location(
        segment, original_start, original_start + length
    )
    return symbol_location, coordinate_location, False


def scan_motif(
    value: SequenceInput,
    motif: str,
    *,
    mode: MotifMode = "exact",
    name: str = "motif",
    strand: Strand | str = Strand.BOTH,
    overlapping: bool = True,
    merge_strands: bool = False,
    max_matches: int = 100_000,
    max_scan_length: int = 1_000_000,
    max_pattern_length: int = 1_000,
    max_scan_cells: int = 50_000_000,
) -> PatternResult[MotifHit]:
    """Scan literal, IUPAC-compatible, or Python-regex motifs on selected strands.

    Regex scanning is bounded by ``max_scan_length`` and ``max_pattern_length``;
    circular origin wrapping is supported only for fixed-length exact/IUPAC motifs.
    """

    sequence, sequence_id = resolve_sequence(value)
    validate_positive_int(max_matches, "max_matches")
    validate_positive_int(max_scan_length, "max_scan_length")
    validate_positive_int(max_pattern_length, "max_pattern_length")
    validate_positive_int(max_scan_cells, "max_scan_cells")
    validate_bool(overlapping, "overlapping")
    validate_bool(merge_strands, "merge_strands")
    if sequence.symbol_length > max_scan_length:
        raise ConfigurationError(
            "Sequence exceeds max_scan_length.",
            context={"symbol_length": sequence.symbol_length, "max_scan_length": max_scan_length},
        )
    if not isinstance(name, str) or not name.strip():
        raise ConfigurationError("Motif name must be non-empty.")
    raw_mode: object = mode
    if not isinstance(raw_mode, str) or raw_mode not in ("exact", "iupac", "regex"):
        raise ConfigurationError("Motif mode must be exact, iupac, or regex.")
    resolved_mode = cast(MotifMode, raw_mode)
    if not isinstance(motif, str) or not motif:
        raise ConfigurationError("motif must be a non-empty string.")
    if len(motif) > max_pattern_length:
        raise ConfigurationError("Motif exceeds max_pattern_length.")
    resolved_motif = motif if resolved_mode == "regex" else validate_iupac_text(motif, "motif")
    regex_work_multiplier = (
        _validate_safe_regex(resolved_motif, sequence.symbol_length)
        if resolved_mode == "regex"
        else 1
    )
    circular = sequence.topology.value == "circular"
    require_ungapped_circular(sequence, "motif scanning")
    if circular and resolved_mode == "regex":
        raise SequenceError(
            "Regex motifs do not have a fixed circular wrap length.",
            code="CIRCULAR_REGEX_UNSUPPORTED",
            hint="Use exact/IUPAC mode or linearize the sequence at a documented origin.",
        )
    selected_strands = coerce_strands(strand)
    estimated_cells = (
        sequence.symbol_length * len(resolved_motif) * len(selected_strands) * regex_work_multiplier
    )
    if estimated_cells > max_scan_cells:
        raise ConfigurationError(
            "Motif scan exceeds max_scan_cells.",
            context={"estimated_cells": estimated_cells, "max_scan_cells": max_scan_cells},
        )

    hits: list[MotifHit] = []
    truncated = False
    seen: dict[tuple[object, object, str], int] = {}
    for item in segments(sequence):
        for resolved_strand in selected_strands:
            oriented = (
                item.text
                if resolved_strand is Strand.FORWARD
                else reverse_complement_text(item.text)
            )
            iterator: Iterator[tuple[int, int, str]]
            if resolved_mode == "regex":
                iterator = _regex_matches(oriented, resolved_motif, overlapping=overlapping)
            else:
                iterator = _fixed_matches(
                    oriented,
                    resolved_motif,
                    mode=resolved_mode,
                    overlapping=overlapping,
                    circular=circular,
                )
            for start, end, matched in iterator:
                symbol_location, coordinate_location, wrapped = _mapped_location(
                    local_start=start,
                    local_end=end,
                    segment_text_length=len(item.text),
                    segment=item,
                    strand=resolved_strand,
                    circular=circular,
                )
                key = (symbol_location, coordinate_location, matched)
                if merge_strands and key in seen:
                    old_index = seen[key]
                    previous = hits[old_index]
                    hits[old_index] = MotifHit(
                        motif_name=previous.motif_name,
                        matched_sequence=previous.matched_sequence,
                        strand=Strand.BOTH,
                        symbol_location=previous.symbol_location,
                        coordinate_location=previous.coordinate_location,
                        score=None,
                        threshold=None,
                        wraps_origin=previous.wraps_origin,
                    )
                    continue
                if len(hits) >= max_matches:
                    truncated = True
                    break
                seen[key] = len(hits)
                hits.append(
                    MotifHit(
                        motif_name=name,
                        matched_sequence=matched,
                        strand=resolved_strand,
                        symbol_location=symbol_location,
                        coordinate_location=coordinate_location,
                        score=None,
                        threshold=None,
                        wraps_origin=wrapped,
                    )
                )
            if truncated:
                break
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="motif_scan",
        method=resolved_mode,
        algorithm_version="1.0",
        parameters={
            "motif": resolved_motif,
            "motif_name": name,
            "strand": str(strand),
            "overlapping": overlapping,
            "merge_strands": merge_strands,
            "max_scan_length": max_scan_length,
            "max_pattern_length": max_pattern_length,
            "estimated_scan_cells": estimated_cells,
            "max_scan_cells": max_scan_cells,
            "target_iupac_rule": "set-intersection" if resolved_mode == "iupac" else "literal",
            "regex_safety": "restricted-dna-subset" if resolved_mode == "regex" else None,
            "regex_work_multiplier": regex_work_multiplier if resolved_mode == "regex" else None,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(reimplementation=resolved_mode == "regex"),
    )
    return cast(PatternResult[MotifHit], result)


def _background(values: Mapping[str, float] | None) -> dict[str, float]:
    if values is None:
        return {base: 0.25 for base in "ACGT"}
    if set(values) != set("ACGT"):
        raise ConfigurationError("PWM background must contain exactly A, C, G, and T.")
    resolved: dict[str, float] = {}
    for base in "ACGT":
        value = values[base]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0
        ):
            raise ConfigurationError("PWM background probabilities must be finite and positive.")
        resolved[base] = float(value)
    total = math.fsum(resolved.values())
    return {base: resolved[base] / total for base in "ACGT"}


def _score_pwm(
    text: str,
    pwm: PWM,
    background: Mapping[str, float],
    pseudocount: float,
) -> float | None:
    if any(base not in "ACGT" for base in text):
        return None
    scores: list[float] = []
    for index, base in enumerate(text):
        column_total = math.fsum(pwm.row(item)[index] + pseudocount for item in "ACGT")
        probability = (pwm.row(base)[index] + pseudocount) / column_total
        if probability <= 0:
            return None
        scores.append(math.log2(probability / background[base]))
    return math.fsum(scores)


def scan_pwm(
    value: SequenceInput,
    pwm: PWM,
    *,
    threshold: float,
    background: Mapping[str, float] | None = None,
    pseudocount: float = 0.0,
    strand: Strand | str = Strand.BOTH,
    max_matches: int = 100_000,
    max_scan_length: int = 1_000_000,
    max_pwm_length: int = 10_000,
    max_score_cells: int = 50_000_000,
) -> PatternResult[MotifHit]:
    """Scan a user-provided PWM with explicit log2-odds threshold and background."""

    sequence, sequence_id = resolve_sequence(value)
    if not isinstance(pwm, PWM):
        raise ConfigurationError("pwm must be a PWM object.")
    validate_finite_score(threshold, "threshold")
    validate_finite_score(pseudocount, "pseudocount")
    if threshold is None:
        raise ConfigurationError("threshold must be a finite number.")
    if pseudocount is None:
        raise ConfigurationError("pseudocount must be a finite non-negative number.")
    if pseudocount < 0:
        raise ConfigurationError("pseudocount must be non-negative.")
    validate_positive_int(max_matches, "max_matches")
    validate_positive_int(max_scan_length, "max_scan_length")
    validate_positive_int(max_pwm_length, "max_pwm_length")
    validate_positive_int(max_score_cells, "max_score_cells")
    if sequence.symbol_length > max_scan_length:
        raise ConfigurationError("Sequence exceeds max_scan_length.")
    if pwm.length > max_pwm_length:
        raise ConfigurationError("PWM exceeds max_pwm_length.")
    resolved_background = _background(background)
    circular = sequence.topology.value == "circular"
    require_ungapped_circular(sequence, "PWM scanning")
    if circular and pwm.length > sequence.symbol_length:
        raise ConfigurationError("A circular PWM cannot be longer than the sequence.")
    selected_strands = coerce_strands(strand)
    possible_windows = sum(
        len(item.text) if circular else max(0, len(item.text) - pwm.length + 1)
        for item in segments(sequence)
    ) * len(selected_strands)
    if possible_windows * pwm.length > max_score_cells:
        raise ConfigurationError(
            "PWM scan exceeds max_score_cells.",
            context={
                "possible_windows": possible_windows,
                "pwm_length": pwm.length,
                "max_score_cells": max_score_cells,
            },
        )

    hits: list[MotifHit] = []
    truncated = False
    for item in segments(sequence):
        for resolved_strand in selected_strands:
            oriented = (
                item.text
                if resolved_strand is Strand.FORWARD
                else reverse_complement_text(item.text)
            )
            searchable = (
                oriented + oriented[: pwm.length - 1] if circular and pwm.length > 1 else oriented
            )
            limit = len(oriented) if circular else len(oriented) - pwm.length + 1
            for start in range(max(0, limit)):
                candidate = searchable[start : start + pwm.length]
                score = _score_pwm(candidate, pwm, resolved_background, float(pseudocount))
                if score is None or score < threshold:
                    continue
                if len(hits) >= max_matches:
                    truncated = True
                    break
                symbol_location, coordinate_location, wrapped = _mapped_location(
                    local_start=start,
                    local_end=start + pwm.length,
                    segment_text_length=len(item.text),
                    segment=item,
                    strand=resolved_strand,
                    circular=circular,
                )
                hits.append(
                    MotifHit(
                        motif_name=pwm.name,
                        matched_sequence=candidate,
                        strand=resolved_strand,
                        symbol_location=symbol_location,
                        coordinate_location=coordinate_location,
                        score=score,
                        threshold=float(threshold),
                        wraps_origin=wrapped,
                    )
                )
            if truncated:
                break
        if truncated:
            break
    result = build_result(
        sequence,
        sequence_id,
        name="pwm_scan",
        method="log2_odds",
        algorithm_version="1.0",
        parameters={
            "pwm_name": pwm.name,
            "pwm_length": pwm.length,
            "threshold": float(threshold),
            "background": resolved_background,
            "pseudocount": float(pseudocount),
            "strand": str(strand),
            "ambiguous_target_policy": "skip-window",
            "max_scan_length": max_scan_length,
            "max_pwm_length": max_pwm_length,
            "score_cells": possible_windows * pwm.length,
            "max_score_cells": max_score_cells,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="FIMO-style PWM log-odds scoring boundary",
            reference_version="DNAKit-1.0",
        ),
    )
    return cast(PatternResult[MotifHit], result)


def scan_promoter_motifs(
    value: SequenceInput,
    *,
    motifs: Mapping[str, str] | None = None,
    strand: Strand | str = Strand.BOTH,
    max_matches: int = 100_000,
    max_scan_length: int = 1_000_000,
    max_scan_cells: int = 50_000_000,
    max_motifs: int = 10_000,
) -> PatternResult[MotifHit]:
    """Scan predefined/user consensus motifs without predicting promoter activity."""

    sequence, sequence_id = resolve_sequence(value)
    if motifs is not None and not isinstance(motifs, Mapping):
        raise ConfigurationError("motifs must be a mapping from names to IUPAC definitions.")
    definitions = dict(PROMOTER_MOTIFS if motifs is None else motifs)
    if not definitions:
        raise ConfigurationError("At least one promoter motif definition is required.")
    if any(not isinstance(motif_name, str) or not motif_name.strip() for motif_name in definitions):
        raise ConfigurationError("Promoter motif names must be non-empty strings.")
    validate_positive_int(max_matches, "max_matches")
    validate_positive_int(max_scan_length, "max_scan_length")
    validate_positive_int(max_scan_cells, "max_scan_cells")
    validate_positive_int(max_motifs, "max_motifs")
    if len(definitions) > max_motifs:
        raise ConfigurationError("Promoter motif panel exceeds max_motifs.")
    definitions = {
        motif_name: validate_iupac_text(motif, f"motifs[{motif_name!r}]")
        for motif_name, motif in definitions.items()
    }
    selected_strands = coerce_strands(strand)
    estimated_total_cells = (
        sequence.symbol_length
        * len(selected_strands)
        * sum(len(motif) for motif in definitions.values())
    )
    if estimated_total_cells > max_scan_cells:
        raise ConfigurationError(
            "Promoter motif panel exceeds max_scan_cells.",
            context={
                "estimated_total_cells": estimated_total_cells,
                "max_scan_cells": max_scan_cells,
            },
        )
    hits: list[MotifHit] = []
    truncated = False
    for motif_name in sorted(definitions):
        remaining = max_matches - len(hits)
        if remaining <= 0:
            truncated = True
            break
        partial = scan_motif(
            sequence,
            definitions[motif_name],
            mode="iupac",
            name=motif_name,
            strand=strand,
            max_matches=remaining,
            max_scan_length=max_scan_length,
            max_scan_cells=max_scan_cells,
        )
        hits.extend(partial.hits)
        truncated = truncated or partial.truncated
    result = build_result(
        sequence,
        sequence_id,
        name="promoter_motif_scan",
        method="iupac_consensus_candidates",
        algorithm_version="1.0",
        parameters={
            "motifs": definitions,
            "catalog": "DNAKit built-in boundary catalog" if motifs is None else "user-provided",
            "catalog_version": PROMOTER_CATALOG_VERSION if motifs is None else None,
            "activity_prediction": False,
            "strand": str(strand),
            "max_scan_length": max_scan_length,
            "estimated_total_scan_cells": estimated_total_cells,
            "max_scan_cells": max_scan_cells,
            "max_motifs": max_motifs,
        },
        hits=hits,
        max_matches=max_matches,
        truncated=truncated,
        provenance=pattern_provenance(
            reimplementation=True,
            reference_name="consensus promoter motif definitions",
            reference_version=PROMOTER_CATALOG_VERSION if motifs is None else "user-provided",
        ),
    )
    return cast(PatternResult[MotifHit], result)


def scan_tf_pwm(
    value: SequenceInput,
    tf_name: str,
    pwm: PWM,
    *,
    threshold: float,
    background: Mapping[str, float] | None = None,
    pseudocount: float = 0.0,
    strand: Strand | str = Strand.BOTH,
    max_matches: int = 100_000,
    max_scan_length: int = 1_000_000,
    max_pwm_length: int = 10_000,
    max_score_cells: int = 50_000_000,
) -> PatternResult[MotifHit]:
    """Scan a user-supplied TF PWM; no binding-strength or database claim is made."""

    if not isinstance(tf_name, str) or not tf_name.strip():
        raise ConfigurationError("tf_name must be non-empty.")
    if not isinstance(pwm, PWM):
        raise ConfigurationError("pwm must be a PWM object.")
    renamed = PWM(tf_name, {base: pwm.row(base) for base in "ACGT"})
    result = scan_pwm(
        value,
        renamed,
        threshold=threshold,
        background=background,
        pseudocount=pseudocount,
        strand=strand,
        max_matches=max_matches,
        max_scan_length=max_scan_length,
        max_pwm_length=max_pwm_length,
        max_score_cells=max_score_cells,
    )
    parameters = dict(result.parameters)
    parameters.update({"tf_name": tf_name, "binding_strength_prediction": False})
    return PatternResult(
        name="tf_motif_scan",
        method=result.method,
        algorithm_version=result.algorithm_version,
        sequence_id=result.sequence_id,
        parameters=type(result.parameters)(parameters),
        hits=result.hits,
        inspected_symbol_count=result.inspected_symbol_count,
        gap_count=result.gap_count,
        unknown_gap_count=result.unknown_gap_count,
        max_matches=result.max_matches,
        truncated=result.truncated,
        coordinate_system=result.coordinate_system,
        gap_policy=result.gap_policy,
        topology=result.topology,
        provenance=result.provenance,
        issues=result.issues,
    )


__all__ = [
    "MAX_PWM_LENGTH",
    "PROMOTER_CATALOG_VERSION",
    "PROMOTER_MOTIFS",
    "PWM",
    "scan_motif",
    "scan_promoter_motifs",
    "scan_pwm",
    "scan_tf_pwm",
]
