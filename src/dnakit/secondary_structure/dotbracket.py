"""Backend-neutral dot-bracket parsing and probability-derived metrics."""

from __future__ import annotations

import math
from collections.abc import Iterable
from itertools import islice

from dnakit.core import DNASequence
from dnakit.exceptions import ConfigurationError
from dnakit.thermodynamics._shared import canonical_linear_symbols

from .results import (
    AccessibilityWindow,
    PairProbabilityResult,
    SecondaryBasePair,
    SecondaryStem,
    SecondaryStructureSummary,
)

_OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
_CLOSE_TO_OPEN = {value: key for key, value in _OPEN_TO_CLOSE.items()}
_MAX_STRANDS = 100
_MAX_TOTAL_NT = 10_000


def _strands(values: Iterable[DNASequence]) -> tuple[str, ...]:
    if not isinstance(values, Iterable):
        raise ConfigurationError(
            "strands must be an iterable of DNASequence objects.",
            code="INVALID_SECONDARY_STRUCTURE_STRANDS",
        )
    items = tuple(islice(iter(values), _MAX_STRANDS + 1))
    if not 1 <= len(items) <= _MAX_STRANDS:
        raise ConfigurationError(
            f"strands must contain 1-{_MAX_STRANDS} entries.",
            code="INVALID_SECONDARY_STRUCTURE_STRANDS",
        )
    if any(not isinstance(item, DNASequence) for item in items):
        raise ConfigurationError(
            "Every strand must be a DNASequence.",
            code="INVALID_SECONDARY_STRUCTURE_STRANDS",
        )
    symbols = tuple(
        canonical_linear_symbols(
            item,
            operation="secondary structure",
            min_length=1,
            max_length=_MAX_TOTAL_NT,
        )
        for item in items
    )
    if sum(map(len, symbols)) > _MAX_TOTAL_NT:
        raise ConfigurationError(
            f"Secondary-structure input exceeds {_MAX_TOTAL_NT} total nucleotides.",
            code="SECONDARY_STRUCTURE_SIZE_LIMIT",
        )
    return symbols


def _locations(strands: tuple[str, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(
        (strand_index, local_index)
        for strand_index, sequence in enumerate(strands)
        for local_index in range(len(sequence))
    )


def analyze_dot_bracket(
    strands: Iterable[DNASequence],
    dot_bracket: str,
    *,
    three_prime_window: int = 5,
) -> SecondaryStructureSummary:
    """Parse dot-parens-plus or extended bracket notation into structure metrics."""

    symbols = _strands(strands)
    if not isinstance(dot_bracket, str) or not dot_bracket:
        raise ConfigurationError("dot_bracket must be non-empty text.", code="INVALID_DOT_BRACKET")
    if (
        isinstance(three_prime_window, bool)
        or not isinstance(three_prime_window, int)
        or not 1 <= three_prime_window <= _MAX_TOTAL_NT
    ):
        raise ConfigurationError(
            "three_prime_window must be a positive bounded integer.",
            code="INVALID_DOT_BRACKET_WINDOW",
        )
    parts = dot_bracket.split("+")
    if len(parts) != len(symbols) or any(
        len(part) != len(sequence) for part, sequence in zip(parts, symbols, strict=True)
    ):
        raise ConfigurationError(
            "Dot-bracket strand separators and lengths must match the input strands.",
            code="DOT_BRACKET_LENGTH_MISMATCH",
            context={
                "expected_lengths": tuple(map(len, symbols)),
                "observed_lengths": tuple(map(len, parts)),
            },
        )
    structure = "".join(parts)
    invalid = sorted(set(structure) - {".", *_OPEN_TO_CLOSE, *_CLOSE_TO_OPEN})
    if invalid:
        raise ConfigurationError(
            "Dot-bracket text contains unsupported symbols.",
            code="INVALID_DOT_BRACKET_SYMBOL",
            context={"symbols": tuple(invalid)},
        )
    locations = _locations(symbols)
    stacks: dict[str, list[int]] = {character: [] for character in _OPEN_TO_CLOSE}
    pairs: list[SecondaryBasePair] = []
    for global_index, character in enumerate(structure):
        if character in _OPEN_TO_CLOSE:
            stacks[character].append(global_index)
        elif character in _CLOSE_TO_OPEN:
            opener = _CLOSE_TO_OPEN[character]
            if not stacks[opener]:
                raise ConfigurationError(
                    "Dot-bracket text closes a pair before it opens.",
                    code="UNBALANCED_DOT_BRACKET",
                    context={"global_index": global_index, "symbol": character},
                )
            first = stacks[opener].pop()
            first_strand, first_local = locations[first]
            second_strand, second_local = locations[global_index]
            pairs.append(
                SecondaryBasePair(
                    first_global_index=first,
                    second_global_index=global_index,
                    first_strand_index=first_strand,
                    first_local_index=first_local,
                    second_strand_index=second_strand,
                    second_local_index=second_local,
                    bracket_type=opener + character,
                    inter_strand=first_strand != second_strand,
                )
            )
    unclosed = tuple((symbol, tuple(indices)) for symbol, indices in stacks.items() if indices)
    if unclosed:
        raise ConfigurationError(
            "Dot-bracket text contains unclosed pairs.",
            code="UNBALANCED_DOT_BRACKET",
            context={"unclosed": unclosed},
        )
    ordered_pairs = tuple(sorted(pairs, key=lambda item: item.first_global_index))
    by_indices = {
        (item.first_global_index, item.second_global_index): item for item in ordered_pairs
    }
    remaining = set(by_indices)
    stems: list[SecondaryStem] = []
    while remaining:
        start = min(remaining)
        run: list[SecondaryBasePair] = []
        current = start
        while current in remaining:
            pair = by_indices[current]
            run.append(pair)
            remaining.remove(current)
            current = (current[0] + 1, current[1] - 1)
        stems.append(
            SecondaryStem(
                base_pairs=tuple(run),
                length=len(run),
                inter_strand=all(item.inter_strand for item in run),
            )
        )
    hairpin_loops: list[int] = []
    for stem in stems:
        innermost = stem.base_pairs[-1]
        if innermost.inter_strand:
            continue
        nested_pair_exists = any(
            innermost.first_global_index < other.first_global_index
            and other.second_global_index < innermost.second_global_index
            for other in ordered_pairs
        )
        if not nested_pair_exists:
            hairpin_loops.append(innermost.second_local_index - innermost.first_local_index - 1)
    three_prime_stems = tuple(
        stem
        for stem in stems
        if stem.inter_strand
        and any(
            pair.first_local_index >= len(symbols[pair.first_strand_index]) - three_prime_window
            or pair.second_local_index
            >= len(symbols[pair.second_strand_index]) - three_prime_window
            for pair in stem.base_pairs
        )
    )
    if not ordered_pairs:
        structure_type = "unstructured" if len(symbols) == 1 else "unbound-strands"
    elif len(symbols) == 1:
        structure_type = "hairpin" if hairpin_loops else "intramolecular-fold"
    elif len(symbols) == 2 and symbols[0] == symbols[1]:
        structure_type = "self-dimer"
    elif len(symbols) == 2:
        structure_type = "heterodimer"
    else:
        structure_type = "multi-strand-complex"
    total_nt = sum(map(len, symbols))
    stem_lengths = tuple(stem.length for stem in stems)
    return SecondaryStructureSummary(
        strands_5to3=symbols,
        dot_bracket=dot_bracket,
        base_pairs=ordered_pairs,
        stems=tuple(stems),
        structure_type=structure_type,
        base_pair_count=len(ordered_pairs),
        paired_base_fraction=(2.0 * len(ordered_pairs) / total_nt),
        stem_lengths=stem_lengths,
        hairpin_count=len(hairpin_loops),
        hairpin_loop_lengths=tuple(hairpin_loops),
        max_contiguous_pair_count=max(stem_lengths, default=0),
        three_prime_window=three_prime_window,
        three_prime_dimer=bool(three_prime_stems),
        three_prime_dimer_max_contiguous_pairs=max(
            (stem.length for stem in three_prime_stems), default=0
        ),
        method="deterministic-dot-bracket-parser-v1",
    )


def pair_probability_metrics(
    strands: Iterable[DNASequence],
    probability_matrix: Iterable[Iterable[float]],
    *,
    accessibility_window_size: int = 4,
) -> PairProbabilityResult:
    """Validate a dense NUPACK-style matrix and derive per-base accessibility."""

    symbols = _strands(strands)
    total_nt = sum(map(len, symbols))
    if (
        isinstance(accessibility_window_size, bool)
        or not isinstance(accessibility_window_size, int)
        or accessibility_window_size <= 0
    ):
        raise ConfigurationError(
            "accessibility_window_size must be a positive integer.",
            code="INVALID_ACCESSIBILITY_WINDOW",
        )
    if not isinstance(probability_matrix, Iterable):
        raise ConfigurationError(
            "probability_matrix must be a square numeric matrix.",
            code="INVALID_PAIR_PROBABILITY_MATRIX",
        )
    raw_rows = tuple(islice(iter(probability_matrix), total_nt + 1))
    if len(raw_rows) != total_nt:
        raise ConfigurationError(
            "Pair-probability matrix size does not match total strand length.",
            code="PAIR_PROBABILITY_SIZE_MISMATCH",
        )
    rows: list[tuple[float, ...]] = []
    for row in raw_rows:
        if isinstance(row, (str, bytes)) or not isinstance(row, Iterable):
            raise ConfigurationError(
                "Every pair-probability row must be iterable.",
                code="INVALID_PAIR_PROBABILITY_MATRIX",
            )
        raw_values = tuple(islice(iter(row), total_nt + 1))
        if len(raw_values) != total_nt:
            raise ConfigurationError(
                "Pair-probability matrix must be square.",
                code="INVALID_PAIR_PROBABILITY_MATRIX",
            )
        values: list[float] = []
        for value in raw_values:
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or not 0.0 <= value <= 1.0
            ):
                raise ConfigurationError(
                    "Pair probabilities must be finite values in [0, 1].",
                    code="INVALID_PAIR_PROBABILITY_MATRIX",
                )
            values.append(float(value))
        if not math.isclose(math.fsum(values), 1.0, abs_tol=1e-5):
            raise ConfigurationError(
                "Each dense pair-probability row must sum to one including its diagonal.",
                code="INVALID_PAIR_PROBABILITY_ROW_SUM",
            )
        rows.append(tuple(values))
    for first in range(total_nt):
        for second in range(first + 1, total_nt):
            if not math.isclose(rows[first][second], rows[second][first], abs_tol=1e-5):
                raise ConfigurationError(
                    "Pair-probability matrix must be symmetric.",
                    code="ASYMMETRIC_PAIR_PROBABILITY_MATRIX",
                )
    unpaired = tuple(rows[index][index] for index in range(total_nt))
    pairing = tuple(math.fsum(rows[index]) - rows[index][index] for index in range(total_nt))
    off_diagonal = tuple(
        tuple(0.0 if row_index == column_index else value for column_index, value in enumerate(row))
        for row_index, row in enumerate(rows)
    )
    windows: list[AccessibilityWindow] = []
    offset = 0
    for sequence in symbols:
        if len(sequence) >= accessibility_window_size:
            for local_start in range(len(sequence) - accessibility_window_size + 1):
                start = offset + local_start
                end = start + accessibility_window_size
                windows.append(
                    AccessibilityWindow(
                        start=start,
                        end=end,
                        mean_unpaired_probability=math.fsum(unpaired[start:end])
                        / accessibility_window_size,
                    )
                )
        offset += len(sequence)
    most_accessible = (
        None
        if not windows
        else max(windows, key=lambda item: (item.mean_unpaired_probability, -item.start)).start
    )
    return PairProbabilityResult(
        strands_5to3=symbols,
        pair_probabilities=off_diagonal,
        pairing_probabilities_by_base=pairing,
        unpaired_probabilities_by_base=unpaired,
        accessibility_window_size=accessibility_window_size,
        accessibility_windows=tuple(windows),
        most_accessible_window_start=most_accessible,
        method="nupack-dense-pair-matrix-derived-marginals-v1",
        applicability=(
            "Window accessibility is the arithmetic mean of marginal unpaired probabilities; "
            "it is not the joint probability that every base in the window is unpaired."
        ),
    )


def ensemble_defect_from_probabilities(
    target: SecondaryStructureSummary,
    probabilities: PairProbabilityResult,
) -> float:
    """Return the normalized expected number of incorrectly paired nucleotides."""

    if not isinstance(target, SecondaryStructureSummary) or not isinstance(
        probabilities, PairProbabilityResult
    ):
        raise ConfigurationError(
            "target and probabilities have invalid types.",
            code="INVALID_ENSEMBLE_DEFECT_INPUT",
        )
    total_nt = sum(map(len, target.strands_5to3))
    if probabilities.strands_5to3 != target.strands_5to3:
        raise ConfigurationError(
            "Target and probability results must describe identical strands.",
            code="ENSEMBLE_DEFECT_STRAND_MISMATCH",
        )
    if len(probabilities.unpaired_probabilities_by_base) != total_nt:
        raise ConfigurationError(
            "Target and probability matrix lengths differ.",
            code="ENSEMBLE_DEFECT_LENGTH_MISMATCH",
        )
    target_partner = {index: index for index in range(total_nt)}
    for pair in target.base_pairs:
        target_partner[pair.first_global_index] = pair.second_global_index
        target_partner[pair.second_global_index] = pair.first_global_index
    correct_probabilities = []
    for index, partner in target_partner.items():
        if partner == index:
            correct_probabilities.append(probabilities.unpaired_probabilities_by_base[index])
        else:
            correct_probabilities.append(probabilities.pair_probabilities[index][partner])
    return math.fsum(1.0 - value for value in correct_probabilities) / total_nt


def target_structure_probability(
    target_free_energy_kcal_per_mol: float,
    ensemble_free_energy_kcal_per_mol: float,
    *,
    temperature_celsius: float = 37.0,
) -> float:
    """Calculate P(target) from target and ensemble standard free energies."""

    values = (target_free_energy_kcal_per_mol, ensemble_free_energy_kcal_per_mol)
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
        for value in values
    ):
        raise ConfigurationError(
            "Free energies must be finite numeric data.",
            code="INVALID_TARGET_STRUCTURE_ENERGY",
        )
    if (
        isinstance(temperature_celsius, bool)
        or not isinstance(temperature_celsius, (int, float))
        or not math.isfinite(temperature_celsius)
        or not 0.0 <= temperature_celsius <= 100.0
    ):
        raise ConfigurationError(
            "temperature_celsius must be in [0, 100].",
            code="INVALID_TARGET_STRUCTURE_TEMPERATURE",
        )
    exponent = -(
        float(target_free_energy_kcal_per_mol) - float(ensemble_free_energy_kcal_per_mol)
    ) / (0.00198720425864083 * (float(temperature_celsius) + 273.15))
    return min(1.0, max(0.0, math.exp(min(0.0, exponent))))


__all__ = [
    "analyze_dot_bracket",
    "ensemble_defect_from_probabilities",
    "pair_probability_metrics",
    "target_structure_probability",
]
