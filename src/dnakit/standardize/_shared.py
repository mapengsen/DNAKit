"""Private character-analysis kernels shared by normalization and validation."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping

from dnakit.core.gap import Gap

from .results import (
    AmbiguityReport,
    ProbabilityResolution,
    SymbolCount,
    SymbolOccurrence,
)

CANONICAL_BASES = frozenset("ACGT")
IUPAC_AMBIGUITY = frozenset("RYSWKMBDHVN")
IUPAC_BASES = CANONICAL_BASES | IUPAC_AMBIGUITY
IUPAC_EXPANSION: Mapping[str, frozenset[str]] = {
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}


def ambiguity_report(
    parts: Iterable[str | Gap],
    *,
    denominator_includes_gap: bool = False,
    base_priors: Mapping[str, float] | None = None,
) -> AmbiguityReport:
    """Describe every ambiguity symbol using normalized symbol coordinates."""
    positions: dict[str, list[int]] = defaultdict(list)
    occurrences: list[SymbolOccurrence] = []
    resolutions: list[ProbabilityResolution] = []
    symbol_offset = 0
    gap_total = 0
    denominator_known = True

    for part_index, part in enumerate(parts):
        if isinstance(part, Gap):
            if denominator_includes_gap:
                if part.length is None:
                    denominator_known = False
                else:
                    gap_total += part.length
            continue
        for part_offset, symbol in enumerate(part):
            if symbol in IUPAC_AMBIGUITY:
                positions[symbol].append(symbol_offset)
                occurrences.append(
                    SymbolOccurrence(
                        symbol=symbol,
                        symbol_offset=symbol_offset,
                        part_index=part_index,
                        part_offset=part_offset,
                    )
                )
                if base_priors is not None:
                    resolutions.append(
                        ProbabilityResolution(
                            symbol=symbol,
                            symbol_offset=symbol_offset,
                            probabilities=_probabilities(symbol, base_priors),
                        )
                    )
            symbol_offset += 1

    total = sum(len(item) for item in positions.values())
    denominator = symbol_offset + gap_total if denominator_known else None
    fraction = None if not denominator else total / denominator
    return AmbiguityReport(
        total_count=total,
        denominator=denominator,
        fraction=fraction,
        by_symbol=tuple(
            SymbolCount(symbol=symbol, count=len(items), positions=tuple(items))
            for symbol, items in sorted(positions.items())
        ),
        occurrences=tuple(occurrences),
        probability_resolutions=tuple(resolutions),
    )


def _probabilities(symbol: str, priors: Mapping[str, float]) -> tuple[tuple[str, float], ...]:
    allowed = IUPAC_EXPANSION[symbol]
    total = sum(priors[base] for base in allowed)
    if total == 0.0:
        uniform = 1.0 / len(allowed)
        return tuple((base, uniform) for base in sorted(allowed))
    return tuple((base, priors[base] / total) for base in sorted(allowed))
