"""Small dependency-free sequence summaries for reproducible local analyses."""

from __future__ import annotations

from collections.abc import Iterable
from statistics import mean, median


def sequence_length_gc_summary(sequences: Iterable[str]) -> dict[str, object]:
    """Summarize lengths and GC fraction for an iterable of DNA strings.

    Whitespace is ignored and bases are case-insensitive. ``N`` and other
    IUPAC ambiguity symbols contribute to length but not to GC or AT counts.
    The function consumes the iterable once and returns JSON-serializable data.
    """

    lengths: list[int] = []
    gc_fractions: list[float] = []
    for raw in sequences:
        if not isinstance(raw, str):
            raise TypeError("sequences must contain strings")
        sequence = "".join(raw.split()).upper().replace("U", "T")
        invalid = sorted(set(sequence) - set("ACGTRYSWKMBDHVN"))
        if invalid:
            raise ValueError(f"invalid DNA symbols: {', '.join(invalid)}")
        lengths.append(len(sequence))
        gc_fractions.append((sequence.count("G") + sequence.count("C")) / len(sequence) if sequence else 0.0)

    if not lengths:
        return {
            "count": 0,
            "lengths": [],
            "min_length": 0,
            "max_length": 0,
            "mean_length": 0.0,
            "median_length": 0.0,
            "mean_gc_fraction": 0.0,
        }
    return {
        "count": len(lengths),
        "lengths": lengths,
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": mean(lengths),
        "median_length": median(lengths),
        "mean_gc_fraction": mean(gc_fractions),
    }
