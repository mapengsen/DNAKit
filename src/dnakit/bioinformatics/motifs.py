"""Native PWM/PSSM construction, consensus, score distributions and motif scans."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dnakit.bioinformatics._common import input_path, matrix, positive, require, sequence_text
from dnakit.exceptions import ConfigurationError

_IUPAC = {
    frozenset(bases): symbol
    for symbol, bases in {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": "AG",
        "Y": "CT",
        "S": "CG",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT",
    }.items()
}


@dataclass(frozen=True)
class MotifModel:
    name: str
    counts: Any
    probabilities: Any
    log_odds: Any
    background: tuple[float, ...]

    @property
    def consensus(self) -> str:
        return "".join(
            "ACGT"[int(self.probabilities[:, index].argmax())]
            for index in range(self.counts.shape[1])
        )

    def degenerate_consensus(self, *, minimum_probability: float = 0.2) -> str:
        if not 0 < minimum_probability <= 1:
            raise ConfigurationError("minimum_probability must be in (0,1].")
        result = []
        for index in range(self.counts.shape[1]):
            selected = frozenset(
                base
                for row, base in enumerate("ACGT")
                if self.probabilities[row, index] >= minimum_probability
            )
            result.append(_IUPAC[selected] if selected else self.consensus[index])
        return "".join(result)


def build_motif(
    sites: Iterable[Any] | None = None,
    *,
    counts: Any = None,
    name: str = "motif",
    background: Mapping[str, float] | None = None,
    pseudocount: float = 0.5,
) -> MotifModel:
    """Create a 4 x width model in A,C,G,T row order; pseudocount is per cell."""
    np = require("numpy", "bio-analysis")
    if (sites is None) == (counts is None):
        raise ConfigurationError("Provide exactly one of aligned sites or counts.")
    if not name or not math.isfinite(pseudocount) or pseudocount <= 0:
        raise ConfigurationError("Provide a name and positive finite pseudocount.")
    if sites is not None:
        sequences = tuple(sequence_text(site) for site in sites)
        if (
            not sequences
            or len({len(site) for site in sequences}) != 1
            or any(set(site) - set("ACGT") for site in sequences)
        ):
            raise ConfigurationError("Sites must be equal-length canonical DNA sequences.")
        data = np.asarray(
            [
                [
                    sum(site[index] == base for site in sequences)
                    for index in range(len(sequences[0]))
                ]
                for base in "ACGT"
            ],
            dtype=float,
        )
    else:
        data = matrix(counts).copy()
    if data.shape[0] != 4 or (data < 0).any() or (data.sum(axis=0) == 0).any():
        raise ConfigurationError("Counts need four nonnegative rows and positive column totals.")
    bg = dict(background or dict.fromkeys("ACGT", 0.25))
    if (
        set(bg) != set("ACGT")
        or any(not math.isfinite(p) or p <= 0 for p in bg.values())
        or not math.isclose(sum(bg.values()), 1)
    ):
        raise ConfigurationError(
            "Background must be positive A/C/G/T probabilities summing to one."
        )
    probabilities = (data + pseudocount) / (data.sum(axis=0) + 4 * pseudocount)
    background_tuple = tuple(bg[base] for base in "ACGT")
    scores = np.log2(probabilities / np.asarray(background_tuple)[:, None])
    for array in (data, probabilities, scores):
        array.setflags(write=False)
    return MotifModel(name, data, probabilities, scores, background_tuple)


def fpr_threshold(
    model: MotifModel, fpr: float = 1e-4, *, precision: int = 1000, max_states: int = 200_000
) -> float:
    """Conservative background FPR threshold using a quantized dynamic program.

    Each score is rounded down; adding width/precision to the threshold bounds
    rounding error. This is a model/background FPR, not a genome-wide q-value.
    """
    np = require("numpy", "bio-analysis")
    if not 0 < fpr < 1:
        raise ConfigurationError("fpr must lie in (0,1).")
    positive(precision, "precision")
    positive(max_states, "max_states")
    quantized = np.floor(model.log_odds * precision).astype(int)
    distribution: dict[int, float] = {0: 1.0}
    for index in range(quantized.shape[1]):
        updated: defaultdict[int, float] = defaultdict(float)
        for score, probability in distribution.items():
            for row in range(4):
                updated[score + int(quantized[row, index])] += probability * model.background[row]
        if len(updated) > max_states:
            raise ConfigurationError("PWM score distribution exceeds max_states; lower precision.")
        distribution = dict(updated)
    tail = 0.0
    threshold = max(distribution) + 1
    for score in sorted(distribution, reverse=True):
        if tail + distribution[score] > fpr:
            break
        tail += distribution[score]
        threshold = score
    return float((threshold + quantized.shape[1]) / precision)


def scan_motif_model(
    sequence: Any,
    model: MotifModel,
    *,
    threshold: float,
    both_strands: bool = True,
    circular: bool = False,
) -> tuple[dict[str, Any], ...]:
    """Native log-odds scan; reverse hits report forward-reference half-open coordinates.

    Circular hits may have end > sequence_length to express origin crossing.
    Windows containing ambiguity are skipped, not assigned artificial scores.
    """
    text = sequence_text(sequence)
    width = model.counts.shape[1]
    if not math.isfinite(threshold) or width > len(text):
        raise ConfigurationError("Need a finite threshold and motif no longer than the sequence.")
    orientations = [("+", text)]
    if both_strands:
        orientations.append(
            ("-", text.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1])
        )
    result = []
    for strand, oriented in orientations:
        searchable = oriented + oriented[: width - 1] if circular else oriented
        for start in range(len(oriented) if circular else len(oriented) - width + 1):
            word = searchable[start : start + width]
            if set(word) - set("ACGT"):
                continue
            score = sum(
                float(model.log_odds["ACGT".index(base), index]) for index, base in enumerate(word)
            )
            if score >= threshold:
                left = start if strand == "+" else len(text) - start - width
                if circular:
                    left %= len(text)
                result.append(
                    {
                        "start": left,
                        "end": left + width,
                        "strand": strand,
                        "score": score,
                        "motif": model.name,
                    }
                )
    return tuple(sorted(result, key=lambda hit: (hit["start"], hit["strand"])))


def read_motifs(
    path: str | Path, *, format: str = "jaspar", pseudocount: float = 0.5
) -> tuple[MotifModel, ...]:
    """Use Biopython only as a JASPAR/MEME/TRANSFAC file codec."""
    with input_path(path).open(encoding="utf-8") as handle:
        parsed = require("Bio.motifs", "bio-formats").parse(handle, format)
        return tuple(
            build_motif(
                counts=[motif.counts[base] for base in "ACGT"],
                name=motif.name or f"motif_{index}",
                pseudocount=pseudocount,
            )
            for index, motif in enumerate(parsed)
        )


def dotplot(
    left: Any,
    right: Any | None = None,
    *,
    word_length: int = 10,
    max_matches: int = 1_000_000,
    axes: Any = None,
) -> Any:
    """Exact-word forward and reverse-complement dot plot (native word index)."""
    a = sequence_text(left)
    b = a if right is None else sequence_text(right)
    positive(word_length, "word_length")
    positive(max_matches, "max_matches")
    if word_length > min(len(a), len(b)):
        raise ConfigurationError("word_length exceeds a sequence length.")
    plt = require("matplotlib.pyplot", "bio-plots")
    if axes is None:
        _, axes = plt.subplots()
    index: defaultdict[str, list[int]] = defaultdict(list)
    for start in range(len(a) - word_length + 1):
        word = a[start : start + word_length]
        if not set(word) - set("ACGT"):
            index[word].append(start)
    count = 0
    for reverse in (False, True):
        xs: list[int] = []
        ys: list[int] = []
        for start in range(len(b) - word_length + 1):
            word = b[start : start + word_length]
            if reverse:
                word = word.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            matches = index.get(word, ())
            count += len(matches)
            if count > max_matches:
                raise ConfigurationError("Dot plot exceeds max_matches; increase word_length.")
            xs.extend(matches)
            ys.extend([start] * len(matches))
        axes.scatter(xs, ys, s=2, label="reverse complement" if reverse else "forward")
    axes.set(xlabel="Left position (0-based)", ylabel="Right position (0-based)")
    axes.legend()
    return axes
