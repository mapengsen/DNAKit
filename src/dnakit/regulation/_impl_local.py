"""Native helpers behind the regulation candidate functions (PWM scoring, enrichment tests, in-silico mutagenesis)."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

BASES = "ACGT"
MIN_PROBABILITY = 1e-6  # zero probabilities are floored here so log-odds stay finite (about -17.9 bits for a 0.25 background)
MAX_SEQUENCE = 5_000_000
_COMP = str.maketrans("ACGT", "TGCA")


# ----------------------------------------------------------------------------- option helpers
def opt_str(options: Mapping[str, object], key: str, default: str, allowed: Sequence[str] | None = None) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"Option {key!r} must be a non-empty string.", code="REGULATION_OPTION_INVALID")
    value = value.strip().lower()
    if allowed is not None and value not in allowed:
        raise ConfigurationError(
            f"Option {key!r}={value!r} is not supported.", code="REGULATION_OPTION_INVALID", hint=f"Use one of: {', '.join(allowed)}."
        )
    return value


def opt_int(options: Mapping[str, object], key: str, default: int, minimum: int, maximum: int) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"Option {key!r} must be an integer in [{minimum}, {maximum}].", code="REGULATION_OPTION_INVALID"
        )
    return value


def opt_float(options: Mapping[str, object], key: str, default: float, minimum: float, maximum: float) -> float:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"Option {key!r} must be a number in [{minimum}, {maximum}].", code="REGULATION_OPTION_INVALID"
        )
    return float(value)


def entity(options: Mapping[str, object], default: str, allowed: Sequence[str]) -> str:
    return opt_str(options, "entity", default, allowed)


# ----------------------------------------------------------------------------- sequences
def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def clean_dna(text: object, name: str = "sequence", allow_n: bool = True) -> str:
    if not isinstance(text, str):
        raise InputFormatError(f"{name} must be a DNA string.", code="REGULATION_SEQUENCE_INVALID")
    seq = re.sub(r"\s+", "", text).upper().replace("U", "T")
    allowed = set(BASES + ("N" if allow_n else ""))
    if not seq or set(seq) - allowed:
        raise InputFormatError(f"{name} is empty or has characters outside ACGT{'N' if allow_n else ''}.",
                               code="REGULATION_SEQUENCE_INVALID")
    if len(seq) > MAX_SEQUENCE:
        raise ConfigurationError(f"{name} is longer than {MAX_SEQUENCE} bases.", code="REGULATION_SEQUENCE_TOO_LONG")
    return seq


def _is_file(text: str) -> bool:
    if "\n" in text or "\t" in text or len(text) >= 4096:
        return False
    try:
        return Path(text).is_file()
    except OSError:
        return False


def read_text(value: object) -> str:
    if isinstance(value, Path) or (isinstance(value, str) and _is_file(value)):
        return Path(value).read_text(encoding="utf-8", errors="replace")
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value
    raise InputFormatError("Expected text, bytes or a path.", code="REGULATION_INPUT_INVALID")


def parse_sequences(value: object, name: str = "sequences", max_items: int = 100_000) -> list[tuple[str, str]]:
    """A DNA string, FASTA text/path, a list of strings or records (``sequence``/``id``), or ``{id: sequence}``."""

    if isinstance(value, Mapping):
        if "sequence" in value and isinstance(value["sequence"], str):
            items = [(str(value.get("id", "seq1")), value["sequence"])]
        else:
            items = [(str(k), str(v)) for k, v in value.items()]
    elif isinstance(value, (str, bytes, Path)):
        text = read_text(value)
        if text.lstrip().startswith(">"):
            items, cur, chunks = [], None, []
            for line in text.splitlines():
                if line.startswith(">"):
                    if cur is not None:
                        items.append((cur, "".join(chunks)))
                    cur, chunks = (line[1:].split()[0] if line[1:].strip() else f"seq{len(items) + 1}"), []
                else:
                    chunks.append(line.strip())
            if cur is not None:
                items.append((cur, "".join(chunks)))
        else:
            lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
            items = [(f"seq{i + 1}", ln) for i, ln in enumerate(lines)]
    elif isinstance(value, Sequence):
        items = []
        for i, v in enumerate(value):
            if isinstance(v, Mapping):
                items.append((str(v.get("id", f"seq{i + 1}")), str(v.get("sequence", ""))))
            else:
                items.append((f"seq{i + 1}", str(v)))
    else:
        raise InputFormatError(f"Unsupported {name} input.", code="REGULATION_INPUT_INVALID")
    if not items or len(items) > max_items:
        raise InputFormatError(f"Provide between 1 and {max_items} {name}.", code="REGULATION_INPUT_INVALID")
    return [(n, clean_dna(s, f"{name} {n!r}")) for n, s in items]


# ----------------------------------------------------------------------------- PWM
@dataclass
class Pwm:
    """Position probability matrix (rows = positions, columns A, C, G, T) and its log2-odds matrix."""

    name: str
    probs: list[list[float]]
    background: list[float]
    log_odds: list[list[float]]

    @property
    def length(self) -> int:
        return len(self.probs)

    @property
    def max_score(self) -> float:
        return sum(max(row) for row in self.log_odds)

    @property
    def min_score(self) -> float:
        return sum(min(row) for row in self.log_odds)


def _background(options: Mapping[str, object]) -> list[float]:
    raw = options.get("background")
    if raw is None:
        return [0.25] * 4
    if not isinstance(raw, Mapping) or set(raw) != set(BASES):
        raise ConfigurationError("background must map A, C, G and T to probabilities.", code="REGULATION_OPTION_INVALID")
    values = [float(raw[b]) for b in BASES]
    if any(not math.isfinite(v) or v <= 0 for v in values):
        raise ConfigurationError("background probabilities must be positive.", code="REGULATION_OPTION_INVALID")
    total = sum(values)
    return [v / total for v in values]


def _matrix_rows(value: object, label: str) -> list[list[float]]:
    """Rows (positions x 4) from {A,C,G,T: [..]}, a list of 4-number rows, or a JASPAR-style text block."""

    if isinstance(value, Mapping):
        if all(b in value for b in BASES):
            cols = [list(value[b]) for b in BASES]  # type: ignore[call-overload]
            if len({len(c) for c in cols}) != 1:
                raise InputFormatError(f"{label}: A/C/G/T rows have different lengths.", code="REGULATION_PWM_INVALID")
            return [[float(cols[b][i]) for b in range(4)] for i in range(len(cols[0]))]
        for key in ("matrix", "pwm", "ppm", "pfm", "counts", "probabilities"):
            if key in value:
                return _matrix_rows(value[key], label)
        raise InputFormatError(f"{label}: expected A/C/G/T keys or a matrix.", code="REGULATION_PWM_INVALID")
    if isinstance(value, (str, bytes, Path)):
        text = read_text(value)
        letters: dict[str, list[float]] = {}
        for line in text.splitlines():
            m = re.match(r"^\s*([ACGT])\s*\[?([^\]]*)\]?\s*$", line)
            if m:
                letters[m.group(1)] = [float(x) for x in m.group(2).split()]
        if len(letters) != 4:
            raise InputFormatError(f"{label}: no JASPAR-style A/C/G/T rows found.", code="REGULATION_PWM_INVALID")
        return _matrix_rows(letters, label)
    if isinstance(value, Sequence) and value and all(isinstance(r, Sequence) and not isinstance(r, str) and len(r) == 4 for r in value):
        return [[float(x) for x in r] for r in value]  # type: ignore[union-attr]
    raise InputFormatError(f"{label}: unsupported matrix; use A/C/G/T lists, a Lx4 list or JASPAR text.",
                           code="REGULATION_PWM_INVALID")


def build_pwm(value: object, name: str, options: Mapping[str, object]) -> Pwm:
    """Turn counts or probabilities into a Pwm.

    ``matrix_type``: ``"auto"`` (default: columns that all sum to 1 +- 1e-3 are probabilities, otherwise counts),
    ``"counts"`` or ``"probabilities"``. Counts become ``(count + pseudocount) / (column total + 4 pseudocount)``
    (``pseudocount`` default 0). Probabilities are floored at 1e-6 and renormalised per column.
    """

    rows = _matrix_rows(value, name)
    if not rows or len(rows) > 1_000:
        raise InputFormatError(f"{name}: a PWM needs 1-1000 positions.", code="REGULATION_PWM_INVALID")
    if any(not math.isfinite(x) or x < 0 for r in rows for x in r):
        raise InputFormatError(f"{name}: matrix entries must be finite and non-negative.", code="REGULATION_PWM_INVALID")
    if any(sum(r) <= 0 for r in rows):
        raise InputFormatError(f"{name}: a column sums to zero.", code="REGULATION_PWM_INVALID")
    mtype = opt_str(options, "matrix_type", "auto", ("auto", "counts", "probabilities"))
    if mtype == "auto":
        mtype = "probabilities" if all(abs(sum(r) - 1.0) <= 1e-3 for r in rows) else "counts"
    pseudo = opt_float(options, "pseudocount", 0.0, 0.0, 1e6)
    probs = []
    for r in rows:
        if mtype == "counts":
            total = sum(r) + 4 * pseudo
            p = [(x + pseudo) / total for x in r]
        else:
            p = [x / sum(r) for x in r]
        p = [max(x, MIN_PROBABILITY) for x in p]
        total = sum(p)
        probs.append([x / total for x in p])
    bg = _background(options)
    return Pwm(name, probs, bg, [[math.log2(p[b] / bg[b]) for b in range(4)] for p in probs])


def parse_motifs(options: Mapping[str, object], key: str = "motifs") -> list[Pwm]:
    """``motifs`` (or ``motif``): one matrix, a list of matrices, or ``{name: matrix}``."""

    raw = options.get(key, options.get("motif"))
    if raw is None:
        raise ConfigurationError(f"This operation needs options['{key}'] (a PWM).", code="REGULATION_OPTION_INVALID",
                                 hint="Pass counts/probabilities as {'A': [...], 'C': [...], 'G': [...], 'T': [...]}.")
    named: list[tuple[str, object]]
    single = isinstance(raw, Mapping) and (all(b in raw for b in BASES) or any(k in raw for k in ("matrix", "pwm", "ppm", "pfm", "counts", "probabilities")))
    if single or isinstance(raw, (str, bytes, Path)):
        named = [(str(raw.get("name", "motif1")) if isinstance(raw, Mapping) else "motif1", raw)]
    elif isinstance(raw, Mapping):
        named = [(str(k), v) for k, v in raw.items()]
    elif isinstance(raw, Sequence) and raw and all(isinstance(r, Sequence) and not isinstance(r, str) and len(r) == 4
                                                   and all(isinstance(x, (int, float)) for x in r) for r in raw):
        named = [("motif1", raw)]
    elif isinstance(raw, Sequence):
        named = [(str(v.get("name", f"motif{i + 1}")) if isinstance(v, Mapping) else f"motif{i + 1}", v) for i, v in enumerate(raw)]
    else:
        raise InputFormatError("Unsupported motif input.", code="REGULATION_PWM_INVALID")
    if not named or len(named) > 500:
        raise InputFormatError("Provide between 1 and 500 motifs.", code="REGULATION_PWM_INVALID")
    return [build_pwm(v, n, options) for n, v in named]


def information_content(pwm: Pwm) -> dict[str, object]:
    """Per-column information content sum_b p_b log2(p_b / bg_b) in bits (2 - entropy for a uniform background)."""

    columns = [sum(p[b] * math.log2(p[b] / pwm.background[b]) for b in range(4) if p[b] > 0) for p in pwm.probs]
    return {"per_position": columns, "total_bits": sum(columns)}


def window_scores(seq: str, pwm: Pwm) -> list[tuple[float, float]]:
    """(forward, reverse) log2-odds for every window start; windows with N score -inf; reverse = revcomp window."""

    w = pwm.length
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    out = []
    for s in range(len(seq) - w + 1):
        window = seq[s : s + w]
        if "N" in window:
            out.append((-math.inf, -math.inf))
            continue
        fwd = sum(pwm.log_odds[i][idx[c]] for i, c in enumerate(window))
        rev_window = revcomp(window)
        rev = sum(pwm.log_odds[i][idx[c]] for i, c in enumerate(rev_window))
        out.append((fwd, rev))
    return out


def relative(score: float, pwm: Pwm) -> float:
    span = pwm.max_score - pwm.min_score
    return (score - pwm.min_score) / span if span > 0 else 1.0


def scan(seq: str, pwm: Pwm, threshold: float | None, min_relative: float | None, strand: str = "both") -> list[dict[str, object]]:
    """All windows with score >= threshold (absolute log2-odds) or relative score >= min_relative."""

    hits = []
    for s, (fwd, rev) in enumerate(window_scores(seq, pwm)):
        for st, score in (("+", fwd), ("-", rev)):
            if (strand == "forward" and st == "-") or (strand == "reverse" and st == "+"):
                continue
            if not math.isfinite(score):
                continue
            ok = score >= threshold if threshold is not None else relative(score, pwm) >= (min_relative or 0.0)
            if ok:
                hits.append({"start": s, "end": s + pwm.length, "strand": st, "score": score,
                             "relative_score": relative(score, pwm),
                             "matched": seq[s : s + pwm.length] if st == "+" else revcomp(seq[s : s + pwm.length])})
    hits.sort(key=lambda h: (h["start"], h["strand"]))  # type: ignore[arg-type,return-value]
    return hits


def best_hit(seq: str, pwm: Pwm, strand: str = "both") -> dict[str, object] | None:
    best = None
    for s, (fwd, rev) in enumerate(window_scores(seq, pwm)):
        for st, score in (("+", fwd), ("-", rev)):
            if (strand == "forward" and st == "-") or (strand == "reverse" and st == "+") or not math.isfinite(score):
                continue
            if best is None or score > best["score"]:  # type: ignore[operator]
                best = {"start": s, "end": s + pwm.length, "strand": st, "score": score}
    if best is not None:
        best["relative_score"] = relative(float(best["score"]), pwm)  # type: ignore[arg-type]
    return best


def resolve_scan_thresholds(options: Mapping[str, object]) -> tuple[float | None, float | None]:
    threshold = options.get("threshold")
    if threshold is not None:
        return opt_float(options, "threshold", 0.0, -1e6, 1e6), None
    return None, opt_float(options, "min_relative_score", 0.8, 0.0, 1.0)


def strand_option(options: Mapping[str, object]) -> str:
    return opt_str(options, "strand", "both", ("both", "forward", "reverse"))


# ----------------------------------------------------------------------------- PWM comparison
def compare_pwms(a: Pwm, b: Pwm, min_overlap: int) -> dict[str, object]:
    """Best ungapped alignment of two PWMs by mean per-column Pearson correlation of the A/C/G/T probabilities.

    Columns with zero variance correlate 0. Both orientations of ``b`` are tried; ``offset`` is the shift of b
    relative to a (positive: b starts after a starts). Alignments need at least ``min_overlap`` columns.
    """

    def corr(p: Sequence[float], q: Sequence[float]) -> float:
        mp, mq = sum(p) / 4, sum(q) / 4
        sp = sum((x - mp) ** 2 for x in p)
        sq = sum((x - mq) ** 2 for x in q)
        return sum((x - mp) * (y - mq) for x, y in zip(p, q)) / math.sqrt(sp * sq) if sp > 0 and sq > 0 else 0.0

    best: dict[str, object] | None = None
    rc_probs = [list(reversed(p)) for p in reversed(b.probs)]  # reverse complement: reverse rows and A<->T, C<->G
    for orientation, probs in (("forward", b.probs), ("reverse_complement", rc_probs)):
        for offset in range(-(len(probs) - min_overlap), a.length - min_overlap + 1):
            pairs = [(a.probs[i], probs[i - offset]) for i in range(max(0, offset), min(a.length, offset + len(probs)))]
            if len(pairs) < min_overlap:
                continue
            score = sum(corr(p, q) for p, q in pairs) / len(pairs)
            if best is None or score > best["similarity"] + 1e-12:  # type: ignore[operator]
                best = {"similarity": score, "offset": offset, "orientation": orientation, "overlap": len(pairs)}
    if best is None:
        raise ConfigurationError("The motifs are shorter than min_overlap.", code="REGULATION_OPTION_INVALID")
    return best


# ----------------------------------------------------------------------------- enrichment statistics
def _log_choose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, population: int, successes: int, draws: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(population N, successes K, draws n)."""

    low, high = max(0, draws - (population - successes)), min(draws, successes)
    if k <= low:
        return 1.0
    if k > high:
        return 0.0
    if high - k > 5_000_000:
        raise ConfigurationError("Counts are too large for the exact test.", code="REGULATION_COUNTS_TOO_LARGE")
    denom = _log_choose(population, draws)
    terms = [_log_choose(successes, x) + _log_choose(population - successes, draws - x) - denom for x in range(k, high + 1)]
    top = max(terms)
    return min(1.0, math.exp(top) * sum(math.exp(t - top) for t in terms))


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""

    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    if n - k > 5_000_000:
        raise ConfigurationError("Counts are too large for the exact test.", code="REGULATION_COUNTS_TOO_LARGE")
    logp, logq = math.log(p), math.log1p(-p)
    terms = [_log_choose(n, x) + x * logp + (n - x) * logq for x in range(k, n + 1)]
    top = max(terms)
    return min(1.0, math.exp(top) * sum(math.exp(t - top) for t in terms))


def benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, pvalues[index] * n / rank)
        adjusted[index] = running
    return adjusted


# ----------------------------------------------------------------------------- variant effect / ISM
def best_score(seq: str, pwm: Pwm) -> float | None:
    hit = best_hit(seq, pwm)
    return float(hit["score"]) if hit else None  # type: ignore[arg-type]


def parse_variants(value: object) -> list[dict[str, object]]:
    """Variants as ``{"id", "ref_sequence", "alt_sequence"}`` or ``{"id", "sequence", "position", "ref", "alt"}``
    (0-based position of ``ref`` inside ``sequence``; ref must match the sequence)."""

    if isinstance(value, Mapping) and ("ref_sequence" in value or "sequence" in value):
        value = [value]
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value or len(value) > 10_000:
        raise InputFormatError("Provide 1-10000 variant records.", code="REGULATION_INPUT_INVALID")
    out = []
    for i, v in enumerate(value):
        if not isinstance(v, Mapping):
            raise InputFormatError("Each variant must be a mapping.", code="REGULATION_INPUT_INVALID")
        vid = str(v.get("id", f"variant{i + 1}"))
        if "ref_sequence" in v and "alt_sequence" in v:
            ref, alt = clean_dna(v["ref_sequence"], f"variant {vid!r} ref_sequence"), clean_dna(v["alt_sequence"], f"variant {vid!r} alt_sequence")
        elif {"sequence", "position", "ref", "alt"} <= set(v):
            seq = clean_dna(v["sequence"], f"variant {vid!r} sequence")
            pos, r, a = v["position"], clean_dna(v["ref"], "ref"), clean_dna(v["alt"], "alt")
            if isinstance(pos, bool) or not isinstance(pos, int) or pos < 0 or seq[pos : pos + len(r)] != r:
                raise InputFormatError(f"Variant {vid!r}: the reference allele does not match the sequence at position {pos}.",
                                       code="REGULATION_VARIANT_INVALID")
            ref, alt = seq, seq[:pos] + a + seq[pos + len(r) :]
        else:
            raise InputFormatError(f"Variant {vid!r} needs ref_sequence+alt_sequence or sequence+position+ref+alt.",
                                   code="REGULATION_INPUT_INVALID")
        out.append({"id": vid, "ref_sequence": ref, "alt_sequence": alt})
    return out


def variant_effects(variants: Sequence[Mapping[str, object]], pwms: Sequence[Pwm]) -> list[dict[str, object]]:
    rows = []
    for v in variants:
        per_motif = []
        for pwm in pwms:
            ref_hit, alt_hit = best_hit(str(v["ref_sequence"]), pwm), best_hit(str(v["alt_sequence"]), pwm)
            if ref_hit is None and alt_hit is None:
                continue
            ref_s = float(ref_hit["score"]) if ref_hit else None  # type: ignore[arg-type]
            alt_s = float(alt_hit["score"]) if alt_hit else None  # type: ignore[arg-type]
            per_motif.append({"motif": pwm.name, "ref_score": ref_s, "alt_score": alt_s,
                              "delta": alt_s - ref_s if ref_s is not None and alt_s is not None else None})
        deltas = [abs(m["delta"]) for m in per_motif if m["delta"] is not None]  # type: ignore[arg-type]
        rows.append({"id": v["id"], "motifs": per_motif, "max_abs_delta": max(deltas) if deltas else None})
    return rows


def ism_matrix(seq: str, pwm: Pwm) -> list[dict[str, object]]:
    """Single-base substitution effects on the best PWM score (max over windows and both strands).

    effect(i, b) = best_score(seq with base b at i) - best_score(seq). Unscorable sequences (shorter than the PWM,
    or only N windows) are rejected earlier.
    """

    w, n = pwm.length, len(seq)
    base = [max(f, r) for f, r in window_scores(seq, pwm)]
    reference = max(base)
    out = []
    for i, ref_base in enumerate(seq):
        lo, hi = max(0, i - w + 1), min(i, n - w)  # window starts covering position i
        outside = max([base[s] for s in range(0, lo)] + [base[s] for s in range(hi + 1, len(base))] + [-math.inf])
        effects: dict[str, float | None] = {}
        for b in BASES:
            if b == ref_base:
                continue
            mutated = seq[:i] + b + seq[i + 1 :]
            inner = max(max(window_scores(mutated[s : s + w], pwm)[0]) for s in range(lo, hi + 1))
            value = max(outside, inner)
            effects[b] = value - reference if math.isfinite(value) else None
        out.append({"position": i, "ref": ref_base, "effects": effects})
    return out
