"""Native implementations behind the assembly candidate functions.

Everything here is pure standard library, deterministic and bounded.  The public
wrappers in ``dnakit.assembly`` call these helpers and wrap the plain results in the
DNAKit candidate envelope.  Definitions used (also repeated in the public docstrings):

* ``Nx`` / ``Lx`` / ``NGx`` follow ``dnakit.assembly.local`` (cumulative length of the
  longest-first sorted contigs reaching x percent of the assembly / expected genome).
* ``auN`` = sum(L_i ** 2) / sum(L_i) (area under the Nx curve).
* GC fraction = (G + C) / (A + C + G + T), N and other IUPAC symbols are excluded.
* A gap is a maximal run of ``N``/``n``.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from dnakit.exceptions import ConfigurationError, InputFormatError

from .local import _nx

MAX_FILE_BYTES = 200_000_000
MAX_TOTAL_BASES = 100_000_000
MAX_KMER_BASES = 5_000_000
MAX_GENOMES = 5_000
MAX_GENES = 500_000
MAX_GRAPH_NODES = 2_000

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

Contig = tuple[str, "str | None", int]  # (name, sequence or None, length)


# --------------------------------------------------------------------------- input parsing


def reverse_complement(sequence: str) -> str:
    """Reverse complement of an ACGT string (other symbols are left unchanged)."""

    return sequence.translate(_COMPLEMENT)[::-1]


def read_text_input(value: object, *, what: str = "input") -> str | None:
    """Return text for a ``Path``/existing-file string input, a plain string, or None for non-text."""

    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    path: Path | None = None
    if isinstance(value, Path):
        path = value
    elif (
        isinstance(value, str)
        and "\n" not in value
        and 0 < len(value) < 4096
        and not value.startswith(">")
    ):
        candidate = Path(value)
        try:
            if candidate.is_file():
                path = candidate
        except OSError:
            path = None
    if path is not None:
        if not path.is_file():
            raise InputFormatError(
                f"{what} file does not exist: {path}",
                code="ASSEMBLY_INPUT_FILE",
                hint="Pass an existing file.",
            )
        if path.stat().st_size > MAX_FILE_BYTES:
            raise InputFormatError(
                f"{what} file is larger than {MAX_FILE_BYTES} bytes.",
                code="ASSEMBLY_INPUT_TOO_LARGE",
            )
        return path.read_text(encoding="utf-8", errors="replace")
    if isinstance(value, str):
        return value
    return None


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into ``(name, UPPERCASE sequence)``; the name is the first header word."""

    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts).upper()))
            header = line[1:].split()
            if not header:
                raise InputFormatError(
                    "FASTA header without a record name.",
                    code="ASSEMBLY_FASTA_FORMAT",
                    hint="Give every '>' header an identifier.",
                )
            name = header[0]
            parts = []
        else:
            if name is None:
                raise InputFormatError(
                    "FASTA sequence data before the first '>' header.",
                    code="ASSEMBLY_FASTA_FORMAT",
                    hint="The input must be FASTA text.",
                )
            parts.append(line)
    if name is not None:
        records.append((name, "".join(parts).upper()))
    return records


def _from_items(items: Sequence[object]) -> list[Contig]:
    contigs: list[Contig] = []
    for index, item in enumerate(items):
        default_name = f"contig_{index + 1}"
        if isinstance(item, bool):
            raise InputFormatError(
                f"contig {index} must be a sequence string or a positive length.",
                code="ASSEMBLY_CONTIG_FORMAT",
            )
        if isinstance(item, (int, float)):
            contigs.append((default_name, None, _positive_length(item, index)))
        elif isinstance(item, str):
            contigs.append((default_name, item.strip().upper(), len(item.strip())))
        elif isinstance(item, Mapping):
            name = str(item.get("name", item.get("id", default_name)))
            seq = item.get("sequence", item.get("seq"))
            if isinstance(seq, str):
                seq = seq.strip().upper()
                contigs.append((name, seq, len(seq)))
            else:
                length = item.get("length", item.get("size"))
                contigs.append((name, None, _positive_length(length, index)))
        else:
            raise InputFormatError(
                f"contig {index} has unsupported type {type(item).__name__}.",
                code="ASSEMBLY_CONTIG_FORMAT",
            )
    return contigs


def _positive_length(value: object, index: int) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or (isinstance(value, float) and not value.is_integer())
        or value <= 0
    ):
        raise InputFormatError(
            f"contig length {index} must be a positive integer.", code="ASSEMBLY_CONTIG_FORMAT"
        )
    return int(value)


def parse_contigs(
    value: object, *, require_sequences: bool = False, what: str = "assembly"
) -> list[Contig]:
    """Normalise assembly input into contigs.

    Accepted forms: FASTA text or file path; ``{name: sequence}``; ``{name: length}``;
    ``{"contigs"|"records"|"sequences"|"lengths": ...}``; a list of sequences, lengths or
    ``{"name", "sequence"|"length"}`` mappings.
    """

    contigs: list[Contig]
    text = read_text_input(value, what=what)
    if text is not None:
        stripped = text.lstrip()
        if stripped.startswith(">"):
            contigs = [(name, seq, len(seq)) for name, seq in parse_fasta(text)]
        elif re.fullmatch(r"[A-Za-z*\-\s]+", stripped or "!"):
            seq = "".join(stripped.split()).upper()
            contigs = [("contig_1", seq, len(seq))]
        else:
            raise InputFormatError(
                f"{what} text is neither FASTA nor a bare sequence.",
                code="ASSEMBLY_INPUT_FORMAT",
                hint="Pass FASTA text, a FASTA path, a mapping or a list of sequences/lengths.",
            )
    elif isinstance(value, Mapping):
        for key in ("contigs", "records", "sequences", "lengths", "contig_lengths", "assembly"):
            if key in value:
                return parse_contigs(value[key], require_sequences=require_sequences, what=what)
        if value and all(isinstance(v, str) for v in value.values()):
            contigs = [(str(k), v.strip().upper(), len(v.strip())) for k, v in value.items()]
        elif value and all(
            isinstance(v, (int, float)) and not isinstance(v, bool) for v in value.values()
        ):
            contigs = [
                (str(k), None, _positive_length(v, i)) for i, (k, v) in enumerate(value.items())
            ]
        else:
            raise InputFormatError(
                f"{what} mapping must map contig names to sequences or lengths.",
                code="ASSEMBLY_INPUT_FORMAT",
            )
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        contigs = _from_items(value)
    else:
        raise InputFormatError(
            f"{what} has unsupported type {type(value).__name__}.",
            code="ASSEMBLY_INPUT_FORMAT",
            hint="Pass FASTA text, a path, a mapping or a list.",
        )
    if not contigs:
        raise InputFormatError(f"{what} contains no contigs.", code="ASSEMBLY_EMPTY")
    if sum(length for _, _, length in contigs) > MAX_TOTAL_BASES:
        raise InputFormatError(
            f"{what} is larger than {MAX_TOTAL_BASES} bases.", code="ASSEMBLY_INPUT_TOO_LARGE"
        )
    if require_sequences and any(seq is None for _, seq, _ in contigs):
        raise InputFormatError(
            f"{what} must contain sequences, not only lengths.",
            code="ASSEMBLY_SEQUENCES_REQUIRED",
            hint="Pass FASTA text or a {name: sequence} mapping.",
        )
    return contigs


def options_dict(options: Mapping[str, object] | None) -> dict[str, object]:
    return dict(options or {})


def get_int(
    settings: Mapping[str, object], key: str, default: int | None, *, low: int, high: int
) -> int | None:
    value = settings.get(key, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(
            f"options.{key} must be an integer in [{low}, {high}].",
            code="ASSEMBLY_OPTION",
            hint=f"Set options['{key}'] to an integer in [{low}, {high}].",
        )
    return value


def get_float(
    settings: Mapping[str, object],
    key: str,
    default: float,
    *,
    low: float,
    high: float,
    low_open: bool = False,
) -> float:
    value = settings.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(f"options.{key} must be a number.", code="ASSEMBLY_OPTION")
    ok = (low < value if low_open else low <= value) and value <= high
    if not ok:
        raise ConfigurationError(
            f"options.{key} must be in {'(' if low_open else '['}{low}, {high}].",
            code="ASSEMBLY_OPTION",
        )
    return float(value)


# --------------------------------------------------------------------------- contig statistics


def sequence_composition(sequence: str) -> dict[str, int]:
    counts = Counter(sequence)
    acgt = {base: counts.get(base, 0) for base in "ACGT"}
    n_count = counts.get("N", 0)
    return {**acgt, "N": n_count, "other": len(sequence) - sum(acgt.values()) - n_count}


def gap_runs(sequence: str) -> list[tuple[int, int]]:
    """Maximal runs of N as 0-based half-open intervals."""

    return [(m.start(), m.end()) for m in re.finditer(r"N+", sequence)]


def assembly_stats(contigs: Sequence[Contig], genome_size: int | None = None) -> dict[str, Any]:
    lengths = [length for _, _, length in contigs]
    total = sum(lengths)
    n50, l50 = _nx(lengths, 50.0)
    n90, l90 = _nx(lengths, 90.0)
    stats: dict[str, Any] = {
        "contig_count": len(lengths),
        "total_length": total,
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": total / len(lengths),
        "n50": n50,
        "l50": l50,
        "n90": n90,
        "l90": l90,
        "aun": sum(length * length for length in lengths) / total,
    }
    if genome_size is not None:
        ng50, lg50 = _nx(lengths, 50.0, genome_size)
        stats["ng50"] = ng50
        stats["lg50"] = lg50
    if all(seq is not None for _, seq, _ in contigs):
        comp: Counter[str] = Counter()
        gaps = 0
        for _, seq, _ in contigs:
            comp.update(sequence_composition(seq or ""))  # type: ignore[arg-type]
            gaps += len(gap_runs(seq or ""))
        acgt = sum(comp[b] for b in "ACGT")
        stats["gc_fraction"] = (comp["G"] + comp["C"]) / acgt if acgt else None
        stats["n_bases"] = comp["N"]
        stats["gap_count"] = gaps
        stats["other_symbol_bases"] = comp["other"]
    return stats


def canonical_sequence(sequence: str) -> str:
    return min(sequence, reverse_complement(sequence))


def fasta_text(records: Iterable[tuple[str, str]], width: int = 80) -> str:
    lines: list[str] = []
    for name, seq in records:
        lines.append(f">{name}")
        lines.extend(seq[i : i + width] for i in range(0, len(seq), width))
    return "\n".join(lines) + ("\n" if lines else "")


def check_assembly_input(value: object) -> dict[str, Any]:
    """Validate assembly input and list every problem found."""

    errors: list[str] = []
    try:
        contigs = parse_contigs(value)
    except InputFormatError as exc:
        return {"valid": False, "format": None, "contig_count": 0, "errors": [exc.message]}
    names = [name for name, _, _ in contigs]
    duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicates:
        errors.append(f"duplicate contig names: {', '.join(duplicates[:20])}")
    has_sequences = all(seq is not None for _, seq, _ in contigs)
    invalid: dict[str, list[str]] = {}
    if has_sequences:
        for name, seq, _ in contigs:
            bad = sorted(set(re.sub(r"[ACGTNRYSWKMBDHV]", "", seq or "")))
            if bad:
                invalid[name] = bad
    if invalid:
        first = ", ".join(f"{n}: {''.join(b)}" for n, b in list(invalid.items())[:10])
        errors.append(f"non-IUPAC symbols in {len(invalid)} contig(s) ({first})")
    return {
        "valid": not errors,
        "format": "sequences" if has_sequences else "lengths",
        "contig_count": len(contigs),
        "total_length": sum(length for _, _, length in contigs),
        "duplicate_names": duplicates,
        "invalid_symbols": invalid,
        "errors": errors,
    }


# --------------------------------------------------------------------------- k-mers


def canonical_kmer_counts(sequences: Iterable[str], k: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    for seq in sequences:
        for segment in re.split(r"[^ACGT]+", seq.upper()):
            for i in range(len(segment) - k + 1):
                kmer = segment[i : i + k]
                counts[min(kmer, reverse_complement(kmer))] += 1
    return counts


def completeness(contigs: Sequence[Contig], settings: Mapping[str, object]) -> dict[str, Any]:
    stats = assembly_stats(contigs)
    total = stats["total_length"]
    result: dict[str, Any] = {"assembly_length": total, "contig_count": stats["contig_count"]}
    methods: list[str] = []
    reference_length = settings.get("reference_length")
    reference = settings.get("reference")
    if reference_length is not None:
        if (
            isinstance(reference_length, bool)
            or not isinstance(reference_length, (int, float))
            or reference_length <= 0
        ):
            raise ConfigurationError(
                "options.reference_length must be a positive number.", code="ASSEMBLY_OPTION"
            )
        ratio = total / float(reference_length)
        result.update(
            reference_length=int(reference_length),
            length_ratio=ratio,
            length_completeness=min(ratio, 1.0),
        )
        methods.append("length")
    if reference is not None:
        ref_contigs = parse_contigs(reference, require_sequences=True, what="options.reference")
        if any(seq is None for _, seq, _ in contigs):
            raise InputFormatError(
                "k-mer completeness needs assembly sequences, not only lengths.",
                code="ASSEMBLY_SEQUENCES_REQUIRED",
            )
        k = get_int(settings, "k", 21, low=3, high=63) or 21
        if sum(length for _, _, length in ref_contigs) + total > MAX_KMER_BASES:
            raise InputFormatError(
                f"k-mer completeness is limited to {MAX_KMER_BASES} total bases.",
                code="ASSEMBLY_INPUT_TOO_LARGE",
                hint="Use options.reference_length for a length-only estimate.",
            )
        ref_kmers = set(canonical_kmer_counts((s or "" for _, s, _ in ref_contigs), k))
        asm_kmers = set(canonical_kmer_counts((s or "" for _, s, _ in contigs), k))
        if not ref_kmers:
            raise InputFormatError(
                "reference has no ACGT k-mers of the requested length.",
                code="ASSEMBLY_REFERENCE_EMPTY",
            )
        shared = len(ref_kmers & asm_kmers)
        result.update(
            k=k,
            reference_distinct_kmers=len(ref_kmers),
            shared_kmers=shared,
            kmer_completeness=shared / len(ref_kmers),
        )
        methods.append("kmer")
    if not methods:
        raise ConfigurationError(
            "calculate_assembly_completeness needs options.reference_length or options.reference.",
            code="ASSEMBLY_OPTION",
            hint="Provide the expected genome length or a reference assembly.",
        )
    result["methods"] = methods
    return result


# --------------------------------------------------------------------------- BUSCO


_BUSCO_RE = re.compile(
    r"C:\s*([\d.]+)%\s*\[\s*S:\s*([\d.]+)%\s*,\s*D:\s*([\d.]+)%\s*\]\s*,\s*F:\s*([\d.]+)%\s*,\s*M:\s*([\d.]+)%\s*,\s*n:\s*(\d+)"
)


def busco_completeness(value: object) -> dict[str, Any]:
    counts: dict[str, int]
    if isinstance(value, Mapping):
        aliases = {
            "S": ("S", "single", "single_copy", "complete_single", "complete_single_copy"),
            "D": ("D", "duplicated", "complete_duplicated"),
            "F": ("F", "fragmented"),
            "M": ("M", "missing"),
        }
        counts = {}
        for letter, names in aliases.items():
            found = [value[n] for n in names if n in value]
            if not found:
                raise InputFormatError(
                    f"BUSCO counts need a '{letter}' entry ({'/'.join(names)}).",
                    code="BUSCO_COUNTS",
                    hint="Pass {'single': .., 'duplicated': .., 'fragmented': .., 'missing': ..}.",
                )
            number = found[0]
            if isinstance(number, bool) or not isinstance(number, int) or number < 0:
                raise InputFormatError(
                    f"BUSCO count '{letter}' must be a non-negative integer.", code="BUSCO_COUNTS"
                )
            counts[letter] = number
    else:
        text = read_text_input(value, what="BUSCO summary")
        match = _BUSCO_RE.search(text or "")
        if match is None:
            raise InputFormatError(
                "no BUSCO short-summary line (C:..%[S:..%,D:..%],F:..%,M:..%,n:..) found.",
                code="BUSCO_SUMMARY",
                hint="Pass the short_summary text/file or a mapping of counts.",
            )
        n = int(match.group(6))
        counts = {
            key: round(float(match.group(i)) * n / 100.0)
            for key, i in (("S", 2), ("D", 3), ("F", 4), ("M", 5))
        }
        counts["_n"] = n
    n_total = counts.pop("_n", None) or sum(counts.values())
    if n_total <= 0:
        raise InputFormatError("BUSCO total (n) must be positive.", code="BUSCO_COUNTS")
    complete = counts["S"] + counts["D"]
    pct = lambda c: 100.0 * c / n_total  # noqa: E731
    return {
        "n": n_total,
        "single": counts["S"],
        "duplicated": counts["D"],
        "complete": complete,
        "fragmented": counts["F"],
        "missing": counts["M"],
        "complete_percent": pct(complete),
        "single_percent": pct(counts["S"]),
        "duplicated_percent": pct(counts["D"]),
        "fragmented_percent": pct(counts["F"]),
        "missing_percent": pct(counts["M"]),
        "summary": f"C:{pct(complete):.1f}%[S:{pct(counts['S']):.1f}%,D:{pct(counts['D']):.1f}%],F:{pct(counts['F']):.1f}%,M:{pct(counts['M']):.1f}%,n:{n_total}",
    }


# --------------------------------------------------------------------------- k-mer histogram models


def kmer_histogram(value: object, settings: Mapping[str, object]) -> dict[int, int]:
    """Coverage -> number of distinct k-mers, from a histogram or (bounded) sequences."""

    if isinstance(value, Mapping) and "histogram" in value:
        value = value["histogram"]
    hist: dict[int, int] = {}

    def add(cov: object, count: object) -> None:
        try:
            c = int(float(cov))  # type: ignore[arg-type]
            n = int(float(count))  # type: ignore[arg-type]
        except (TypeError, ValueError):
            raise InputFormatError(
                "histogram entries must be numeric (coverage, count).", code="KMER_HISTOGRAM"
            ) from None
        if c < 1 or n < 0:
            raise InputFormatError(
                "histogram coverage must be >= 1 and counts >= 0.", code="KMER_HISTOGRAM"
            )
        hist[c] = hist.get(c, 0) + n

    text = read_text_input(value, what="k-mer input")
    if isinstance(value, Mapping) and value and all(re.fullmatch(r"\d+", str(k)) for k in value):
        for cov, count in value.items():
            add(cov, count)
    elif (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and value
        and all(
            isinstance(row, Sequence) and not isinstance(row, (str, bytes)) and len(row) == 2
            for row in value
        )
    ):
        for cov, count in value:
            add(cov, count)
    elif text is not None and re.fullmatch(
        r"(?:\s*\d+(?:\.\d+)?[ \t,]+\d+(?:\.\d+)?\s*\n?)+", text
    ):
        for line in text.splitlines():
            if line.strip():
                a, b = re.split(r"[ \t,]+", line.strip())
                add(a, b)
    else:
        k = get_int(settings, "k", 21, low=3, high=63) or 21
        if isinstance(value, Mapping) or text is not None or isinstance(value, Sequence):
            contigs = parse_contigs(value, require_sequences=True, what="k-mer input")
        else:
            raise InputFormatError("unsupported k-mer histogram input.", code="KMER_HISTOGRAM")
        if sum(length for _, _, length in contigs) > MAX_KMER_BASES:
            raise InputFormatError(
                f"counting k-mers is limited to {MAX_KMER_BASES} bases; pass a precomputed histogram instead.",
                code="ASSEMBLY_INPUT_TOO_LARGE",
            )
        counts = canonical_kmer_counts((s or "" for _, s, _ in contigs), k)
        for multiplicity in counts.values():
            hist[multiplicity] = hist.get(multiplicity, 0) + 1
    if not hist or sum(hist.values()) == 0:
        raise InputFormatError("the k-mer histogram is empty.", code="KMER_HISTOGRAM")
    return hist


def histogram_peak(hist: Mapping[int, int], settings: Mapping[str, object]) -> tuple[int, int]:
    """Return ``(min_coverage, peak_coverage)``: user-given, else first valley and highest bin after it."""

    max_cov = max(hist)
    min_cov = get_int(settings, "min_coverage", None, low=1, high=10**9)
    peak = get_int(settings, "peak_coverage", None, low=1, high=10**9)
    if min_cov is None:
        for x in range(min(hist), max_cov):
            if hist.get(x, 0) <= hist.get(x + 1, 0):
                min_cov = x
                break
        if min_cov is None:
            raise InputFormatError(
                "no error/valley cutoff found: the histogram decreases monotonically.",
                code="KMER_HISTOGRAM_PEAK",
                hint="Set options.min_coverage and options.peak_coverage.",
            )
    if peak is None:
        candidates = (
            [x for x in hist if x > min_cov]
            if "min_coverage" not in settings
            else [x for x in hist if x >= min_cov]
        )
        if not candidates:
            raise InputFormatError(
                "no k-mer coverage above the cutoff.",
                code="KMER_HISTOGRAM_PEAK",
                hint="Set options.peak_coverage.",
            )
        peak = max(candidates, key=lambda x: (hist[x], -x))
    return min_cov, peak


def estimate_genome_size(hist: Mapping[int, int], settings: Mapping[str, object]) -> dict[str, Any]:
    min_cov, peak = histogram_peak(hist, settings)
    max_cov = get_int(settings, "max_coverage", None, low=1, high=10**12)
    kept = {x: n for x, n in hist.items() if x >= min_cov and (max_cov is None or x <= max_cov)}
    total_kmers = sum(x * n for x, n in kept.items())
    size = total_kmers / peak
    return {
        "genome_size": size,
        "peak_coverage": peak,
        "min_coverage": min_cov,
        "max_coverage": max_cov,
        "total_kmer_instances": total_kmers,
        "distinct_kmers_used": sum(kept.values()),
        "model": "G = sum_{x>=min_coverage} x * n(x) / peak_coverage",
    }


def estimate_repeat_content(
    hist: Mapping[int, int], settings: Mapping[str, object]
) -> dict[str, Any]:
    min_cov, peak = histogram_peak(hist, settings)
    factor = get_float(settings, "repeat_factor", 1.5, low=1.0, high=100.0, low_open=True)
    kept = {x: n for x, n in hist.items() if x >= min_cov}
    total = sum(x * n for x, n in kept.items())
    repeat = sum(x * n for x, n in kept.items() if x > factor * peak)
    classes: dict[int, int] = {}
    for x, n in kept.items():
        if x > factor * peak:
            copies = max(2, round(x / peak))
            classes[copies] = classes.get(copies, 0) + n
    return {
        "repeat_fraction": repeat / total if total else 0.0,
        "repeat_positions_estimate": repeat / peak,
        "genome_size_estimate": total / peak,
        "repeat_distinct_kmers": sum(classes.values()),
        "copy_number_classes": {str(c): classes[c] for c in sorted(classes)},
        "peak_coverage": peak,
        "min_coverage": min_cov,
        "repeat_factor": factor,
        "definition": "repeat_fraction = sum_{x > repeat_factor*peak} x*n(x) / sum_{x >= min_coverage} x*n(x)",
    }


def estimate_het(hist: Mapping[int, int], settings: Mapping[str, object]) -> dict[str, Any]:
    k = get_int(settings, "k", 21, low=3, high=127) or 21
    min_cov, peak = histogram_peak(hist, settings)
    het_upper = get_float(settings, "het_upper", 0.75, low=0.0, high=1.0, low_open=True)
    hom_upper = get_float(settings, "hom_upper", 1.5, low=1.0, high=100.0, low_open=True)
    het = sum(n for x, n in hist.items() if min_cov <= x <= het_upper * peak)
    hom = sum(n for x, n in hist.items() if het_upper * peak < x <= hom_upper * peak)
    if hom == 0:
        raise InputFormatError(
            "no k-mers in the homozygous window; cannot estimate heterozygosity.",
            code="KMER_HISTOGRAM_PEAK",
            hint="Check options.peak_coverage.",
        )
    ratio = het / hom
    rate = ratio / (k * (2.0 + ratio))
    return {
        "heterozygosity": rate,
        "het_distinct_kmers": het,
        "hom_distinct_kmers": hom,
        "het_to_hom_ratio": ratio,
        "k": k,
        "peak_coverage": peak,
        "min_coverage": min_cov,
        "model": "diploid, isolated SNPs spaced > k apart: het k-mers = 2*k*r*L, hom k-mers = L*(1-k*r); r = R / (k*(2+R)) with R = het/hom",
    }


# --------------------------------------------------------------------------- pangenome matrices


def presence_matrix(
    value: object, *, what: str = "pangenome"
) -> tuple[list[str], list[str], list[list[int]]]:
    """Return ``(genomes, genes, rows)`` where ``rows[g][j]`` is 1 if gene g is present in genome j."""

    if not isinstance(value, Mapping) or not value:
        raise InputFormatError(
            f"{what} must be a non-empty mapping (see the docstring for accepted forms).",
            code="PANGENOME_INPUT",
            hint="Pass {genome: [genes]} or {'genomes': [...], 'genes': [...], 'matrix': [[...]]}.",
        )
    if "matrix" in value and "genomes" in value and "genes" in value:
        genomes = [str(g) for g in value["genomes"]]  # type: ignore[index, union-attr]
        genes = [str(g) for g in value["genes"]]  # type: ignore[index, union-attr]
        matrix = value["matrix"]  # type: ignore[index]
        if not isinstance(matrix, Sequence) or len(matrix) != len(genes):
            raise InputFormatError("matrix needs one row per gene.", code="PANGENOME_INPUT")
        rows: list[list[int]] = []
        for row in matrix:
            if not isinstance(row, Sequence) or len(row) != len(genomes):
                raise InputFormatError(
                    "every matrix row needs one value per genome.", code="PANGENOME_INPUT"
                )
            rows.append([_presence(v) for v in row])
    elif all(isinstance(v, (list, tuple, set, frozenset)) for v in value.values()):
        genomes = [str(g) for g in value]
        gene_sets = [set(str(x) for x in v) for v in value.values()]
        genes = sorted(set().union(*gene_sets))
        rows = [[1 if gene in gs else 0 for gs in gene_sets] for gene in genes]
    elif all(isinstance(v, Mapping) for v in value.values()):
        genes = [str(g) for g in value]
        genomes = []
        for entry in value.values():
            for genome in entry:  # type: ignore[union-attr]
                if str(genome) not in genomes:
                    genomes.append(str(genome))
        rows = [
            [_presence(entry.get(g, entry.get(str(g), 0))) for g in genomes]
            for entry in value.values()
        ]  # type: ignore[union-attr]
    else:
        raise InputFormatError(
            f"unrecognised {what} layout.",
            code="PANGENOME_INPUT",
            hint="Pass {genome: [genes]}, {gene: {genome: 0/1}} or the genomes/genes/matrix form.",
        )
    if len(set(genomes)) != len(genomes) or len(set(genes)) != len(genes):
        raise InputFormatError("genome and gene names must be unique.", code="PANGENOME_INPUT")
    if not genomes or not genes:
        raise InputFormatError("the matrix has no genomes or no genes.", code="PANGENOME_INPUT")
    if len(genomes) > MAX_GENOMES or len(genes) > MAX_GENES:
        raise InputFormatError(
            f"the matrix exceeds {MAX_GENOMES} genomes or {MAX_GENES} genes.",
            code="ASSEMBLY_INPUT_TOO_LARGE",
        )
    return genomes, genes, rows


def _presence(value: object) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (int, float)) and value >= 0:
        return 1 if value > 0 else 0
    raise InputFormatError(
        "presence/absence values must be booleans or non-negative numbers.", code="PANGENOME_INPUT"
    )


def matrix_dict(
    genomes: Sequence[str], genes: Sequence[str], rows: Sequence[Sequence[int]]
) -> dict[str, Any]:
    return {"genomes": list(genomes), "genes": list(genes), "matrix": [list(r) for r in rows]}


def pangenome_partition(
    genomes: Sequence[str],
    genes: Sequence[str],
    rows: Sequence[Sequence[int]],
    core_threshold: float,
    cloud_threshold: float = 0.15,
) -> dict[str, Any]:
    n = len(genomes)
    needed = math.ceil(core_threshold * n - 1e-9)
    core: list[str] = []
    shell: list[str] = []
    cloud: list[str] = []
    fractions: dict[str, float] = {}
    for gene, row in zip(genes, rows, strict=False):
        present = sum(row)
        fractions[gene] = present / n
        if present >= needed:
            core.append(gene)
        elif present / n >= cloud_threshold:
            shell.append(gene)
        else:
            cloud.append(gene)
    return {
        "genome_count": n,
        "total_genes": len(genes),
        "core_threshold": core_threshold,
        "core_min_genomes": needed,
        "core_genes": core,
        "core_count": len(core),
        "accessory_genes": shell + cloud,
        "accessory_count": len(shell) + len(cloud),
        "shell_genes": shell,
        "cloud_genes": cloud,
        "cloud_threshold": cloud_threshold,
        "gene_fractions": fractions,
    }


def core_threshold_option(settings: Mapping[str, object]) -> float:
    return get_float(settings, "core_threshold", 0.99, low=0.0, high=1.0, low_open=True)


def new_gene_deltas(gene_sets: Sequence[set[int]], order: Sequence[int]) -> list[int]:
    seen: set[int] = set()
    deltas: list[int] = []
    for index in order:
        new = gene_sets[index] - seen
        deltas.append(len(new))
        seen |= gene_sets[index]
    return deltas


def fit_heaps(deltas: Sequence[float]) -> dict[str, Any]:
    """Fit delta(N) = kappa * N**(-alpha) by ordinary least squares on log-log, N = 2..G."""

    points = [(math.log(n), math.log(d)) for n, d in enumerate(deltas, start=2) if d > 0]
    if len(points) < 3:
        raise InputFormatError(
            "Heaps' law needs at least three genomes with new genes (>= 4 genomes).",
            code="PANGENOME_TOO_FEW_GENOMES",
            hint="Provide at least four genomes.",
        )
    mx = sum(p[0] for p in points) / len(points)
    my = sum(p[1] for p in points) / len(points)
    sxx = sum((p[0] - mx) ** 2 for p in points)
    slope = sum((p[0] - mx) * (p[1] - my) for p in points) / sxx
    intercept = my - slope * mx
    return {
        "alpha": -slope,
        "kappa": math.exp(intercept),
        "points_used": len(points),
        "points_dropped": len(deltas) - len(points),
    }


def pangenome_openness(
    genomes: Sequence[str],
    genes: Sequence[str],
    rows: Sequence[Sequence[int]],
    permutations: int,
    seed: int,
) -> dict[str, Any]:
    import random

    g = len(genomes)
    if g < 4:
        raise InputFormatError(
            "pangenome openness needs at least four genomes.",
            code="PANGENOME_TOO_FEW_GENOMES",
            hint="Provide at least four genomes.",
        )
    gene_sets: list[set[int]] = [set() for _ in range(g)]
    for gi, row in enumerate(rows):
        for j, v in enumerate(row):
            if v:
                gene_sets[j].add(gi)
    if permutations == 0:
        orders = [list(range(g))]
    else:
        rng = random.Random(seed)
        orders = []
        for _ in range(permutations):
            order = list(range(g))
            rng.shuffle(order)
            orders.append(order)
    totals = [0.0] * g
    for order in orders:
        for i, d in enumerate(new_gene_deltas(gene_sets, order)):
            totals[i] += d
    mean_deltas = [t / len(orders) for t in totals]
    fit = fit_heaps(mean_deltas[1:])
    pan = [sum(mean_deltas[: i + 1]) for i in range(g)]
    return {
        **fit,
        "open": fit["alpha"] < 1.0,
        "genome_count": g,
        "permutations": len(orders),
        "seed": seed if permutations else None,
        "mean_new_genes_by_genome_number": mean_deltas,
        "mean_pangenome_size_by_genome_number": pan,
        "model": "Heaps' law: mean number of new genes found when the N-th genome is added = kappa * N^(-alpha), fitted by OLS on log-log for N=2..G; open if alpha < 1, closed if alpha >= 1 (Tettelin et al. 2008; Rasko et al. 2008)",
    }


# --------------------------------------------------------------------------- repeats, telomeres


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def repeats_in_assembly(
    contigs: Sequence[Contig], k: int, min_copies: int, min_length: int
) -> dict[str, Any]:
    if sum(length for _, _, length in contigs) > MAX_KMER_BASES:
        raise InputFormatError(
            f"k-mer repeat detection is limited to {MAX_KMER_BASES} bases.",
            code="ASSEMBLY_INPUT_TOO_LARGE",
        )
    counts = canonical_kmer_counts((s or "" for _, s, _ in contigs), k)
    per_contig: dict[str, list[dict[str, int]]] = {}
    repeat_bases = 0
    total = 0
    for name, seq, length in contigs:
        total += length
        seq = seq or ""
        spans: list[tuple[int, int]] = []
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if (
                set(kmer) <= set("ACGT")
                and counts[min(kmer, reverse_complement(kmer))] >= min_copies
            ):
                spans.append((i, i + k))
        intervals = [(a, b) for a, b in merge_intervals(spans) if b - a >= min_length]
        per_contig[name] = [{"start": a, "end": b, "length": b - a} for a, b in intervals]
        repeat_bases += sum(b - a for a, b in intervals)
    return {
        "k": k,
        "min_copies": min_copies,
        "min_length": min_length,
        "repeat_bases": repeat_bases,
        "total_bases": total,
        "repeat_fraction": repeat_bases / total if total else 0.0,
        "repeats": per_contig,
        "definition": "positions covered by canonical k-mers that occur >= min_copies times in the assembly; overlapping/adjacent k-mers merged; 0-based half-open intervals",
    }


def _longest_run(runs: list[tuple[int, int, int]]) -> tuple[int, int, int] | None:
    """Return the longest (start, end, copies) run, the leftmost on ties, or None."""

    return max(runs, key=lambda run: (run[1] - run[0], -run[0])) if runs else None


def telomeres(
    contigs: Sequence[Contig], motif: str, min_repeats: int, max_distance: int
) -> dict[str, Any]:
    if not re.fullmatch(r"[ACGT]{2,12}", motif):
        raise ConfigurationError(
            "options.motif must be 2-12 bases of ACGT.", code="ASSEMBLY_OPTION"
        )
    rotations = {motif[i:] + motif[:i] for i in range(len(motif))}
    rc = reverse_complement(motif)
    rotations |= {rc[i:] + rc[:i] for i in range(len(rc))}
    units = sorted(rotations)
    results: dict[str, dict[str, Any]] = {}
    for name, seq, length in contigs:
        seq = seq or ""
        runs: list[tuple[int, int, int]] = []
        for unit in units:
            for m in re.finditer(f"(?:{unit}){{{min_repeats},}}", seq):
                runs.append((m.start(), m.end(), (m.end() - m.start()) // len(unit)))
        left = [r for r in runs if r[0] <= max_distance and r[0] <= len(seq) - r[1]]
        right = [r for r in runs if len(seq) - r[1] <= max_distance and len(seq) - r[1] < r[0]]
        entry: dict[str, Any] = {"length": length}
        for side, rs in (("left", left), ("right", right)):
            best = _longest_run(rs)
            entry[side] = (
                None if best is None else {"start": best[0], "end": best[1], "copies": best[2]}
            )
        entry["telomeres_at_ends"] = int(entry["left"] is not None) + int(
            entry["right"] is not None
        )
        results[name] = entry
    return {
        "motif": motif,
        "min_repeats": min_repeats,
        "max_distance": max_distance,
        "contigs": results,
        "contigs_with_telomere": sum(1 for e in results.values() if e["telomeres_at_ends"]),
        "contigs_with_both_ends": sum(1 for e in results.values() if e["telomeres_at_ends"] == 2),
        "definition": "a telomere end is a perfect tandem run of >= min_repeats copies of the motif or its reverse complement (any rotation) lying within max_distance bases of the nearer contig end (ties go to the left end)",
    }


# --------------------------------------------------------------------------- GFA, de Bruijn graphs


def parse_gfa(text: str) -> dict[str, Any]:
    segments: dict[str, dict[str, Any]] = {}
    links: list[tuple[str, str, str, str, str]] = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        kind = fields[0]
        if kind == "S":
            if len(fields) < 3:
                raise InputFormatError(
                    f"GFA line {lineno}: S line needs a name and sequence.", code="GFA_FORMAT"
                )
            name, seq = fields[1], fields[2]
            tags = {f[:2]: f[5:] for f in fields[3:] if len(f) > 5 and f[2] == ":"}
            length = len(seq) if seq != "*" else int(tags["LN"]) if "LN" in tags else None
            depth = None
            for tag in ("DP", "dp", "KC"):
                if tag in tags:
                    try:
                        depth = float(tags[tag])
                    except ValueError:
                        depth = None
                    break
            segments[name] = {
                "sequence": None if seq == "*" else seq.upper(),
                "length": length,
                "depth": depth,
            }
        elif kind == "L":
            if len(fields) < 6 or fields[2] not in "+-" or fields[4] not in "+-":
                raise InputFormatError(f"GFA line {lineno}: malformed L line.", code="GFA_FORMAT")
            links.append((fields[1], fields[2], fields[3], fields[4], fields[5]))
    if not segments:
        raise InputFormatError(
            "no GFA segments (S lines) found.",
            code="GFA_FORMAT",
            hint="Pass GFA1 text or a path to a .gfa file.",
        )
    if len(segments) > 1_000_000:
        raise InputFormatError("GFA has too many segments.", code="ASSEMBLY_INPUT_TOO_LARGE")
    for a, _, b, _, _ in links:
        if a not in segments or b not in segments:
            raise InputFormatError(
                f"link refers to unknown segment {a if a not in segments else b}.",
                code="GFA_FORMAT",
            )
    return {"segments": segments, "links": links}


def load_gfa(value: object) -> dict[str, Any]:
    text = read_text_input(value, what="GFA")
    if text is None:
        raise InputFormatError("GFA input must be text or a path.", code="GFA_FORMAT")
    return parse_gfa(text)


def _link_key(a: str, ao: str, b: str, bo: str) -> tuple[str, str, str, str]:
    flip = {"+": "-", "-": "+"}
    forward = (a, ao, b, bo)
    reverse = (b, flip[bo], a, flip[ao])
    return min(forward, reverse)


def compare_gfa(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    use_sequences = all(s["sequence"] for g in (left, right) for s in g["segments"].values())

    def translate(graph: dict[str, Any]) -> dict[str, tuple[str, bool]]:
        """Segment name -> (matching key, whether the stored sequence is the reverse of the key's orientation)."""

        out: dict[str, tuple[str, bool]] = {}
        for name, seg in graph["segments"].items():
            if use_sequences:
                canon = canonical_sequence(seg["sequence"])
                out[name] = (canon, canon != seg["sequence"])
            else:
                out[name] = (name, False)
        return out

    flip = {"+": "-", "-": "+"}

    def link_set(
        graph: dict[str, Any], trans: dict[str, tuple[str, bool]]
    ) -> set[tuple[str, str, str, str]]:
        keys = set()
        for a, ao, b, bo, _ in graph["links"]:
            (ka, fa), (kb, fb) = trans[a], trans[b]
            ao2, bo2 = (flip[ao] if fa else ao), (flip[bo] if fb else bo)
            keys.add(_link_key(ka, ao2, kb, bo2))
        return keys

    ltrans, rtrans = translate(left), translate(right)
    lset = {key for key, _ in ltrans.values()}
    rset = {key for key, _ in rtrans.values()}
    llinks, rlinks = link_set(left, ltrans), link_set(right, rtrans)

    def jaccard(a: set[Any], b: set[Any]) -> float | None:
        union = a | b
        return len(a & b) / len(union) if union else None

    def total(graph: dict[str, Any]) -> int:
        return sum(s["length"] or 0 for s in graph["segments"].values())

    return {
        "compared_by": "canonical segment sequence" if use_sequences else "segment name",
        "left": {
            "segments": len(left["segments"]),
            "links": len(left["links"]),
            "total_segment_length": total(left),
        },
        "right": {
            "segments": len(right["segments"]),
            "links": len(right["links"]),
            "total_segment_length": total(right),
        },
        "shared_segments": len(lset & rset),
        "segment_jaccard": jaccard(lset, rset),
        "left_only_segments": len(lset - rset),
        "right_only_segments": len(rset - lset),
        "shared_links": len(llinks & rlinks),
        "link_jaccard": jaccard(llinks, rlinks),
        "identical": lset == rset and llinks == rlinks,
    }


def graph_svg(graph: dict[str, Any], size: int = 640) -> str:
    names = sorted(graph["segments"])
    if len(names) > MAX_GRAPH_NODES:
        raise InputFormatError(
            f"plot_graph is limited to {MAX_GRAPH_NODES} segments.", code="ASSEMBLY_INPUT_TOO_LARGE"
        )
    centre = size / 2
    radius = size / 2 - 40
    pos: dict[str, tuple[float, float]] = {}
    for i, name in enumerate(names):
        angle = 2 * math.pi * i / len(names)
        pos[name] = (centre + radius * math.cos(angle), centre + radius * math.sin(angle))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" viewBox="0 0 {size} {size}">'
    ]
    for a, _, b, _, _ in graph["links"]:
        (x1, y1), (x2, y2) = pos[a], pos[b]
        parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#888" stroke-width="1"/>'
        )
    for name in names:
        x, y = pos[name]
        length = graph["segments"][name]["length"] or 1
        r = 4 + 2 * math.log10(length + 1)
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="#4c78a8"><title>{escape(name)} ({length} bp)</title></circle>'
        )
        if len(names) <= 60:
            parts.append(
                f'<text x="{x + r + 2:.1f}" y="{y + 3:.1f}" font-size="9">{escape(name)}</text>'
            )
    parts.append("</svg>")
    return "\n".join(parts)


def de_bruijn_unitigs(contigs: Sequence[Contig], k: int) -> dict[str, Any]:
    """Compacted (forward-strand) de Bruijn graph of the input sequences."""

    if sum(length for _, _, length in contigs) > 2_000_000:
        raise InputFormatError(
            "build_graph_genome is limited to 2,000,000 input bases.",
            code="ASSEMBLY_INPUT_TOO_LARGE",
        )
    succ: dict[str, set[str]] = {}
    pred: dict[str, set[str]] = {}
    sources: dict[str, set[str]] = {}
    for name, seq, _ in contigs:
        for segment in re.split(r"[^ACGT]+", (seq or "").upper()):
            previous: str | None = None
            for i in range(len(segment) - k + 1):
                kmer = segment[i : i + k]
                succ.setdefault(kmer, set())
                pred.setdefault(kmer, set())
                sources.setdefault(kmer, set()).add(name)
                if previous is not None:
                    succ[previous].add(kmer)
                    pred[kmer].add(previous)
                previous = kmer
    if not succ:
        raise InputFormatError(
            f"no ACGT k-mers of length {k} in the input.",
            code="ASSEMBLY_NO_KMERS",
            hint="Use a smaller options.k.",
        )

    def is_start(kmer: str) -> bool:
        if len(pred[kmer]) != 1:
            return True
        (p,) = pred[kmer]
        return len(succ[p]) != 1

    unitigs: list[list[str]] = []
    owner: dict[str, int] = {}
    for kmer in sorted(succ):
        if is_start(kmer):
            path = [kmer]
            while len(succ[path[-1]]) == 1:
                (nxt,) = succ[path[-1]]
                if len(pred[nxt]) != 1 or nxt in owner or nxt == kmer:
                    break
                path.append(nxt)
            for node in path:
                owner[node] = len(unitigs)
            unitigs.append(path)
    for kmer in sorted(succ):  # isolated cycles
        if kmer not in owner:
            path = [kmer]
            owner[kmer] = len(unitigs)
            cursor = kmer
            while True:
                (nxt,) = succ[cursor]
                if nxt in owner:
                    break
                owner[nxt] = len(unitigs)
                path.append(nxt)
                cursor = nxt
            unitigs.append(path)
    records = []
    for i, path in enumerate(unitigs):
        seq = path[0] + "".join(node[-1] for node in path[1:])
        colors = sorted(set().union(*(sources[n] for n in path)))
        records.append(
            {
                "name": f"utg{i + 1}",
                "sequence": seq,
                "length": len(seq),
                "kmers": len(path),
                "sources": colors,
            }
        )
    links: list[tuple[str, str]] = []
    for i, path in enumerate(unitigs):
        for nxt in sorted(succ[path[-1]]):
            links.append((records[i]["name"], records[owner[nxt]]["name"]))
    gfa = ["H\tVN:Z:1.0"]
    gfa += [f"S\t{r['name']}\t{r['sequence']}\tLN:i:{r['length']}" for r in records]
    gfa += [f"L\t{a}\t+\t{b}\t+\t{k - 1}M" for a, b in links]
    return {
        "k": k,
        "kmer_count": len(succ),
        "unitig_count": len(records),
        "link_count": len(links),
        "unitigs": records,
        "links": [{"from": a, "to": b, "overlap": k - 1} for a, b in links],
        "gfa": "\n".join(gfa) + "\n",
        "note": "forward-strand (non-canonical) de Bruijn graph; nodes are k-mers, edges are observed (k-1)-overlaps; unitigs are maximal non-branching paths",
    }


# --------------------------------------------------------------------------- scaffolding, filtering, parsing


def scaffold_layout(contigs: Sequence[Contig], layout: object, default_gap: int) -> dict[str, Any]:
    by_name = {name: seq for name, seq, _ in contigs}
    if len(by_name) != len(contigs):
        raise InputFormatError(
            "contig names must be unique to scaffold.", code="ASSEMBLY_DUPLICATE_NAMES"
        )
    if not isinstance(layout, Sequence) or isinstance(layout, (str, bytes)) or not layout:
        raise ConfigurationError(
            "options.layout must be a non-empty list of scaffold definitions.",
            code="ASSEMBLY_OPTION",
            hint="Use [{'name': 'scaf1', 'contigs': [{'contig': 'c1', 'orientation': '+', 'gap': 100}, 'c2']}].",
        )
    used: set[str] = set()
    scaffolds: list[dict[str, Any]] = []
    agp: list[dict[str, Any]] = []
    for sindex, definition in enumerate(layout):
        if not isinstance(definition, Mapping) or "contigs" not in definition:
            raise ConfigurationError(
                f"layout entry {sindex} must be a mapping with a 'contigs' list.",
                code="ASSEMBLY_OPTION",
            )
        sname = str(definition.get("name", f"scaffold_{sindex + 1}"))
        pieces: list[str] = []
        position = 0
        members = definition["contigs"]
        for mindex, member in enumerate(members):  # type: ignore[arg-type]
            entry = member if isinstance(member, Mapping) else {"contig": member}
            cname = str(entry.get("contig", entry.get("name")))
            orient = str(entry.get("orientation", "+"))
            if cname not in by_name or by_name[cname] is None:
                raise ConfigurationError(
                    f"layout refers to unknown contig {cname!r}.", code="ASSEMBLY_OPTION"
                )
            if cname in used:
                raise ConfigurationError(
                    f"contig {cname!r} is used more than once in the layout.",
                    code="ASSEMBLY_OPTION",
                )
            if orient not in ("+", "-"):
                raise ConfigurationError("orientation must be '+' or '-'.", code="ASSEMBLY_OPTION")
            used.add(cname)
            seq = by_name[cname] or ""
            seq = seq if orient == "+" else reverse_complement(seq)
            agp.append(
                {
                    "scaffold": sname,
                    "start": position,
                    "end": position + len(seq),
                    "type": "contig",
                    "contig": cname,
                    "orientation": orient,
                }
            )
            pieces.append(seq)
            position += len(seq)
            if mindex < len(members) - 1:  # type: ignore[arg-type]
                gap = entry.get("gap", default_gap)
                if isinstance(gap, bool) or not isinstance(gap, int) or gap < 0:
                    raise ConfigurationError(
                        "gap must be a non-negative integer.", code="ASSEMBLY_OPTION"
                    )
                if gap:
                    agp.append(
                        {
                            "scaffold": sname,
                            "start": position,
                            "end": position + gap,
                            "type": "gap",
                            "contig": None,
                            "orientation": None,
                        }
                    )
                    pieces.append("N" * gap)
                    position += gap
        scaffolds.append(
            {
                "name": sname,
                "sequence": "".join(pieces),
                "length": position,
                "contig_count": len(members),
            }
        )  # type: ignore[arg-type]
    unplaced = [name for name, _, _ in contigs if name not in used]
    return {"scaffolds": scaffolds, "agp": agp, "unplaced_contigs": unplaced}


_SPADES_RE = re.compile(r"NODE_(\d+)_length_(\d+)_cov_([\d.]+)")
_MEGAHIT_RE = re.compile(r"flag=(\d+)\s+multi=([\d.]+)\s+len=(\d+)")


def parse_assembler_output(value: object) -> dict[str, Any]:
    text = read_text_input(value, what="assembler output")
    if text is None:
        raise InputFormatError(
            "assembler output must be FASTA/GFA text or a path.", code="ASSEMBLY_INPUT_FORMAT"
        )
    first = next((line for line in text.splitlines() if line.strip()), "")
    if first.startswith(">"):
        records = []
        name = None
        header = ""
        seq_parts: list[str] = []
        entries: list[tuple[str, str, str]] = []
        for line in text.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    entries.append((name, header, "".join(seq_parts).upper()))
                header = line[1:]
                name = header.split()[0] if header.split() else ""
                seq_parts = []
            elif line:
                seq_parts.append(line)
        if name is not None:
            entries.append((name, header, "".join(seq_parts).upper()))
        convention = None
        for name, header, seq in entries:
            rec: dict[str, Any] = {"name": name, "length": len(seq)}
            acgt = sum(seq.count(b) for b in "ACGT")
            rec["gc_fraction"] = (seq.count("G") + seq.count("C")) / acgt if acgt else None
            m = _SPADES_RE.search(header)
            if m:
                convention = "spades"
                rec.update(declared_length=int(m.group(2)), coverage=float(m.group(3)))
            m = _MEGAHIT_RE.search(header)
            if m:
                convention = "megahit"
                rec.update(declared_length=int(m.group(3)), multiplicity=float(m.group(2)))
            records.append(rec)
        if not records:
            raise InputFormatError("no FASTA records found.", code="ASSEMBLY_FASTA_FORMAT")
        contigs = [(r["name"], None, r["length"]) for r in records if r["length"] > 0]
        return {
            "format": "fasta",
            "header_convention": convention,
            "records": records,
            "stats": assembly_stats(contigs) if contigs else None,
        }
    graph = parse_gfa(text)
    records = [
        {"name": n, "length": s["length"], "depth": s["depth"]}
        for n, s in graph["segments"].items()
    ]
    lengths = [(r["name"], None, r["length"]) for r in records if r["length"]]
    return {
        "format": "gfa",
        "header_convention": None,
        "records": records,
        "link_count": len(graph["links"]),
        "stats": assembly_stats(lengths) if lengths else None,
    }


_LEVELS = ("CRITICAL", "FATAL", "ERROR", "WARNING", "WARN", "INFO", "DEBUG")
_LEVEL_RE = re.compile(r"\b(" + "|".join(_LEVELS) + r")\b", re.IGNORECASE)


def parse_log_text(text: str, max_lines: int = 20) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    errors: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    total = 0
    for number, line in enumerate(text.splitlines(), start=1):
        total += 1
        match = _LEVEL_RE.search(line)
        if not match:
            continue
        level = match.group(1).upper()
        level = {"WARN": "WARNING", "FATAL": "CRITICAL"}.get(level, level)
        counts[level] += 1
        entry = {"line": number, "text": line.strip()[:500]}
        if level in ("ERROR", "CRITICAL") and len(errors) < max_lines:
            errors.append(entry)
        elif level == "WARNING" and len(warnings) < max_lines:
            warnings.append(entry)
    return {
        "line_count": total,
        "level_counts": {
            level: counts.get(level, 0)
            for level in ("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG")
        },
        "has_errors": bool(counts.get("ERROR") or counts.get("CRITICAL")),
        "errors": errors,
        "warnings": warnings,
        "definition": "a line's level is the first word ERROR/CRITICAL/FATAL/WARNING/WARN/INFO/DEBUG found in it (case-insensitive)",
    }
