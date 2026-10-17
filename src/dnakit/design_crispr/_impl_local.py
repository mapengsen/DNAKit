"""Native helpers behind the design_crispr candidate functions (guides, primers, fragments, HDR geometry)."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_GUIDES = 100_000
MAX_SEQUENCE = 10_000_000
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]",
    "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


# ----------------------------------------------------------------------------- option helpers
def opt_str(options: Mapping[str, object], key: str, default: str, allowed: Sequence[str] | None = None) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"Option {key!r} must be a non-empty string.", code="CRISPR_OPTION_INVALID")
    if allowed is not None and value.strip().lower() not in [a.lower() for a in allowed]:
        raise ConfigurationError(
            f"Option {key!r}={value!r} is not supported.", code="CRISPR_OPTION_INVALID", hint=f"Use one of: {', '.join(allowed)}."
        )
    return value.strip() if allowed is None else next(a for a in allowed if a.lower() == value.strip().lower())


def opt_int(options: Mapping[str, object], key: str, default: int, minimum: int, maximum: int) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"Option {key!r} must be an integer in [{minimum}, {maximum}].", code="CRISPR_OPTION_INVALID"
        )
    return value


def opt_float(options: Mapping[str, object], key: str, default: float, minimum: float, maximum: float) -> float:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"Option {key!r} must be a number in [{minimum}, {maximum}].", code="CRISPR_OPTION_INVALID"
        )
    return float(value)


def opt_bool(options: Mapping[str, object], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(f"Option {key!r} must be a boolean.", code="CRISPR_OPTION_INVALID")
    return value


# ----------------------------------------------------------------------------- sequences
def clean_dna(text: object, name: str = "sequence", allow_iupac: bool = False, max_length: int = MAX_SEQUENCE) -> str:
    if not isinstance(text, str):
        raise InputFormatError(f"{name} must be a DNA string.", code="CRISPR_SEQUENCE_INVALID")
    seq = re.sub(r"\s+", "", text).upper().replace("U", "T")
    allowed = set(_IUPAC) if allow_iupac else set("ACGT")
    if not seq or set(seq) - allowed:
        raise InputFormatError(
            f"{name} is empty or contains characters outside {'IUPAC DNA' if allow_iupac else 'ACGT'}.",
            code="CRISPR_SEQUENCE_INVALID",
        )
    if len(seq) > max_length:
        raise ConfigurationError(f"{name} is longer than {max_length} bases.", code="CRISPR_SEQUENCE_TOO_LONG")
    return seq


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _is_file(text: str) -> bool:
    if "\n" in text or "\t" in text or len(text) >= 4096:
        return False
    try:
        return Path(text).is_file()
    except OSError:  # e.g. file name too long: it is text, not a path
        return False


def read_text(value: object) -> str:
    if isinstance(value, Path) or (isinstance(value, str) and _is_file(value)):
        return Path(value).read_text(encoding="utf-8", errors="replace")
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value
    raise InputFormatError("Expected text, bytes or a path.", code="CRISPR_INPUT_INVALID")


def parse_fasta_or_plain(text: str) -> list[tuple[str, str]]:
    """FASTA records, or one record per non-empty line when no '>' header is present."""

    if text.lstrip().startswith(">"):
        items: list[tuple[str, str]] = []
        name, chunks = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    items.append((name, "".join(chunks)))
                name, chunks = (line[1:].split()[0] if line[1:].strip() else f"seq{len(items) + 1}"), []
            else:
                chunks.append(line.strip())
        if name is not None:
            items.append((name, "".join(chunks)))
        return items
    return [(f"seq{i + 1}", line.strip()) for i, line in enumerate(text.splitlines()) if line.strip()]


def parse_sequences(value: object, name: str = "sequences", max_items: int = 10_000, allow_iupac: bool = False) -> list[tuple[str, str]]:
    """A string/FASTA/path, a list of strings, or a mapping id -> sequence; returns cleaned (id, sequence)."""

    if isinstance(value, Mapping):
        if set(value) <= {"sequence", "id", "name"} and "sequence" in value:
            items = [(str(value.get("id", value.get("name", "seq1"))), str(value["sequence"]))]
        else:
            items = [(str(k), str(v)) for k, v in value.items()]
    elif isinstance(value, (str, bytes, Path)):
        items = parse_fasta_or_plain(read_text(value))
    elif isinstance(value, Sequence):
        items = []
        for i, v in enumerate(value):
            if isinstance(v, Mapping):
                items.append((str(v.get("id", f"seq{i + 1}")), str(v.get("sequence", ""))))
            else:
                items.append((f"seq{i + 1}", str(v)))
    else:
        raise InputFormatError(f"Unsupported {name} input.", code="CRISPR_INPUT_INVALID")
    if not items or len(items) > max_items:
        raise InputFormatError(f"Provide between 1 and {max_items} {name}.", code="CRISPR_INPUT_INVALID")
    return [(n, clean_dna(s, f"{name} {n!r}", allow_iupac)) for n, s in items]


def parse_guides(value: object, max_items: int = MAX_GUIDES, strict: bool = True) -> list[dict[str, object]]:
    """Guides as strings, FASTA, a mapping id -> guide, or dict records (``guide``/``sequence`` plus extra fields)."""

    records: list[dict[str, object]] = []
    if isinstance(value, Mapping):
        if "guide" in value or "sequence" in value:
            value = [value]
        elif "guides" in value:
            return parse_guides(value["guides"], max_items, strict)
        else:
            value = [{"id": k, "guide": v} for k, v in value.items()]
    if isinstance(value, (str, bytes, Path)):
        value = [{"id": n, "guide": s} for n, s in parse_fasta_or_plain(read_text(value))]
    if not isinstance(value, Sequence):
        raise InputFormatError("Unsupported guide input.", code="CRISPR_INPUT_INVALID")
    if not value or len(value) > max_items:
        raise InputFormatError(f"Provide between 1 and {max_items} guides.", code="CRISPR_INPUT_INVALID")
    for i, item in enumerate(value):
        if isinstance(item, Mapping):
            raw = item.get("guide", item.get("sequence"))
            record = {k: v for k, v in item.items() if k not in ("guide", "sequence")}
            record.setdefault("id", f"guide{i + 1}")
        elif isinstance(item, str):
            raw, record = item, {"id": f"guide{i + 1}"}
        else:
            raise InputFormatError("Each guide must be a string or a mapping.", code="CRISPR_INPUT_INVALID")
        if not isinstance(raw, str):
            raise InputFormatError("Guide records need a 'guide' or 'sequence' string.", code="CRISPR_INPUT_INVALID")
        if strict:
            record["sequence"] = clean_dna(raw, f"guide {record['id']!r}")
        else:
            record["sequence"] = re.sub(r"\s+", "", raw).upper().replace("U", "T")
        records.append(record)
    return records


# ----------------------------------------------------------------------------- guide properties
def gc_fraction(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)


def longest_run(seq: str, bases: str = "ACGT") -> int:
    best = run = 0
    previous = ""
    for ch in seq:
        run = run + 1 if ch == previous else 1
        previous = ch
        if ch in bases:
            best = max(best, run)
    return best


def motif_positions(seq: str, motif: str, both_strands: bool = True) -> list[dict[str, object]]:
    """Occurrences of an IUPAC motif; reverse-strand hits are found by scanning the reverse complement."""

    out: list[dict[str, object]] = []
    for strand, pattern in (("+", motif),) + ((("-", revcomp(motif)),) if both_strands and revcomp(motif) != motif else ()):
        regex = re.compile("(?=(" + "".join(_IUPAC[c] for c in pattern) + "))")
        out.extend({"start": m.start(), "strand": strand} for m in regex.finditer(seq))
    return sorted(out, key=lambda o: (o["start"], o["strand"]))  # type: ignore[arg-type,return-value]


def resolve_motifs(options: Mapping[str, object]) -> dict[str, str]:
    """``exclude_motifs`` (IUPAC strings) and ``avoid_enzymes`` (names from the built-in restriction catalog)."""

    motifs: dict[str, str] = {}
    raw = options.get("exclude_motifs", ())
    if isinstance(raw, str) or not isinstance(raw, Sequence):
        raise ConfigurationError("exclude_motifs must be a list of IUPAC strings.", code="CRISPR_OPTION_INVALID")
    for m in raw:
        motifs[str(m).upper()] = clean_dna(str(m), "exclude motif", allow_iupac=True)
    enzymes = options.get("avoid_enzymes", ())
    if isinstance(enzymes, str) or not isinstance(enzymes, Sequence):
        raise ConfigurationError("avoid_enzymes must be a list of enzyme names.", code="CRISPR_OPTION_INVALID")
    if enzymes:
        from dnakit.patterns import BUILTIN_RESTRICTION_ENZYMES

        for name in enzymes:
            if name not in BUILTIN_RESTRICTION_ENZYMES:
                raise ConfigurationError(
                    f"Unknown enzyme {name!r}.", code="CRISPR_OPTION_INVALID",
                    hint=f"Built-in: {', '.join(sorted(BUILTIN_RESTRICTION_ENZYMES))}; or use exclude_motifs.",
                )
            motifs[name] = BUILTIN_RESTRICTION_ENZYMES[name].recognition_sequence
    return motifs


def guide_checks(seq: str, *, min_gc: float, max_gc: float, max_homopolymer: int, avoid_poly_t: bool,
                 motifs: Mapping[str, str]) -> dict[str, object]:
    """Compute GC, homopolymer, poly-T and excluded-motif properties plus the list of failed rules."""

    gc = gc_fraction(seq)
    run = longest_run(seq)
    poly_t = "TTTT" in seq
    hits = [name for name, motif in motifs.items() if motif_positions(seq, motif)]
    failed = []
    if gc < min_gc:
        failed.append("gc_below_min")
    if gc > max_gc:
        failed.append("gc_above_max")
    if run > max_homopolymer:
        failed.append("homopolymer_run")
    if avoid_poly_t and poly_t:
        failed.append("poly_t")
    if hits:
        failed.append("excluded_motif:" + ",".join(hits))
    return {"gc": gc, "longest_homopolymer": run, "poly_t": poly_t, "excluded_motifs": hits, "failed": failed}


def proxy_score(seq: str) -> float:
    """HEURISTIC on-target proxy in [0, 100] (not a trained model; no experimental validation).

    Start at 100 and subtract rule-of-thumb penalties: GC outside 40-70 % costs 15 (outside 20-80 % costs 40
    instead); a TTTT run (RNA Pol III terminator) costs 40; any A/C/G homopolymer of length >= 4 costs 15.
    """

    gc = gc_fraction(seq)
    score = 100.0
    if gc < 0.20 or gc > 0.80:
        score -= 40
    elif gc < 0.40 or gc > 0.70:
        score -= 15
    if "TTTT" in seq:
        score -= 40
    if longest_run(seq, "ACG") >= 4:
        score -= 15
    return max(score, 0.0)


# ----------------------------------------------------------------------------- PAM / guide enumeration
PAM_RULES = {"SpCas9": ("NGG", "3prime", 20), "SaCas9": ("NNGRRT", "3prime", 21), "AsCas12a": ("TTTV", "5prime", 20)}


def resolve_nuclease(options: Mapping[str, object]) -> tuple[str, str, str, int]:
    """Return (name, pam, pam_side, guide_length) from ``nuclease`` or a custom ``pam``/``pam_side``/``guide_length``."""

    if "pam" in options:
        pam = clean_dna(str(options["pam"]), "pam", allow_iupac=True, max_length=20)
        side = opt_str(options, "pam_side", "3prime", ("3prime", "5prime"))
        return str(options.get("nuclease", "custom")), pam, side, opt_int(options, "guide_length", 20, 10, 40)
    name = opt_str(options, "nuclease", "SpCas9", tuple(PAM_RULES))
    pam, side, length = PAM_RULES[name]
    return name, pam, side, opt_int(options, "guide_length", length, 10, 40)


def cut_site(nuclease: str, side: str, strand: str, start: int, end: int) -> int | None:
    """Forward-coordinate cut boundary 3 bp from the PAM for 3'-PAM Cas9-type nucleases (SpCas9, SaCas9); else None."""

    if nuclease not in ("SpCas9", "SaCas9") or side != "3prime":
        return None
    return end - 3 if strand == "forward" else start + 3


def enumerate_guides(sequence: str, options: Mapping[str, object]) -> list[dict[str, object]]:
    """PAM-adjacent guide candidates on both strands (delegates to ``dnakit.patterns.scan_pam_candidates``)."""

    from dnakit.core import DNASequence
    from dnakit.patterns import PAMRule, scan_pam_candidates

    name, pam, side, length = resolve_nuclease(options)
    strand = opt_str(options, "strand", "both", ("both", "forward", "reverse"))
    rule = PAMRule(name, pam, side, length)  # type: ignore[arg-type]
    result = scan_pam_candidates(
        DNASequence(sequence), rule, guide_length=length, strand=strand, max_matches=MAX_GUIDES, verbose=False
    )
    out = []
    for hit in result.hits:
        st = hit.strand.value
        start, end = hit.guide_coordinate_location.start, hit.guide_coordinate_location.end
        out.append({
            "sequence": hit.guide_sequence, "pam": hit.pam_sequence, "strand": st, "start": start, "end": end,
            "pam_start": hit.pam_coordinate_location.start, "cut_site": cut_site(name, side, st, start, end),
        })
    out.sort(key=lambda g: (g["start"], g["strand"]))  # type: ignore[arg-type,return-value]
    return out


# ----------------------------------------------------------------------------- off-target search
def hamming_scan(guide: str, reference: str, max_mismatches: int, pam: str | None, pam_side: str) -> list[dict[str, object]]:
    """Exact ungapped Hamming search of ``guide`` on both strands of ``reference``.

    With ``pam`` (IUPAC) a site only counts when the PAM adjacent to the protospacer matches it. Coordinates are
    0-based half-open on the forward strand of ``reference``; mismatch positions are 1-based from the guide 5' end.
    """

    length = len(guide)
    plen = len(pam) if pam else 0
    pam_re = re.compile("".join(_IUPAC[c] for c in pam)) if pam else None
    hits: list[dict[str, object]] = []
    n = len(reference)
    for strand, text in (("+", reference), ("-", revcomp(reference))):
        for i in range(0, n - length + 1):
            window = text[i : i + length]
            mismatches = []
            for j in range(length):
                if window[j] != guide[j]:
                    mismatches.append(j + 1)
                    if len(mismatches) > max_mismatches:
                        break
            if len(mismatches) > max_mismatches:
                continue
            if pam_re is not None:
                if pam_side == "3prime":
                    pam_text = text[i + length : i + length + plen]
                else:
                    pam_text = text[i - plen : i] if i - plen >= 0 else ""
                if len(pam_text) != plen or not pam_re.fullmatch(pam_text):
                    continue
            start = i if strand == "+" else n - (i + length)
            hits.append({"strand": strand, "start": start, "end": start + length, "mismatches": len(mismatches),
                         "mismatch_positions": mismatches, "site": window})
    hits.sort(key=lambda h: (h["mismatches"], h["start"], h["strand"]))  # type: ignore[arg-type,return-value]
    return hits


# ----------------------------------------------------------------------------- primers
def longest_common_substring(a: str, b: str) -> tuple[int, int, int]:
    """Return (length, start_in_a, start_in_b) of the longest common substring (O(len(a)*len(b)))."""

    best = (0, 0, 0)
    previous = [0] * (len(b) + 1)
    for i in range(1, len(a) + 1):
        current = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                current[j] = previous[j - 1] + 1
                if current[j] > best[0]:
                    best = (current[j], i - current[j], j - current[j])
        previous = current
    return best


def three_prime_run(a: str, b: str) -> int:
    """Longest L such that the last L bases of ``a`` are fully complementary (antiparallel) to a stretch of ``b``."""

    rc_b = revcomp(b)
    for length in range(min(len(a), len(b)), 0, -1):
        if a[-length:] in rc_b:
            return length
    return 0


def find_hairpins(seq: str, min_stem: int, min_loop: int) -> list[dict[str, object]]:
    """Maximal perfect Watson-Crick stems with a loop of at least ``min_loop`` unpaired bases.

    A stem is a run of nested pairs (seq[i+k], seq[j-k]) with A-T / C-G complementarity. Only maximal stems are
    reported: the outer pair cannot be extended outwards, and the stem is extended inwards while the loop stays
    >= ``min_loop``. No mismatches, bulges or G-U pairs; no energy model.
    """

    n = len(seq)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}

    def pairs(x: int, y: int) -> bool:
        return comp.get(seq[x]) == seq[y]

    found: list[dict[str, object]] = []
    for i in range(n):
        for j in range(i + 1, n):
            if not pairs(i, j) or (i > 0 and j < n - 1 and pairs(i - 1, j + 1)):
                continue
            length = 1
            while pairs(i + length, j - length) and (j - i + 1 - 2 * (length + 1)) >= min_loop:
                length += 1
            loop = j - i + 1 - 2 * length
            if length >= min_stem and loop >= min_loop:
                right_start = j - length + 1
                found.append({
                    "stem_length": length, "left_start": i, "right_start": right_start, "loop_length": loop,
                    "stem_5prime": seq[i : i + length], "stem_3prime": seq[right_start : j + 1],
                    "involves_3prime_end": j == n - 1,
                })
    found.sort(key=lambda h: (-h["stem_length"], h["left_start"], h["right_start"]))  # type: ignore[operator]
    return found


# ----------------------------------------------------------------------------- SVG
def svg_guide_map(length: int, guides: Sequence[Mapping[str, object]]) -> str:
    width, pad = 800, 20
    scale = (width - 2 * pad) / max(length, 1)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="90">',
             f'<line x1="{pad}" y1="45" x2="{width - pad}" y2="45" stroke="#333" stroke-width="2"/>',
             f'<text x="{pad}" y="85" font-size="10">0</text><text x="{width - pad - 30}" y="85" font-size="10">{length}</text>']
    for g in guides:
        x = pad + float(g["start"]) * scale  # type: ignore[arg-type]
        w = max((float(g["end"]) - float(g["start"])) * scale, 1.0)  # type: ignore[arg-type]
        up = g.get("strand") in ("forward", "+")
        y = 25 if up else 50
        color = "#4e79a7" if up else "#e15759"
        parts.append(
            f'<rect class="guide" data-start="{g["start"]}" data-end="{g["end"]}" data-strand="{g.get("strand")}" '
            f'x="{x:.2f}" y="{y}" width="{w:.2f}" height="15" fill="{color}"/>'
        )
    parts.append("</svg>")
    return "".join(parts)


# ----------------------------------------------------------------------------- off-target drivers
def scan_offtargets(guides: Sequence[Mapping[str, object]], refs: Sequence[tuple[str, str]], max_mismatches: int,
                    pam: str | None, side: str) -> list[dict[str, object]]:
    total = sum(len(s) for _, s in refs)
    if len(guides) * total * 2 > 50_000_000:
        raise ConfigurationError(
            "Off-target search is too large.", code="CRISPR_SEARCH_TOO_LARGE",
            hint="Use fewer guides or shorter references (limit: guides x reference bases x 2 <= 50 million).",
        )
    out = []
    for g in guides:
        hits: list[dict[str, object]] = []
        for rid, ref in refs:
            hits.extend({"reference": rid, **h} for h in hamming_scan(str(g["sequence"]), ref, max_mismatches, pam, side))
        counts: dict[int, int] = {}
        for h in hits:
            counts[int(h["mismatches"])] = counts.get(int(h["mismatches"]), 0) + 1  # type: ignore[call-overload]
        out.append({"id": g["id"], "sequence": g["sequence"], "hit_count": len(hits),
                    "counts_by_mismatch": {str(k): counts[k] for k in sorted(counts)}, "hits": hits})
    return out


def add_offtargets(guides: Sequence[dict[str, object]], refs: Sequence[tuple[str, str]], options: Mapping[str, object]) -> None:
    """Annotate guide dicts with exact_sites / near_sites / off_target_sites (see ``design``)."""

    if not guides:
        return
    _, pam, side, _ = resolve_nuclease(options)
    max_mm = opt_int(options, "max_mismatches", 3, 0, 6)
    for g, scanned in zip(guides, scan_offtargets(guides, refs, max_mm, pam, side)):
        exact = int(scanned["counts_by_mismatch"].get("0", 0))  # type: ignore[attr-defined]
        near = int(scanned["hit_count"]) - exact  # type: ignore[call-overload]
        g.update({"exact_sites": exact, "near_sites": near, "off_target_sites": near + max(exact - 1, 0)})


# ----------------------------------------------------------------------------- primers (Tm, pair evaluation)
MAX_PRIMER = 200
WALLACE_MAX_AUTO = 13
NN_MAX = 60


def clean_primer(text: object, name: str, rna: bool = False) -> str:
    """Upper-case an unambiguous primer; ``U`` is accepted (as ``T``) only when ``rna`` is true."""

    if not isinstance(text, str):
        raise InputFormatError(f"{name} must be a DNA string.", code="CRISPR_SEQUENCE_INVALID")
    seq = re.sub(r"\s+", "", text).upper()
    if rna:
        seq = seq.replace("U", "T")
    bad = sorted(set(seq) - set("ACGT"))
    if not seq or bad:
        what = "is empty" if not seq else f"contains {''.join(bad)!r}"
        hint = "Only A, C, G and T are accepted: Tm is undefined for ambiguous (IUPAC/N) bases."
        if "U" in bad:
            hint = "Set options['rna'] = true to read U as T."
        raise InputFormatError(f"{name} {what}.", code="CRISPR_SEQUENCE_INVALID", hint=hint)
    if len(seq) > MAX_PRIMER:
        raise ConfigurationError(f"{name} is longer than {MAX_PRIMER} bases.", code="CRISPR_SEQUENCE_TOO_LONG")
    return seq


def parse_primers(value: object, rna: bool = False, max_items: int = 1_000) -> list[tuple[str, str]]:
    """A primer string, ``{"id", "sequence"}``, ``{name: sequence}`` or a list of strings/``{"id", "sequence"}``."""

    if isinstance(value, str):
        raw = [("primer1", value)]
    elif isinstance(value, Mapping):
        if "sequence" in value and set(value) <= {"sequence", "id", "name"}:
            raw = [(str(value.get("id", value.get("name", "primer1"))), value["sequence"])]
        else:
            raw = [(str(k), v) for k, v in value.items()]
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        raw = []
        for i, item in enumerate(value):
            if isinstance(item, Mapping):
                raw.append((str(item.get("id", item.get("name", f"primer{i + 1}"))), item.get("sequence")))
            else:
                raw.append((f"primer{i + 1}", item))
    else:
        raise InputFormatError("Provide a primer string, a list or a mapping.", code="CRISPR_INPUT_INVALID")
    if not raw or len(raw) > max_items:
        raise InputFormatError(f"Provide between 1 and {max_items} primers.", code="CRISPR_INPUT_INVALID")
    return [(n, clean_primer(s, f"primer {n!r}", rna)) for n, s in raw]


def tm_conditions(options: Mapping[str, object]) -> dict[str, float]:
    conditions = {
        "sodium_molar": opt_float(options, "sodium_molar", 0.05, 0.0, 2.0),
        "potassium_molar": opt_float(options, "potassium_molar", 0.0, 0.0, 2.0),
        "strand_concentration_molar": opt_float(options, "strand_concentration_molar", 250e-9, 1e-12, 1.0),
    }
    if not 0.0 < conditions["sodium_molar"] + conditions["potassium_molar"] <= 2.0:
        raise ConfigurationError("sodium_molar + potassium_molar must be in (0, 2] mol/L.", code="CRISPR_OPTION_INVALID")
    return conditions


def primer_tm(seq: str, method: str, conditions: Mapping[str, float]) -> tuple[float, str]:
    """Tm in Celsius and the method actually used (``wallace`` or ``nearest_neighbor``)."""

    used = method if method != "auto" else ("wallace" if len(seq) <= WALLACE_MAX_AUTO else "nearest_neighbor")
    if used == "wallace":
        return float(2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))), used
    if not 2 <= len(seq) <= NN_MAX:
        raise ConfigurationError(
            f"The nearest-neighbour Tm is defined for 2-{NN_MAX} nt; the primer has {len(seq)}.",
            code="CRISPR_PRIMER_LENGTH", hint="Use options['method'] = 'wallace' or a shorter primer.",
        )
    from dnakit.core import DNASequence
    from dnakit.thermodynamics import ThermodynamicConditions, melting_temperature

    result = melting_temperature(
        DNASequence(seq), method="nearest_neighbor", conditions=ThermodynamicConditions(**conditions), verbose=False
    )
    return round(float(result.tm_celsius), 2), used


def self_complementary_run(seq: str) -> int:
    """Longest contiguous stretch of ``seq`` whose reverse complement also occurs in ``seq`` (self-dimer type)."""

    return longest_common_substring(seq, revcomp(seq))[0]


def gc_clamp(seq: str, span: int = 5) -> int:
    tail = seq[-span:]
    return tail.count("G") + tail.count("C")


# ----------------------------------------------------------------------------- DNA data storage (toy 2 bits/base code)
_BASE_OF_BITS = "ACGT"
MAX_STORAGE_BYTES = 1_000_000


def bytes_to_dna(data: bytes) -> str:
    """Four bases per byte, most significant bit pair first: A=00, C=01, G=10, T=11."""

    table = ["".join(_BASE_OF_BITS[(b >> s) & 3] for s in (6, 4, 2, 0)) for b in range(256)]
    return "".join(table[b] for b in data)


def dna_to_bytes(dna: str) -> bytes:
    """Inverse of :func:`bytes_to_dna`; the length must be a multiple of 4 and the bases ACGT."""

    seq = re.sub(r"\s+", "", dna).upper()
    if len(seq) % 4 or set(seq) - set(_BASE_OF_BITS):
        raise InputFormatError("DNA to decode must be ACGT with a length that is a multiple of 4.", code="CRISPR_SEQUENCE_INVALID")
    code = {b: i for i, b in enumerate(_BASE_OF_BITS)}
    return bytes(
        (code[seq[i]] << 6) | (code[seq[i + 1]] << 4) | (code[seq[i + 2]] << 2) | code[seq[i + 3]] for i in range(0, len(seq), 4)
    )
