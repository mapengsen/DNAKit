"""Native replication-strand-asymmetry and DNA-repair-signature computations.

All functions are deterministic, bounded and standard-library only. Coordinates are 0-based,
half-open. GC-skew based origin/terminus positions are heuristics (see the public docstrings).
"""

from __future__ import annotations

import html
import itertools
import math
import random
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._candidate_runtime import candidate_ok
from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_SEQUENCE = 20_000_000
_SEQ_LINE = re.compile(r"^[A-Za-z*.\-]*$")
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def run_op(
    operation: str, fn: Any, value: object, options: Mapping[str, object] | None
) -> dict[str, object]:
    if options is not None and not isinstance(options, Mapping):
        raise ConfigurationError(
            "options must be a mapping or None.",
            code="INVALID_OPTION",
            hint="Pass options={'key': value}.",
        )
    opts = dict(options or {})
    warnings: list[str] = []
    result = fn(value, opts, warnings)
    return candidate_ok(
        "replication_repair", operation, value, result, options=opts, warnings=warnings
    )


# --------------------------------------------------------------------------- option / input helpers


def opt_int(
    options: Mapping[str, object], key: str, default: int, *, minimum: int, maximum: int
) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(
            f"option {key!r} must be an integer in [{minimum}, {maximum}]; got {value!r}.",
            code="INVALID_OPTION",
        )
    return value


def opt_num(
    options: Mapping[str, object], key: str, default: float, *, minimum: float, maximum: float
) -> float:
    value = options.get(key, default)
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or not minimum <= value <= maximum
    ):
        raise ConfigurationError(
            f"option {key!r} must be a number in [{minimum}, {maximum}]; got {value!r}.",
            code="INVALID_OPTION",
        )
    return float(value)


def opt_bool(options: Mapping[str, object], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(
            f"option {key!r} must be a boolean; got {value!r}.", code="INVALID_OPTION"
        )
    return value


def opt_choice(
    options: Mapping[str, object], key: str, default: str, allowed: Sequence[str]
) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(
            f"option {key!r} must be one of {', '.join(allowed)}; got {value!r}.",
            code="INVALID_OPTION",
        )
    return value


def dna_text(value: object, name: str = "sequence") -> str:
    """Upper-case DNA from a plain string, FASTA text, a file path or a ``{'sequence': ...}`` mapping."""

    if isinstance(value, Mapping):
        for key in ("sequence", "seq"):
            if key in value:
                return dna_text(value[key], name)
        raise InputFormatError(f"{name}: mapping has no 'sequence' key.", code="INVALID_INPUT")
    if isinstance(value, Path):
        value = value.read_text(encoding="utf-8")
    if not isinstance(value, str):
        raise InputFormatError(
            f"{name}: expected a DNA string, FASTA text, a path or {{'sequence': ...}}; got {type(value).__name__}.",
            code="INVALID_INPUT",
        )
    text = value
    if "\n" not in text and 0 < len(text) < 4096 and not text.startswith(">"):
        try:
            candidate = Path(text)
            if candidate.is_file():
                text = candidate.read_text(encoding="utf-8")
        except OSError:
            pass
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if lines and lines[0].startswith(">"):
        body = []
        for line in lines[1:]:
            if line.startswith(">"):
                raise InputFormatError(
                    f"{name}: multi-record FASTA is not supported here.",
                    code="INVALID_INPUT",
                    hint="Pass one record at a time.",
                )
            body.append(line)
        lines = body
    seq = "".join(lines).upper().replace("U", "T")
    if not seq:
        raise InputFormatError(
            f"{name} is empty.", code="EMPTY_INPUT", hint="Provide a DNA sequence."
        )
    if not re.fullmatch(r"[ACGTNRYKMSWBDHV]+", seq):
        raise InputFormatError(
            f"{name} contains characters that are not IUPAC DNA letters.", code="INVALID_SEQUENCE"
        )
    if len(seq) > MAX_SEQUENCE:
        raise InputFormatError(
            f"{name} is longer than {MAX_SEQUENCE} bases.", code="INPUT_TOO_LARGE"
        )
    return seq


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _round(x: float | None, digits: int = 6) -> float | None:
    return None if x is None else round(x, digits)


def write_svg(options: Mapping[str, object], svg: str) -> str | None:
    target = options.get("output_path")
    if target is None:
        return None
    if not isinstance(target, (str, Path)) or Path(target).suffix.lower() != ".svg":
        raise ConfigurationError(
            "output_path must be a path ending in .svg.", code="INVALID_OPTION"
        )
    path = Path(target)
    if not path.parent.is_dir():
        raise ConfigurationError(f"directory {path.parent} does not exist.", code="INVALID_OPTION")
    if path.exists() and not opt_bool(options, "overwrite", False):
        raise ConfigurationError(
            f"{path} already exists.",
            code="OUTPUT_EXISTS",
            hint="Pass overwrite=True to replace it.",
        )
    path.write_text(svg, encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------- GC skew


def skew_windows(seq: str, window: int, step: int, kind: str) -> list[dict[str, Any]]:
    """Per-window skew: GC skew (G-C)/(G+C) or AT skew (A-T)/(A+T) over full windows only."""

    a, b = ("G", "C") if kind == "GC" else ("A", "T")
    n = len(seq)
    # prefix counts make each window O(1)
    pa = [0] * (n + 1)
    pb = [0] * (n + 1)
    ca = cb = 0
    for i, ch in enumerate(seq):
        if ch == a:
            ca += 1
        elif ch == b:
            cb += 1
        pa[i + 1], pb[i + 1] = ca, cb
    rows = []
    cumulative = 0.0
    for start in range(0, n - window + 1, step):
        end = start + window
        x, y = pa[end] - pa[start], pb[end] - pb[start]
        skew = (x - y) / (x + y) if x + y else 0.0
        cumulative += skew
        rows.append(
            {
                "start": start,
                "end": end,
                "skew": skew,
                "cumulative_skew": cumulative,
                "count_" + a: x,
                "count_" + b: y,
            }
        )
    return rows


def _skew_params(options: Mapping[str, object], n: int) -> tuple[int, int, str]:
    window = opt_int(options, "window", min(1000, n), minimum=2, maximum=10_000_000)
    step = opt_int(options, "step", window, minimum=1, maximum=10_000_000)
    kind = opt_choice(options, "skew_type", "GC", ("GC", "AT"))
    if window > n:
        raise ConfigurationError(
            f"window={window} is longer than the sequence ({n} bases).",
            code="INVALID_OPTION",
            hint="Use a smaller window.",
        )
    return window, step, kind


def _ori_ter(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    low = min(rows, key=lambda r: (r["cumulative_skew"], r["start"]))
    high = max(rows, key=lambda r: (r["cumulative_skew"], -r["start"]))
    return {
        "origin_estimate": (low["start"] + low["end"]) // 2,
        "origin_window": [low["start"], low["end"]],
        "terminus_estimate": (high["start"] + high["end"]) // 2,
        "terminus_window": [high["start"], high["end"]],
        "cumulative_skew_range": _round(high["cumulative_skew"] - low["cumulative_skew"]),
    }


def op_calculate(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    seq = dna_text(value)
    window, step, kind = _skew_params(options, len(seq))
    rows = skew_windows(seq, window, step, kind)
    a, b = ("G", "C") if kind == "GC" else ("A", "T")
    total_a, total_b = seq.count(a), seq.count(b)
    result: dict[str, Any] = {
        "skew_type": kind,
        "length": len(seq),
        "window": window,
        "step": step,
        "windows": len(rows),
        "overall_skew": _round((total_a - total_b) / (total_a + total_b))
        if total_a + total_b
        else None,
        **_ori_ter(rows),
        "heuristic": "Origin = centre of the window with the minimum cumulative skew, terminus = maximum (Lobry 1996; Grigoriev 1998). "
        "Valid only for circular bacterial-like chromosomes with a leading/lagging-strand composition bias; the positions are "
        "estimates with window-size resolution, not experimental origins.",
        "definition": f"skew = ({a}-{b})/({a}+{b}) per window; cumulative skew = running sum of the window skews.",
    }
    if len(rows) < 4:
        warnings.append("fewer than 4 windows: the origin/terminus estimate is not meaningful.")
    limit = opt_int(options, "max_windows_returned", 20_000, minimum=1, maximum=1_000_000)
    if opt_bool(options, "include_windows", True):
        result["window_table"] = [
            {k: _round(v) if isinstance(v, float) else v for k, v in r.items()}
            for r in rows[:limit]
        ]
        if len(rows) > limit:
            warnings.append(
                f"window_table truncated to {limit} of {len(rows)} windows (max_windows_returned)."
            )
    return result


# --------------------------------------------------------------------------- DnaA boxes


def dnaa_boxes(seq: str, max_mismatches: int) -> list[dict[str, Any]]:
    """Occurrences of the DnaA-box consensus TTATNCACA on both strands (Hamming <= max_mismatches)."""

    consensus = "TTATNCACA"
    out = []
    n = len(seq)
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for i in range(n - 8):
            mism = 0
            for x, c in zip(s[i : i + 9], consensus, strict=False):
                if c != "N" and x != c:
                    mism += 1
                    if mism > max_mismatches:
                        break
            if mism <= max_mismatches:
                start = i if strand == "+" else n - i - 9
                out.append(
                    {
                        "start": start,
                        "end": start + 9,
                        "strand": strand,
                        "mismatches": mism,
                        "sequence": s[i : i + 9],
                    }
                )
    return sorted(out, key=lambda r: (r["start"], r["strand"]))


def op_predict(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    seq = dna_text(value)
    mm = opt_int(options, "max_mismatches", 0, minimum=0, maximum=2)
    window = opt_int(options, "window", 500, minimum=20, maximum=100_000)
    top_n = opt_int(options, "top_n", 3, minimum=1, maximum=50)
    boxes = dnaa_boxes(seq, mm)
    n = len(seq)
    candidates = []
    if boxes:
        starts = [b["start"] for b in boxes]
        best_for: dict[int, int] = {}
        for anchor in starts:
            lo, hi = anchor, anchor + window
            count = sum(1 for s in starts if lo <= s and s + 9 <= hi)
            best_for[anchor] = count
        used: list[tuple[int, int]] = []
        for anchor, count in sorted(best_for.items(), key=lambda kv: (-kv[1], kv[0])):
            if any(anchor < u1 and u0 < anchor + window for u0, u1 in used):
                continue
            used.append((anchor, anchor + window))
            members = [b for b in boxes if anchor <= b["start"] and b["end"] <= anchor + window]
            candidates.append(
                {
                    "start": anchor,
                    "end": min(n, anchor + window),
                    "dnaa_boxes": count,
                    "box_positions": [(b["start"], b["strand"]) for b in members],
                }
            )
            if len(candidates) >= top_n:
                break
    result: dict[str, Any] = {
        "length": n,
        "dnaa_box_consensus": "TTATNCACA",
        "max_mismatches": mm,
        "boxes_found": len(boxes),
        "candidate_regions": candidates,
        "heuristic": "oriC is predicted as the non-overlapping windows (default 500 bp) with the most DnaA boxes "
        "(TTATNCACA on either strand, Hamming <= max_mismatches). DnaA-box density alone is a weak, "
        "bacterial-specific predictor; compare with the GC-skew minimum below and with experimental data.",
    }
    if n >= 16:
        window_s = max(2, min(1000, n // 4))
        rows = skew_windows(seq, window_s, window_s, "GC")
        ori = _ori_ter(rows)
        result["gc_skew_origin_estimate"] = ori["origin_estimate"]
        if candidates:
            for c in candidates:
                centre = (c["start"] + c["end"]) // 2
                d = abs(centre - ori["origin_estimate"])
                c["circular_distance_to_gc_skew_origin"] = min(d, n - d)
    if not candidates:
        warnings.append("no DnaA box found; no oriC candidate can be proposed.")
    return result


# --------------------------------------------------------------------------- repair motifs


def _primitive(unit: str) -> bool:
    k = len(unit)
    return not any(k % d == 0 and unit == unit[:d] * (k // d) for d in range(1, k))


_DEFAULT_MIN_COPIES = {1: 8, 2: 5, 3: 4, 4: 3}


def op_detect(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    seq = dna_text(value)
    max_results = opt_int(options, "max_results", 1000, minimum=1, maximum=100_000)
    raw = options.get("min_copies")
    min_copies = dict(_DEFAULT_MIN_COPIES)
    if raw is not None:
        if not isinstance(raw, Mapping):
            raise ConfigurationError(
                "min_copies must be a mapping unit_length -> copies.", code="INVALID_OPTION"
            )
        for k, v in raw.items():
            if isinstance(v, bool) or not isinstance(v, int) or v < 2 or int(k) < 1 or int(k) > 6:
                raise ConfigurationError(
                    "min_copies entries need unit length 1-6 and at least 2 copies.",
                    code="INVALID_OPTION",
                )
            min_copies[int(k)] = v
    gatc = [m.start() for m in re.finditer(r"(?=GATC)", seq)]
    tracts: list[dict[str, Any]] = []
    for k, m in sorted(min_copies.items()):
        found = []
        for match in re.finditer(rf"(?=(([ACGT]{{{k}}})\2{{{m - 1},}}))", seq):
            unit = match.group(2)
            if not _primitive(unit):
                continue
            text = match.group(1)
            start = match.start()
            found.append(
                {
                    "start": start,
                    "end": start + len(text),
                    "unit": unit,
                    "copies": len(text) // k,
                    "length": len(text),
                }
            )
        # the same repeat is seen in every reading frame: keep the longest tract of each overlapping group
        kept: list[dict[str, Any]] = []
        for tract in sorted(found, key=lambda t: (-t["length"], t["start"])):
            if all(tract["end"] <= o["start"] or tract["start"] >= o["end"] for o in kept):
                kept.append(tract)
        tracts.extend(kept)
    tracts.sort(key=lambda t: (t["start"], t["unit"]))
    n = len(seq)
    return {
        "length": n,
        "gatc_sites": {
            "count": len(gatc),
            "density_per_kb": _round(len(gatc) / n * 1000),
            "expected_random_per_kb": _round(1000 / 256, 4),
            "positions": gatc[:max_results],
            "note": "GATC is the Dam-methylation site used by E. coli MutHLS to tell the parental from the new strand; it is "
            "palindromic, so each site is one position on both strands.",
        },
        "microsatellite_tracts": {
            "count": len(tracts),
            "min_copies": {str(k): v for k, v in sorted(min_copies.items())},
            "tracts": tracts[:max_results],
            "note": "Perfect repeats of 1-4 bp units are the loci that mismatch-repair-deficient cells expand or contract "
            "(microsatellite instability); presence does not imply instability.",
        },
        "heuristic": "Exact string scan for the GATC motif and perfect mono- to tetranucleotide tandem repeats at or above the "
        "minimum copy numbers; no statistical model of instability is applied.",
    }


# --------------------------------------------------------------------------- annotate (strand asymmetry)


def _features(value: object) -> list[dict[str, Any]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise InputFormatError(
            "'features' must be a non-empty list of {id, start, end, strand} mappings.",
            code="INVALID_INPUT",
        )
    out = []
    for i, f in enumerate(value):
        if not isinstance(f, Mapping) or "start" not in f or "end" not in f or "strand" not in f:
            raise InputFormatError(
                f"feature {i + 1} needs start, end and strand.", code="INVALID_INPUT"
            )
        start, end, strand = f["start"], f["end"], f["strand"]
        if (
            isinstance(start, bool)
            or isinstance(end, bool)
            or not isinstance(start, int)
            or not isinstance(end, int)
            or start < 0
            or end <= start
        ):
            raise InputFormatError(
                f"feature {i + 1}: need integers 0 <= start < end.", code="INVALID_INPUT"
            )
        if strand not in ("+", "-", 1, -1):
            raise InputFormatError(
                f"feature {i + 1}: strand must be '+' or '-'.", code="INVALID_INPUT"
            )
        out.append(
            {
                "id": str(f.get("id", f.get("name", f"feature{i + 1}"))),
                "start": start,
                "end": end,
                "strand": "+" if strand in ("+", 1) else "-",
            }
        )
    return out


def op_annotate(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    if not isinstance(value, Mapping) or "features" not in value:
        raise InputFormatError(
            "input must be a mapping with 'features' and either 'length', 'sequence' or both origin/terminus.",
            code="INVALID_INPUT",
            hint="Example: {'length': 1000000, 'origin': 0, 'terminus': 500000, 'features': [{'id': 'g1', 'start': 10, 'end': 900, 'strand': '+'}]}",
        )
    features = _features(value["features"])
    seq = dna_text(value["sequence"]) if "sequence" in value else None
    length = value.get("length", len(seq) if seq else None)
    if isinstance(length, bool) or not isinstance(length, int) or length < 2:
        raise InputFormatError(
            "provide 'length' (integer >= 2) or 'sequence'.", code="INVALID_INPUT"
        )
    origin = options.get("origin", value.get("origin"))
    terminus = options.get("terminus", value.get("terminus"))
    source = "supplied"
    if origin is None:
        if seq is None:
            raise InputFormatError(
                "origin is not given and no sequence is available to estimate it.",
                code="INVALID_INPUT",
                hint="Provide 'origin' or 'sequence'.",
            )
        window = max(2, min(1000, length // 4))
        est = _ori_ter(skew_windows(seq, window, window, "GC"))
        origin, source = est["origin_estimate"], "cumulative GC-skew estimate (heuristic)"
        if terminus is None:
            terminus = est["terminus_estimate"]
    if terminus is None:
        terminus = (origin + length // 2) % length
    for name, pos in (("origin", origin), ("terminus", terminus)):
        if isinstance(pos, bool) or not isinstance(pos, int) or not 0 <= pos < length:
            raise InputFormatError(
                f"{name} must be an integer in [0, {length}).", code="INVALID_INPUT"
            )
    if origin == terminus:
        raise InputFormatError("origin and terminus must differ.", code="INVALID_INPUT")
    right_span = (
        terminus - origin
    ) % length  # clockwise (increasing coordinate) origin -> terminus
    rows = []
    for f in features:
        if f["end"] > length:
            raise InputFormatError(
                f"feature {f['id']} ends beyond the chromosome length {length}.",
                code="INVALID_INPUT",
            )
        mid = (f["start"] + f["end"]) // 2
        offset = (mid - origin) % length
        right = offset < right_span
        replichore = "right" if right else "left"
        leading = (f["strand"] == "+") == right
        distance = offset if right else (length - offset) % length
        rows.append(
            {
                "id": f["id"],
                "start": f["start"],
                "end": f["end"],
                "strand": f["strand"],
                "replichore": replichore,
                "orientation": "co-directional (leading strand)"
                if leading
                else "head-on (lagging strand)",
                "leading_strand": leading,
                "distance_from_origin_along_replichore": distance,
            }
        )
    n_lead = sum(1 for r in rows if r["leading_strand"])
    return {
        "length": length,
        "origin": origin,
        "terminus": terminus,
        "origin_source": source,
        "right_replichore_length": right_span,
        "left_replichore_length": length - right_span,
        "features": rows,
        "leading_strand_features": n_lead,
        "lagging_strand_features": len(rows) - n_lead,
        "leading_strand_fraction": _round(n_lead / len(rows)),
        "definition": "Circular chromosome; the right replichore runs origin -> terminus in increasing coordinates, the left replichore "
        "terminus -> origin. A feature (position = midpoint) is on the leading strand, i.e. transcribed co-directionally with the "
        "replication fork, when it is on the + strand in the right replichore or on the - strand in the left replichore.",
    }


# --------------------------------------------------------------------------- compare


def _pearson(x: Sequence[float], y: Sequence[float]) -> float | None:
    n = len(x)
    mx, my = math.fsum(x) / n, math.fsum(y) / n
    sxx = math.fsum((a - mx) ** 2 for a in x)
    syy = math.fsum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return None
    return math.fsum((a - mx) * (b - my) for a, b in zip(x, y, strict=False)) / math.sqrt(sxx * syy)


def _resample(values: Sequence[float], points: int) -> list[float]:
    n = len(values)
    out = []
    for i in range(points):
        pos = i * (n - 1) / (points - 1) if points > 1 else 0.0
        lo = math.floor(pos)
        hi = min(lo + 1, n - 1)
        out.append(values[lo] + (values[hi] - values[lo]) * (pos - lo))
    return out


def op_compare(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    if isinstance(value, Mapping) and "a" in value and "b" in value:
        pair = [value["a"], value["b"]]
        names = ["a", "b"]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == 2:
        pair = list(value)
        names = ["first", "second"]
    else:
        raise InputFormatError(
            "input must be {'a': seq, 'b': seq} or a list of exactly two sequences.",
            code="INVALID_INPUT",
        )
    kind = opt_choice(options, "skew_type", "GC", ("GC", "AT"))
    points = opt_int(options, "points", 100, minimum=10, maximum=10_000)
    profiles = []
    summary = {}
    for name, item in zip(names, pair, strict=False):
        seq = dna_text(item, name)
        window = opt_int(
            options, "window", max(2, min(1000, len(seq) // 20)), minimum=2, maximum=10_000_000
        )
        if window > len(seq):
            raise ConfigurationError(
                f"window={window} exceeds the length of sequence {name!r}.", code="INVALID_OPTION"
            )
        rows = skew_windows(seq, window, window, kind)
        if len(rows) < 2:
            raise InputFormatError(
                f"sequence {name!r} yields fewer than 2 windows.", code="INPUT_TOO_SHORT"
            )
        est = _ori_ter(rows)
        profiles.append(_resample([r["cumulative_skew"] for r in rows], points))
        summary[name] = {
            "length": len(seq),
            "windows": len(rows),
            "origin_estimate": est["origin_estimate"],
            "origin_relative_position": _round(est["origin_estimate"] / len(seq)),
            "terminus_estimate": est["terminus_estimate"],
            "terminus_relative_position": _round(est["terminus_estimate"] / len(seq)),
        }
    r = _pearson(profiles[0], profiles[1])
    ori_diff = abs(
        summary[names[0]]["origin_relative_position"]
        - summary[names[1]]["origin_relative_position"]
    )
    return {
        "skew_type": kind,
        "sequences": summary,
        "cumulative_skew_correlation": _round(r),
        "origin_relative_position_difference": _round(min(ori_diff, 1 - ori_diff)),
        "definition": f"Cumulative {kind} skew of each sequence is resampled to {points} points over relative position and compared with Pearson "
        "correlation (+1 = identical shape, -1 = inverted skew such as a complementary or oppositely oriented chromosome). "
        "Origin positions are the heuristic cumulative-skew minima.",
    }


# --------------------------------------------------------------------------- replication timing from coverage


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    mid = n // 2
    return ordered[mid] if n % 2 else (ordered[mid - 1] + ordered[mid]) / 2


def _depths(value: object, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise InputFormatError(
            f"'{name}' must be a non-empty list of window depths.", code="INVALID_INPUT"
        )
    out = []
    for v in value:
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0:
            raise InputFormatError(
                f"'{name}' contains a negative or non-numeric depth.", code="INVALID_INPUT"
            )
        out.append(float(v))
    if len(out) > 5_000_000:
        raise InputFormatError(f"'{name}' has more than 5,000,000 windows.", code="INPUT_TOO_LARGE")
    return out


def op_map(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    if not isinstance(value, Mapping) or "s_phase" not in value or "g1" not in value:
        raise InputFormatError(
            "input must be {'s_phase': [depth per window], 'g1': [depth per window]}.",
            code="INVALID_INPUT",
            hint="Windows must be the same fixed-size genomic windows in the same order for both libraries.",
        )
    s = _depths(value["s_phase"], "s_phase")
    g = _depths(value["g1"], "g1")
    if len(s) != len(g):
        raise InputFormatError(
            f"s_phase has {len(s)} windows but g1 has {len(g)}.", code="WINDOW_MISMATCH"
        )
    pseudo = opt_num(options, "pseudocount", 0.0, minimum=0.0, maximum=1e6)
    smooth = opt_int(options, "smoothing_windows", 1, minimum=1, maximum=10_001)
    if smooth % 2 == 0:
        raise ConfigurationError("smoothing_windows must be odd.", code="INVALID_OPTION")
    window_size = opt_int(options, "window_size", 1, minimum=1, maximum=10_000_000_000)
    ms, mg = _median(s), _median(g)
    if ms <= 0 or mg <= 0:
        raise InputFormatError(
            "the median depth of s_phase and g1 must be positive for median normalisation.",
            code="INVALID_INPUT",
        )
    ratios = []
    for i, (a, b) in enumerate(zip(s, g, strict=False)):
        num, den = a / ms + pseudo, b / mg + pseudo
        if num <= 0 or den <= 0:
            raise InputFormatError(
                f"window {i} has zero normalised depth; log2 ratio undefined.",
                code="ZERO_DEPTH",
                hint="Set options['pseudocount'] (for example 0.5) or drop empty windows.",
            )
        ratios.append(math.log2(num / den))
    if smooth > 1:
        half = smooth // 2
        ratios = [_median(ratios[max(0, i - half) : i + half + 1]) for i in range(len(ratios))]
    centre = opt_num(options, "center", 0.0, minimum=-100.0, maximum=100.0)
    windows = [
        {
            "start": i * window_size,
            "end": (i + 1) * window_size,
            "log2_s_over_g1": _round(r),
            "timing": "early" if r > centre else "late" if r < centre else "intermediate",
        }
        for i, r in enumerate(ratios)
    ]
    peaks = [
        i
        for i in range(1, len(ratios) - 1)
        if ratios[i] > ratios[i - 1] and ratios[i] >= ratios[i + 1] and ratios[i] > centre
    ]
    early = sum(1 for r in ratios if r > centre)
    return {
        "windows": len(ratios),
        "pseudocount": pseudo,
        "smoothing_windows": smooth,
        "early_fraction": _round(early / len(ratios)),
        "late_fraction": _round(sum(1 for r in ratios if r < centre) / len(ratios)),
        "early_peaks": [
            {"window": i, "start": i * window_size, "log2_s_over_g1": _round(ratios[i])}
            for i in peaks[:1000]
        ],
        "window_table": windows[
            : opt_int(options, "max_windows_returned", 20_000, minimum=1, maximum=1_000_000)
        ],
        "definition": "Replication timing score = log2( (S-phase depth / median S depth) / (G1 depth / median G1 depth) ) per window "
        "(Koren et al. 2012 style): higher = replicated earlier in S phase, because early-replicating regions have more copies in S-phase "
        "cells. 'early' / 'late' means above / below `center` (default 0, the median-normalised ratio of 1). Optional running-median smoothing.",
    }


# --------------------------------------------------------------------------- replication simulation


def _origins(options: Mapping[str, object], length: int) -> list[tuple[float, float]]:
    raw = options.get("origins")
    out: list[tuple[float, float]] = []
    if raw is None:
        n = opt_int(options, "n_origins", 1, minimum=1, maximum=10_000)
        seed = opt_int(options, "seed", 0, minimum=0, maximum=2**32 - 1)
        rng = random.Random(seed)
        out = [(float(rng.randrange(length)), 0.0) for _ in range(n)]
        return out
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)) or not raw:
        raise ConfigurationError(
            "origins must be a non-empty list of positions or {position, firing_time} mappings.",
            code="INVALID_OPTION",
        )
    for item in raw:
        if isinstance(item, Mapping):
            pos, t = item.get("position"), item.get("firing_time", 0.0)
        else:
            pos, t = item, 0.0
        if isinstance(pos, bool) or not isinstance(pos, (int, float)) or not 0 <= pos < length:
            raise ConfigurationError(
                f"origin position {pos!r} must be in [0, {length}).", code="INVALID_OPTION"
            )
        if isinstance(t, bool) or not isinstance(t, (int, float)) or t < 0:
            raise ConfigurationError(
                f"origin firing_time {t!r} must be a non-negative number.", code="INVALID_OPTION"
            )
        out.append((float(pos), float(t)))
    return out


def op_simulate(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    if isinstance(value, Mapping) and "length" in value:
        value = value["length"]
    if isinstance(value, int) and not isinstance(value, bool):
        length = value
        if length < 2 or length > 5_000_000_000:
            raise InputFormatError(
                "genome length must be an integer in [2, 5e9].", code="INVALID_INPUT"
            )
    else:
        length = len(dna_text(value))
    speed = opt_num(options, "fork_speed", 1000.0, minimum=1e-9, maximum=1e12)  # bp per time unit
    circular = opt_bool(options, "circular", False)
    window = opt_int(options, "window", max(1, length // 100), minimum=1, maximum=length)
    origins = _origins(options, length)

    def timing(x: float) -> float:
        best = math.inf
        for pos, t in origins:
            d = abs(x - pos)
            if circular:
                d = min(d, length - d)
            best = min(best, t + d / speed)
        return best

    n_windows = -(-length // window)
    if n_windows > 200_000:
        raise ConfigurationError(
            "window yields more than 200000 windows.",
            code="INVALID_OPTION",
            hint="Increase options['window'].",
        )
    centres = [min(length - 1, i * window + window // 2) for i in range(n_windows)]
    times = [timing(c) for c in centres]
    sorted_origins = sorted(origins)
    terminations = []
    ring = sorted_origins + (
        [(sorted_origins[0][0] + length, sorted_origins[0][1])] if circular else []
    )
    for (p0, t0), (p1, t1) in itertools.pairwise(ring):
        meet = (p0 + p1) / 2 + speed * (t1 - t0) / 2
        if p0 < meet < p1:
            terminations.append(
                {
                    "position": _round(meet % length if circular else meet, 3),
                    "time": _round(t0 + (meet - p0) / speed, 6),
                }
            )
    duration = max(times)
    steps = opt_int(options, "time_points", 5, minimum=2, maximum=1000)
    curve = []
    for k in range(steps + 1):
        t = duration * k / steps
        curve.append(
            {
                "time": _round(t, 6),
                "fraction_replicated": _round(sum(1 for x in times if x <= t) / len(times)),
            }
        )
    return {
        "length": length,
        "circular": circular,
        "fork_speed_bp_per_time": speed,
        "origins": [{"position": p, "firing_time": t} for p, t in sorted_origins],
        "s_phase_duration": _round(duration, 6),
        "mean_replication_time": _round(math.fsum(times) / len(times), 6),
        "termination_sites": terminations,
        "window": window,
        "timing_profile": [
            {"start": i * window, "end": min(length, (i + 1) * window), "time": _round(t, 6)}
            for i, t in enumerate(times)
        ],
        "replication_curve": curve,
        "definition": "Deterministic bidirectional-fork model: all origins fire at their firing_time and forks move at fork_speed in both "
        "directions; the replication time of position x is min over origins of firing_time + |x - origin| / fork_speed (passive "
        "replication by a neighbour's fork included). Termination sites are where converging forks meet. No stochastic firing, "
        "fork stalling or origin licensing is modelled; random origins (n_origins, seed) are uniform positions fired at time 0.",
    }


# --------------------------------------------------------------------------- recombination map


def _map_distance(r: float, fn: str) -> float:
    """Map distance in Morgans from a recombination fraction (Haldane, Kosambi) or r itself ('none')."""

    if fn == "none":
        return r
    if not 0 <= r < 0.5:
        raise InputFormatError(
            f"recombination fraction {r} must be in [0, 0.5) for the {fn} map function.",
            code="INVALID_INPUT",
        )
    if fn == "haldane":
        return -0.5 * math.log(1 - 2 * r)
    return 0.25 * math.log((1 + 2 * r) / (1 - 2 * r))


def op_build_recombination_map(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    markers = value["markers"] if isinstance(value, Mapping) and "markers" in value else value
    if not isinstance(markers, Sequence) or isinstance(markers, (str, bytes)) or len(markers) < 2:
        raise InputFormatError(
            "input must be a list (or {'markers': [...]}) of at least two markers ordered by position.",
            code="INVALID_INPUT",
            hint="Each marker: {'id', 'position' (bp), and 'cm' OR 'recombination_fraction' OR 'recombinants'+'meioses' relative to the previous marker}.",
        )
    fn = opt_choice(options, "map_function", "haldane", ("haldane", "kosambi", "none"))
    rows = []
    cumulative_cm = 0.0
    previous = None
    for i, m in enumerate(markers):
        if (
            not isinstance(m, Mapping)
            or not isinstance(m.get("position"), int)
            or isinstance(m.get("position"), bool)
            or m["position"] < 0
        ):
            raise InputFormatError(
                f"marker {i + 1} needs a non-negative integer 'position'.", code="INVALID_INPUT"
            )
        pos = m["position"]
        name = str(m.get("id", f"marker{i + 1}"))
        if previous is None:
            if "cm" in m:
                cumulative_cm = float(m["cm"])
            previous = pos
            rows.append({"id": name, "position": pos, "cumulative_cm": _round(cumulative_cm)})
            continue
        if pos <= previous:
            raise InputFormatError(
                f"marker {name!r}: positions must be strictly increasing.", code="INVALID_INPUT"
            )
        if "cm" in m:
            cm = float(m["cm"])
            if cm < cumulative_cm:
                raise InputFormatError(
                    f"marker {name!r}: cM must not decrease.", code="INVALID_INPUT"
                )
            delta_cm, r = cm - cumulative_cm, None
        else:
            if "recombination_fraction" in m:
                r = float(m["recombination_fraction"])
            elif "recombinants" in m and "meioses" in m and m["meioses"]:
                r = m["recombinants"] / m["meioses"]
            else:
                raise InputFormatError(
                    f"marker {name!r} needs 'cm', 'recombination_fraction' or 'recombinants'+'meioses'.",
                    code="INVALID_INPUT",
                )
            if not 0 <= r <= 0.5:
                raise InputFormatError(
                    f"marker {name!r}: recombination fraction {r} must be in [0, 0.5].",
                    code="INVALID_INPUT",
                )
            delta_cm = 100.0 * _map_distance(r, fn)
        cumulative_cm += delta_cm
        span = pos - previous
        rows.append(
            {
                "id": name,
                "position": pos,
                "interval_start": previous,
                "interval_end": pos,
                "recombination_fraction": _round(r) if r is not None else None,
                "delta_cm": _round(delta_cm),
                "cumulative_cm": _round(cumulative_cm),
                "rate_cm_per_mb": _round(delta_cm / (span / 1e6)),
            }
        )
        previous = pos
    intervals = [r for r in rows if "delta_cm" in r]
    total_cm = sum(r["delta_cm"] for r in intervals)
    span_total = rows[-1]["position"] - rows[0]["position"]
    return {
        "map_function": fn
        if any(
            "recombination_fraction" in r and r["recombination_fraction"] is not None for r in rows
        )
        else "not used (cM given)",
        "markers": len(rows),
        "total_cm": _round(total_cm),
        "span_bp": span_total,
        "mean_rate_cm_per_mb": _round(total_cm / (span_total / 1e6)),
        "map": rows,
        "definition": "Interval rate = delta cM / (interval length in Mb). Distances from recombination fractions r use the Haldane "
        "(d = -1/2 ln(1-2r)) or Kosambi (d = 1/4 ln((1+2r)/(1-2r))) map function (Morgans x 100 = cM); r may be given directly or as "
        "recombinants/meioses. Positions are physical bp coordinates.",
    }


# --------------------------------------------------------------------------- repair pathway


def op_classify_repair_pathway(
    value: object, options: Mapping[str, object], warnings: list[str]
) -> dict[str, Any]:
    if not isinstance(value, Mapping) or "reference" not in value or "edited" not in value:
        raise InputFormatError(
            "input must be {'reference': <unedited sequence>, 'edited': <repaired sequence>[, 'donor': <HDR template product>]}.",
            code="INVALID_INPUT",
        )
    ref = dna_text(value["reference"], "reference")
    edited = dna_text(value["edited"], "edited")
    donor = dna_text(value["donor"], "donor") if value.get("donor") else None
    min_mh = opt_int(options, "min_microhomology", 2, minimum=1, maximum=20)
    ssa_mh = opt_int(options, "ssa_min_homology", 25, minimum=min_mh + 1, maximum=10_000)
    max_ins = opt_int(options, "max_nhej_insertion", 20, minimum=1, maximum=1000)
    if ref == edited:
        return {
            "pathway": "no_edit",
            "reason": "edited sequence equals the reference.",
            "rule": _RULE,
        }
    if donor is not None and edited == donor:
        return {
            "pathway": "HDR",
            "reason": "the repaired sequence equals the supplied donor/template product (template-directed repair).",
            "rule": _RULE,
        }
    lr, le = len(ref), len(edited)
    p = 0
    while p < min(lr, le) and ref[p] == edited[p]:
        p += 1
    s = 0
    while s < min(lr, le) and ref[lr - 1 - s] == edited[le - 1 - s]:
        s += 1
    sp = min(s, min(lr, le) - p)  # suffix that does not overlap the matched prefix
    core_ref, core_edit = ref[p : lr - sp], edited[p : le - sp]
    out: dict[str, Any] = {
        "reference_length": lr,
        "edited_length": le,
        "deleted_bases": len(core_ref),
        "inserted_bases": len(core_edit),
        "first_difference": p,
    }
    if lr == le and len(core_ref) == len(core_edit):
        out.update(
            pathway="undetermined",
            reason="same-length substitution(s); no end-joining deletion/insertion signature and no donor match.",
        )
    elif core_edit == "" and core_ref:
        d = lr - le
        mh = max(0, p + s - le)  # bases by which the deletion can slide = homology at the junction
        mh = min(mh, d)
        out["microhomology_length"] = mh
        if mh >= ssa_mh:
            out.update(
                pathway="SSA",
                reason=f"deletion of {d} bp flanked by {mh} bp of identical sequence (>= {ssa_mh}): single-strand annealing signature.",
            )
        elif mh >= min_mh:
            out.update(
                pathway="MMEJ",
                reason=f"deletion of {d} bp with {mh} bp microhomology at the junction (>= {min_mh}): microhomology-mediated end joining signature.",
            )
        else:
            out.update(
                pathway="NHEJ",
                reason=f"deletion of {d} bp with {mh} bp (< {min_mh}) junction microhomology: classical end-joining signature.",
            )
    elif core_ref == "" and core_edit:
        n = le - lr
        if n <= max_ins:
            out.update(
                pathway="NHEJ",
                reason=f"insertion of {n} bp without deletion (<= {max_ins}): fill-in/ligation signature of classical end joining.",
            )
        else:
            out.update(
                pathway="undetermined",
                reason=f"insertion of {n} bp exceeds {max_ins}; could be template-directed or an integration event.",
            )
    else:
        out.update(
            pathway="undetermined",
            reason=f"complex indel ({len(core_ref)} bp replaced by {len(core_edit)} bp); no single-pathway signature applies.",
        )
    out["rule"] = _RULE
    return out


_RULE = (
    "Heuristic indel-signature rule, not a measurement of the pathway used: edited == donor -> HDR; pure deletion whose two "
    "junctions share >= ssa_min_homology (25) identical bp -> SSA; >= min_microhomology (2) bp -> MMEJ; fewer -> NHEJ; pure insertion "
    "<= max_nhej_insertion (20) bp -> NHEJ; everything else (substitutions, complex indels, long insertions) -> undetermined. "
    "Junction homology = how far the deletion can slide: prefix_match + suffix_match - len(edited)."
)


# --------------------------------------------------------------------------- plot


def op_plot(value: object, options: Mapping[str, object], warnings: list[str]) -> dict[str, Any]:
    seq = dna_text(value)
    window, step, kind = _skew_params(options, len(seq))
    rows = skew_windows(seq, window, step, kind)
    est = _ori_ter(rows)
    width, height, left, right, top, bottom = 760, 380, 60, 20, 30, 45
    pw, ph = width - left - right, height - top - bottom
    values = [r["cumulative_skew"] for r in rows]
    vmin, vmax = min(values), max(values)
    span = (vmax - vmin) or 1.0
    n = len(seq)

    def px(pos: float) -> float:
        return left + pw * pos / n

    def py(v: float) -> float:
        return top + ph * (1 - (v - vmin) / span)

    pts = " ".join(
        f"{px((r['start'] + r['end']) / 2):.1f},{py(r['cumulative_skew']):.1f}" for r in rows
    )
    title = f"Cumulative {kind} skew"
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="sans-serif" font-size="11">',
        f"<title>{html.escape(title)}</title>",
        f'<rect width="{width}" height="{height}" fill="white"/>',
        f'<text x="{width / 2}" y="18" text-anchor="middle" font-size="13">{html.escape(title)} (window {window}, step {step})</text>',
        f'<rect x="{left}" y="{top}" width="{pw}" height="{ph}" fill="none" stroke="#444"/>',
        f'<polyline fill="none" stroke="#1f5fbf" stroke-width="1.5" points="{pts}"/>',
        f'<line x1="{px(est["origin_estimate"]):.1f}" x2="{px(est["origin_estimate"]):.1f}" y1="{top}" y2="{top + ph}" stroke="#2a9d3f" stroke-dasharray="4 3"/>',
        f'<text x="{px(est["origin_estimate"]) + 4:.1f}" y="{top + 12}" fill="#2a9d3f">origin (min) {est["origin_estimate"]}</text>',
        f'<line x1="{px(est["terminus_estimate"]):.1f}" x2="{px(est["terminus_estimate"]):.1f}" y1="{top}" y2="{top + ph}" stroke="#c0392b" stroke-dasharray="4 3"/>',
        f'<text x="{px(est["terminus_estimate"]) + 4:.1f}" y="{top + 26}" fill="#c0392b">terminus (max) {est["terminus_estimate"]}</text>',
        f'<text x="{left + pw / 2}" y="{height - 8}" text-anchor="middle">Position (bp)</text>',
        f'<text x="14" y="{top + ph / 2}" transform="rotate(-90 14 {top + ph / 2})" text-anchor="middle">cumulative skew</text>',
        f'<text x="{left - 5}" y="{top + 4}" text-anchor="end">{vmax:.2f}</text>',
        f'<text x="{left - 5}" y="{top + ph}" text-anchor="end">{vmin:.2f}</text>',
        "</svg>",
    ]
    text = "\n".join(svg)
    return {
        "svg": text,
        "output_path": write_svg(options, text),
        "skew_type": kind,
        "windows": len(rows),
        "origin_estimate": est["origin_estimate"],
        "terminus_estimate": est["terminus_estimate"],
        "note": "Origin/terminus markers are cumulative-skew extrema (heuristic).",
    }
