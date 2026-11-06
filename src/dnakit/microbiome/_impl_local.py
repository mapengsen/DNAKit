"""Native implementations behind the microbiome candidate functions.

Everything here is deterministic standard-library code over small in-memory abundance tables.
Optional scientific libraries (scipy) are imported lazily and only where a function says so.
"""

from __future__ import annotations

import math
import random
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from dnakit.exceptions import ConfigurationError, InputFormatError

MAX_SAMPLES = 5_000
MAX_TAXA = 200_000
MAX_CELLS = 20_000_000
MAX_TEXT_BYTES = 50_000_000
MAX_READS_PER_SAMPLE = 5_000_000


# ----------------------------------------------------------------------------- option helpers
def opt_str(
    options: Mapping[str, object], key: str, default: str, allowed: Sequence[str] | None = None
) -> str:
    value = options.get(key, default)
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(
            f"Option {key!r} must be a non-empty string.", code="MICROBIOME_OPTION_INVALID"
        )
    value = value.strip().lower()
    if allowed is not None and value not in allowed:
        raise ConfigurationError(
            f"Option {key!r}={value!r} is not supported.",
            code="MICROBIOME_OPTION_INVALID",
            hint=f"Use one of: {', '.join(allowed)}.",
        )
    return value


def opt_float(
    options: Mapping[str, object],
    key: str,
    default: float,
    minimum: float = -math.inf,
    maximum: float = math.inf,
) -> float:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ConfigurationError(
            f"Option {key!r} must be a finite number.", code="MICROBIOME_OPTION_INVALID"
        )
    if not minimum <= value <= maximum:
        raise ConfigurationError(
            f"Option {key!r}={value} is outside [{minimum}, {maximum}].",
            code="MICROBIOME_OPTION_INVALID",
        )
    return float(value)


def opt_int(
    options: Mapping[str, object], key: str, default: int, minimum: int, maximum: int
) -> int:
    value = options.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigurationError(
            f"Option {key!r} must be an integer.", code="MICROBIOME_OPTION_INVALID"
        )
    if not minimum <= value <= maximum:
        raise ConfigurationError(
            f"Option {key!r}={value} is outside [{minimum}, {maximum}].",
            code="MICROBIOME_OPTION_INVALID",
        )
    return value


def opt_bool(options: Mapping[str, object], key: str, default: bool) -> bool:
    value = options.get(key, default)
    if not isinstance(value, bool):
        raise ConfigurationError(
            f"Option {key!r} must be a boolean.", code="MICROBIOME_OPTION_INVALID"
        )
    return value


def entity(options: Mapping[str, object], default: str, allowed: Sequence[str]) -> str:
    return opt_str(options, "entity", default, allowed)


# ----------------------------------------------------------------------------- abundance tables
@dataclass
class Table:
    """Samples x taxa count matrix (rows are samples)."""

    samples: list[str]
    taxa: list[str]
    counts: list[list[float]]

    def as_dict(self) -> dict[str, object]:
        return {"samples": self.samples, "taxa": self.taxa, "counts": self.counts}


def _number(value: object, where: str) -> float:
    if isinstance(value, bool):
        raise InputFormatError(
            "Abundance values must be numbers.", code="MICROBIOME_TABLE_INVALID", hint=where
        )
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise InputFormatError(
            f"Non-numeric abundance value {value!r}.", code="MICROBIOME_TABLE_INVALID", hint=where
        ) from exc
    if not math.isfinite(number) or number < 0:
        raise InputFormatError(
            f"Abundance values must be finite and non-negative, got {value!r}.",
            code="MICROBIOME_TABLE_INVALID",
            hint=where,
        )
    return number


def _check_shape(
    samples: Sequence[str], taxa: Sequence[str], counts: Sequence[Sequence[float]]
) -> None:
    if not samples or not taxa:
        raise InputFormatError("The abundance table is empty.", code="MICROBIOME_TABLE_EMPTY")
    if len(samples) > MAX_SAMPLES or len(taxa) > MAX_TAXA or len(samples) * len(taxa) > MAX_CELLS:
        raise ConfigurationError(
            "The abundance table exceeds the size limits.",
            code="MICROBIOME_TABLE_TOO_LARGE",
            hint=f"At most {MAX_SAMPLES} samples, {MAX_TAXA} taxa and {MAX_CELLS} cells.",
        )
    if len(set(samples)) != len(samples) or len(set(taxa)) != len(taxa):
        raise InputFormatError(
            "Sample and taxon identifiers must be unique.", code="MICROBIOME_TABLE_INVALID"
        )
    if len(counts) != len(samples) or any(len(row) != len(taxa) for row in counts):
        raise InputFormatError(
            "Count matrix shape does not match samples x taxa.", code="MICROBIOME_TABLE_INVALID"
        )


def _looks_like_path(text: str) -> bool:
    if "\n" in text or "\t" in text or len(text) >= 4096:
        return False
    try:
        return Path(text).is_file()
    except OSError:  # e.g. file name too long: it is text, not a path
        return False


def read_text_input(value: object) -> str:
    """Return text from a str/bytes/Path input; a str naming an existing file is read from disk."""

    if isinstance(value, Path) or (isinstance(value, str) and _looks_like_path(value)):
        path = Path(value)
        if path.stat().st_size > MAX_TEXT_BYTES:
            raise ConfigurationError("Input file is too large.", code="MICROBIOME_INPUT_TOO_LARGE")
        return path.read_text(encoding="utf-8", errors="replace")
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        return value
    raise InputFormatError("Expected text, bytes or a path.", code="MICROBIOME_INPUT_INVALID")


def parse_delimited_table(text: str, orientation: str = "taxa_rows") -> Table:
    """Parse a QIIME/BIOM-TSV style table. Lines starting with '#' other than the header are comments."""

    lines = [line for line in text.splitlines() if line.strip()]
    header_index = None
    for index, line in enumerate(lines):
        if (
            line.startswith("#")
            and not line.startswith("#OTU")
            and not line.startswith("#Feature")
            and not line.startswith("#ID")
        ):
            continue
        header_index = index
        break
    if header_index is None:
        raise InputFormatError("No table header found.", code="MICROBIOME_TABLE_EMPTY")
    sep = "\t" if "\t" in lines[header_index] else ","
    header = [cell.strip() for cell in lines[header_index].lstrip("#").split(sep)]
    if len(header) < 2:
        raise InputFormatError(
            "The header needs an identifier column and at least one data column.",
            code="MICROBIOME_TABLE_INVALID",
        )
    labels: list[str] = []
    rows: list[list[float]] = []
    for line in lines[header_index + 1 :]:
        if line.startswith("#"):
            continue
        cells = [cell.strip() for cell in line.split(sep)]
        if len(cells) != len(header):
            raise InputFormatError(
                f"Row {cells[0]!r} has {len(cells)} fields, expected {len(header)}.",
                code="MICROBIOME_TABLE_INVALID",
            )
        labels.append(cells[0])
        rows.append([_number(cell, f"row {cells[0]!r}") for cell in cells[1:]])
    columns = header[1:]
    if orientation == "taxa_rows":  # rows are taxa, columns are samples
        counts = [[rows[t][s] for t in range(len(rows))] for s in range(len(columns))]
        table = Table(columns, labels, counts)
    else:
        table = Table(labels, columns, rows)
    _check_shape(table.samples, table.taxa, table.counts)
    return table


def parse_table(value: object, options: Mapping[str, object]) -> Table:
    """Accept {sample: {taxon: count}}, {samples, taxa, counts}, a list of count rows, or table text/path."""

    if isinstance(value, Mapping):
        if "counts" in value:
            counts = [[_number(c, "counts") for c in row] for row in value["counts"]]  # type: ignore[union-attr]
            samples = [
                str(s) for s in value.get("samples", [f"sample{i + 1}" for i in range(len(counts))])
            ]  # type: ignore[union-attr]
            width = len(counts[0]) if counts else 0
            taxa = [str(t) for t in value.get("taxa", [f"taxon{i + 1}" for i in range(width)])]  # type: ignore[union-attr]
            table = Table(samples, taxa, counts)
        else:
            samples = [str(s) for s in value]
            taxa: list[str] = []
            seen: set[str] = set()
            for row in value.values():
                if not isinstance(row, Mapping):
                    raise InputFormatError(
                        "A sample-to-taxon mapping needs mapping values.",
                        code="MICROBIOME_TABLE_INVALID",
                    )
                for name in row:
                    if str(name) not in seen:
                        seen.add(str(name))
                        taxa.append(str(name))
            counts = [
                [_number(row.get(t, row.get(str(t), 0)), f"sample {s!r}") for t in taxa]
                for s, row in zip(samples, value.values(), strict=False)
            ]
            table = Table(samples, taxa, counts)
    elif isinstance(value, (str, bytes, Path)):
        orientation = opt_str(options, "orientation", "taxa_rows", ("taxa_rows", "samples_rows"))
        return parse_delimited_table(read_text_input(value), orientation)
    elif isinstance(value, Sequence):
        counts = [[_number(c, "counts") for c in row] for row in value]
        width = len(counts[0]) if counts else 0
        samples = [
            str(s) for s in options.get("samples", [f"sample{i + 1}" for i in range(len(counts))])
        ]  # type: ignore[union-attr]
        taxa = [str(t) for t in options.get("taxa", [f"taxon{i + 1}" for i in range(width)])]  # type: ignore[union-attr]
        table = Table(samples, taxa, counts)
    else:
        raise InputFormatError(
            "Unsupported abundance table input.", code="MICROBIOME_TABLE_INVALID"
        )
    _check_shape(table.samples, table.taxa, table.counts)
    return table


def require_depth(table: Table) -> None:
    if any(sum(row) <= 0 for row in table.counts):
        raise ConfigurationError(
            "Every sample needs a positive total abundance.",
            code="MICROBIOME_EMPTY_SAMPLE",
            hint="Remove empty samples first (see microbiome.qc).",
        )


# ----------------------------------------------------------------------------- alpha diversity
def _is_int(x: float) -> bool:
    return abs(x - round(x)) < 1e-9


def alpha_metrics(row: Sequence[float], log_base: float = math.e) -> dict[str, float | None]:
    """Textbook per-sample richness/diversity metrics.

    observed = #taxa with count>0; shannon H = -sum p ln p; simpson (Gini-Simpson) = 1 - sum p^2;
    inverse_simpson = 1 / sum p^2; pielou = H / ln(observed) (None when observed <= 1);
    chao1 (bias-corrected) = S_obs + F1 (F1 - 1) / (2 (F2 + 1)), goods_coverage = 1 - F1 / N
    with F1/F2 the numbers of singleton/doubleton taxa (integer counts only, else None).
    """

    total = sum(row)
    present = [x for x in row if x > 0]
    observed = len(present)
    p = [x / total for x in present]
    shannon = -sum(v * math.log(v) for v in p) / math.log(log_base)
    sum_sq = sum(v * v for v in p)
    integer = all(_is_int(x) for x in present)
    f1 = sum(1 for x in present if round(x) == 1) if integer else 0
    f2 = sum(1 for x in present if round(x) == 2) if integer else 0
    return {
        "observed": float(observed),
        "shannon": shannon,
        "simpson": 1.0 - sum_sq,
        "inverse_simpson": 1.0 / sum_sq,
        "pielou_evenness": shannon / (math.log(observed) / math.log(log_base))
        if observed > 1
        else None,
        "chao1": observed + f1 * (f1 - 1) / (2.0 * (f2 + 1)) if integer else None,
        "goods_coverage": 1.0 - f1 / total if integer else None,
    }


# ----------------------------------------------------------------------------- beta diversity
def _clr(row: Sequence[float], pseudocount: float) -> list[float]:
    logs = [math.log(x + pseudocount) for x in row]
    mean = sum(logs) / len(logs)
    return [v - mean for v in logs]


def beta_distance(
    a: Sequence[float], b: Sequence[float], metric: str, pseudocount: float = 1.0
) -> float:
    if metric == "braycurtis":
        denom = sum(a) + sum(b)
        return sum(abs(x - y) for x, y in zip(a, b, strict=False)) / denom if denom else 0.0
    if metric == "jaccard":
        union = sum(1 for x, y in zip(a, b, strict=False) if x > 0 or y > 0)
        inter = sum(1 for x, y in zip(a, b, strict=False) if x > 0 and y > 0)
        return 1.0 - inter / union if union else 0.0
    if metric == "aitchison":
        ca, cb = _clr(a, pseudocount), _clr(b, pseudocount)
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(ca, cb, strict=False)))
    if metric == "euclidean":
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=False)))
    raise ConfigurationError("Unknown beta-diversity metric.", code="MICROBIOME_OPTION_INVALID")


def distance_matrix(table: Table, metric: str, pseudocount: float = 1.0) -> list[list[float]]:
    n = len(table.samples)
    out = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            out[i][j] = out[j][i] = beta_distance(
                table.counts[i], table.counts[j], metric, pseudocount
            )
    return out


# ----------------------------------------------------------------------------- rarefaction
def rarefy(table: Table, depth: int, seed: int, drop_short: bool = True) -> tuple[Table, list[str]]:
    """Subsample every sample to ``depth`` reads without replacement (random.Random(seed) per sample)."""

    rng = random.Random(seed)
    samples: list[str] = []
    rows: list[list[float]] = []
    dropped: list[str] = []
    for name, row in zip(table.samples, table.counts, strict=False):
        if not all(_is_int(x) for x in row):
            raise ConfigurationError(
                "Rarefaction needs integer counts.", code="MICROBIOME_TABLE_INVALID"
            )
        total = int(round(sum(row)))
        if total > MAX_READS_PER_SAMPLE:
            raise ConfigurationError(
                "Sample depth exceeds the rarefaction limit.", code="MICROBIOME_TABLE_TOO_LARGE"
            )
        if total < depth:
            if drop_short:
                dropped.append(name)
                continue
            raise ConfigurationError(
                f"Sample {name!r} has {total} reads, fewer than depth {depth}.",
                code="MICROBIOME_DEPTH_TOO_HIGH",
            )
        population = [i for i, x in enumerate(row) for _ in range(int(round(x)))]
        picked = rng.sample(population, depth)
        new = [0.0] * len(row)
        for index in picked:
            new[index] += 1
        samples.append(name)
        rows.append(new)
    return Table(samples, list(table.taxa), rows), dropped


# ----------------------------------------------------------------------------- taxonomy
_PREFIX_RANK = {
    "d": "domain",
    "k": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
_POSITIONAL = ("domain", "phylum", "class", "order", "family", "genus", "species")
_RANK_ORDER = _POSITIONAL
_PLACEHOLDER = {
    "",
    "unassigned",
    "unclassified",
    "unknown",
    "uncultured",
    "ambiguous",
    "incertae sedis",
    "na",
    "none",
}


def parse_taxonomy(text: str) -> dict[str, object]:
    """Split a lineage ('k__Bacteria; p__Firmicutes', 'd__..;p__..' or unprefixed 'A;B;C') into ranks."""

    parts = [p.strip() for p in re.split(r"[;|]", text) if p.strip() != ""]
    ranks: dict[str, str | None] = {rank: None for rank in _RANK_ORDER}
    prefixed = [m for m in (re.match(r"^([A-Za-z])__(.*)$", p) for p in parts) if m]
    if parts and len(prefixed) == len(parts):
        for match in prefixed:
            rank = _PREFIX_RANK.get(match.group(1).lower())
            if rank is None:
                continue
            name = match.group(2).strip()
            ranks[rank] = None if name.lower() in _PLACEHOLDER else name
    else:
        if len(parts) > len(_POSITIONAL):
            raise InputFormatError(
                f"Lineage has {len(parts)} unprefixed levels; at most {len(_POSITIONAL)} are supported.",
                code="MICROBIOME_TAXONOMY_INVALID",
                hint="Use rank prefixes such as 'g__Escherichia'.",
            )
        for rank, part in zip(_POSITIONAL, parts, strict=False):
            ranks[rank] = None if part.lower() in _PLACEHOLDER else part
    assigned = [(rank, ranks[rank]) for rank in _RANK_ORDER if ranks[rank]]
    return {
        "ranks": ranks,
        "lowest_rank": assigned[-1][0] if assigned else None,
        "lowest_taxon": assigned[-1][1] if assigned else None,
    }


# ----------------------------------------------------------------------------- statistics helpers
def pearson(x: Sequence[float], y: Sequence[float]) -> float | None:
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y, strict=False)) / math.sqrt(sxx * syy)


def average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = mean_rank
        i = j + 1
    return ranks


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


# ----------------------------------------------------------------------------- k-mer ANI estimate
_COMP = str.maketrans("ACGT", "TGCA")


def kmer_set(sequence: str, k: int, canonical: bool) -> set[str]:
    seq = sequence.upper()
    out: set[str] = set()
    for i in range(len(seq) - k + 1):
        word = seq[i : i + k]
        if set(word) - set("ACGT"):
            continue
        if canonical:
            rc = word.translate(_COMP)[::-1]
            word = min(word, rc)
        out.add(word)
    return out


def kmer_ani(a: set[str], b: set[str], k: int) -> dict[str, float | None]:
    """Jaccard/containment of two k-mer sets and the Mash-formula distance D = -ln(2j/(1+j))/k."""

    inter = len(a & b)
    union = len(a | b)
    jaccard = inter / union if union else 0.0
    containment_a = inter / len(a) if a else 0.0
    containment_b = inter / len(b) if b else 0.0
    if jaccard > 0:
        distance = -math.log(2 * jaccard / (1 + jaccard)) / k
        distance = max(distance, 0.0)
        ani = 1.0 - distance
    else:
        distance, ani = None, None
    return {
        "shared_kmers": float(inter),
        "jaccard": jaccard,
        "containment_a_in_b": containment_a,
        "containment_b_in_a": containment_b,
        "mash_distance": distance,
        "ani_estimate": ani,
    }


def parse_sequences(
    value: object, minimum: int = 1, maximum: int = 200, max_bases: int = 100_000_000
) -> list[tuple[str, str]]:
    """Accept a list of strings, an id->sequence mapping, or FASTA text/path; returns (id, ACGTN sequence)."""

    items: list[tuple[str, str]] = []
    if isinstance(value, Mapping):
        items = [(str(k), str(v)) for k, v in value.items()]
    elif isinstance(value, (str, bytes, Path)):
        text = read_text_input(value)
        if text.lstrip().startswith(">"):
            name, chunks = None, []
            for line in text.splitlines():
                if line.startswith(">"):
                    if name is not None:
                        items.append((name, "".join(chunks)))
                    name, chunks = (
                        (line[1:].split()[0] if line[1:].strip() else f"seq{len(items) + 1}"),
                        [],
                    )
                else:
                    chunks.append(line.strip())
            if name is not None:
                items.append((name, "".join(chunks)))
        else:
            items = [
                (f"seq{i + 1}", line.strip())
                for i, line in enumerate(text.splitlines())
                if line.strip()
            ]
    elif isinstance(value, Sequence):
        items = [(f"seq{i + 1}", str(v)) for i, v in enumerate(value)]
    else:
        raise InputFormatError("Unsupported sequence input.", code="MICROBIOME_INPUT_INVALID")
    cleaned: list[tuple[str, str]] = []
    for name, seq in items:
        seq = re.sub(r"\s+", "", seq).upper().replace("U", "T")
        if not seq or set(seq) - set("ACGTRYSWKMBDHVN"):
            raise InputFormatError(
                f"Sequence {name!r} is empty or has non-DNA characters.",
                code="MICROBIOME_SEQUENCE_INVALID",
            )
        cleaned.append((name, seq))
    if len(cleaned) < minimum or len(cleaned) > maximum:
        raise InputFormatError(
            f"Expected between {minimum} and {maximum} sequences, got {len(cleaned)}.",
            code="MICROBIOME_INPUT_INVALID",
        )
    if sum(len(s) for _, s in cleaned) > max_bases:
        raise ConfigurationError("Input sequences are too long.", code="MICROBIOME_INPUT_TOO_LARGE")
    return cleaned


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


# ----------------------------------------------------------------------------- SVG
def svg_stacked_bars(table: Table, top_n: int) -> str:
    """Stacked relative-abundance bars (one per sample), top_n taxa by total abundance plus 'Other'."""

    palette = [
        "#4e79a7",
        "#f28e2b",
        "#e15759",
        "#76b7b2",
        "#59a14f",
        "#edc948",
        "#b07aa1",
        "#ff9da7",
        "#9c755f",
        "#bab0ac",
    ]
    totals = [
        sum(table.counts[s][t] for s in range(len(table.samples))) for t in range(len(table.taxa))
    ]
    order = sorted(range(len(table.taxa)), key=lambda t: (-totals[t], table.taxa[t]))
    top = order[:top_n]
    bar_w, gap, height, left, bottom = 28, 8, 200, 40, 60
    width = left + len(table.samples) * (bar_w + gap) + 200
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height + bottom + 20}">'
    ]
    for s, name in enumerate(table.samples):
        row = table.counts[s]
        total = sum(row) or 1.0
        x = left + s * (bar_w + gap)
        y = height
        values = [
            (table.taxa[t], row[t] / total, palette[i % len(palette)]) for i, t in enumerate(top)
        ]
        other = 1.0 - sum(v for _, v, _ in values)
        if other > 1e-12:
            values.append(("Other", other, "#cccccc"))
        for taxon, fraction, color in values:
            if fraction <= 0:
                continue
            h = fraction * height
            y -= h
            parts.append(
                f'<rect class="bar" data-sample="{_esc(name)}" data-taxon="{_esc(taxon)}" data-fraction="{fraction:.6f}" '
                f'x="{x}" y="{y:.2f}" width="{bar_w}" height="{h:.2f}" fill="{color}"/>'
            )
        parts.append(f'<text x="{x}" y="{height + 14}" font-size="9">{_esc(name)}</text>')
    lx = left + len(table.samples) * (bar_w + gap) + 10
    for i, t in enumerate(top):
        parts.append(
            f'<rect x="{lx}" y="{10 + 14 * i}" width="10" height="10" fill="{palette[i % len(palette)]}"/>'
        )
        parts.append(
            f'<text x="{lx + 14}" y="{19 + 14 * i}" font-size="10">{_esc(table.taxa[t])}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def _esc(text: str) -> str:
    return (
        text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")
    )
