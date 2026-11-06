"""QC report assembly and SVG figures for the native sequencing-QC functions."""

from __future__ import annotations

import html
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit.exceptions import ConfigurationError, InputFormatError

from ._impl_metrics import (
    MAX_PHRED,
    TENX_EXTRA_ADAPTERS,
    Warnings,
    _adapter_content,
    _basecalling,
    _gc_bias,
    _keys,
    _length_stats,
    _platform,
    _qscore_summary,
    _qual_hist,
    _quality,
    _yield,
    adapters_option,
    dup_levels,
    op_barcode_quality_tenx,
    op_duplication_tenx,
    op_yield_tenx,
    tenx_regions,
)
from ._impl_reads import (
    Rd,
    hist_percentile,
    opt_bool,
    opt_choice,
    opt_int,
    parse_pair,
    round6,
)

# --------------------------------------------------------------------------- output


def write_text(options: Mapping[str, object], text: str, suffixes: Sequence[str]) -> str | None:
    """Write ``text`` to ``options['output_path']`` (refusing to overwrite unless ``overwrite`` is true)."""

    target = options.get("output_path")
    if target is None:
        return None
    if not isinstance(target, (str, Path)) or not str(target) or "\0" in str(target):
        raise ConfigurationError("output_path must be a non-empty path.", code="INVALID_OPTION")
    path = Path(target)
    if path.suffix.lower() not in suffixes:
        raise ConfigurationError(
            f"output_path must end with one of {', '.join(suffixes)}; got {path.suffix!r}.",
            code="INVALID_OPTION",
        )
    if not path.parent.is_dir():
        raise ConfigurationError(f"directory {path.parent} does not exist.", code="INVALID_OPTION")
    if path.exists() and not opt_bool(options, "overwrite", False):
        raise ConfigurationError(
            f"{path} already exists.",
            code="OUTPUT_EXISTS",
            hint="Pass overwrite=True to replace it.",
        )
    path.write_text(text, encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------- flags (FastQC thresholds)


def _flag(value: float | None, warn: float, fail: float, *, higher_is_worse: bool) -> str:
    if value is None:
        return "not_available"
    if higher_is_worse:
        return "fail" if value > fail else "warn" if value > warn else "pass"
    return "fail" if value < fail else "warn" if value < warn else "pass"


def _flags(
    quality: dict[str, Any] | None,
    n_by_position: Sequence[float],
    adapter_max: float | None,
    remaining_percent: float,
    gc_deviation: float | None,
) -> dict[str, str]:
    flags: dict[str, str] = {}
    if quality is not None:
        rows = quality["per_position"]
        min_q1 = min(r["q1"] for r in rows)
        min_med = min(r["median"] for r in rows)
        base = (
            "fail"
            if min_q1 < 5 or min_med < 20
            else "warn"
            if min_q1 < 10 or min_med < 25
            else "pass"
        )
        flags["per_base_quality"] = base
        mode = quality["per_read_mean_quality"]["mode"]
        flags["per_sequence_quality"] = _flag(mode, 27, 20, higher_is_worse=False)
    else:
        flags["per_base_quality"] = flags["per_sequence_quality"] = "not_available"
    flags["per_base_n_content"] = _flag(
        max(n_by_position) * 100 if n_by_position else 0.0, 5, 20, higher_is_worse=True
    )
    flags["adapter_content"] = (
        "not_available"
        if adapter_max is None
        else "fail"
        if adapter_max >= 10
        else "warn"
        if adapter_max >= 5
        else "pass"
    )
    flags["sequence_duplication"] = _flag(100.0 - remaining_percent, 20, 50, higher_is_worse=True)
    flags["gc_content_shape"] = _flag(gc_deviation, 15, 30, higher_is_worse=True)
    return flags


def _compact_read_set(
    reads: list[Rd], adapters: Mapping[str, str], options: Mapping[str, object]
) -> dict[str, Any]:
    has_q = any(r.qual for r in reads)
    quality = _quality(reads) if has_q else None
    gc = _gc_bias(reads, None)
    adapter = _adapter_content(reads, adapters, 12, 0)
    keys = [r.seq for r in reads]
    dup = dup_levels(Counter(keys))
    maxlen = max(len(r.seq) for r in reads)
    n_pos = [
        sum(1 for r in reads if len(r.seq) > i and r.seq[i] == "N")
        / max(1, sum(1 for r in reads if len(r.seq) > i))
        for i in range(maxlen)
    ]
    length = _length_stats(reads, options)
    length.pop("histogram")
    section: dict[str, Any] = {
        "yield": _yield(reads),
        "read_length": length,
        "platform": {k: v for k, v in _platform(reads).items() if k != "definition"},
        "gc_content": {
            k: gc[k]
            for k in (
                "mean_gc_percent",
                "stdev_gc_percent",
                "mode_gc_percent",
                "deviation_from_normal_percent",
            )
        },
        "adapter_content": {
            "fraction_reads_with_any_adapter": adapter["fraction_reads_with_any_adapter"],
            "max_cumulative_percent": {
                n: a["max_cumulative_percent"] for n, a in adapter["adapters"].items()
            },
        },
        "duplication": {
            k: dup[k] for k in ("distinct_sequences", "percent_remaining_after_deduplication")
        },
        "max_n_percent_at_a_position": round6(max(n_pos) * 100),
    }
    if quality is not None:
        section["qscore"] = _qscore_summary(reads, "probability")
        medians = [r["median"] for r in quality["per_position"]]
        q1s = [r["q1"] for r in quality["per_position"]]
        section["per_base_quality"] = {
            "positions": len(medians),
            "min_median": min(medians),
            "min_lower_quartile": min(q1s),
        }
        section["per_read_mean_quality"] = {
            k: quality["per_read_mean_quality"][k] for k in ("mean", "median", "mode")
        }
    max_adapter = max(a["max_cumulative_percent"] for a in adapter["adapters"].values())
    section["flags"] = _flags(
        quality,
        n_pos,
        max_adapter,
        dup["percent_remaining_after_deduplication"],
        gc["deviation_from_normal_percent"],
    )
    return section


def _report_text(title: str, sections: Mapping[str, Any]) -> str:
    lines = [title, "=" * len(title)]
    for name, section in sections.items():
        lines.append(f"\n[{name}]")
        y = section.get("yield", {})
        lines.append(
            f"reads: {y.get('reads')}  bases: {y.get('bases')}  mean length: {y.get('mean_read_length')}"
        )
        if y.get("q30_fraction") is not None:
            lines.append(f"Q20 fraction: {y['q20_fraction']}  Q30 fraction: {y['q30_fraction']}")
        gc = section.get("gc_content", {})
        lines.append(
            f"mean GC%: {gc.get('mean_gc_percent')}  platform: {section.get('platform', {}).get('platform')}"
        )
        for flag, state in section.get("flags", {}).items():
            lines.append(f"  {flag}: {state}")
    return "\n".join(lines) + "\n"


def op_report(value: object, options: Mapping[str, object], warnings: Warnings) -> dict[str, Any]:
    first, second = parse_pair(value, options)
    adapters = adapters_option(options)
    sections: dict[str, Any] = (
        {"read1": _compact_read_set(first, adapters, options)}
        if second is not None
        else {"reads": _compact_read_set(first, adapters, options)}
    )
    result: dict[str, Any] = {"paired": second is not None}
    if second is not None:
        sections["read2"] = _compact_read_set(second, adapters, options)
        keys, _ = _keys(value, options)
        pair_dup = dup_levels(Counter(keys))
        result["pair_duplication"] = {
            k: pair_dup[k] for k in ("distinct_sequences", "percent_remaining_after_deduplication")
        }
    result["sections"] = sections
    result["flag_definitions"] = (
        "FastQC thresholds: per-base quality warn if lower quartile < 10 or median < 25 (fail < 5 / < 20) at any position; "
        "per-sequence quality warn if modal read mean Q < 27 (fail < 20); N content warn > 5% (fail > 20%) at any position; "
        "adapter content warn >= 5% (fail >= 10%) cumulative; duplication warn if non-unique reads > 20% (fail > 50%); "
        "GC shape warn if deviation from the fitted normal > 15% (fail > 30%)."
    )
    text = _report_text("DNAKit sequencing QC report", sections)
    result["report_text"] = text
    result["output_path"] = write_text(options, text, (".txt", ".md"))
    return result


def op_report_tenx(
    value: object, options: Mapping[str, object], warnings: Warnings
) -> dict[str, Any]:
    regions, meta = tenx_regions(value, options)
    adapters = adapters_option(options, TENX_EXTRA_ADAPTERS)
    sections: dict[str, Any] = {}
    for name in ("barcode", "umi", "rna_read"):
        reads = regions[name]
        has_q = any(r.qual for r in reads)
        entry: dict[str, Any] = {"yield": _yield(reads)}
        if has_q:
            entry["basecalling"] = _basecalling(reads, 20.0)
        sections[name] = entry
    rna = regions["rna_read"]
    has_q = any(r.qual for r in rna)
    sections["rna_read"]["compact"] = _compact_read_set(rna, adapters, options)
    result: dict[str, Any] = {
        **meta,
        "sections": sections,
        "barcode_and_umi": op_barcode_quality_tenx(value, options, warnings),
        "duplication": {
            k: v
            for k, v in op_duplication_tenx(value, options, warnings).items()
            if k in ("duplicate_rate", "distinct_sequences", "definition")
        },
        "yield": {
            k: v for k, v in op_yield_tenx(value, options, warnings).items() if k != "regions"
        },
    }
    q30 = {n: s.get("basecalling", {}).get("q30_fraction") for n, s in sections.items()}
    lines = [
        "DNAKit 10x Chromium QC report",
        "=============================",
        f"chemistry {meta['chemistry']}: {meta['barcode_length']} bp barcode + {meta['umi_length']} bp UMI; read pairs: {meta['read_pairs']}",
        "Q30 fraction (barcode / UMI / RNA read): "
        + " / ".join(str(q30[n]) for n in ("barcode", "umi", "rna_read")),
        f"duplicate rate (proxy): {result['duplication']['duplicate_rate']}",
        f"distinct barcodes: {result['yield']['distinct_barcodes']}",
    ]
    if has_q:
        for flag, state in sections["rna_read"]["compact"]["flags"].items():
            lines.append(f"  RNA read {flag}: {state}")
    text = "\n".join(lines) + "\n"
    result["report_text"] = text
    result["output_path"] = write_text(options, text, (".txt", ".md"))
    return result


# --------------------------------------------------------------------------- SVG


def _svg_open(width: int, height: int, title: str) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        'font-family="sans-serif" font-size="11">',
        f"<title>{html.escape(title)}</title>",
        f'<rect width="{width}" height="{height}" fill="white"/>',
        f'<text x="{width / 2}" y="16" text-anchor="middle" font-size="13">{html.escape(title)}</text>',
    ]


def _bin_hist(hist: list[list[int]], width: int) -> list[tuple[str, list[int]]]:
    out = []
    for start in range(0, len(hist), width):
        merged = [0] * (MAX_PHRED + 1)
        chunk = hist[start : start + width]
        for h in chunk:
            for q, c in enumerate(h):
                merged[q] += c
        label = str(start + 1) if len(chunk) == 1 else f"{start + 1}-{start + len(chunk)}"
        out.append((label, merged))
    return out


def op_plot_per_base_quality(
    value: object, options: Mapping[str, object], warnings: Warnings
) -> dict[str, Any]:
    first, second = parse_pair(value, options)
    mate = opt_choice(options, "mate", "read1", ("read1", "read2"))
    if mate == "read2" and second is None:
        raise InputFormatError("mate='read2' needs a read pair.", code="PAIR_REQUIRED")
    reads = second if mate == "read2" and second is not None else first
    hist = _qual_hist(reads)
    if not hist:
        raise InputFormatError(
            "per-base quality needs FASTQ reads with Phred qualities.", code="QUALITY_REQUIRED"
        )
    max_columns = opt_int(options, "max_columns", 150, minimum=10, maximum=1000)
    group = -(-len(hist) // max_columns)
    bins = _bin_hist(hist, group)
    top = max(41, max(max((q for q, c in enumerate(h) if c), default=0) for _, h in bins) + 1)
    left, right, tp, bottom = 50, 20, 28, 50
    plot_w = max(300, min(900, 6 * len(bins)))
    plot_h = 300
    w, h_ = left + plot_w + right, tp + plot_h + bottom

    def y(q: float) -> float:
        return tp + plot_h * (1 - q / top)

    step = plot_w / len(bins)
    svg = _svg_open(w, h_, f"Quality scores per base ({mate})")
    for lo, hi, color in ((0, 20, "#f8d7da"), (20, 28, "#fff3cd"), (28, top, "#d4edda")):
        svg.append(
            f'<rect x="{left}" y="{y(min(hi, top))}" width="{plot_w}" height="{y(lo) - y(min(hi, top))}" fill="{color}"/>'
        )
    svg.append(
        f'<rect x="{left}" y="{tp}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#444"/>'
    )
    for q in range(0, int(top) + 1, 10):
        svg.append(f'<line x1="{left - 4}" x2="{left}" y1="{y(q)}" y2="{y(q)}" stroke="#444"/>')
        svg.append(f'<text x="{left - 7}" y="{y(q) + 4}" text-anchor="end">{q}</text>')
    means = []
    for i, (_label, counts) in enumerate(bins):
        n = sum(counts)
        if n == 0:
            continue
        cx = left + step * (i + 0.5)
        p10, q1, med, q3, p90 = (
            hist_percentile(counts, f) or 0.0 for f in (0.10, 0.25, 0.5, 0.75, 0.90)
        )
        means.append((cx, sum(q * c for q, c in enumerate(counts)) / n))
        bw = max(1.0, step * 0.7)
        svg.append(f'<line x1="{cx}" x2="{cx}" y1="{y(p90)}" y2="{y(p10)}" stroke="#333"/>')
        svg.append(
            f'<rect x="{cx - bw / 2}" y="{y(q3)}" width="{bw}" height="{max(0.5, y(q1) - y(q3))}" fill="#ffd966" stroke="#333" stroke-width="0.5"/>'
        )
        svg.append(
            f'<line x1="{cx - bw / 2}" x2="{cx + bw / 2}" y1="{y(med)}" y2="{y(med)}" stroke="#c00" stroke-width="1.5"/>'
        )
    if len(means) > 1:
        svg.append(
            '<polyline fill="none" stroke="#1f5fbf" stroke-width="1.2" points="'
            + " ".join(f"{x:.1f},{yy:.1f}" for x, yy in [(x, y(m)) for x, m in means])
            + '"/>'
        )
    for i in range(0, len(bins), max(1, len(bins) // 12)):
        svg.append(
            f'<text x="{left + step * (i + 0.5)}" y="{tp + plot_h + 14}" text-anchor="middle" font-size="9">{bins[i][0]}</text>'
        )
    svg.append(
        f'<text x="{left + plot_w / 2}" y="{h_ - 8}" text-anchor="middle">Position in read (bp){f", {group} positions per box" if group > 1 else ""}</text>'
    )
    svg.append(
        f'<text x="12" y="{tp + plot_h / 2}" transform="rotate(-90 12 {tp + plot_h / 2})" text-anchor="middle">Phred quality</text>'
    )
    svg.append("</svg>")
    text = "\n".join(svg)
    path = write_text(options, text, (".svg",))
    return {
        "svg": text,
        "output_path": path,
        "positions": len(hist),
        "positions_per_box": group,
        "legend": "box = lower to upper quartile, red line = median, whiskers = 10th-90th percentile, blue line = mean; "
        "background: green Q>=28, yellow 20-28, red <20.",
    }


def _bar_chart(title: str, rows: Sequence[tuple[str, float | None]]) -> str:
    left, plot_w, row_h = 230, 360, 24
    height = 40 + row_h * len(rows) + 20
    width = left + plot_w + 70
    svg = _svg_open(width, height, title)
    svg.append(f'<line x1="{left}" x2="{left}" y1="28" y2="{height - 18}" stroke="#444"/>')
    for i, (label, v) in enumerate(rows):
        yy = 36 + i * row_h
        svg.append(
            f'<text x="{left - 8}" y="{yy + 13}" text-anchor="end">{html.escape(label)}</text>'
        )
        if v is None:
            svg.append(f'<text x="{left + 6}" y="{yy + 13}" fill="#666">n/a</text>')
            continue
        bar = max(0.0, min(100.0, v)) / 100 * plot_w
        svg.append(f'<rect x="{left}" y="{yy}" width="{bar:.1f}" height="16" fill="#3b7dd8"/>')
        svg.append(f'<text x="{left + bar + 5:.1f}" y="{yy + 13}">{v:.2f}%</text>')
    svg.append(
        f'<text x="{left + plot_w / 2}" y="{height - 4}" text-anchor="middle" fill="#555">percent (0-100)</text>'
    )
    svg.append("</svg>")
    return "\n".join(svg)


def _percent(x: float | None) -> float | None:
    return None if x is None else 100.0 * x


def _summary_rows(
    reads: list[Rd], prefix: str, adapters: Mapping[str, str]
) -> list[tuple[str, float | None]]:
    y = _yield(reads)
    adapter = _adapter_content(reads, adapters, 12, 0)
    dup = dup_levels(Counter(r.seq for r in reads))
    n_pct = 100.0 * y["n_bases"] / y["bases"] if y["bases"] else None
    return [
        (f"{prefix}bases >= Q20", _percent(y["q20_fraction"])),
        (f"{prefix}bases >= Q30", _percent(y["q30_fraction"])),
        (f"{prefix}GC content", _percent(y["gc_fraction"])),
        (f"{prefix}N bases", n_pct),
        (f"{prefix}duplicated reads", 100.0 - dup["percent_remaining_after_deduplication"]),
        (f"{prefix}reads with adapter", _percent(adapter["fraction_reads_with_any_adapter"])),
    ]


def op_plot_qc_summary(
    value: object, options: Mapping[str, object], warnings: Warnings
) -> dict[str, Any]:
    first, second = parse_pair(value, options)
    adapters = adapters_option(options)
    rows = _summary_rows(first, "read1: " if second is not None else "", adapters)
    if second is not None:
        rows += _summary_rows(second, "read2: ", adapters)
    svg = _bar_chart("Sequencing QC summary", rows)
    return {
        "svg": svg,
        "output_path": write_text(options, svg, (".svg",)),
        "metrics": {label: round6(v) for label, v in rows},
    }


def op_plot_qc_summary_tenx(
    value: object, options: Mapping[str, object], warnings: Warnings
) -> dict[str, Any]:
    regions, meta = tenx_regions(value, options)
    adapters = adapters_option(options, TENX_EXTRA_ADAPTERS)
    rows: list[tuple[str, float | None]] = []
    for name, label in (
        ("barcode", "Q30 bases in barcode"),
        ("umi", "Q30 bases in UMI"),
        ("rna_read", "Q30 bases in RNA read"),
    ):
        rows.append((label, _percent(_yield(regions[name])["q30_fraction"])))
    adapter = _adapter_content(regions["rna_read"], adapters, 12, 0)
    rows.append(
        ("RNA reads with adapter/TSO/polyA", _percent(adapter["fraction_reads_with_any_adapter"]))
    )
    dup = op_duplication_tenx(value, options, [])
    rows.append(("duplicate rate (proxy)", _percent(dup["duplicate_rate"])))
    rows.append(("GC content of RNA read", _percent(_yield(regions["rna_read"])["gc_fraction"])))
    whitelist_rows = op_yield_tenx(value, options, [])
    if "valid_barcode_fraction" in whitelist_rows:
        rows.append(("valid barcodes", _percent(whitelist_rows["valid_barcode_fraction"])))
    svg = _bar_chart(f"10x Chromium QC summary ({meta['chemistry']})", rows)
    return {
        **meta,
        "svg": svg,
        "output_path": write_text(options, svg, (".svg",)),
        "metrics": {label: round6(v) for label, v in rows},
    }
