#!/usr/bin/env python3
"""Manhattan + QQ + QC summary for a GWAS sumstats file.

Auto-detects column names. Subsamples low-significance variants for speed.
Writes PNG + PDF for both plots, a QC summary JSON, and a sidecar manifest.

Usage:
    manhattan_qq.py --in <sumstats> --out-prefix <prefix> [flags]

See ../SKILL.md for full options and behavior.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_dist


CHROM_CANDIDATES = ["hm_chrom", "chromosome", "chrom", "CHR", "chr"]
POS_CANDIDATES = ["hm_pos", "base_pair_location", "position", "pos", "BP", "bp"]
P_CANDIDATES = ["p_value", "P", "pval", "p", "PVAL", "P_VAL"]
RSID_CANDIDATES = ["hm_rsid", "rsid", "SNP", "variant_id"]


def detect_columns(header: list[str], chrom_col, pos_col, p_col, rsid_col):
    def pick(explicit, candidates):
        if explicit is not None:
            if explicit not in header:
                raise SystemExit(
                    f"explicit column {explicit!r} not in header: {header}"
                )
            return explicit
        for c in candidates:
            if c in header:
                return c
        return None

    chrom = pick(chrom_col, CHROM_CANDIDATES)
    pos = pick(pos_col, POS_CANDIDATES)
    p = pick(p_col, P_CANDIDATES)
    rsid = pick(rsid_col, RSID_CANDIDATES)
    if not (chrom and pos and p):
        raise SystemExit(
            f"could not auto-detect required columns from header={header}. "
            f"chrom={chrom}, pos={pos}, p={p}. "
            f"Pass --chrom-col / --pos-col / --p-col explicitly."
        )
    return chrom, pos, p, rsid


def read_sumstats(path: Path, args) -> tuple[pd.DataFrame, dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        header_line = f.readline().rstrip("\n").rstrip("\r")
    sep = "\t" if "\t" in header_line else ","
    header = header_line.split(sep)
    chrom_col, pos_col, p_col, rsid_col = detect_columns(
        header, args.chrom_col, args.pos_col, args.p_col, args.rsid_col
    )
    usecols = [chrom_col, pos_col, p_col]
    if rsid_col:
        usecols.append(rsid_col)
    df = pd.read_csv(
        path,
        sep=sep,
        usecols=usecols,
        dtype={chrom_col: "string", pos_col: "Int64", p_col: "float64"},
        na_values=["NA", "", ".", "NaN", "nan"],
        low_memory=False,
    )
    df = df.rename(
        columns={chrom_col: "chrom", pos_col: "pos", p_col: "p_value"}
        | ({rsid_col: "rsid"} if rsid_col else {})
    )
    df = df.dropna(subset=["chrom", "pos", "p_value"])
    df = df[(df["p_value"] > 0) & (df["p_value"] <= 1)]
    autosomes = {str(i) for i in range(1, 23)}
    allowed = set(autosomes)
    if args.include_x:
        allowed |= {"X", "23"}
    df = df[df["chrom"].astype(str).str.lstrip("chr").isin(allowed)].copy()
    df["chrom_str"] = df["chrom"].astype(str).str.lstrip("chr").replace({"X": "23"})
    df["chrom_i"] = df["chrom_str"].astype(int)
    df["pos"] = df["pos"].astype(int)
    cols_meta = {
        "chrom_col": chrom_col,
        "pos_col": pos_col,
        "p_col": p_col,
        "rsid_col": rsid_col,
        "sep": sep,
    }
    return df, cols_meta


def lead_prune(df: pd.DataFrame, p_thresh: float, window_bp: int) -> pd.DataFrame:
    gw = df[df["p_value"] < p_thresh].sort_values("p_value").reset_index(drop=True)
    if gw.empty:
        return gw
    leads = []
    masked = np.zeros(len(gw), dtype=bool)
    for i in range(len(gw)):
        if masked[i]:
            continue
        leads.append(i)
        c = gw.at[i, "chrom_i"]
        p = gw.at[i, "pos"]
        same_chr = gw["chrom_i"].values == c
        near = np.abs(gw["pos"].values - p) < window_bp
        masked |= same_chr & near
    return gw.iloc[leads].reset_index(drop=True)


def plot_manhattan(df: pd.DataFrame, args, n_gw: int, out_png: Path, out_pdf: Path) -> None:
    top = df[df["p_value"] < 1e-3]
    rest = df[df["p_value"] >= 1e-3]
    if len(rest) > 0:
        rest = rest.sample(frac=args.subsample_frac, random_state=42)
    plot_df = (
        pd.concat([top, rest])
        .sort_values(["chrom_i", "pos"])
        .reset_index(drop=True)
    )
    chrom_sizes = plot_df.groupby("chrom_i")["pos"].max()
    chroms_present = sorted(chrom_sizes.index.tolist())
    offset = {chroms_present[0]: 0}
    cum = 0
    for c in chroms_present:
        cum_prev = cum
        cum += int(chrom_sizes[c])
        offset[c] = cum_prev
    next_after = cum
    plot_df["x"] = plot_df["chrom_i"].map(offset) + plot_df["pos"]
    plot_df["neglog10p"] = -np.log10(plot_df["p_value"])
    tick_pos = []
    tick_labels = []
    for c in chroms_present:
        tick_pos.append(offset[c] + int(chrom_sizes[c]) / 2)
        tick_labels.append("X" if c == 23 else str(c))
    fig, ax = plt.subplots(figsize=(14, 4.5))
    colors = ["#2b6cb0", "#a0aec0"]
    for i, c in enumerate(chroms_present):
        sub = plot_df[plot_df["chrom_i"] == c]
        if sub.empty:
            continue
        ax.scatter(sub["x"], sub["neglog10p"], s=2, c=colors[i % 2], rasterized=True)
    ax.axhline(-math.log10(args.gws_threshold), color="red", lw=0.8, ls="--", label=f"GWS {args.gws_threshold:.0e}")
    ax.axhline(-math.log10(args.suggestive_threshold), color="orange", lw=0.6, ls=":", label=f"Sugg {args.suggestive_threshold:.0e}")
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_xlim(0, next_after)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_ylim(0, max(plot_df["neglog10p"].max() * 1.05, 10))
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_qq(p: np.ndarray, out_png: Path, out_pdf: Path) -> None:
    p_sorted = np.sort(p)
    n = len(p_sorted)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(p_sorted)
    # Always keep significant tail; log-thin the bulk.
    keep_mask = obs > 1
    idx_lo = np.where(~keep_mask)[0]
    if len(idx_lo) > 0:
        log_idx = np.unique(
            np.logspace(0, np.log10(len(idx_lo)), num=3000).astype(int) - 1
        )
        log_idx = log_idx[log_idx < len(idx_lo)]
        keep = np.zeros_like(keep_mask)
        keep[idx_lo[log_idx]] = True
        keep |= keep_mask
    else:
        keep = keep_mask
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(exp[keep], obs[keep], s=4, c="#2b6cb0", rasterized=True)
    exp_max = float(exp.max())
    obs_max = float(obs.max())
    # Diagonal is the null y=x, bounded to the expected-p range so it does
    # not stretch out and visually pinch the bulk into a thin vertical stripe
    # when there are extreme-significance outliers (e.g. MHC).
    ax.plot([0, exp_max], [0, exp_max], color="red", lw=0.8)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax.set_xlim(0, exp_max * 1.05)
    ax.set_ylim(0, obs_max * 1.05)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_dataset_okg_id(args, in_path: Path) -> Optional[str]:
    if args.dataset_okg_id:
        return args.dataset_okg_id
    # Look for a sibling .fetch.json from the gwas-fetch skill.
    sidecar = in_path.with_suffix(in_path.suffix + ".fetch.json")
    if not sidecar.exists():
        # also handle the case where the file is .h.tsv.gz and sidecar is alongside
        sidecar = Path(str(in_path) + ".fetch.json")
    if sidecar.exists():
        try:
            data = json.loads(sidecar.read_text())
            okg_ids = data.get("okg_node_ids") or {}
            return okg_ids.get("dataset")
        except Exception:
            return None
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Manhattan + QQ + QC for GWAS sumstats")
    ap.add_argument("--in", dest="in_path", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--chrom-col")
    ap.add_argument("--pos-col")
    ap.add_argument("--p-col")
    ap.add_argument("--rsid-col")
    ap.add_argument("--gws-threshold", type=float, default=5e-8)
    ap.add_argument("--suggestive-threshold", type=float, default=1e-5)
    ap.add_argument("--lead-window-kb", type=int, default=500)
    ap.add_argument("--top-loci", type=int, default=20)
    ap.add_argument("--include-x", action="store_true")
    ap.add_argument("--subsample-frac", type=float, default=0.01)
    ap.add_argument("--okg-repo", type=Path, default=None)
    ap.add_argument("--dataset-okg-id", default=None)
    args = ap.parse_args()

    in_path = args.in_path.resolve()
    out_prefix = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    print(f"reading {in_path}", file=sys.stderr)
    df, cols_meta = read_sumstats(in_path, args)
    n_total = len(df)
    if n_total < 100_000:
        print(
            f"error: only {n_total:,} valid variants after filtering — "
            f"Manhattan/QQ are uninformative below 100k.",
            file=sys.stderr,
        )
        return 2

    chi2 = chi2_dist.isf(df["p_value"].values, df=1)
    chi2 = chi2[np.isfinite(chi2)]
    lambda_gc = float(np.median(chi2) / 0.4549)

    n_gw = int((df["p_value"] < args.gws_threshold).sum())
    n_sugg = int((df["p_value"] < args.suggestive_threshold).sum())

    leads = lead_prune(df, args.gws_threshold, args.lead_window_kb * 1000)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    manhattan_png = Path(str(out_prefix) + ".manhattan.png")
    manhattan_pdf = Path(str(out_prefix) + ".manhattan.pdf")
    qq_png = Path(str(out_prefix) + ".qq.png")
    qq_pdf = Path(str(out_prefix) + ".qq.pdf")
    summary_path = Path(str(out_prefix) + ".qc_summary.json")
    manifest_path = Path(str(out_prefix) + ".manhattan_qq.json")

    print(f"rendering Manhattan -> {manhattan_png}", file=sys.stderr)
    plot_manhattan(df, args, n_gw, manhattan_png, manhattan_pdf)
    print(f"rendering QQ -> {qq_png}", file=sys.stderr)
    plot_qq(df["p_value"].values, qq_png, qq_pdf)

    top_loci_records = []
    for _, row in leads.head(args.top_loci).iterrows():
        rec = {
            "chrom": int(row["chrom_i"]),
            "pos": int(row["pos"]),
            "p_value": float(row["p_value"]),
        }
        if "rsid" in df.columns and "rsid" in row.index:
            rsid = row.get("rsid")
            if pd.notna(rsid):
                rec["rsid"] = str(rsid)
        top_loci_records.append(rec)

    summary = {
        "n_variants_total": int(n_total),
        "n_genome_wide_5e8": n_gw,
        "n_suggestive_1e5": n_sugg,
        "lambda_gc": lambda_gc,
        f"n_independent_loci_{args.lead_window_kb}kb": int(len(leads)),
        "top_loci": top_loci_records,
    }
    summary_path.write_text(json.dumps(summary, indent=2))

    okg_id = resolve_dataset_okg_id(args, in_path)

    manifest = {
        "skill": "manhattan-qq",
        "version": "0.1",
        "input": str(in_path),
        "input_sha256": sha256(in_path),
        "input_columns_detected": cols_meta,
        "outputs": {
            "manhattan_png": str(manhattan_png),
            "manhattan_pdf": str(manhattan_pdf),
            "qq_png": str(qq_png),
            "qq_pdf": str(qq_pdf),
            "qc_summary": str(summary_path),
        },
        "thresholds": {
            "gws": args.gws_threshold,
            "suggestive": args.suggestive_threshold,
        },
        "lead_window_kb": args.lead_window_kb,
        "include_x": bool(args.include_x),
        "subsample_frac": args.subsample_frac,
        "qc": summary,
        "okg_node_ids": ({"dataset": okg_id} if okg_id else {}),
        "captured_at": datetime.now(timezone.utc).isoformat(),
    }
    if lambda_gc > 2:
        manifest["warnings"] = [
            f"lambda_gc={lambda_gc:.3f} > 2.0 — likely population stratification or other systematic inflation. "
            f"Run `ldsc h2` and inspect the LDSC intercept before downstream analyses."
        ]
    manifest_path.write_text(json.dumps(manifest, indent=2))

    print(json.dumps({
        "n_variants_total": n_total,
        "n_genome_wide_5e8": n_gw,
        "n_suggestive_1e5": n_sugg,
        "lambda_gc": round(lambda_gc, 4),
        "n_independent_loci": int(len(leads)),
        "top_loci_preview": top_loci_records[:5],
        "outputs": [str(manhattan_png), str(qq_png), str(summary_path), str(manifest_path)],
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
