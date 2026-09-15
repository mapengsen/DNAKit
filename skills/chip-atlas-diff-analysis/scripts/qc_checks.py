"""
Automated QC pass for ChIP-Atlas diff-analysis results.

Each check returns one or more warnings of the form:

    {
        "severity": "major" | "minor",
        "issue":    "<short headline>",
        "details":  "<full sentence(s) with numbers>",
        "recommendation": "<action item>",
        ...
    }

Checks
------
* sex-chromosome confounds (all regions + significant-only)
* non-standard contigs
* implausibly small regions (<10 bp)
* chrM enrichment in ChIP/ATAC/DNase-seq modes
* dense significant-region clusters
* MA-plot asymmetry (low-count direction skew)
* unusually small median region size
* low sample size (n<=3 major, n<5 minor)

Inputs : parsed DataFrame, group sizes, analysis_type, optional genome
Outputs: list of warning dicts
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable, List, Optional

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas numpy")

try:
    from scripts.filter_regions import get_standard_chroms
except ImportError:
    from filter_regions import get_standard_chroms


def run_qc_checks(
    frame: pd.DataFrame,
    n_samples_a: int,
    n_samples_b: int,
    analysis_type: str = "diffbind",
    genome: Optional[str] = None,
) -> List[dict]:
    """Return the list of QC warnings for `frame`."""
    warnings: List[dict] = []
    if frame.empty:
        return warnings

    print("\n   Running QC checks…")

    warnings.extend(_check_sex_chrom_all(frame))
    warnings.extend(_check_sex_chrom_significant(frame))
    warnings.extend(_check_nonstandard_contigs(frame, genome=genome))
    warnings.extend(_check_small_regions(frame))
    warnings.extend(_check_chrM(frame, analysis_type))
    warnings.extend(_check_genomic_clusters(frame))
    warnings.extend(_check_ma_asymmetry(frame))
    warnings.extend(_check_region_size_distribution(frame, analysis_type))
    warnings.extend(_check_sample_size(n_samples_a, n_samples_b))

    warnings = _dedupe_sex_chrom_pairs(warnings)
    _print_summary(warnings)
    return warnings


def _dedupe_sex_chrom_pairs(warnings: List[dict]) -> List[dict]:
    """If significant-only and all-region warnings overlap, keep the specific one."""
    significant_chroms = {
        w["issue"].split(" ")[0]
        for w in warnings
        if "significant regions" in w.get("issue", "")
    }
    if not significant_chroms:
        return warnings
    return [
        w for w in warnings
        if not (
            "sex chromosome confound" in w.get("issue", "")
            and "significant regions" not in w.get("issue", "")
            and w["issue"].split(" ")[0] in significant_chroms
        )
    ]


def _print_summary(warnings: List[dict]) -> None:
    if not warnings:
        print("   ✓ No QC issues detected")
        return
    print(f"   ⚠ {len(warnings)} QC warning(s):")
    for warning in warnings:
        marker = "[MAJOR]" if warning["severity"] == "major" else "[minor]"
        print(f"   {marker} {warning['issue']}")
        print(f"      {warning['details']}")
        print(f"      -> {warning['recommendation']}")


def _check_sex_chrom_all(frame: pd.DataFrame) -> List[dict]:
    findings: List[dict] = []
    for chrom in ("chrY", "chrX"):
        subset = frame[frame["chrom"] == chrom]
        if len(subset) < 3:
            continue

        a = int((subset["direction"] == "A_enriched").sum())
        b = int((subset["direction"] == "B_enriched").sum())
        if max(a, b) / len(subset) < 0.8:
            continue
        if "mean_count_a" not in subset.columns or "mean_count_b" not in subset.columns:
            continue

        depleted_col = "mean_count_b" if a >= b else "mean_count_a"
        enriched_col = "mean_count_a" if a >= b else "mean_count_b"
        dep_mean = subset[depleted_col].mean()
        enr_mean = subset[enriched_col].mean()
        if enr_mean <= 0 or dep_mean / enr_mean >= 0.1:
            continue

        n_sig = int(subset.get("significant", pd.Series(False, index=subset.index)).sum())
        enriched_group = "A" if a >= b else "B"
        findings.append({
            "severity": "major",
            "issue": f"{chrom} sex chromosome confound",
            "details": (
                f"{len(subset)} {chrom} regions point overwhelmingly toward group {enriched_group} "
                f"(mean count {enr_mean:.1f} vs {dep_mean:.1f} in the depleted group); "
                f"{n_sig} reach FDR<0.05. This pattern is the signature of a sex mismatch between "
                "groups (e.g. male vs female lines), not genuine differential binding."
            ),
            "recommendation": (
                f"Exclude {chrom} regions from biological interpretation; report them as probable "
                "sex-chromosome artifacts."
            ),
        })
    return findings


def _check_sex_chrom_significant(frame: pd.DataFrame) -> List[dict]:
    findings: List[dict] = []
    if "significant" not in frame.columns:
        return findings

    for chrom in ("chrY", "chrX"):
        sig_subset = frame[(frame["chrom"] == chrom) & frame["significant"]]
        if len(sig_subset) < 2:
            continue

        a = int((sig_subset["direction"] == "A_enriched").sum())
        b = int((sig_subset["direction"] == "B_enriched").sum())
        if max(a, b) / len(sig_subset) < 0.67:
            continue
        if "mean_count_a" not in sig_subset.columns or "mean_count_b" not in sig_subset.columns:
            continue

        depleted_col = "mean_count_b" if a >= b else "mean_count_a"
        enriched_col = "mean_count_a" if a >= b else "mean_count_b"
        dominant_direction = "A_enriched" if a >= b else "B_enriched"
        dominant = sig_subset[sig_subset["direction"] == dominant_direction]
        dep_mean = dominant[depleted_col].mean()
        enr_mean = dominant[enriched_col].mean()
        if enr_mean <= 0 or dep_mean / enr_mean >= 0.1:
            continue

        enriched_group = "A" if a >= b else "B"
        findings.append({
            "severity": "major",
            "issue": f"{chrom} sex chromosome confound (significant regions)",
            "details": (
                f"{max(a, b)}/{len(sig_subset)} significant {chrom} hits favor group {enriched_group} "
                f"with mean counts {enr_mean:.1f} vs {dep_mean:.1f} in the other group. "
                "These top results are almost certainly sex-chromosome artifacts."
            ),
            "recommendation": (
                f"Filter {chrom} regions before drawing biological conclusions."
            ),
        })
    return findings


def _check_nonstandard_contigs(frame: pd.DataFrame, genome: Optional[str]) -> List[dict]:
    standard = get_standard_chroms(genome)
    mask = ~frame["chrom"].isin(standard)
    if not mask.any():
        return []

    n = int(mask.sum())
    pct = n / len(frame) * 100
    contigs = int(frame.loc[mask, "chrom"].nunique())
    n_sig = int((mask & frame.get("significant", False)).sum())
    return [{
        "severity": "minor",
        "issue": "Non-standard contigs present",
        "details": (
            f"{n} of {len(frame)} regions ({pct:.1f}%) sit on {contigs} non-standard contigs "
            f"(random scaffolds, alt loci); {n_sig} reach FDR<0.05. These are commonly mapping "
            "artifacts."
        ),
        "recommendation": (
            "Default QC drops them. Pre-QC rows live in `diff_regions_unfiltered.csv`."
        ),
    }]


def _check_small_regions(frame: pd.DataFrame) -> List[dict]:
    if "region_size" not in frame.columns:
        return []
    mask = frame["region_size"] < 10
    if not mask.any():
        return []
    n = int(mask.sum())
    n_sig = int((mask & frame.get("significant", False)).sum())
    smallest = int(frame.loc[mask, "region_size"].min())
    return [{
        "severity": "minor",
        "issue": "Implausibly small regions",
        "details": (
            f"{n} regions are <10 bp (smallest: {smallest} bp), {n_sig} of them significant. "
            "Below this size most calls are edge artifacts; TF motifs are ~20 bp and histone "
            "marks 200 bp+."
        ),
        "recommendation": "Default QC drops <10 bp; original data is in `diff_regions_unfiltered.csv`.",
    }]


def _check_chrM(frame: pd.DataFrame, analysis_type: str) -> List[dict]:
    if analysis_type == "dmr":
        return []
    subset = frame[frame["chrom"] == "chrM"]
    if len(subset) < 3:
        return []
    n_sig = int(subset.get("significant", pd.Series(False, index=subset.index)).sum())
    if n_sig == 0:
        return []
    total_sig = int(frame.get("significant", pd.Series(False, index=frame.index)).sum())
    pct = n_sig / total_sig * 100 if total_sig else 0
    severity = "major" if n_sig >= 5 or pct >= 10 else "minor"
    return [{
        "severity": severity,
        "issue": f"chrM enrichment ({n_sig} significant regions)",
        "details": (
            f"{n_sig} of {total_sig} significant regions ({pct:.0f}%) are on chrM. "
            "Mitochondrial DNA lacks histones, so significant chrM peaks in ChIP/ATAC/DNase data "
            "almost always reflect non-specific antibody binding or library-prep contamination."
        ),
        "recommendation": (
            "Exclude chrM regions from biological interpretation; consider removing chrM entirely."
        ),
    }]


def _check_genomic_clusters(
    frame: pd.DataFrame,
    min_regions: int = 10,
    window_bp: int = 2_000_000,
) -> List[dict]:
    findings: List[dict] = []
    if "significant" not in frame.columns:
        return findings

    sig = frame[frame["significant"]].copy()
    if len(sig) < min_regions:
        return findings

    for chrom, chrom_sig in sig.groupby("chrom"):
        chrom_sig = chrom_sig.sort_values("chromStart").reset_index(drop=True)
        starts = chrom_sig["chromStart"].astype(int).to_numpy()
        ends = chrom_sig["chromEnd"].astype(int).to_numpy()
        if len(starts) < min_regions:
            continue

        best = (0, 0, 0)  # (count, start, end)
        j = 0
        for i in range(len(starts)):
            while j < len(starts) and starts[j] - starts[i] <= window_bp:
                j += 1
            count = j - i
            if count > best[0]:
                best = (count, int(starts[i]), int(ends[j - 1]))

        if best[0] < min_regions:
            continue

        count, win_start, win_end = best
        window_size = win_end - win_start
        cluster = chrom_sig[
            (chrom_sig["chromStart"] >= win_start)
            & (chrom_sig["chromStart"] <= win_start + window_bp)
        ]
        n_a = int((cluster["direction"] == "A_enriched").sum())
        n_b = int((cluster["direction"] == "B_enriched").sum())
        dominant = "A" if n_a >= n_b else "B"

        depleted_col = "mean_count_b" if dominant == "A" else "mean_count_a"
        zero_note = ""
        if depleted_col in cluster.columns:
            mean_depleted = cluster[depleted_col].mean()
            if mean_depleted < 1.0:
                zero_note = (
                    f" The depleted group's mean count is near zero ({mean_depleted:.2f}), "
                    "consistent with a copy-number deletion or amplification rather than "
                    "differential binding."
                )

        findings.append({
            "severity": "minor",
            "issue": f"{chrom} genomic cluster ({count} significant regions in {window_size / 1e6:.1f} Mb)",
            "details": (
                f"{count} significant regions cluster within {chrom}:{win_start:,}-{win_end:,} "
                f"({window_size / 1e6:.1f} Mb). {max(n_a, n_b) / count * 100:.0f}% favor group {dominant}."
                f"{zero_note}"
            ),
            "recommendation": (
                f"Inspect {chrom}:{win_start:,}-{win_end:,} in a genome browser. "
                "If a CNV explains the cluster, exclude it from interpretation."
            ),
            "flagged_regions": [
                (row["chrom"], int(row["chromStart"]), int(row["chromEnd"]))
                for _, row in cluster.iterrows()
            ],
        })
    return findings


def _check_ma_asymmetry(frame: pd.DataFrame) -> List[dict]:
    needed = {"significant", "mean_count_a", "mean_count_b", "direction"}
    if not needed.issubset(frame.columns):
        return []

    sig = frame[frame["significant"]].copy()
    if len(sig) < 5:
        return []

    def log2_mean(row):
        return np.log2(max((row["mean_count_a"] + row["mean_count_b"]) / 2, 0.1))

    sig_a = sig[sig["direction"] == "A_enriched"]
    sig_b = sig[sig["direction"] == "B_enriched"]
    if len(sig_a) < 3 or len(sig_b) < 3:
        return []

    median_a = sig_a.apply(log2_mean, axis=1).median()
    median_b = sig_b.apply(log2_mean, axis=1).median()
    diff = abs(median_a - median_b)
    if diff < 1.0:
        return []

    if median_a < median_b:
        low_group, low_median, high_median, low_n = "A", median_a, median_b, len(sig_a)
    else:
        low_group, low_median, high_median, low_n = "B", median_b, median_a, len(sig_b)

    return [{
        "severity": "minor",
        "issue": f"MA plot asymmetry (group {low_group} hits at low counts)",
        "details": (
            f"Significant hits enriched in group {low_group} have median log2 mean count "
            f"{low_median:.1f} vs {high_median:.1f} on the other side (difference: {diff:.1f}). "
            f"That direction's {low_n} hits cluster at low counts — typical of noise-driven or "
            "normalization-driven differential signal."
        ),
        "recommendation": (
            f"Inspect the MA plot. Focus interpretation on higher-count significant regions when "
            f"comparing group {low_group} hits to group {('B' if low_group == 'A' else 'A')}."
        ),
    }]


def _check_region_size_distribution(frame: pd.DataFrame, analysis_type: str) -> List[dict]:
    if analysis_type == "dmr":
        return []
    if "region_size" not in frame.columns or "significant" not in frame.columns:
        return []

    sig = frame[frame["significant"]]
    if len(sig) < 5:
        return []

    median_size = float(sig["region_size"].median())
    if median_size >= 100:
        return []

    return [{
        "severity": "minor",
        "issue": f"Unusually small significant regions (median {median_size:.0f} bp)",
        "details": (
            f"Median significant region size is {median_size:.0f} bp — below the typical 100–500 bp "
            "TF footprint and 200 bp+ histone marks. Peak fragmentation may be inflating the region "
            "count."
        ),
        "recommendation": (
            "Browse a few representative regions in IGV to check whether sub-peaks sit inside "
            "larger consensus peaks."
        ),
    }]


def _check_sample_size(n_a: int, n_b: int) -> List[dict]:
    min_n = min(n_a, n_b)
    if min_n <= 3:
        return [{
            "severity": "major",
            "issue": f"Very low sample size (n={n_a} vs n={n_b})",
            "details": (
                f"edgeR's dispersion estimate is unstable with n≤3 per group. FDR values should be "
                "treated as approximate rankings rather than calibrated false-discovery rates."
            ),
            "recommendation": (
                "Trust large |logFC| and very small q-values; validate downstream findings with "
                "independent replicates or orthogonal assays."
            ),
        }]
    if min_n < 5:
        return [{
            "severity": "minor",
            "issue": f"Low sample size (n={n_a} vs n={n_b})",
            "details": (
                "With n<5 per group, edgeR FDR estimates can be noisy for moderate effect sizes."
            ),
            "recommendation": "Lean on large effects; validate marginal hits.",
        }]
    return []


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", required=True, help="Parsed diff regions CSV")
    parser.add_argument("--n-a", type=int, required=True)
    parser.add_argument("--n-b", type=int, required=True)
    parser.add_argument("--analysis-type", default="diffbind", choices=["diffbind", "dmr"])
    parser.add_argument("--genome", default=None)
    args = parser.parse_args(argv)

    frame = pd.read_csv(args.csv)
    warnings = run_qc_checks(frame, args.n_a, args.n_b, args.analysis_type, args.genome)
    print(f"✓ {len(warnings)} QC warning(s) found")
    return 0


if __name__ == "__main__":
    sys.exit(main())
