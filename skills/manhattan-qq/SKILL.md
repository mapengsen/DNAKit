---
name: manhattan-qq
description: Produce a Manhattan plot, QQ plot, and a QC summary (genomic inflation
  λ_GC, count of genome-wide-significant and suggestive hits, list of independent
  lead SNPs) from a GWAS summary statistics file. Auto-detects common column-name
  conventions (`SNP/rsid/variant_id`, `CHR/chromosome/hm_chrom`, `BP/position/base_pair_location`,
  `P/p_value`). Subsamples low-significance variants to keep plots fast (~seconds)
  on full HM3-scale or 1KG-scale sumstats. Writes a sidecar `.manhattan_qq.json` manifest
  citing the input dataset's OKG node when `$OKG_REPO` is set. Use when a user asks
  for a Manhattan plot, QQ plot, genomic-inflation check, or a quick QC pass on a
  GWAS sumstats file.
license: MIT
metadata:
  author: stephendorn
  version: '0.1'
  compatibility: Requires Python 3.9+, numpy, pandas, scipy, matplotlib. No external
    services or data downloads.
---

## DNAKit execution

Read [DNAKIT.md](DNAKIT.md) first for actual runtime, resource paths and validation status.


## What this skill does

Given a GWAS sumstats TSV (gzip-OK), produces:

1. **Manhattan plot** — alternating-color autosome scatter, with horizontal lines at 5e-8 (genome-wide) and 1e-5 (suggestive). All variants with p < 1e-3 are plotted in full; the long, uninformative low-significance tail is downsampled (default 1%) so the plot renders in seconds even on 10M-variant sumstats.
2. **QQ plot** — observed vs. expected −log₁₀(p) under the null, with a log-spaced thinning of the bulk so the figure stays small. The plot title reports the genomic-inflation factor λ_GC (= median(χ²) / 0.4549, where χ² is recovered from the p-values).
3. **QC summary JSON** with:
   - `n_variants_total`
   - `n_genome_wide_5e8`, `n_suggestive_1e5`
   - `lambda_gc`
   - `n_independent_loci_<window>kb` — count of independent lead SNPs after pruning by ±window-kb on each chromosome
   - `top_loci`: chrom / pos / p_value / rsid (when available) for the top N (default 20) lead SNPs

4. **Sidecar manifest** `<out_prefix>.manhattan_qq.json` recording inputs, OKG citation (when resolvable), the QC summary, and SHA-256 of the input file.

The skill is **plot + QC only**. It does not filter, lift, or modify the sumstats. For build alignment use the `liftover` skill; for h² / partitioned-h² use the `ldsc` skill.

## Inputs you may need to elicit

- **`--in <path>`** (required) — sumstats TSV or TSV.gz. The skill auto-detects these column conventions:
  - **chrom**: `hm_chrom`, `chromosome`, `chrom`, `CHR`, `chr`
  - **pos**: `hm_pos`, `base_pair_location`, `position`, `pos`, `BP`, `bp`
  - **p**: `p_value`, `P`, `pval`, `p`, `PVAL`, `P_VAL`
  - **rsid** *(optional)*: `hm_rsid`, `rsid`, `SNP`, `variant_id`
  - Pass `--chrom-col / --pos-col / --p-col / --rsid-col` to override.
- **`--out-prefix <prefix>`** (required) — `<prefix>.manhattan.{png,pdf}`, `<prefix>.qq.{png,pdf}`, `<prefix>.qc_summary.json`, `<prefix>.manhattan_qq.json`.
- **`--gws-threshold <float>`** — default `5e-8`.
- **`--suggestive-threshold <float>`** — default `1e-5`.
- **`--lead-window-kb <int>`** — default `500`. Distance for the simple lead-SNP pruning.
- **`--top-loci <int>`** — default `20`. How many leads to record in `top_loci`.
- **`--include-x`** — include chrX (default: autosomes only). chrY / chrMT are always excluded.
- **`--subsample-frac <float>`** — fraction of p ≥ 1e-3 variants to draw in the Manhattan (default `0.01`). All p < 1e-3 are always drawn.
- **`--okg-repo <path>`** — override `$OKG_REPO`.
- **`--dataset-okg-id <node_id>`** — explicit OKG `dataset_metadata` node ID to record in the manifest (e.g. `dataset:gwas_catalog:gcst004131`). If omitted, the skill looks for a sibling `<in>.fetch.json` (from `gwas-fetch`) and lifts the dataset node from there.

## How to execute

1. Run: `python3 scripts/manhattan_qq.py --in <sumstats> --out-prefix <prefix>` with any overrides.
2. The script:
   a. Reads only the columns it needs (chrom, pos, p, optional rsid) — never loads the full sumstats into memory.
   b. Drops non-autosomal rows (unless `--include-x`), NA p-values, and pathological p=0 rows.
   c. Computes λ_GC from the recovered χ².
   d. Lead-prunes (greedy, sort by p ascending, mask ±window-kb around each accepted lead, repeat).
   e. Renders Manhattan + QQ as PNG (150 dpi) and PDF.
   f. Writes the QC summary JSON and the sidecar manifest.
3. Report back:
   - `n_variants_total`, `n_genome_wide_5e8`, `lambda_gc`, `n_independent_loci`
   - Top 5 loci (chrom:pos, rsid if available, p-value)
   - Paths to the four output files
4. **Suggested next step**: if `n_genome_wide_5e8 > 0` and λ_GC is reasonable (0.95–1.15), chain into the `ldsc` skill for h². If λ_GC > 1.20, flag possible stratification / inflation issues before downstream work.

## Refusal triggers

- Fewer than 100,000 valid variants after filtering — Manhattan / QQ are uninformative below that scale.
- λ_GC > 2 — heavy inflation; do not silently produce a plot that misleads. Flag, write the summary, but warn that the GWAS likely needs reprocessing (LDSC intercept will catch this too).
- No column matching any of the autodetect patterns AND no explicit `--*-col` override — the script lists the detected columns and asks the user to specify.

## Examples

```bash
# 1. Vanilla run on a GWAS Catalog harmonised file
python3 scripts/manhattan_qq.py \
    --in raw/28067908-GCST004131-EFO_0003767.h.tsv.gz \
    --out-prefix qc/ibd_gcst004131

# 2. Override column names for a non-standard sumstats schema
python3 scripts/manhattan_qq.py \
    --in mydata.tsv --out-prefix qc/mydata \
    --chrom-col CHR --pos-col POS --p-col PVAL --rsid-col rsID

# 3. Include chrX and tighten the lead-SNP window to 250 kb
python3 scripts/manhattan_qq.py \
    --in sumstats.tsv.gz --out-prefix qc/run --include-x --lead-window-kb 250

# 4. With OKG provenance — picks up the dataset_metadata node from gwas-fetch sidecar
OKG_REPO=~/Lab/KG/okg python3 scripts/manhattan_qq.py \
    --in raw/28067908-GCST004131-EFO_0003767.h.tsv.gz \
    --out-prefix qc/ibd_gcst004131
```

## Sidecar manifest schema

`<out_prefix>.manhattan_qq.json`:

```json
{
  "skill": "manhattan-qq",
  "version": "0.1",
  "input": "raw/28067908-GCST004131-EFO_0003767.h.tsv.gz",
  "input_sha256": "...",
  "outputs": {
    "manhattan_png": "qc/ibd_gcst004131.manhattan.png",
    "manhattan_pdf": "qc/ibd_gcst004131.manhattan.pdf",
    "qq_png": "qc/ibd_gcst004131.qq.png",
    "qq_pdf": "qc/ibd_gcst004131.qq.pdf",
    "qc_summary": "qc/ibd_gcst004131.qc_summary.json"
  },
  "thresholds": {"gws": 5e-8, "suggestive": 1e-5},
  "lead_window_kb": 500,
  "include_x": false,
  "qc": {
    "n_variants_total": 9633631,
    "n_genome_wide_5e8": 12345,
    "n_suggestive_1e5": 23456,
    "lambda_gc": 1.142,
    "n_independent_loci_500kb": 89,
    "top_loci": [{"chrom": 16, "pos": 50745926, "rsid": "rs2066845", "p_value": 1.2e-180}, ...]
  },
  "okg_node_ids": {"dataset": "dataset:gwas_catalog:gcst004131"},
  "captured_at": "2026-05-21T..."
}
```

## Notes and edge cases

- **No method/software OKG nodes.** Manhattan/QQ is not a method; matplotlib is generic. The manifest only cites the input dataset (when known). Don't fabricate a "method:manhattan_plot" node.
- **λ_GC** assumes a 1-df χ² distribution. It will be biased for GWAS that report case/control z-stats with non-standard scaling; use LDSC's intercept (`ldsc h2`) as the more rigorous inflation diagnostic. This skill's λ_GC is a quick-look only.
- **Lead-SNP pruning** is a naïve greedy ±window-kb approach (NOT LD-based clumping). For real clumping use PLINK `--clump` or the `finemap region` step. This is fine for "how many independent loci roughly" but should not be cited as a final hit list.
- **Subsampling**: variants with p ≥ 1e-3 are subsampled (default 1%) for the Manhattan plot only — the QC summary numbers are computed on the full sumstats.
- **Memory**: only chrom/pos/p (+ optional rsid) are read; even 20M-variant files stay well under 1 GB resident.
