# Output File Formats

The ChIP-Atlas Diff Analysis pipeline produces three main artefacts that
end up in `raw_results/`: a custom 8-column BED, a standard BED9-with-GFF3
companion for IGV, and an IGV session XML. The skill's parser adds derived
columns and saves CSV copies for downstream consumption.

## Main result — `wabi_result.bed`

A tab-separated file with **no header** and the following eight columns:

| # | Column | Type | Meaning |
|---|--------|------|---------|
| 1 | `chrom` | str | UCSC-style chromosome (e.g. `chr1`, `chrX`) |
| 2 | `chromStart` | int | 0-based start coordinate |
| 3 | `chromEnd` | int | half-open end coordinate |
| 4 | `counts_a` | str | Comma-separated TMM-normalized counts, one per group-A experiment |
| 5 | `counts_b` | str | Comma-separated TMM-normalized counts, one per group-B experiment |
| 6 | `logFC` | float | log2 fold change (positive = enriched in A, negative = enriched in B) |
| 7 | `pvalue` | float | edgeR quasi-likelihood F-test p-value (or metilene region-level p in DMR mode) |
| 8 | `qvalue` | float | Benjamini-Hochberg corrected q-value (FDR) |

Example row:

```
chr6	36651800	36652050	48.2,31.4,36.0	3.6,3.9,4.2	3.17	2.05e-11	7.4e-09
```

Interpretation: chr6:36,651,800–36,652,050 has mean TMM-normalized counts
of ~38.5 in group A vs ~3.9 in group B, log2FC = 3.17, FDR = 7.4 × 10⁻⁹.

### Derived columns added by `parse_diff_bed.py`

| Column | Type | How it's computed |
|--------|------|-------------------|
| `region_size` | int | `chromEnd − chromStart` |
| `direction` | str | `A_enriched` if `logFC > 0`, else `B_enriched` |
| `significant` | bool | `qvalue < 0.05` |
| `mean_count_a` | float | mean of `counts_a` values |
| `mean_count_b` | float | mean of `counts_b` values |
| `score` | int | `min(1000, -log10(qvalue) * 100)`; useful when re-exporting to BED |

The parser sorts ascending by qvalue, breaking ties by `|logFC|` then
position, so the first rows are always the strongest hits.

## IGV companion — `wabi_result.igv.bed`

This file is **BED9** with URL-encoded GFF3 attributes in the `name`
field. It is meant for visualization in IGV (and is loaded automatically
by `wabi_result.xml`).

| # | Column | Notes |
|---|--------|------|
| 1 | `chrom` | UCSC-style |
| 2 | `chromStart` | 0-based |
| 3 | `chromEnd` | half-open |
| 4 | `name` | URL-encoded GFF3 — typically `LogFC=…;P-value=…;Q-value=…` |
| 5 | `score` | 0–1000, derived from −log10(Q) capped at 1000 |
| 6 | `strand` | always `.` |
| 7 | `thickStart` | equals `chromStart` |
| 8 | `thickEnd` | equals `chromEnd` |
| 9 | `itemRgb` | color = direction (see below) |

### Color encoding

| Color | RGB | Direction |
|-------|-----|-----------|
| Orange | `222,131,68` | Enriched in group A (positive logFC) |
| Blue | `106,153,208` | Enriched in group B (negative logFC) |

The skill's parser can fall back to this file if the 8-column main BED
is missing — it decodes the GFF3 attributes back into `logFC`, `pvalue`,
and `qvalue` columns.

## IGV session — `wabi_result.xml`

A standard IGV session file. Load it via **File → Open Session…** in IGV
to get both the BigWig coverage tracks (for all submitted experiments)
and the differential BED, color-coded by direction. This is usually the
fastest way to eyeball top hits.

## Region-size heuristics

Different feature classes produce characteristic region widths. Use this
as a sanity check when the QC report flags unusually small regions.

| Width | Typical biology |
|-------|----------------|
| 100–500 bp | Transcription factor binding sites (motif + flank) |
| 500–2,000 bp | Narrow histone marks (H3K4me3), enhancer elements |
| 2,000–10,000 bp | Broad histone marks (H3K27me3, H3K36me3) |
| > 10,000 bp | Super-enhancers, large chromatin domains |
| < 50 bp | Often peak-edge artifacts; the parser drops <10 bp by default |

## CSV exports produced by the skill

`export_all.py` writes filtered slices of the parsed table:

- `diff_regions_all.csv` — default QC-filtered rows (non-standard contigs
  and regions <10 bp removed).
- `diff_regions_significant.csv` — `qvalue < 0.05` subset.
- `diff_regions_top50.csv` — top 50 by Q-value, with a `nearest_gene`
  column when UCSC annotation succeeded.
- `diff_regions_unfiltered.csv` — original pre-QC rows; written only when
  QC actually removed something.

The CSVs preserve column order and types from the parsed DataFrame, so
they round-trip through `pandas.read_csv` cleanly.
