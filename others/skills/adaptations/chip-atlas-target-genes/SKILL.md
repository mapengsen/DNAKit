---
id: "skill_47d8c3659b0c4e74b3fce7061d7a8bf7"
name: chip-atlas-target-genes
display-name: "ChIP-Atlas Target Genes"
description: "For a single transcription factor or chromatin regulator, download ChIP-Atlas's pre-computed per-gene MACS2 binding scores across every public ChIP-seq experiment, rank the targets, optionally subset to specific cell types, and integrate the pre-embedded STRING interaction scores. Use when you have one TF in mind and want a ranked list of its candidate target genes."
category: "epigenomics"
keywords: "ChIP-Atlas, target genes, MACS2, transcription factor, ChIP-seq, regulator targets, STRING, epigenomics"
version: "1.0"
last-updated: "2026-05-20"
---

## DNAKit execution

Use `dnakit skills doctor chip-atlas-target-genes` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# ChIP-Atlas Target Genes

Ask: "*what genes does this transcription factor bind?*" — and let the
ChIP-Atlas data server hand you a pre-computed wide-format TSV with the
MACS2 binding score for every gene, across every public ChIP-seq
experiment for the chosen factor. STRING protein-interaction scores are
already embedded in the same file (no separate STRING API call).

## When to Use This Skill

Use this skill when you have:

- **One specific TF (or chromatin regulator)** and want a ranked list of
  its candidate targets.
- **A cell-type-specific question** — re-subset the per-experiment scores
  to your cell line / tissue of interest.
- **A starting point for regulatory network construction** — top targets
  become the nodes one hop away from the factor.
- **STRING evidence integration** — a column of pre-embedded STRING
  scores is included; no separate STRING query is needed.

Don't use this skill for:

- ❌ **"Which TFs bind my gene list?"** — that's `chip-atlas-peak-enrichment`.
- ❌ **Comparing binding between two conditions** — use
  `chip-atlas-diff-analysis`.
- ❌ **Histone modifications** — the target-gene endpoint only serves
  non-histone antigens. Use the peak-enrichment workflow for histone
  marks.
- ❌ **Genomes outside the supported ten** (hg38, hg19, mm10, mm9, rn6,
  dm6, dm3, ce11, ce10, sacCer3).
- ❌ **Offline analysis** — requires HTTPS to the ChIP-Atlas data server.

**Core idea.** A static TSV at
`https://chip-atlas.dbcls.jp/data/{genome}/target/{Protein}.{distance}.tsv`
contains one row per gene and columns: `Target_genes`,
`{Protein}|Average` (mean MACS2 score across experiments), per-experiment
columns `{SRX_ID}|{cell_type}`, and `STRING`. The skill downloads the
TSV, parses it into a gene-level summary plus a wide-format
experiment matrix, optionally re-subsets by cell type (recomputing the
average), and writes ranked CSVs, plots, and a markdown report.

## Quick Start

```python
from scripts.load_example_query import load_example_query
from scripts.run_target_genes_workflow import run_target_genes_workflow
from scripts.generate_all_plots import generate_all_plots
from scripts.export_all import export_all

query = load_example_query("tp53")

results = run_target_genes_workflow(
    protein=query["protein"],
    genome=query["genome"],
    distance=query["distance"],
    top_n=500,
    output_dir="target_genes_results",
)
generate_all_plots(results, output_dir="target_genes_results", top_n=25)
export_all(results, output_dir="target_genes_results")
```

Runtime for `tp53`: ~10–30 seconds (the TSV is ~13 MB and downloads
quickly; parsing dominates).

## Installation

| Package | Version | License | Commercial use | Install |
|---------|---------|---------|---------------|---------|
| pandas | ≥1.3 | BSD-3-Clause | Permitted | `pip install pandas` |
| numpy | ≥1.20 | BSD-3-Clause | Permitted | `pip install numpy` |
| requests | ≥2.25 | Apache-2.0 | Permitted | `pip install requests` |
| matplotlib | ≥3.4 | PSF-based | Permitted | `pip install matplotlib` |
| seaborn | ≥0.12 | BSD-3-Clause | Permitted | `pip install seaborn` |
| plotnine | ≥0.12 | MIT | Permitted | `pip install plotnine` |
| plotnine-prism | ≥0.1 | MIT | Permitted | `pip install plotnine-prism` |

One-liner:

```bash
pip install "pandas>=1.3" "numpy>=1.20" "requests>=2.25" "matplotlib>=3.4" "seaborn>=0.12" "plotnine>=0.12" "plotnine-prism>=0.1"
```

System requirements: outbound HTTPS to `https://chip-atlas.dbcls.jp` to
fetch the TSV. Larger TFs (TP53, CTCF, MYC) can be 10–20 MB and require
~2 minutes of wall time on a slow connection.

## Inputs

**Required**

- `protein` — case-sensitive gene symbol of the factor (e.g. `"TP53"`,
  `"CTCF"`, `"MYC"`). Histone marks (`"H3K4me3"`, etc.) return 404 —
  they are not available via this endpoint.

**Optional**

- `genome` — one of `hg38` (default), `hg19`, `mm10`, `mm9`, `rn6`,
  `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3`.
- `distance` — TSS window in kilobases: `1`, `5` (default), or `10`.
  Larger windows include more distal binding events.
- `top_n` — keep only the top N genes by `avg_score` (default 500).
- `min_score` — drop genes whose mean MACS2 score is below this
  threshold (default 0, i.e. no filter).
- `cell_types` — list of cell-type names to subset to. **Re-computes**
  `avg_score`, `binding_rate`, `max_score`, and `num_bound` using only
  the matching experiments.
- `min_string_score` — drop genes with STRING score below this value.
- `min_binding_rate` — drop genes bound in fewer than this fraction of
  experiments (0.0–1.0).

## Outputs

All written to `output_dir`:

| File | Format | Purpose / consumer |
|------|--------|--------------------|
| `analysis_object.pkl` | Python pickle | Full results dict for downstream use |
| `target_genes_all.csv` | CSV | All retained target genes (gene, avg_score, max_score, binding_rate, num_bound, num_experiments, string_score, colocated_group) |
| `target_genes_top50.csv` | CSV | Top 50 by `avg_score` |
| `target_genes_with_string.csv` | CSV | Genes with `string_score > 0` (only when present) |
| `experiment_scores_top50.csv` | CSV | Wide-format per-experiment scores for the top 50 |
| `target_genes_top_targets.{png,svg}` | image | Horizontal bar plot of top targets, coloured by binding rate |
| `target_genes_score_distribution.{png,svg}` | image | log10-scaled distribution of mean scores |
| `target_genes_heatmap.{png,svg}` | image | Clustered heatmap (top genes × top experiments) |
| `target_genes_string_vs_binding.{png,svg}` | image | STRING vs binding scatter (falls back to binding-rate scatter when STRING is absent) |
| `summary_report.md` | markdown | Human-readable analysis report with experiment-composition table and caveats |

The pickle contains `target_genes`, `experiment_data`, `cell_types`,
`protein`, `parameters`, and `metadata` (including `colocated_info`).

## Clarification Questions

1. **⚠️ CRITICAL: ASK THIS FIRST.** Which protein/TF are we querying?
   - Symbol is **case-sensitive** (`TP53`, not `tp53`).
   - Histone marks (H3K4me3, H3K27ac, etc.) are unavailable here.
   - If unavailable, offer a bundled example: `tp53` (large dataset,
     ~395 experiments, ~16k genes), `e2f1` (cell cycle), or `myc`
     (moderate).
2. Which **genome assembly**? Default `hg38`.
3. Which **TSS distance window**? `±5 kb` (default) for broad regulatory
   scope, `±1 kb` for proximal-only, `±10 kb` for distal capture.
4. **Cell-type filter?** Leave at "all" by default; ask the user if they
   want to narrow to a specific cell line (e.g. `["MCF-7"]` or
   `["K-562", "GM12878"]`).
5. **Top-N cap?** 500 is plenty for most downstream needs; rarely is
   keeping >1000 rows useful.

When the user picks an example dataset, do NOT re-ask about genome,
distance, or filters — they are pre-set.

## Standard Workflow

1. **Load query** — `load_example_query(name)` or
   `load_user_query(protein, genome, distance)`. The user-query helper
   validates the genome and distance values and reminds you about
   case-sensitivity.
2. **Run workflow** with `run_target_genes_workflow(...)`. It validates
   antigen availability (HTTP HEAD), downloads the TSV, parses it into
   summary + experiment DataFrames, applies filters, and detects
   co-located genes (multiple genes sharing identical score vectors
   because they sit within the same TSS window).
3. **Plots** with `generate_all_plots(results, output_dir, top_n=25)`.
   Four PNG/SVG plots: top targets, score distribution, heatmap,
   STRING vs binding scatter.
4. **Export** with `export_all(results, output_dir)` — writes the
   pickle, CSVs, and markdown summary.

## When Scripts Fail

1. **Fix & retry (~90%).** Wrong case in the protein name? Connection
   wobble? Update and rerun.
2. **Edit the script (~5%).** Tweak filters or output names; keep the
   contract intact so `export_all.py` can still consume the dict.
3. **Read for inspiration (~4%).** Adapt the approach but match the
   contract (column names, CSV file names) used by downstream scripts.
4. **Write from scratch (~1%).** Explain why.

## Decision Points

- **Average scores include zeros.** A gene with score 1000 in 5 of 395
  experiments has an average around 12. The summary report exposes
  `max_score` and `binding_rate` so you can spot strong cell-type-
  specific binding. Lead the user-facing narrative with the average
  ranking but cite `max_score` for outlier cell-specific hits.
- **STRING ≠ ground truth.** STRING coverage is sparse and biased toward
  well-studied interactions. A STRING score of 0 does **not** mean "not
  a target" — even canonical TP53 targets like BBC3/PUMA score 0.
- **Co-located genes.** Multiple genes that share identical score vectors
  fall within the same TSS window (gene pairs at the same locus). The
  `colocated_group` column flags them. For pathway enrichment, consider
  collapsing groups to avoid double-counting loci.
- **Cell-type re-subsetting.** Setting `cell_types=["MCF-7"]` recomputes
  averages, so the rankings reflect MCF-7 binding only. This is the
  cleanest way to remove well-studied-cell-line bias.
- **Threshold heuristic.** Score ≥ 500 ≈ MACS2 Q ≤ 1e-50 per individual
  experiment (the average is diluted by zeros — see caveat above).

## Common Issues

| Symptom | Likely cause | What to do |
|---------|--------------|------------|
| HTTP 404 for the protein | Wrong case (`tp53` vs `TP53`), histone mark requested, or factor without ChIP-Atlas data | The validator prints the exact URL and suggested fixes. See `references/target_genes_data_format.md`. |
| Download timeout | Slow connection or large TF (TP53 is ~13 MB) | Increase `request_timeout` in `download_target_genes` or retry. |
| Memory error on parse | Very wide TSV (hundreds of experiments) | Set `top_n` smaller; subset with `cell_types`. |
| All STRING scores are zero | Factor missing from ChIP-Atlas's pre-embedded STRING table | Normal for less-studied TFs. Binding scores remain valid; the workflow falls back to a binding-rate scatter for plot 4. |
| Empty result after filters | Combination too strict | Lower `min_score`, remove `cell_types` filter, raise `top_n`. |
| Many co-located gene groups | Multiple genes in the same TSS window | Expected at `±5 kb`; consider `±1 kb` for proximal-only or collapse groups via the `colocated_group` column before enrichment. |
| `Could not find Average column` | Server returned an unexpected schema | Confirm the URL works in a browser; report the schema mismatch. |
| Heatmap clustering fails | Too few genes / experiments retained | Lower the per-plot `n_genes` / `n_experiments` inputs. |
| SVG export warnings | Optional `cairosvg` missing | PNGs always produced. |
| Re-subsetting yields very few experiments | Cell type spelled differently in headers | The summary lists available cell types; check spelling case-insensitively. |

## Best Practices

- **Quote ranks from the data, not memory.** When listing top targets,
  copy the gene names from `target_genes_top50.csv` verbatim — do not
  re-insert "famous" TP53 targets that the data did not rank highly.
- **Mark biological annotations as background.** "CDKN1A — known cell-
  cycle arrest effector — ranks #1" is fine. Don't let known biology
  blur into a data-derived claim.
- **Always report the experiment composition.** The summary table shows
  top cell types and the share of experiments they contribute. Mention
  that average scores are biased toward heavily represented lines.
- **Mention co-located groups explicitly.** The independent-loci count
  in the summary is the right number when discussing "how many targets"
  the factor has.
- **Cite ChIP-Atlas.** Reference Zou et al. 2024 (and Oki et al. 2018)
  when summarising.

## Suggested Next Steps

1. **Run peak enrichment** with the top targets to identify co-regulating
   factors (`chip-atlas-peak-enrichment`).
2. **Cell-type-specific rerun** — re-run with `cell_types=[...]` matching
   the user's system to remove cell-line bias.
3. **Pathway / GO enrichment** on the top target set to interpret the
   regulatory program.
4. **Cross-reference with expression** — pair the targets with RNA-seq
   data to flag the ones that change at the transcript level.
5. **Validate canonical targets are present** — sanity-check that
   well-characterised targets show up in the top N at the expected rank.

## Related Skills

- **chip-atlas-peak-enrichment** — the reverse direction: given a gene
  list, which factors are enriched?
- **chip-atlas-diff-analysis** — differential binding between two
  experiment groups (e.g. before/after treatment).
- **grn-pyscenic** — gene regulatory networks from single-cell data;
  good downstream consumer of a target list.

## References

Reference files under `references/`:

- `target_genes_data_format.md` — the URL pattern, the TSV column
  conventions, and the rules for HEAD-checking antigen availability.
- `macs2_binding_scores.md` — what the binding scores mean, how they
  are computed, and why averages dilute under zero-binding experiments.
- `string_scores.md` — interpretation of the embedded STRING column
  and the "STRING=0 does not mean not a target" caveat.

Primary literature:

- Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore
  chromosome architecture together with large-scale regulome data.
  *Nucleic Acids Research* 52(D1):D45–D53.
  [doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)
- Zou Z. *et al.* (2022). ChIP-Atlas 2021 update. *Nucleic Acids Research*
  50(D1):D67–D72.
  [doi:10.1093/nar/gkab933](https://doi.org/10.1093/nar/gkab933)
- Oki S. *et al.* (2018). ChIP-Atlas: a data-mining suite powered by full
  integration of public ChIP-seq data. *EMBO Reports* 19(12):e46255.
  [doi:10.15252/embr.201846255](https://doi.org/10.15252/embr.201846255)
- Szklarczyk D. *et al.* (2023). The STRING database in 2023.
  *Nucleic Acids Research* 51(D1):D638–D646.
  [doi:10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000)
- Zhang Y. *et al.* (2008). Model-based analysis of ChIP-seq (MACS).
  *Genome Biology* 9(9):R137.
  [doi:10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)

ChIP-Atlas portal: <https://chip-atlas.org>
