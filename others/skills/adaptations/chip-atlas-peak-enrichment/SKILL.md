---
id: "skill_cda3e9616b9a4f38a226358a40de96a1"
name: chip-atlas-peak-enrichment
display-name: "ChIP-Atlas Peak Enrichment"
description: "Submit a gene list to the ChIP-Atlas Enrichment Analysis API and recover the transcription factors, chromatin modifiers, and histone marks whose public ChIP-seq peaks cluster near those gene promoters. Use when you have a candidate gene set and want a ranked, FDR-controlled list of regulators across 430,000+ public ChIP-seq experiments."
category: "epigenomics"
keywords: "ChIP-Atlas, peak enrichment, transcription factors, ChIP-seq, Fisher exact, regulator discovery, gene list, epigenomics"
version: "1.0"
last-updated: "2026-05-20"
---

## DNAKit execution

Use `dnakit skills doctor chip-atlas-peak-enrichment` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# ChIP-Atlas Peak Enrichment

Ask: "*which factors bind near these genes?*" — and let the ChIP-Atlas
server search the union of public ChIP-seq, ATAC-seq, and histone-mark
experiments for you. The API performs a per-experiment Fisher's exact test
of overlap between your gene promoters and that experiment's MACS2 peaks,
then reports BH-corrected q-values, fold enrichment, and gene-level
overlap counts.

## When to Use This Skill

Use this skill when you have:

- **A differential expression gene list** and want to know which TFs likely
  drive it.
- **A pathway or signature** and want to validate whether known regulators
  show ChIP-seq enrichment.
- **A list of co-regulated candidates** and want to discover novel
  regulators (chromatin remodelers, histone modifications) without
  downloading ChIP-seq raw data.
- **Cell-type-specific questions** — the API allows filtering to one of
  ~20 cell-class categories (Blood, Neural, Breast, etc.).

Don't use this skill for:

- ❌ **Finding the targets of a single TF.** Use `chip-atlas-target-genes`,
  which goes the other direction.
- ❌ **Differential binding between two conditions.** Use
  `chip-atlas-diff-analysis`.
- ❌ **Single-gene lookups.** The web UI at chip-atlas.org is faster for one
  gene.
- ❌ **Offline / air-gapped analysis.** Requires HTTPS to the ChIP-Atlas
  WABI host and (optionally) Ensembl for coordinate verification.
- ❌ **Genomes outside the supported ten** (hg38, hg19, mm10, mm9, rn6,
  dm6, dm3, ce11, ce10, sacCer3).

**Core idea.** You hand the API a list of gene symbols. The server builds
TSS windows (default ±5 kb), counts how many of your genes have at least
one MACS2 peak from each public experiment, runs Fisher's exact against
the RefSeq background (~18,600 genes for hg38), and returns one row per
experiment with `experiment_id`, `antigen`, `cell_type`, fold enrichment,
p-value, and q-value (Benjamini-Hochberg across all tested experiments).

## Quick Start

```python
from scripts.load_example_data import load_example_data
from scripts.run_enrichment_workflow import run_enrichment_workflow
from scripts.generate_all_plots import generate_all_plots
from scripts.export_all import export_all

data = load_example_data("tp53_targets")

results = run_enrichment_workflow(
    gene_list=data["genes"],
    genome="hg38",
    antigen_class="TFs and others",
    cell_class="All cell types",
    threshold=50,
    output_dir="enrichment_results",
)
generate_all_plots(results, output_dir="enrichment_results", top_n=15)
export_all(results, output_dir="enrichment_results")
```

Runtime: ~1–3 minutes for the 5-gene example.

## Installation

| Package | Version | License | Commercial use | Install |
|---------|---------|---------|---------------|---------|
| pandas | ≥1.3 | BSD-3-Clause | Permitted | `pip install pandas` |
| numpy | ≥1.20 | BSD-3-Clause | Permitted | `pip install numpy` |
| requests | ≥2.25 | Apache-2.0 | Permitted | `pip install requests` |
| matplotlib | ≥3.4 | PSF-based | Permitted | `pip install matplotlib` |
| plotnine | ≥0.10 | MIT | Permitted | `pip install plotnine` |
| plotnine-prism | ≥0.3 | MIT | Permitted | `pip install plotnine-prism` |

One-liner:

```bash
pip install "pandas>=1.3" "numpy>=1.20" "requests>=2.25" "matplotlib>=3.4" "plotnine>=0.10" "plotnine-prism>=0.3"
```

System requirements: outbound HTTPS to `https://dtn1.ddbj.nig.ac.jp` (the
ChIP-Atlas WABI host) and, optionally, `https://rest.ensembl.org` for
gene-to-region coordinate lookups (not required for the enrichment result
itself).

## Inputs

**Required**

- `gene_list` — list of HGNC (human) / MGI (mouse) / equivalent gene
  symbols. Must contain at least 3 genes; 5–100 genes is the sweet spot.
  Sources: a Python list, a one-per-line text file, a CSV with a `gene`
  column, or a DESeq2/edgeR/limma results table loaded by
  `load_user_data`.

**Optional**

- `genome` — one of `hg38` (default), `hg19`, `mm10`, `mm9`, `rn6`, `dm6`,
  `dm3`, `ce11`, `ce10`, `sacCer3`.
- `antigen_class` — `"TFs and others"` (default), `"Histone"`,
  `"RNA polymerase"`, `"ATAC-Seq"`, `"DNase-seq"`, `"Bisulfite-Seq"`,
  `"Input control"`, or `"Annotation tracks"`.
- `cell_class` — `"All cell types"` (default) or one of the cell-class
  categories listed in `references/api_format.md`.
- `threshold` — MACS2 stringency: `50` (default, ~p<1e-5), `100`
  (~p<1e-10), `200` (~p<1e-20), or `500` (~p<1e-50). Higher = fewer but
  more confident peaks.
- `distance_up`, `distance_down` — TSS window in bp (default 5000 each).
- `antigen_filter`, `cell_type_filter`, `min_fold_enrichment`,
  `max_qvalue` — client-side post-filters applied after the API returns.

## Outputs

All written to `output_dir`:

| File | Format | Purpose / consumer |
|------|--------|--------------------|
| `analysis_object.pkl` | Python pickle | Full results dict for downstream use |
| `enrichment_results_all.csv` | CSV | Every experiment the API returned (experiment_id, antigen, cell_type, fold_enrichment, p_value, q_value, significant, overlap_rate, regions_with_overlaps, total_regions, total_peaks) |
| `enrichment_results_significant.csv` | CSV | Subset with `q_value < 0.05` |
| `enrichment_results_top20.csv` | CSV | Top 20 by q-value with ≥2 gene overlaps required |
| `failed_genes.txt` | text | Genes that Ensembl could not resolve (omitted if empty) |
| `chipatlas_enrichment.{png,svg}` | image | 4-panel summary: top factors, p-value histogram, overlap vs fold-enrichment scatter, volcano plot |
| `summary_report.md` | markdown | Per-experiment and aggregated-by-factor tables, overlap-rate stats, caveats |

The pickle contains: `enrichment_results`, `input_genes`,
`input_regions`, `failed_genes`, `metadata`, `api_total_regions`, and the
`parameters` used.

## Clarification Questions

1. **⚠️ CRITICAL: ASK THIS FIRST.** What gene list are we analysing?
   - Format options: Python list, plain text (one gene per line), CSV
     with a `gene` column, or a DE results file from DESeq2/edgeR/limma.
   - Minimum 3 genes; 5–100 is recommended.
   - If unavailable, offer a bundled example: `tp53_targets` (5 genes,
     fast), `immune_response` (20 cytokine/TF genes), or
     `oncogenic_signature` (50-gene panel covering proliferation,
     apoptosis, and immune pathways).
2. Which **genome assembly**? Default `hg38`. Full set in the Inputs
   section.
3. Which **antigen class** matters? `"TFs and others"` is the default
   for "which TFs are enriched?"; pick `"Histone"` for histone marks,
   `"ATAC-Seq"` for accessibility.
4. Restrict to a **cell class** (Blood, Neural, Breast, …) or leave at
   `"All cell types"`? Restricting cuts runtime and helps when the user
   has a system-of-interest.
5. **Peak threshold?** `50` is the balanced default. Raise to `100` or
   `200` only when probing whether the signal survives a stricter cutoff.

## Standard Workflow

1. **Load genes** — `load_example_data(name)` or `load_user_data(path)`.
   Both functions surface skipped / failed gene entries so the user
   knows about retired symbols.
2. **Run enrichment** with `run_enrichment_workflow(...)`. The function
   submits to ChIP-Atlas, polls (10-minute default timeout, 15-second
   poll interval), retrieves the TSV result, parses it into the standard
   contract (raw p/q values, fold enrichment, overlap fractions), and
   optionally queries Ensembl to populate `input_regions`. Ensembl
   failures are non-fatal — they only affect the independent coordinate
   verification, not the enrichment numbers.
3. **Generate plots** with `generate_all_plots(results, output_dir,
   top_n=15)`. The 4-panel figure renders to both PNG and SVG; SVG falls
   back gracefully if the optional dependency is missing.
4. **Export** with `export_all(results, output_dir)`. CSVs, the pickle,
   the failed-gene text file (if any), and the markdown report all land
   together.

## When Scripts Fail

1. **Fix & retry (~90%).** Network blip? Wait and rerun. Missing package?
   Install and rerun.
2. **Edit the script (~5%).** Add a filter, tweak a default — keep the
   contract (column names, output filenames) intact for downstream
   compatibility.
3. **Read for inspiration (~4%).** Adapt the approach but follow the same
   contract so `export_all.py` can still consume the result dict.
4. **Write from scratch (~1%).** Last resort; explain why.

## Decision Points

- **Q-value, not fold enrichment, is the primary ranking.** Sentinel
  enrichment values (≥100,000×) happen when the RefSeq background has
  zero peaks; they always come with low overlap counts and non-significant
  q-values. Top-20 ranking requires `regions_with_overlaps >= 2` to
  exclude these.
- **Small gene sets (≤10 genes) are exploratory.** A single gene
  contributes ≥10% of the overlap rate. Lead any narrative with the
  exploratory framing; cite the caveat every time you discuss a 2–10×
  fold-enrichment factor.
- **Aggregated factor table vs per-experiment table.** Multiple
  experiments for the same TF are not independent biological signals.
  Cite the aggregated table when reporting "factor counts"; cite the
  per-experiment table for raw observations.
- **Use Median FE (Sig), not Median FE (All), to gauge enrichment
  strength.** A factor with 137 experiments but only 5 significant looks
  weakly enriched in the "all" median; the significant subset can still
  show strong fold enrichment.
- **Threshold sensitivity check.** Rerunning at threshold=100 (and
  comparing top factors) is the cheapest robustness check. Offer this
  proactively.

## Common Issues

| Symptom | Likely cause | What to do |
|---------|--------------|------------|
| HTTP 400 from submit | Empty `antigenClass` or `cellClass`, unrecognised genome | Default `cellClass` is `"All cell types"`. See `references/api_format.md`. |
| API timeout after 10 minutes | "All cell types" with many experiments | Increase `timeout` on `poll_job_status` or narrow `cell_class`. |
| Ensembl timeout / 0 mappings | Public Ensembl is rate-limited or down | The workflow retries once after 60s. If it still fails, enrichment is unaffected; report this to the user as "no independent coordinate verification." |
| Gene reported as failed | Retired HGNC symbol (e.g. IL8 → CXCL8) | The converter retries with a curated alias map. Check HGNC for the current symbol when a gene still fails. |
| Fold enrichment ≥100,000× | Zero RefSeq background overlap (sentinel) | Filtered out of the top-20 by requiring ≥2 overlaps. Do not cite as a real enrichment. |
| No significant enrichments | Truly no shared regulators, or wrong threshold | Try threshold 50, widen TSS window, or check that the genes really share a biological context. |
| Ensembl mapped 0/N while API analysed N | The two systems use different gene databases (Ensembl vs RefSeq) | Enrichment results stand. Explain the discrepancy in the user-facing summary; don't dismiss it as "optional." |
| API analysed fewer regions than genes submitted | Gene not in RefSeq, duplicate promoter, or unrecognised symbol | Surface the count from `api_total_regions`. The API does not say which gene was dropped, so do not speculate. |
| SVG export warnings | Optional `cairosvg` missing | PNGs are always produced. |
| `Median FE (All)` looks tiny while a few experiments are very enriched | Diluted by non-significant experiments | Quote `Median FE (Sig)` for that factor instead. |
| Many factors flagged significant at small N | Each gene contributes a huge share of overlap | Lead with the exploratory framing; note in the summary that a single gene's inclusion/exclusion can eliminate signal. |

## Best Practices

- **Cite Q-values verbatim.** Numbers in the report and CSVs are
  computed; do not paraphrase or round when describing them.
- **Distinguish experiments from factors.** Six MYC ChIP-seq experiments
  in the top 20 are one factor's signal seen six times, not six
  independent discoveries.
- **Acknowledge availability bias.** Well-studied TFs (TP53, RELA,
  CTCF) have hundreds of experiments and naturally show up in any broad
  enrichment. Note their experiment count from the aggregated table.
- **Flag external biology explicitly.** If you write "BRD4 reads
  acetylated histones at active promoters," tag it as background
  knowledge — not a conclusion from this analysis.
- **Total significant factors.** The header in `summary_report.md` gives
  the full count; never report "10 factors are significant" when more
  exist.

## Suggested Next Steps

1. **Threshold sensitivity** — rerun at `threshold=100` (or 200) and
   compare the top factors. If they survive, the call is robust.
2. **Cell-type-specific rerun** — set `cell_class` to the system that
   matches the experiment to remove well-studied-cell-line bias.
3. **Validate with target-gene direction** — for top factors, run
   `chip-atlas-target-genes` on each TF and check whether your input
   genes appear in *its* top targets.
4. **Motif analysis** — confirm direct binding by looking for the top
   factor's motif in the promoter regions (HOMER, MEME).
5. **Differential binding** — `chip-atlas-diff-analysis` if you have two
   conditions of the top factor's ChIP-seq available publicly.

## Related Skills

- **chip-atlas-target-genes** — reverse-direction sibling: given a TF,
  what are its targets?
- **chip-atlas-diff-analysis** — sibling for differential binding between
  two experiment groups.
- **functional-enrichment-from-degs** — pair regulator discovery with
  pathway / GO enrichment for a richer picture of the gene list.

## References

Reference files under `references/`:

- `api_format.md` — POST parameters, supported genomes, antigen and
  cell-class vocabularies, response columns.
- `enrichment_statistics.md` — Fisher's exact contingency table,
  Benjamini-Hochberg correction, sentinel fold-enrichment values, data
  availability bias.
- `peak_thresholds.md` — what `threshold` means and how to choose it.
- `promoter_definitions.md` — TSS window mechanics and the optional
  Ensembl coordinate verification step.

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
- Zhang Y. *et al.* (2008). Model-based analysis of ChIP-seq (MACS).
  *Genome Biology* 9(9):R137.
  [doi:10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)

ChIP-Atlas portal: <https://chip-atlas.org>
