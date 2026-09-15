---
id: "skill_67d00809f8d640a690a94a34376d19c5"
name: chip-atlas-diff-analysis
display-name: "ChIP-Atlas Differential Region Analysis"
description: "Submit two groups of public ChIP-seq, ATAC-seq, DNase-seq, or Bisulfite-seq experiments to the ChIP-Atlas Diff Analysis API and recover differential peak regions (DPR via edgeR) or differentially methylated regions (DMR via metilene). Use when you have public experiment accessions and want a server-side comparison without downloading raw data."
category: "epigenomics"
keywords: "ChIP-Atlas, differential binding, edgeR, metilene, DPR, DMR, ChIP-seq, ATAC-seq, DNase-seq, bisulfite, epigenomics, peak comparison"
version: "1.0"
last-updated: "2026-05-20"
---

## DNAKit execution

Use `dnakit skills doctor chip-atlas-diff-analysis` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# ChIP-Atlas Differential Region Analysis

Compare two groups of public assay accessions and recover the genomic regions
that differ between them. The ChIP-Atlas server runs an edgeR-based pipeline
on pre-computed BigWig coverage (for ChIP/ATAC/DNase-seq) or runs `metilene`
on aggregated methylation tracks (for Bisulfite-seq). You submit a pair of
experiment ID lists, poll the job, and download a BED archive — no raw FASTQ
or BAM handling required.

## When to Use This Skill

Reach for this skill when you need to:

- **Contrast two conditions** captured in public ChIP-seq, ATAC-seq, or
  DNase-seq experiments (e.g. treated vs vehicle, tumor vs normal).
- **Find DMRs between two methylation panels** without setting up a
  bisulfite alignment pipeline.
- **Compare two cell types or labs** using only their accession lists.
- **Benchmark your own experiments** against well-curated public sets at the
  region level, with FDR-corrected effect sizes.
- **Side-step raw data downloads** — ChIP-Atlas already aligned, called
  peaks, and built BigWig coverage. You hand it accession IDs.

Don't use this skill for:

- ❌ **Single-experiment peak calling** — call MACS2/MACS3 yourself; this
  API requires ≥2 experiments per group.
- ❌ **Discovering which factors bind your gene list.** That is the question
  `chip-atlas-peak-enrichment` answers.
- ❌ **Building a TF-target network for one factor.** Use
  `chip-atlas-target-genes` instead.
- ❌ **Air-gapped or offline analysis** — the workflow needs network access
  to https://chip-atlas.org and (optionally) UCSC for gene annotations.
- ❌ **Genomes other than the ten listed below.** Only hg38, hg19, mm10,
  mm9, rn6, dm6, dm3, ce11, ce10, sacCer3 are supported by the API.

**Core idea.** You hand the API two newline-separated lists of SRX/ERX/DRX
or GSM accessions. The server segments the genome by the union of called
peaks in your experiments, builds a region × experiment count matrix, runs
an edgeR quasi-likelihood F-test (TMM normalization, negative-binomial GLM),
and returns a BED archive. For methylation, `metilene` runs circular-binary-
segmentation plus a Mann-Whitney U test per CpG. The returned BED is custom
8-column: `chrom start end counts_A counts_B logFC pvalue qvalue`.

## Quick Start

```python
from scripts.load_example_data import load_example_data
from scripts.run_diff_workflow import run_diff_workflow
from scripts.generate_all_plots import generate_all_plots
from scripts.export_all import export_all

data = load_example_data("tp53_k562_dmso_vs_daunorubicin")

results = run_diff_workflow(
    experiments_a=data["experiments_a"],
    experiments_b=data["experiments_b"],
    genome=data["genome"],
    analysis_type=data["analysis_type"],
    description_a=data["description_a"],
    description_b=data["description_b"],
    design_caveats=data.get("design_caveats"),
    output_dir="diff_results",
)
generate_all_plots(results, output_dir="diff_results")
export_all(results, output_dir="diff_results")
```

Runtime for the example: ~5–10 minutes (mostly waiting on the ChIP-Atlas
batch job).

## Installation

| Package | Version | License | Commercial use | Install |
|---------|---------|---------|---------------|---------|
| pandas | ≥1.3 | BSD-3-Clause | Permitted | `pip install pandas` |
| numpy | ≥1.20 | BSD-3-Clause | Permitted | `pip install numpy` |
| requests | ≥2.25 | Apache-2.0 | Permitted | `pip install requests` |
| plotnine | ≥0.10 | MIT | Permitted | `pip install plotnine` |
| plotnine-prism | ≥0.2 | MIT | Permitted | `pip install plotnine-prism` |

One-liner:

```bash
pip install "pandas>=1.3" "numpy>=1.20" "requests>=2.25" "plotnine>=0.10" "plotnine-prism>=0.2"
```

System requirements: outbound HTTPS to `https://dtn1.ddbj.nig.ac.jp` (the
ChIP-Atlas WABI host) and, optionally, `https://api.genome.ucsc.edu` for
nearest-gene annotations. The API job can take up to 15 minutes for large
groups; plan timeouts accordingly.

## Inputs

**Required**

- `experiments_a`, `experiments_b` — lists of public experiment accessions.
  Each list must contain **at least 2** IDs; 3+ per group is recommended
  for stable dispersion estimates. Accepted formats: SRA (`SRX`, `ERX`,
  `DRX`) and GEO (`GSM`). Plain text files, CSVs with an ID column, or
  Python lists are all loadable through `load_user_data` /
  `load_experiment_groups`.

**Optional**

- `genome` — one of `hg38` (default), `hg19`, `mm10`, `mm9`, `rn6`, `dm6`,
  `dm3`, `ce11`, `ce10`, `sacCer3`.
- `analysis_type` — `diffbind` (default; edgeR DPR for ChIP/ATAC/DNase-seq)
  or `dmr` (metilene DMR for Bisulfite-seq).
- `description_a`, `description_b` — display labels (e.g. "Treated",
  "Control"). Default to `group_A` / `group_B`.
- `design_caveats` — list of strings recording known design quirks
  (different genotypes, batches, etc.) that get surfaced near the top of
  `summary_report.md` so downstream consumers don't miss them.
- `min_score`, `min_region_size`, `direction_filter` — client-side
  post-filters applied after default QC. Use sparingly; the defaults are
  designed to be honest about what the API returned.

## Outputs

All written to `output_dir`:

| File | Format | Purpose / consumer |
|------|--------|--------------------|
| `analysis_object.pkl` | Python pickle | Full results dict for any downstream script |
| `diff_regions_all.csv` | CSV | Default QC-filtered table (chr/start/end, counts_a, counts_b, logFC, pvalue, qvalue, region_size, direction, significant, mean_count_a, mean_count_b, score) |
| `diff_regions_significant.csv` | CSV | Subset with qvalue < 0.05 |
| `diff_regions_top50.csv` | CSV | Top 50 by Q-value; includes `nearest_gene` column when UCSC lookup succeeds |
| `diff_regions_unfiltered.csv` | CSV | Pre-QC regions (kept only when QC removed rows) — useful for sanity checks |
| `raw_results/` | folder | Verbatim ChIP-Atlas output: `wabi_result.bed`, `wabi_result.igv.bed`, `wabi_result.xml`, `wabi_result.log`, `README.md` |
| `volcano_plot.{png,svg}` | image | logFC vs -log10(Q-value) |
| `chromosome_distribution.{png,svg}` | image | Significant-region counts per standard chromosome, split by direction; sex-chromosome confounds are flagged with an asterisk |
| `region_size_distribution.{png,svg}` | image | Log-scaled region-size histogram |
| `ma_plot.{png,svg}` | image | Mean log2-count vs logFC |
| `summary_report.md` | markdown | Human-readable summary with QC warnings, top-10 ranking, statistical caveats |

The `analysis_object.pkl` dict contains: `diff_regions`,
`diff_regions_unfiltered`, `qc_warnings`, `experiments_a`,
`experiments_b`, `raw_files`, `log_content`, `parameters`, plus tallied
counts. Reload with `pickle.load(open("analysis_object.pkl","rb"))`.

The IGV session file (`raw_results/wabi_result.xml`) loads the BigWig
tracks and the colored differential regions in a single click — orange =
group A enriched, blue = group B enriched.

## Clarification Questions

1. **⚠️ CRITICAL: ASK THIS FIRST.** What two groups of experiment IDs are
   you comparing?
   - Each group needs at least 2 SRX/ERX/DRX/GSM accessions; 3+ per group
     is preferred for trustworthy FDRs.
   - If the user does not have IDs ready, offer one of the bundled examples:
     `tp53_k562_dmso_vs_daunorubicin` (TP53 ChIP-seq in K-562 cells,
     DMSO vs daunorubicin, n=8/8) or `tp53_molm13_vs_k562` (TP53 ChIP-seq,
     MOLM-13 vs K-562, n=4/4 — demonstrates the sex-chromosome QC check).
2. Which **genome assembly** should the analysis run against? Default is
   `hg38`. The full set is hg19, mm10, mm9, rn6, dm6, dm3, ce11, ce10, and
   sacCer3.
3. Which **analysis type**? `diffbind` (default; for ChIP/ATAC/DNase-seq)
   or `dmr` (for Bisulfite-seq).
4. What **labels** should the two groups carry in plots and the report?
   (Auto-populated for the example data.)
5. Are the experiments within each group **true biological replicates**, or
   do they differ (multiple genotypes, tissues, time points)? Capture any
   heterogeneity as `design_caveats` — it will surface prominently in the
   report.

## Standard Workflow

1. **Load the experiment lists.** Either `load_example_data(name)` or
   `load_user_data(path)` / `load_experiment_groups(...)`. Both validate
   that IDs match the `^(SRX|ERX|DRX|GSM)\d+$` pattern and warn (rather
   than crash) on unusual IDs.
2. **Run the diff workflow** with `run_diff_workflow(...)`. It submits the
   job, polls until completion (15-min timeout, 15-second poll interval),
   extracts the ZIP, parses the 8-column BED, runs QC checks on the
   unfiltered table, then applies default QC filtering (drop non-standard
   contigs and regions < 10 bp).
3. **Generate plots** with `generate_all_plots(results, output_dir)`. PNG
   always succeeds; SVG export degrades gracefully if the optional
   dependency is missing.
4. **Export everything** with `export_all(results, output_dir)`. This is
   where the pickle, the CSV tables, the optional UCSC nearest-gene
   annotations, and the markdown report are written.

## When Scripts Fail

If a step raises, follow this hierarchy before reaching for inline code:

1. **Fix and retry (~90% of the time).** Missing package? `pip install …`
   and rerun. Transient network failure? Wait a minute and retry.
2. **Edit the script (~5%).** The scripts are intentionally short. If a
   parameter is missing or you need to add a filter, modify the script and
   note the change in a brief comment.
3. **Read the script for guidance (~4%).** Adapt the approach but keep the
   contract (column names, output paths) consistent with what
   `export_all.py` expects.
4. **Write from scratch (~1%).** Last resort. Explain in your response why
   the script could not be used.

## Decision Points

- **Sample size.** With `n ≤ 3` per group, edgeR's dispersion estimate is
  unstable. The QC system flags this as a major warning; treat FDRs as
  approximate rankings, not strict false-discovery rates.
- **chrY / chrX one-directional regions.** Frequently a sex mismatch
  between the two cell lines. The QC layer flags this when >80% of regions
  on the chromosome point the same way *and* the depleted group has near-
  zero counts. Treat flagged regions as artifacts.
- **chrM regions in `diffbind` mode.** Mitochondrial DNA lacks histones;
  significant chrM peaks in ChIP-seq/ATAC-seq are almost always
  contamination. The QC check is skipped in `dmr` mode where mitochondrial
  methylation differences can be real.
- **Default QC filters.** Non-standard contigs and regions < 10 bp are
  dropped from `diff_regions_all.csv`; the original table is preserved in
  `diff_regions_unfiltered.csv`. Override only if you know why.
- **Dense significant-region clusters.** ≥10 significant regions within a
  2 Mb window on one chromosome often reflect a copy-number amplification
  in one group rather than differential binding. Investigate in IGV.

## Common Issues

| Symptom | Likely cause | What to do |
|---------|--------------|------------|
| HTTP 400 from the submit call | Invalid genome, analysis_type, or malformed IDs | Re-check accession format and that `genome` is one of the supported assemblies. See `references/diff_api_format.md`. |
| Job exceeds the 15-minute poll timeout | Large groups or transient server load | Increase `timeout` on `poll_job_status`, or reduce experiments per group. |
| `Could not parse request ID from API response` | The API returned an unexpected body | The endpoint emits text — not JSON — for diff analysis. Inspect the raw response (printed on failure) for an error message. |
| `No BED file found in ZIP results` | The job finished but produced no usable output | Verify both groups have valid experiments for the chosen genome; some accessions are only available in older builds. |
| Empty differential table | Groups are too similar, or peak calls don't intersect | Confirm the experiments capture the conditions you think they do. |
| All chrY/chrX regions one-directional | Likely sex mismatch between the groups | Exclude sex chromosomes from interpretation; the QC system flags this. |
| Dense significant cluster on one chromosome | Possible CNV/structural variant | Open the region in IGV; QC system flags clusters of ≥10 regions in <2 Mb. |
| `mean significant region size < 100 bp` (ChIP/ATAC) | Peak fragmentation in the differential call | Check representative regions in IGV; small sub-peaks may sit inside larger consensus peaks. |
| UCSC API errors during gene annotation | Network or rate limit | Annotations are optional. The script aborts gene lookup gracefully after three failures; statistical results are unaffected. |
| SVG export warnings | Optional dependency missing | PNGs are always written. Install `cairosvg` if SVGs are required. |
| `Invalid ZIP data` from result fetch | Job is still finishing despite the status check | Wait a few seconds and call `retrieve_zip_results` again with the same `request_id`. |
| `MA plot asymmetry` warning | Significant hits in one direction skew toward low-count regions | Trust the warning's quantitative threshold (>1 log2 median difference). Visual impressions below that threshold are not meaningful. |

## Best Practices

- **State the design caveats up front.** edgeR averages over within-group
  heterogeneity. If group A is "8 different TP53 mutants in K-562",
  highlight that in the design caveats so the reader knows the analysis
  is the *average* mutant-vs-WT effect, not a clean replicate test.
- **Quote numbers from the CSVs, not from memory.** When narrating top
  regions, copy logFC / Q-value / coordinates from `diff_regions_top50.csv`
  or `summary_report.md` verbatim.
- **Rank strictly by Q-value.** Don't skip unannotated regions when listing
  the top differential set. Lacking a nearby gene does not diminish
  statistical significance.
- **Distinguish data from biology.** "A differential region near PCNA" is a
  proximity observation, not a claim of direct regulation.
- **Validate before publishing.** Orthogonal evidence (RNA-seq, perturbation,
  motif analysis) is required for any biological claim downstream.

## Suggested Next Steps

1. **Motif enrichment** with HOMER or MEME-ChIP on the significant region
   set — useful when the differential signal is dominated by TF binding.
2. **Region-to-gene linking** via `chip-atlas-target-genes` (when one of
   the cell types has a well-studied TF) or GREAT for gene ontology of
   distal regions.
3. **Cross-reference differential expression** — pair this with
   `bulk-rnaseq-counts-to-de-deseq2` or `bulk-rnaseq-differential-expression`
   to link regulatory shifts to transcriptional changes.
4. **Visualise in IGV.** Open `raw_results/wabi_result.xml`; the colored
   `.igv.bed` is also loadable standalone.

## Related Skills

- **chip-atlas-peak-enrichment** — sibling skill; same API host. Asks the
  reverse question: which factors enrich near *my* gene set?
- **chip-atlas-target-genes** — sibling skill; pre-computed TF→target
  rankings for a single factor.
- **bulk-rnaseq-counts-to-de-deseq2** — pair differential binding with
  differential expression for the same contrast.

## References

Reference files under `references/`:

- `diff_api_format.md` — full API parameter table, response format,
  experiment ID rules.
- `statistical_methods.md` — how edgeR and metilene pipelines are
  configured server-side, with citations.
- `output_format.md` — column-by-column description of the BED output,
  including the IGV BED9+GFF3 visualization file.

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
- Robinson M.D., McCarthy D.J., Smyth G.K. (2010). edgeR: a Bioconductor
  package for differential expression analysis of digital gene expression
  data. *Bioinformatics* 26(1):139–140.
  [doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)
- Jühling F. *et al.* (2016). metilene: fast and sensitive calling of
  differentially methylated regions from bisulfite sequencing data.
  *Genome Research* 26(2):256–262.
  [doi:10.1101/gr.196394.115](https://doi.org/10.1101/gr.196394.115)

ChIP-Atlas portal: <https://chip-atlas.org>
