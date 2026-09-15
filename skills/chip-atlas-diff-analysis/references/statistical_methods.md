# Statistical Methods Behind the Diff API

This reference explains what happens server-side once the API accepts a
job. Knowing the model assumptions helps decide which findings to trust
and how to communicate uncertainty.

## Differential Peak Regions — edgeR

The `diffbind` mode (used for ChIP-seq, ATAC-seq, and DNase-seq) chains a
fixed pipeline on the server. None of the steps are user-configurable.

1. **Track retrieval.** BigWig coverage and peak BED files are pulled from
   the ChIP-Atlas back-end for the submitted accessions.
2. **bedGraph conversion.** BigWigs are flattened into bedGraph so the
   read-coverage profile can be reconstructed.
3. **Read-count reconstruction.** Total mapped read counts (recorded in
   each track's metadata) are used to scale the coverage profile back to
   raw read counts.
4. **Region universe.** The genome is segmented by the union of peak
   intervals across all submitted experiments — every region that any
   submitted dataset called as a peak becomes a row.
5. **Count matrix.** Read counts per region per experiment are aggregated
   into an `m × n` matrix (`m` regions, `n` experiments across both groups).
6. **Negative-binomial GLM (edgeR).** The matrix is fed to edgeR's
   quasi-likelihood pipeline. Library sizes are normalized with TMM,
   dispersion is estimated empirically, and a QL F-test is applied to
   each region for the two-group contrast.
7. **FDR control.** Per-region p-values are corrected with
   Benjamini-Hochberg; the corrected q-value lands in column 8 of the
   output BED.

### Modeling assumptions

- Within-group experiments are treated as biological replicates. Any
  within-group heterogeneity (genotype differences, donor variability,
  batch effects) contributes to the dispersion estimate as if it were
  technical noise.
- The TMM normalization assumes most regions are *not* differentially
  bound between the two groups. Genome-wide shifts (e.g. global
  histone-mark depletion) violate this assumption and inflate the type-I
  error of the surviving comparisons.
- With small `n` (≤3 per group) the dispersion estimate is unstable. FDR
  values should be treated as approximate rankings rather than calibrated
  error rates.

### Relationship to DiffBind

The pipeline parallels the popular `DiffBind` R package but skips the
BAM-handling step: ChIP-Atlas has already produced BigWig coverage, so
the API just reconstructs counts from those tracks. The trade-off is
losing the ability to recall reads from a custom region set; the region
universe is fixed to the union of called peaks.

### Citation

Robinson M.D., McCarthy D.J., Smyth G.K. (2010). edgeR: a Bioconductor
package for differential expression analysis of digital gene expression
data. *Bioinformatics* 26(1):139–140.
[doi:10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

## Differentially Methylated Regions — metilene

The `dmr` mode runs an entirely different pipeline tailored to
bisulfite-seq.

1. **bedGraph conversion.** Methylation BigWigs are converted to bedGraph.
2. **Per-CpG aggregation.** The `metilene_input.pl` helper assembles a
   matrix of per-CpG β-values across all submitted experiments.
3. **Segmentation + testing.** `metilene` walks the genome with circular
   binary segmentation, identifying contiguous stretches of CpGs whose
   methylation differs between the two groups. At each CpG, a two-sample
   Mann-Whitney U test compares group A vs B.
4. **Multiple-testing correction.** A region-level q-value is reported.

### Modeling notes

- metilene is non-parametric, so it tolerates non-normal β-value
  distributions better than t-test-based DMR callers.
- mtDNA methylation can be biologically meaningful, so chrM is *not*
  excluded automatically in DMR mode (it is in `diffbind` mode).
- The DMR endpoint expects bisulfite-seq experiments specifically; submitting
  ChIP/ATAC accessions with `analysis_type="dmr"` will produce nothing
  useful.

### Citation

Jühling F., Kretzmer H., Bernhart S.H., Otto C., Stadler P.F., Hoffmann S.
(2016). metilene: fast and sensitive calling of differentially methylated
regions from bisulfite sequencing data. *Genome Research* 26(2):256–262.
[doi:10.1101/gr.196394.115](https://doi.org/10.1101/gr.196394.115)

## Output schema (recap)

Both pipelines emit an 8-column BED:

```
chrom  chromStart  chromEnd  counts_A  counts_B  logFC  pvalue  qvalue
```

Counts are comma-separated, one value per submitted experiment.
`logFC > 0` means enrichment in group A. The full column documentation
lives in `output_format.md`.

## When to override the defaults

The API exposes no tuning knobs for the statistical model. If you need to
change normalization, dispersion estimation, or call DMRs with a
different segmentation, you must download the raw tracks yourself and
run edgeR or metilene locally. That is a much larger lift than this
skill — use it only when the default pipeline is genuinely insufficient
(rare).
