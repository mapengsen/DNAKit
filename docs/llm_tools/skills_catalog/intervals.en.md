# Genomic intervals

**14** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`deeptools-ngs-analysis`](#deeptools-ngs-analysis) | sciagent | NGS CLI for ChIP/RNA/ATAC-seq. BAM→bigWig with RPGC/CPM/RPKM, sample correlation/PCA, heatmaps/profiles around features, fingerprints. For alignment use STAR/B… |
| [`geniml`](#geniml) | kdense | Use Geniml for audited local genomic-interval workflows: validate BED and universe contracts, plan Region2Vec or scEmbed runs, inspect model/tokenizer compatib… |
| [`genome-intervals-bed-file-basics`](#genome-intervals-bed-file-basics) | bioskills | Handles BED-format genomic intervals (BED3 through BED12, narrowPeak/broadPeak) and the coordinate-system substrate the whole interval category rests on, with… |
| [`genome-intervals-bedgraph-handling`](#genome-intervals-bedgraph-handling) | bioskills | Generates, normalizes, and converts bedGraph signal tracks (4-column chrom/start/end/value, 0-based half-open) with bedtools genomecov, deepTools bamCoverage/b… |
| [`genome-intervals-bigwig-tracks`](#genome-intervals-bigwig-tracks) | bioskills | Reads, queries, and writes bigWig indexed binary signal tracks (coverage, fold-change, conservation, methylation-rate) with pyBigWig (Python) and the UCSC Kent… |
| [`genome-intervals-coverage-analysis`](#genome-intervals-coverage-analysis) | bioskills | Computes and interprets sequencing read depth and coverage over a genome, windows, or target regions with mosdepth (windowed depth, cumulative distribution, --… |
| [`genome-intervals-gtf-gff-handling`](#genome-intervals-gtf-gff-handling) | bioskills | Parses, queries, converts, and extracts from GTF and GFF3 gene-model annotation files - walking the gene/transcript/exon/CDS hierarchy with gffutils (queryable… |
| [`genome-intervals-interval-arithmetic`](#genome-intervals-interval-arithmetic) | bioskills | Performs set operations on genomic intervals - intersect (-wa/-wb/-wo/-wao/-loj/-c/-v/-u), subtract (-A), merge (-d, -c/-o), complement, cluster, multiinter, u… |
| [`genome-intervals-overlap-significance`](#genome-intervals-overlap-significance) | bioskills | Tests whether two genomic interval sets overlap (colocalize) more than expected by chance using a permutation test against a structured-genome null model. Cove… |
| [`genome-intervals-proximity-operations`](#genome-intervals-proximity-operations) | bioskills | Performs proximity operations on genomic intervals with bedtools (closest, window, flank, slop) and pybedtools - nearest-feature queries with signed/strand-awa… |
| [`genomic-coordinates`](#genomic-coordinates) | kdense | Convert genomic intervals between coordinate conventions, normalise and compare variant representations, and detect assembly or contig-naming mismatches before… |
| [`importing-tracks`](#importing-tracks) | misha | Use when importing external data (WIG / bigWig / BED / TSV / BAM / 2D contacts / array matrices) into a misha track DB, when picking which `gtrack.*import*` /… |
| [`liftover`](#liftover) | statgen-skills | KG-aware genomic liftover for GWAS summary statistics. Resolves source/target genome build (hg18/hg19/hg38) via the statgen-analysis OKG when available, fallin… |
| [`liftover-coordinates`](#liftover-coordinates) | encode-toolkit | Convert genomic coordinates between assembly versions (GRCh37/hg19 to GRCh38/hg38, mm9 to mm10). Guides UCSC liftOver for BED files, CrossMap for VCF/bigWig, a… |

## `deeptools-ngs-analysis` {#deeptools-ngs-analysis}

NGS CLI for ChIP/RNA/ATAC-seq. BAM→bigWig with RPGC/CPM/RPKM, sample correlation/PCA, heatmaps/profiles around features, fingerprints. For alignment use STAR/BWA; for peak calling use MACS2.

- **Source:** [jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/interval-ops/deeptools-ngs-analysis)
- **License:** `CC-BY-4.0`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info deeptools-ngs-analysis
dnakit skills install deeptools-ngs-analysis --agent claude-code   # or --agent codex / cursor
```

## `geniml` {#geniml}

Use Geniml for audited local genomic-interval workflows: validate BED and universe contracts, plan Region2Vec or scEmbed runs, inspect model/tokenizer compatibility, and assess consensus universes.

- **Source:** [K-Dense Inc.; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/geniml)
- **License:** `MIT`
- **Also known as:** `alterlab-geniml`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info geniml
dnakit skills install geniml --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-bed-file-basics` {#genome-intervals-bed-file-basics}

Handles BED-format genomic intervals (BED3 through BED12, narrowPeak/broadPeak) and the coordinate-system substrate the whole interval category rests on, with bedtools (CLI) and pybedtools/pyranges/pandas (Python). Covers the 0-based half-open vs 1-based-closed convention boundary and the start-1/end-unchanged conversion, the silent failures (chrom-name mismatch, CRLF, lexicographic-vs-version sort under -sorted), genome/chrom.sizes generation, sorting contracts, BED12 block invariants, validation, makewindows, cross-assembly liftover (liftOver/CrossMap), and BED&lt;-&gt;VCF/BAM/FASTA conversion. Use when reading, creating, validating, sorting, lifting between genome builds, or converting interval files, preparing inputs for bedtools/tabix/bigBed, or debugging an off-by-one or empty-overlap result.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/bed-file-basics)
- **License:** `MIT`
- **Also known as:** `bio-genome-intervals-bed-file-basics`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-bed-file-basics
dnakit skills install genome-intervals-bed-file-basics --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-bedgraph-handling` {#genome-intervals-bedgraph-handling}

Generates, normalizes, and converts bedGraph signal tracks (4-column chrom/start/end/value, 0-based half-open) with bedtools genomecov, deepTools bamCoverage/bamCompare/bigwigCompare, bedtools unionbedg, and UCSC bedGraphToBigWig. Covers why a raw coverage bedGraph is not comparable across samples until normalized, the CPM/RPKM/BPM/RPGC normalization menu and the conserved-total assumption that makes them wrong under a global perturbation, the strict sorted-non-overlapping-chrom.sizes bedGraphToBigWig contract that silently corrupts a bigWig, effective-genome-size selection, and bin-size aliasing. Use when building or normalizing a coverage/signal track from a BAM, comparing tracks across samples or conditions, converting bedGraph to a browser-ready bigWig, or diagnosing a track that looks plausible but reports wrong heights.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/bedgraph-handling)
- **License:** `MIT`
- **Also known as:** `bio-genome-intervals-bedgraph-handling`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-bedgraph-handling
dnakit skills install genome-intervals-bedgraph-handling --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-bigwig-tracks` {#genome-intervals-bigwig-tracks}

Reads, queries, and writes bigWig indexed binary signal tracks (coverage, fold-change, conservation, methylation-rate) with pyBigWig (Python) and the UCSC Kent tools (bedGraphToBigWig, bigWigToBedGraph, bigWigInfo, bigWigSummary, bigWigAverageOverBed) and deepTools (multiBigwigSummary, computeMatrix, bigwigCompare). Covers the central trap that a wide query returns a precomputed zoom-level summary (by default the mean, which annihilates narrow peaks) not per-base data, when exact=True/values() is mandatory, the NaN-not-zero gap-handling fork, choosing mean vs max vs sum vs coverage by biological question, and the sorted-bedGraph plus chrom.sizes build requirement. Use when extracting signal at regions, computing mean signal per gene/peak, building a browser track from bedGraph, comparing tracks, or building TSS/gene-body metaprofiles.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/bigwig-tracks)
- **License:** `MIT`
- **Also known as:** `bio-genome-intervals-bigwig-tracks`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-bigwig-tracks
dnakit skills install genome-intervals-bigwig-tracks --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-coverage-analysis` {#genome-intervals-coverage-analysis}

Computes and interprets sequencing read depth and coverage over a genome, windows, or target regions with mosdepth (windowed depth, cumulative distribution, --quantize callable BEDs), bedtools genomecov/coverage (bedGraph tracks, per-target stats), samtools depth/coverage (per-base depth, per-contig depth+breadth). Covers the breadth-vs-mean distinction, the cumulative-coverage curve, evenness (CV/Fano/fold-80/Gini), what each tool silently counts (duplicates, secondary/supplementary, MAPQ, read span vs fragment, mate-overlap), the samtools-depth 8000-cap version trap, and the bedtools coverage -a/-b orientation flip. Use when assessing sequencing adequacy, building coverage tracks, computing breadth at a depth threshold, defining callable regions, or QCing target-capture uniformity.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/coverage-analysis)
- **License:** `MIT`
- **Also known as:** `bio-genome-intervals-coverage-analysis`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-coverage-analysis
dnakit skills install genome-intervals-coverage-analysis --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-gtf-gff-handling` {#genome-intervals-gtf-gff-handling}

Parses, queries, converts, and extracts from GTF and GFF3 gene-model annotation files - walking the gene/transcript/exon/CDS hierarchy with gffutils (queryable SQLite DB), converting formats and extracting transcript/CDS/protein FASTA with gffread, slurping to dataframes with gtfparse/pyranges, and sanitizing malformed files with AGAT. Covers the 1-based-inclusive vs 0-based BED coordinate conversion (start-1 only), deriving implicit features (introns/UTRs/TSS), phase-not-frame, the stop-codon-in-or-out-of-CDS convention, and the chr1-vs-1 seqid and gene-ID-version mismatches that silently produce all-zero count matrices and dropped joins. Use when extracting features or sequences from an annotation, converting GTF&lt;-&gt;GFF3 or GTF-&gt;BED, traversing the gene tree, or diagnosing a coordinate/provenance mismatch upstream of counting or DE.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/gtf-gff-handling)
- **License:** `MIT`
- **Also known as:** `bio-genome-intervals-gtf-gff-handling`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-gtf-gff-handling
dnakit skills install genome-intervals-gtf-gff-handling --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-interval-arithmetic` {#genome-intervals-interval-arithmetic}

Performs set operations on genomic intervals - intersect (-wa/-wb/-wo/-wao/-loj/-c/-v/-u), subtract (-A), merge (-d, -c/-o), complement, cluster, multiinter, unionbedg, map, and groupby - with bedtools (CLI) and pybedtools/pyranges/bioframe (Python). Covers the sorted-input contract and the -sorted chromosome-order footgun, reciprocal/fractional overlap (-f/-F/-r/-e) and the A-vs-B asymmetry, -split for spliced/BED12/BAM features, and jaccard/fisher as mechanics only. Use when finding overlapping or unique regions between BED/peak/feature files, building consensus peaksets, removing blacklisted regions, transferring annotation values onto intervals, or computing interval-set similarity; route overlap-significance testing to overlap-significance.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/interval-arithmetic)
- **License:** `MIT`
- **Also known as:** `alterlab-gtars`, `bedtools`, `bedtools-genomic-intervals`, `bio-genome-intervals-interval-arithmetic`, `gtars`, `polars-bio`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-interval-arithmetic
dnakit skills install genome-intervals-interval-arithmetic --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-overlap-significance` {#genome-intervals-overlap-significance}

Tests whether two genomic interval sets overlap (colocalize) more than expected by chance using a permutation test against a structured-genome null model. Covers bedtools fisher (analytic 2x2 screen), bedtools shuffle + jaccard permutation, GAT (isochore/GC-conditioned simulation with FDR), regioneR (flexible permutation, randomizeRegions vs circularRandomizeRegions, localZScore), LOLA (universe-relative Fisher against a region database), and GREAT/rGREAT (regulatory-domain binomial + hypergeometric for ontology-from-regions). Stresses the universe/background choice, matched background, blacklist exclusion, and multiple-testing control. Use when asking whether peaks/regions are enriched at enhancers/TFBS/features, scoring region-set colocalization or region-set enrichment, comparing CNV/SV concordance, or turning an overlap count into a defensible p-value.

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/overlap-significance)
- **License:** `MIT`
- **Also known as:** `bio-genome-intervals-overlap-significance`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-overlap-significance
dnakit skills install genome-intervals-overlap-significance --agent claude-code   # or --agent codex / cursor
```

## `genome-intervals-proximity-operations` {#genome-intervals-proximity-operations}

Performs proximity operations on genomic intervals with bedtools (closest, window, flank, slop) and pybedtools - nearest-feature queries with signed/strand-aware distance, fixed-radius window searches, strand-aware promoter construction, and interval extension. Covers the closest -d/-D a/b/ref/-t/-k/-io/-iu/-id flags, the -D ref strand sign-flip, silent chromosome-end clipping in slop/flank, -t all tie double-counting, and the critical distinction between a geometry answer (nearest TSS) and a biology answer (which gene an element regulates). Use when assigning peaks or variants to genes, defining promoters from a gene model, building distance-to-TSS distributions, finding features within a window, or extending intervals - and when deciding whether nearest-gene is a fair prior (GWAS locus) or a trap (distal enhancer).

- **Source:** [GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/genome-intervals/proximity-operations)
- **License:** `MIT`
- **Also known as:** `bio-genome-intervals-proximity-operations`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genome-intervals-proximity-operations
dnakit skills install genome-intervals-proximity-operations --agent claude-code   # or --agent codex / cursor
```

## `genomic-coordinates` {#genomic-coordinates}

Convert genomic intervals between coordinate conventions, normalise and compare variant representations, and detect assembly or contig-naming mismatches before they corrupt an analysis. Use whenever coordinates cross a format, tool, or assembly boundary - converting between BED, GFF/GTF, VCF, SAM/BAM, WIG, PSL, genePred, Picard interval_list, or region strings; reconciling 0-based half-open with 1-based inclusive; left-aligning or trimming indels; checking whether two variant records describe the same change; mapping genomic to transcript, CDS, or protein positions; auditing a BED/GTF/VCF for convention violations; or diagnosing GRCh37 vs hg19 vs GRCh38 vs T2T, chr-prefix, and liftover problems. Triggers include "off by one", "0-based", "1-based", "half-open", "coordinate system", "left-align", "normalize variant", "bcftools norm", "chr prefix", "wrong genome build", "liftover", "REF mismatch", and "HGVS".

- **Source:** [K-Dense Inc.; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/genomic-coordinates)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info genomic-coordinates
dnakit skills install genomic-coordinates --agent claude-code   # or --agent codex / cursor
```

## `importing-tracks` {#importing-tracks}

Use when importing external data (WIG / bigWig / BED / TSV / BAM / 2D contacts / array matrices) into a misha track DB, when picking which `gtrack.*import*` / `gtrack.create_*` / `gtrack.liftover` / `gtrack.copy` variant to call, or when a track import fails opaquely or produces an all-zero / all-NaN track.

- **Source:** [tanaylab/misha contributors @ `b13d4e2`](https://github.com/tanaylab/misha/tree/b13d4e2593d1883acd43efb30f9bf1132891fc3b/agent-guides/skills/importing-tracks)
- **License:** `MIT`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info importing-tracks
dnakit skills install importing-tracks --agent claude-code   # or --agent codex / cursor
```

## `liftover` {#liftover}

KG-aware genomic liftover for GWAS summary statistics. Resolves source/target genome build (hg18/hg19/hg38) via the statgen-analysis OKG when available, falling back to the GWAS Catalog REST API. Refuses to lift unknown-provenance data and emits a coverage-gap proposal stub. Supports `--gcst (accession)` to fetch+lift directly from the GWAS Catalog. Use when a user asks to lift a sumstats file, asks whether a lift is needed for a GWAS/LD pair, or asks to fetch a GWAS Catalog study and align its build to an LD panel.

- **Source:** [Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/liftover)
- **License:** `MIT`
- **Runs with:** Python modules `pyliftover`; commands `python3`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info liftover
dnakit skills install liftover --agent claude-code   # or --agent codex / cursor
```

## `liftover-coordinates` {#liftover-coordinates}

Convert genomic coordinates between assembly versions (GRCh37/hg19 to GRCh38/hg38, mm9 to mm10). Guides UCSC liftOver for BED files, CrossMap for VCF/bigWig, and handles unmapped regions with provenance logging.

- **Source:** [ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/liftover-coordinates)
- **License:** `AGPL-3.0-only`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info liftover-coordinates
dnakit skills install liftover-coordinates --agent claude-code   # or --agent codex / cursor
```
