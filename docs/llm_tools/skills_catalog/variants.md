# 变异检测与处理

本类共 **33** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`consensus-sequences`](#consensus-sequences) | bioskills | Generate consensus FASTA sequences by applying VCF variants onto a reference with bcftools consensus, or build viral/amplicon consensus with iVar. Use when rec… |
| [`copy-number-allele-specific-copy-number`](#copy-number-allele-specific-copy-number) | bioskills | Infer integer allele-specific copy number, tumor purity, and ploidy from tumor sequencing by jointly modeling read depth (logR) and B-allele frequency (BAF) wi… |
| [`copy-number-cnv-annotation`](#copy-number-cnv-annotation) | bioskills | Annotate copy number variant segments with overlapping genes, dosage-sensitivity scores, cancer driver databases, population frequencies, and clinical-variant… |
| [`copy-number-cnv-visualization`](#copy-number-cnv-visualization) | bioskills | Visualize copy number profiles, segments, allele-specific tracks, and cohort patterns from CNVkit, GATK, ASCAT, FACETS, Sequenza, and other callers. Covers gen… |
| [`copy-number-cnvkit-analysis`](#copy-number-cnvkit-analysis) | bioskills | Detect somatic and germline copy number variants from targeted, exome, and whole-genome sequencing with CNVkit, a read-depth caller that combines on-target and… |
| [`copy-number-copy-ratio-segmentation`](#copy-number-copy-ratio-segmentation) | bioskills | Normalize read-depth copy-ratio profiles and segment them into copy-number regions using circular binary segmentation (CBS, DNAcopy), hidden Markov models, Haa… |
| [`copy-number-focal-amplification-ecdna`](#copy-number-focal-amplification-ecdna) | bioskills | Resolve the architecture of focal oncogene amplifications — extrachromosomal DNA (ecDNA), breakage-fusion-bridge (BFB) cycles, homogeneously staining regions (… |
| [`copy-number-gatk-cnv`](#copy-number-gatk-cnv) | bioskills | Call copy number variants with the GATK best-practices workflows — the somatic CNV pipeline (CollectReadCounts, DenoiseReadCounts with tangent normalization, M… |
| [`copy-number-germline-cnv-interpretation`](#copy-number-germline-cnv-interpretation) | bioskills | Classify constitutional (germline) copy number variants for clinical reporting using the 2019 ACMG/ClinGen technical standards points-based framework, with Cla… |
| [`copy-number-hrd-scoring`](#copy-number-hrd-scoring) | bioskills | Quantify homologous recombination deficiency (HRD) from tumor copy number using the three genomic-scar metrics — loss of heterozygosity (LOH), large-scale stat… |
| [`copy-number-recurrent-cnv`](#copy-number-recurrent-cnv) | bioskills | Identify recurrent and driver copy number alterations across a tumor cohort with GISTIC2 (G-score, Ziggurat deconstruction, focal vs broad/arm-level analysis,… |
| [`copy-number-subclonal-copy-number`](#copy-number-subclonal-copy-number) | bioskills | Resolve subclonal copy number, whole-genome doubling, and copy-number tumor evolution from bulk sequencing with Battenberg, TITAN, and MEDICC2. Covers clonal v… |
| [`data-visualization-lollipop-protein-maps`](#data-visualization-lollipop-protein-maps) | bioskills | Plot per-gene mutation distributions on a protein-domain map (lollipop / needle plots) showing mutation position, recurrence count, and variant classification… |
| [`gatk-variant-calling`](#gatk-variant-calling) | bioskills | Call germline SNPs and indels with GATK HaplotypeCaller and the GVCF joint-genotyping workflow. Covers the local-reassembly + PairHMM mechanism (why HC beats p… |
| [`gigwa-import-and-qc`](#gigwa-import-and-qc) | gigwa-mcp | Import a genotype dataset into a Gigwa server and run standard quality control. Use this when a user wants to load a DArTseq SNP/Silico xlsx report or a VCF (.… |
| [`gigwa-qc-triage`](#gigwa-qc-triage) | gigwa-mcp | Run the full quality-control suite on a Gigwa variant set and give a go/no-go verdict for downstream analysis. Use this when a user wants to triage or vet the… |
| [`outlier-splicing-detection`](#outlier-splicing-detection) | user335-bioskills | Detects aberrant splicing in single rare-disease patients vs a control panel using FRASER 2.0 (Bioconductor; Beta-binomial autoencoder on Intron Jaccard Index,… |
| [`single-cell-cnv-inference`](#single-cell-cnv-inference) | user335-bioskills | Infer large-scale copy-number alterations from tumor single-cell or single-nucleus RNA-seq to separate malignant from normal cells and call subclones, using in… |
| [`splice-variant-prediction`](#splice-variant-prediction) | bioskills | Predicts whether a DNA variant alters mRNA splicing using sequence-based deep-learning tools — SpliceAI (10kb context dilated CNN, clinical default), Pangolin… |
| [`tiledbvcf`](#tiledbvcf) | kdense | Efficient storage and retrieval of genomic variant data using TileDB. Scalable VCF/BCF ingestion, incremental sample addition, compressed storage, parallel que… |
| [`variant-analysis`](#variant-analysis) | tooluniverse | VCF and variant analysis — parsing, annotation, classification (synonymous, missense, frameshift, stop_gained), VAF filtering, coding vs non-coding categorizat… |
| [`variant-annotate`](#variant-annotate) | statgen-skills | Per-rsID variant annotation against gnomAD (v4, GRCh38) and GTEx (v8). For a list of rsIDs, fetches the most-severe coding consequence + gene + per-population… |
| [`variant-annotation`](#variant-annotation) | bioskills | Annotates VCF variants with functional consequences, population frequencies, and pathogenicity scores using bcftools annotate/csq, Ensembl VEP, SnpEff, and ANN… |
| [`variant-calling`](#variant-calling) | bioskills | Call germline SNPs and indels from a BAM/CRAM with bcftools mpileup and call, and select the right calling engine for the job. Use when generating a VCF from a… |
| [`variant-calling-deepvariant`](#variant-calling-deepvariant) | bioskills | Calls germline SNPs and indels with Google DeepVariant, which reframes variant calling as CNN image classification over multi-channel pileup tensors. Covers pl… |
| [`variant-calling-filtering-best-practices`](#variant-calling-filtering-best-practices) | bioskills | Filters germline and somatic variant callsets at the site and genotype level with GATK VQSR (VQSLOD, truth-sensitivity tranches), VETS/ScoreVariantAnnotations,… |
| [`variant-calling-joint-calling`](#variant-calling-joint-calling) | bioskills | Joint genotype a cohort of per-sample gVCFs with GATK (HaplotypeCaller -ERC GVCF -&gt; GenomicsDBImport or CombineGVCFs -&gt; GenotypeGVCFs) or GLnexus for DeepVaria… |
| [`variant-calling-structural-variant-calling`](#variant-calling-structural-variant-calling) | bioskills | Call structural variants (&gt;=50 bp deletions, insertions, inversions, duplications, translocations) from short- or long-read data by reconstructing four orthogo… |
| [`variant-normalization`](#variant-normalization) | bioskills | Left-align and trim indels to parsimonious canonical form, decompose MNPs (atomize), and split multiallelic variants with bcftools norm. Use when comparing var… |
| [`variant-predictor-dms-validation`](#variant-predictor-dms-validation) | tooluniverse | Validate a variant-effect predictor (AlphaMissense, ESM-C SAE, ESM logits, EVE, conservation scores, or any per-variant numeric score) against experimental dee… |
| [`vcf-basics`](#vcf-basics) | bioskills | View, query, and interpret VCF/BCF variant files with bcftools and cyvcf2. Use when inspecting variants, extracting fields with query format strings, convertin… |
| [`vcf-manipulation`](#vcf-manipulation) | bioskills | Combine, split, sort, intersect, and subset VCF/BCF files with bcftools merge, concat, isec, sort, view, and reheader. Use when merging different samples into… |
| [`vcf-statistics`](#vcf-statistics) | bioskills | Compute and interpret VCF quality-control metrics (Ti/Tv, het/hom, novel/known, missingness, HWE, contamination, relatedness) with bcftools stats, vcftools, pl… |

## `consensus-sequences` {#consensus-sequences}

Generate consensus FASTA sequences by applying VCF variants onto a reference with bcftools consensus, or build viral/amplicon consensus with iVar. Use when reconstructing a sample-specific reference or haplotype, deciding -H haplotype vs IUPAC vs all-ALT projection, masking no-coverage sites so a consensus does not manufacture false reference calls, or setting iVar min-depth/min-frequency policy for surveillance genomes.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/consensus-sequences)
- **许可证：**`MIT`
- **运行依赖：**命令 `bcftools`
- **别名：**`bio-consensus-sequences`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info consensus-sequences
dnakit skills install consensus-sequences --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-allele-specific-copy-number` {#copy-number-allele-specific-copy-number}

Infer integer allele-specific copy number, tumor purity, and ploidy from tumor sequencing by jointly modeling read depth (logR) and B-allele frequency (BAF) with ASCAT, Sequenza, FACETS, PURPLE, and PureCN (tumor-only). Covers the purity-ploidy identifiability problem, the diploid-baseline (dipLogR) anchor, major/minor copy number, loss of heterozygosity, sunrise/contour fit diagnostics, and reconciliation of conflicting fits. Use when tumor analysis needs absolute copy number rather than relative log2, when estimating purity and ploidy, calling LOH or copy-neutral LOH, resolving whole-genome doubling, running tumor-only allele-specific calling, or choosing among ASCAT, Sequenza, FACETS, and PureCN.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/allele-specific-copy-number)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-allele-specific-copy-number`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-allele-specific-copy-number
dnakit skills install copy-number-allele-specific-copy-number --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-cnv-annotation` {#copy-number-cnv-annotation}

Annotate copy number variant segments with overlapping genes, dosage-sensitivity scores, cancer driver databases, population frequencies, and clinical-variant content. Covers bedtools/pybedtools interval intersection, AnnotSV comprehensive annotation and ranking, ClinGen haploinsufficiency/triplosensitivity scoring, gnomAD-SV/DGV frequency filtering, COSMIC Cancer Gene Census, and ClinVar overlap. Use when interpreting which genes a CNV affects, distinguishing the driver gene of a focal event from passengers, filtering against population CNVs, separating whole-gene from partial-gene overlap, or preparing CNVs for clinical classification.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/cnv-annotation)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-cnv-annotation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-cnv-annotation
dnakit skills install copy-number-cnv-annotation --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-cnv-visualization` {#copy-number-cnv-visualization}

Visualize copy number profiles, segments, allele-specific tracks, and cohort patterns from CNVkit, GATK, ASCAT, FACETS, Sequenza, and other callers. Covers genome-wide and per-chromosome log2 scatter plots, B-allele-frequency/minor-allele-fraction tracks, ideograms, cohort heatmaps, circos views, and caller-native plots. Use when creating publication CNV figures, choosing which plot answers a given question, diagnosing a wrong diploid baseline visually, displaying loss of heterozygosity, or deciding what depth-only plots cannot reveal.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/cnv-visualization)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-cnv-visualization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-cnv-visualization
dnakit skills install copy-number-cnv-visualization --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-cnvkit-analysis` {#copy-number-cnvkit-analysis}

Detect somatic and germline copy number variants from targeted, exome, and whole-genome sequencing with CNVkit, a read-depth caller that combines on-target and off-target (antitarget) coverage. Covers panel-of-normals construction, flat-reference tumor-only calling, hybrid/amplicon/WGS modes, CBS vs HMM segmentation selection, purity-aware integer calling, and reconciliation against GATK and allele-specific callers. Use when calling CNVs from hybrid-capture panels or exomes, deciding whether CNVkit (depth-only) is the right tool versus an allele-specific caller, building a panel of normals, diagnosing flat-reference false positives, or interpreting log2 ratios into copy-number states.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/cnvkit-analysis)
- **许可证：**`MIT`
- **运行依赖：**命令 `cnvkit`
- **别名：**`bio-copy-number-cnvkit-analysis`, `cnvkit-copy-number`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-cnvkit-analysis
dnakit skills install copy-number-cnvkit-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-copy-ratio-segmentation` {#copy-number-copy-ratio-segmentation}

Normalize read-depth copy-ratio profiles and segment them into copy-number regions using circular binary segmentation (CBS, DNAcopy), hidden Markov models, HaarSeg, and fused-lasso methods. Covers GC-content, mappability, and replication-timing (wave-artifact) bias correction, panel-of-normals/PCA denoising, diploid-baseline centering, and algorithm selection by sequencing depth and event size. Use when choosing a segmentation algorithm, correcting depth bias, diagnosing oversegmentation or a mis-centered baseline, tuning CBS or HMM parameters, or understanding why a downstream CNV caller produced fragmented or shifted segments.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/copy-ratio-segmentation)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-copy-ratio-segmentation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-copy-ratio-segmentation
dnakit skills install copy-number-copy-ratio-segmentation --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-focal-amplification-ecdna` {#copy-number-focal-amplification-ecdna}

Resolve the architecture of focal oncogene amplifications — extrachromosomal DNA (ecDNA), breakage-fusion-bridge (BFB) cycles, homogeneously staining regions (HSR), and linear amplification — from whole-genome sequencing with AmpliconArchitect, the AmpliconSuite pipeline, and AmpliconClassifier. Covers copy-number seed selection, breakpoint-graph reconstruction, balanced-flow optimization, ecDNA classification, and the limits of depth-only amplification calls. Use when a focal amplification needs structural characterization, when distinguishing ecDNA from chromosomal amplification, suspecting ecDNA-driven oncogene amplification or therapy resistance, or selecting copy-number seeds for amplicon reconstruction.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/focal-amplification-ecdna)
- **许可证：**`MIT`
- **运行依赖：**命令 `AmpliconArchitect`
- **别名：**`bio-copy-number-focal-amplification-ecdna`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-focal-amplification-ecdna
dnakit skills install copy-number-focal-amplification-ecdna --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-gatk-cnv` {#copy-number-gatk-cnv}

Call copy number variants with the GATK best-practices workflows — the somatic CNV pipeline (CollectReadCounts, DenoiseReadCounts with tangent normalization, ModelSegments, CallCopyRatioSegments) and the germline GATK-gCNV pipeline (DetermineGermlineContigPloidy, GermlineCNVCaller cohort/case mode, PostprocessGermlineCNVCalls). Covers panel-of-normals construction, AnnotateIntervals/FilterIntervals, allelic-count integration, and QS-based filtering. Use when integrating CNV calling into a GATK variant pipeline, calling rare germline CNVs from an exome cohort, deciding between the somatic and germline GATK workflows, or diagnosing why tangent normalization removed a real event or why gCNV output has low precision.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/gatk-cnv)
- **许可证：**`MIT`
- **运行依赖：**命令 `gatk`
- **别名：**`bio-copy-number-gatk-cnv`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-gatk-cnv
dnakit skills install copy-number-gatk-cnv --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-germline-cnv-interpretation` {#copy-number-germline-cnv-interpretation}

Classify constitutional (germline) copy number variants for clinical reporting using the 2019 ACMG/ClinGen technical standards points-based framework, with ClassifyCNV and AnnotSV for semi-automated scoring. Covers the separate copy-number-loss and copy-number-gain rubrics, the five-tier classification, ClinGen haploinsufficiency/triplosensitivity and dosage-sensitive regions, de novo and segregation evidence, and population-frequency benign evidence. Use when assigning pathogenic/likely-pathogenic/VUS/likely-benign/benign to a constitutional CNV, scoring a CNV against ACMG/ClinGen criteria, or distinguishing the automatable evidence from the case-specific evidence requiring manual input.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/germline-cnv-interpretation)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-germline-cnv-interpretation`, `classify-cnv`, `cnv-acmg-classifier`, `tooluniverse-structural-variant-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-germline-cnv-interpretation
dnakit skills install copy-number-germline-cnv-interpretation --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-hrd-scoring` {#copy-number-hrd-scoring}

Quantify homologous recombination deficiency (HRD) from tumor copy number using the three genomic-scar metrics — loss of heterozygosity (LOH), large-scale state transitions (LST), and telomeric allelic imbalance (TAI) — with scarHRD, and via the whole-genome HRDetect and CHORD models. Covers the genomic instability score, the PARP-inhibitor clinical context, whole-genome-doubling correction, and the scar-versus-state distinction. Use when computing an HRD score for PARP-inhibitor eligibility, deriving LOH/LST/TAI scars from allele-specific copy number, deciding between scar-based and mutational-signature HRD methods, or interpreting an HRD result in a BRCA-reverted or low-purity tumor.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/hrd-scoring)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-hrd-scoring`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-hrd-scoring
dnakit skills install copy-number-hrd-scoring --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-recurrent-cnv` {#copy-number-recurrent-cnv}

Identify recurrent and driver copy number alterations across a tumor cohort with GISTIC2 (G-score, Ziggurat deconstruction, focal vs broad/arm-level analysis, q-values from permutation) and quantify copy-number signatures with the Steele 2022 COSMIC framework and the Drews 2022 CINSignatures framework. Covers driver-gene localization from recurrence peaks, distinguishing focal drivers from arm-level passengers, and the caller-sensitivity caveats of copy-number signatures. Use when finding recurrently amplified or deleted regions in a cohort, localizing driver genes, separating focal from broad events, running GISTIC2, or extracting copy-number mutational signatures.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/recurrent-cnv)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-recurrent-cnv`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-recurrent-cnv
dnakit skills install copy-number-recurrent-cnv --agent claude-code   # 或 --agent codex / cursor
```

## `copy-number-subclonal-copy-number` {#copy-number-subclonal-copy-number}

Resolve subclonal copy number, whole-genome doubling, and copy-number tumor evolution from bulk sequencing with Battenberg, TITAN, and MEDICC2. Covers clonal versus subclonal copy-number states, haplotype phasing for subclonal resolution, cancer cell fraction, whole-genome-doubling detection and timing relative to mutations, mirrored subclonal allelic imbalance, and copy-number phylogenies. Use when a tumor is heterogeneous and bulk data shows non-integer copy number, when calling subclonal CNAs, detecting or timing whole-genome doubling, reconstructing copy-number evolution, or deciding between Battenberg and TITAN.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/copy-number/subclonal-copy-number)
- **许可证：**`MIT`
- **别名：**`bio-copy-number-subclonal-copy-number`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info copy-number-subclonal-copy-number
dnakit skills install copy-number-subclonal-copy-number --agent claude-code   # 或 --agent codex / cursor
```

## `data-visualization-lollipop-protein-maps` {#data-visualization-lollipop-protein-maps}

Plot per-gene mutation distributions on a protein-domain map (lollipop / needle plots) showing mutation position, recurrence count, and variant classification with maftools, g3-lollipop, trackViewer, and ProteinPaint. Use when visualizing recurrent mutation hotspots on a single gene's protein, marking domain boundaries from UniProt/Pfam, comparing missense vs truncating distributions, or contrasting two cohorts on the same lollipop.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/data-visualization/lollipop-protein-maps)
- **许可证：**`MIT`
- **别名：**`bio-data-visualization-lollipop-protein-maps`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info data-visualization-lollipop-protein-maps
dnakit skills install data-visualization-lollipop-protein-maps --agent claude-code   # 或 --agent codex / cursor
```

## `gatk-variant-calling` {#gatk-variant-calling}

Call germline SNPs and indels with GATK HaplotypeCaller and the GVCF joint-genotyping workflow. Covers the local-reassembly + PairHMM mechanism (why HC beats pileup callers on indels), the -ERC GVCF reference-confidence model and &lt;NON_REF&gt; allele, BQSR-vs-DRAGSTR and --dragen-mode error modeling, allele-specific (AS_) annotations, and edge cases (ploidy, Mutect2 mitochondria mode, sex chromosomes/PAR, contamination gating). Use when deciding whether to use HaplotypeCaller vs a pileup or DRAGEN caller, whether BQSR still earns its place, whether to call per-sample GVCFs for a cohort, or how to handle non-diploid, mitochondrial, sex-chromosome, or contaminated samples. Not for post-calling filtering depth (see variant-calling/filtering-best-practices) or cohort joint-genotyping scaling (see variant-calling/joint-calling).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/gatk-variant-calling)
- **许可证：**`MIT`
- **运行依赖：**命令 `gatk`
- **别名：**`bio-gatk-variant-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info gatk-variant-calling
dnakit skills install gatk-variant-calling --agent claude-code   # 或 --agent codex / cursor
```

## `gigwa-import-and-qc` {#gigwa-import-and-qc}

Import a genotype dataset into a Gigwa server and run standard quality control. Use this when a user wants to load a DArTseq SNP/Silico xlsx report or a VCF (.vcf/.vcf.gz) into a Gigwa database, project and run, then verify the import is clean. Covers optional genome-anchoring of DArTseq tag sequences to a reference FASTA, confirming variant and sample counts, and the full QC sweep: per-sample and per-marker call rate, heterozygosity outliers, MAF and missingness filtering, duplicate or clonal accessions, plus an import-quality audit that catches genotype-encoding artifacts. Trigger on requests like 'import this DArTseq report into Gigwa', 'load a VCF and QC it', or 'ingest genotypes and tell me if the import looks clean'. Keywords: Gigwa, DArTseq, VCF import, genotyping, call rate, MAF, heterozygosity, duplicate accessions, import audit.

- **来源：**[Copyright 2026 Anestis Gkanogiannis; Gigwa-MCP contributors @ `d6ccdd1`](https://github.com/gkanogiannis/Gigwa-MCP/tree/d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf/skills/gigwa-import-and-qc)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `gigwa_mcp`; 命令 `gigwa-mcp`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info gigwa-import-and-qc
dnakit skills install gigwa-import-and-qc --agent claude-code   # 或 --agent codex / cursor
```

## `gigwa-qc-triage` {#gigwa-qc-triage}

Run the full quality-control suite on a Gigwa variant set and give a go/no-go verdict for downstream analysis. Use this when a user wants to triage or vet the quality of an already-imported run (identified by its variantSetDbId) before trusting diversity or association results. Runs call rate (worst samples and markers), heterozygosity outliers that flag contamination or off-types, duplicate/clonal accession detection, a MAF/missingness filter preview, and an import-quality audit for genotype-encoding artifacts, then summarises the key numbers into a clear go/no-go recommendation. Trigger on 'QC this Gigwa run', 'is this dataset good enough to analyse', 'check call rates and flag bad samples/markers', or 'triage data quality'. Keywords: Gigwa, quality control, QC, call rate, heterozygosity, duplicate accessions, MAF filter, missingness, import audit, go/no-go.

- **来源：**[Copyright 2026 Anestis Gkanogiannis; Gigwa-MCP contributors @ `d6ccdd1`](https://github.com/gkanogiannis/Gigwa-MCP/tree/d6ccdd1a5eb19b8a0cc224b96ad459868d9c4aaf/skills/gigwa-qc-triage)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `gigwa_mcp`; 命令 `gigwa-mcp`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info gigwa-qc-triage
dnakit skills install gigwa-qc-triage --agent claude-code   # 或 --agent codex / cursor
```

## `outlier-splicing-detection` {#outlier-splicing-detection}

Detects aberrant splicing in single rare-disease patients vs a control panel using FRASER 2.0 (Bioconductor; Beta-binomial autoencoder on Intron Jaccard Index, default delta cutoff 0.1, q hyperparameter), OUTRIDER (gene-level outlier expression via autoencoder denoising), LeafcutterMD (Dirichlet-multinomial outlier mode of LeafCutter for annotation-free junctions), and DROP (Snakemake pipeline integrating FRASER2 + OUTRIDER + monoallelic expression for clinical diagnostics). The statistical model is fundamentally different from differential splicing — single-sample-vs-cohort outlier detection rather than two-group comparison. Standard tool in EU rare-disease (Solve-RD) and NIH UDN programs. Use when applying RNA-seq to undiagnosed Mendelian disease, validating predicted splice variants in clinical samples, or detecting cryptic splicing in disease tissue.

- **来源：**[GPTomics/bioSkills contributors; see retained LICENSE @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alternative-splicing/outlier-splicing-detection)
- **许可证：**`MIT`
- **运行依赖：**命令 `Rscript`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info outlier-splicing-detection
dnakit skills install outlier-splicing-detection --agent claude-code   # 或 --agent codex / cursor
```

## `single-cell-cnv-inference` {#single-cell-cnv-inference}

Infer large-scale copy-number alterations from tumor single-cell or single-nucleus RNA-seq to separate malignant from normal cells and call subclones, using inferCNV, copyKAT, Numbat, and SCEVAN. Use when separating malignant from normal cells in a tumor scRNA-seq dataset, inferring chromosome-arm CNVs or aneuploidy from expression, calling tumor subclones from single cells, choosing a CNV-inference method (reference-based vs reference-free, expression-only vs allele-aware), or deciding which cells are tumor before downstream analysis.

- **来源：**[GPTomics/bioSkills contributors; see retained LICENSE @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/single-cell/cnv-inference)
- **许可证：**`MIT`
- **运行依赖：**命令 `Rscript`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info single-cell-cnv-inference
dnakit skills install single-cell-cnv-inference --agent claude-code   # 或 --agent codex / cursor
```

## `splice-variant-prediction` {#splice-variant-prediction}

Predicts whether a DNA variant alters mRNA splicing using sequence-based deep-learning tools — SpliceAI (10kb context dilated CNN, clinical default), Pangolin (multi-tissue), MMSplice (modular per-region CNN with calibrated ΔPSI), SpliceTransformer/TrASPr (tissue-aware transformers), SpliceVault (empirical 300K-RNA lookup of likely mis-splicing outcomes), CADD-Splice (composite score). Applies the ClinGen SVI 2023 framework for ACMG/AMP variant interpretation (PVS1, PP3, BP4 evidence codes), HGVS splicing nomenclature (c.123+1G&gt;A, c.123-3T&gt;G, r.spl?), extended-window scoring for deep-intronic pseudoexons, tissue-specific predictions, branchpoint variant detection (BPHunter, LaBranchoR), and splice-switching ASO design. Use when interpreting splice impact of clinical variants, prioritizing VUS, identifying deep-intronic pathogenic variants, or designing ASOs.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/alternative-splicing/splice-variant-prediction)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `SpliceAI`, `dnakit`
- **别名：**`bio-splice-variant-prediction`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info splice-variant-prediction
dnakit skills install splice-variant-prediction --agent claude-code   # 或 --agent codex / cursor
```

## `tiledbvcf` {#tiledbvcf}

Efficient storage and retrieval of genomic variant data using TileDB. Scalable VCF/BCF ingestion, incremental sample addition, compressed storage, parallel queries, and export capabilities for population genomics.

- **来源：**[Jeremy Leipzig; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/tiledbvcf)
- **许可证：**`MIT`
- **别名：**`alterlab-tiledbvcf`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info tiledbvcf
dnakit skills install tiledbvcf --agent claude-code   # 或 --agent codex / cursor
```

## `variant-analysis` {#variant-analysis}

VCF and variant analysis — parsing, annotation, classification (synonymous, missense, frameshift, stop_gained), VAF filtering, coding vs non-coding categorization, multi-condition variant comparison. Use for VCF parsing, variant fraction calculations (denominator = coding subset only, NOT all variants), and per-sample mutation profiling.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-variant-analysis)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `dnakit`, `pandas`
- **别名：**`tooluniverse-variant-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-analysis
dnakit skills install variant-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `variant-annotate` {#variant-annotate}

Per-rsID variant annotation against gnomAD (v4, GRCh38) and GTEx (v8). For a list of rsIDs, fetches the most-severe coding consequence + gene + per-population allele frequency from gnomAD, plus the top single-tissue eQTL and sQTL from GTEx. Designed for the post-fine-mapping "what is this SNP doing" lookup that researchers normally do by hand at gnomad.broadinstitute.org and gtexportal.org. Caches per-rsID JSON under ~/.cache/variant-annotate/ so repeated runs are free. Use when a user has a credible set, a GWAS top-hit list, or any rsID list and wants coding/regulatory context attached.

- **来源：**[Copyright 2026 Stephen Dorn; statgen-skills contributors @ `dc5f0c6`](https://github.com/svdorn/statgen-skills/tree/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/variant-annotate)
- **许可证：**`MIT`
- **运行依赖：**命令 `python3`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-annotate
dnakit skills install variant-annotate --agent claude-code   # 或 --agent codex / cursor
```

## `variant-annotation` {#variant-annotation}

Annotates VCF variants with functional consequences, population frequencies, and pathogenicity scores using bcftools annotate/csq, Ensembl VEP, SnpEff, and ANNOVAR. Use when deciding which annotation engine and version to pin, which transcript set to report on (RefSeq vs Ensembl vs MANE Select/Plus Clinical, and why VEP --pick is dangerous clinically), how to reconcile HGVS 3'-shifting with VCF left-alignment, which consequence plus NMD status governs PVS1 eligibility, which single calibrated predictor to use for PP3/BP4 (REVEL, AlphaMissense, CADD, SpliceAI deltas), or how to read gnomAD v2/v3/v4 grpmax filtering allele frequency instead of one global AF cutoff. Not for ACMG combining rules or final classification (see variant-calling/clinical-interpretation).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/variant-annotation)
- **许可证：**`MIT`
- **别名：**`bio-variant-annotation`, `ensembl-vep-variant-annotator`, `genetic-variant-annotation`, `snpeff-variant-annotation`, `tooluniverse-variant-functional-annotation`, `vcf-annotator`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-annotation
dnakit skills install variant-annotation --agent claude-code   # 或 --agent codex / cursor
```

## `variant-calling` {#variant-calling}

Call germline SNPs and indels from a BAM/CRAM with bcftools mpileup and call, and select the right calling engine for the job. Use when generating a VCF from aligned reads, choosing between bcftools, GATK HaplotypeCaller, DeepVariant, and DRAGEN, setting ploidy for haploid/organelle/polyploid/sex-chromosome calling, or deciding whether pileup-based calling is good enough versus a local-reassembly caller for indels and difficult regions. Not for cohort joint genotyping (see variant-calling/joint-calling), GATK-specific workflows (see variant-calling/gatk-variant-calling), deep-learning calling (see variant-calling/deepvariant), or somatic/low-VAF detection.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/variant-calling)
- **许可证：**`MIT`
- **运行依赖：**命令 `bcftools`
- **别名：**`bio-variant-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-calling
dnakit skills install variant-calling --agent claude-code   # 或 --agent codex / cursor
```

## `variant-calling-deepvariant` {#variant-calling-deepvariant}

Calls germline SNPs and indels with Google DeepVariant, which reframes variant calling as CNN image classification over multi-channel pileup tensors. Covers platform-specific model selection (WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA), one-shot run_deepvariant vs the three-stage make_examples/call_variants/postprocess_variants pipeline, GPU acceleration of call_variants, DeepTrio for family/trio and de-novo calling, and joint genotyping of gVCFs with GLnexus (not GenotypeGVCFs). Use when deciding DeepVariant vs GATK vs DRAGEN, picking the right --model_type for a sequencing platform, avoiding post-hoc GATK hard filters or BQSR that degrade CNN calls, calling de-novo variants in a trio, merging a DeepVariant cohort, or weighing GIAB-trained benchmark accuracy before clinical deployment.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/deepvariant)
- **许可证：**`MIT`
- **运行依赖：**命令 `DeepVariant`
- **别名：**`bio-variant-calling-deepvariant`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-calling-deepvariant
dnakit skills install variant-calling-deepvariant --agent claude-code   # 或 --agent codex / cursor
```

## `variant-calling-filtering-best-practices` {#variant-calling-filtering-best-practices}

Filters germline and somatic variant callsets at the site and genotype level with GATK VQSR (VQSLOD, truth-sensitivity tranches), VETS/ScoreVariantAnnotations, NVScoreVariants, hard filters with per-annotation thresholds, and bcftools/cyvcf2 expressions, plus Ti/Tv-based QC. Use when deciding between VQSR, hard filtering, and ML recalibration by cohort size and platform, setting SNP vs indel thresholds, replicating the missing-annotation-passes rule so hom-alt sites survive, applying genotype-level GQ/DP filters, or validating filter impact. Not for VCF normalization (see variant-calling/variant-normalization) or summary statistics (see variant-calling/vcf-statistics).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/filtering-best-practices)
- **许可证：**`MIT`
- **别名：**`bio-variant-calling-filtering-best-practices`, `vcf-variant-filtering`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-calling-filtering-best-practices
dnakit skills install variant-calling-filtering-best-practices --agent claude-code   # 或 --agent codex / cursor
```

## `variant-calling-joint-calling` {#variant-calling-joint-calling}

Joint genotype a cohort of per-sample gVCFs with GATK (HaplotypeCaller -ERC GVCF -&gt; GenomicsDBImport or CombineGVCFs -&gt; GenotypeGVCFs) or GLnexus for DeepVariant gVCFs, producing a squared-off sample-by-site genotype matrix. Use when deciding between joint genotyping and merging single-sample callsets (never bcftools merge as absent==hom-ref), choosing GenomicsDBImport vs CombineGVCFs by cohort size and memory, solving the N+1 problem so a new sample does not force re-calling everyone, understanding cohort rescue of low-coverage het sites, handling the spanning-deletion star allele and GQ/PL recomputation at the joint step, scaling to biobank cohorts by interval sharding, or picking DeepVariant+GLnexus over the GATK path on throughput. Not for single-sample calling (see variant-calling/gatk-variant-calling) or VQSR/hard-filter mechanism (see variant-calling/filtering-best-practices).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/joint-calling)
- **许可证：**`MIT`
- **运行依赖：**命令 `GATK`
- **别名：**`bio-variant-calling-joint-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-calling-joint-calling
dnakit skills install variant-calling-joint-calling --agent claude-code   # 或 --agent codex / cursor
```

## `variant-calling-structural-variant-calling` {#variant-calling-structural-variant-calling}

Call structural variants (&gt;=50 bp deletions, insertions, inversions, duplications, translocations) from short- or long-read data by reconstructing four orthogonal signals (discordant pairs, split reads via the SA tag, read depth, local assembly). Covers Manta, DELLY, LUMPY/smoove, GRIDSS2, SvABA for short reads and Sniffles2, cuteSV, pbsv, dipcall/PAV for long reads, each mapped to the signals it fuses and the blind spots that follow. Use when choosing an SV caller from its signal set and failure modes, decoding the SVLEN-sign / symbolic-vs-BND / CIPOS VCF representation minefield, force-genotyping a cohort matrix instead of unioning discovery VCFs, merging populations with sequence-aware Truvari vs position-only SURVIVOR, parameterizing a Truvari benchmark, or deciding when short-read insertion recall forces a switch to long reads. Not for pure copy-number dosage (see copy-number/cnvkit-analysis).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/structural-variant-calling)
- **许可证：**`MIT`
- **运行依赖：**命令 `manta`
- **别名：**`bio-variant-calling-structural-variant-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-calling-structural-variant-calling
dnakit skills install variant-calling-structural-variant-calling --agent claude-code   # 或 --agent codex / cursor
```

## `variant-normalization` {#variant-normalization}

Left-align and trim indels to parsimonious canonical form, decompose MNPs (atomize), and split multiallelic variants with bcftools norm. Use when comparing variants across callers or cohorts, preparing a VCF for database annotation or ClinVar/dbSNP matching, merging VCFs, reconciling vt-vs-bcftools representation discordance, or resolving the VCF-left-align vs HGVS-3'-rule clash.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/variant-normalization)
- **许可证：**`MIT`
- **运行依赖：**命令 `bcftools`
- **别名：**`bio-variant-normalization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-normalization
dnakit skills install variant-normalization --agent claude-code   # 或 --agent codex / cursor
```

## `variant-predictor-dms-validation` {#variant-predictor-dms-validation}

Validate a variant-effect predictor (AlphaMissense, ESM-C SAE, ESM logits, EVE, conservation scores, or any per-variant numeric score) against experimental deep mutational scanning (DMS) data. Computes per-variant predictor scores, splits variants into neutral vs disruptive groups by DMS effect, runs a Mann-Whitney U test on the predictor scores, and sweeps the stratification thresholds for robustness. Use when you need to know whether a predictor's scores track real functional disruption on a specific protein.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-variant-predictor-dms-validation)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-variant-predictor-dms-validation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-predictor-dms-validation
dnakit skills install variant-predictor-dms-validation --agent claude-code   # 或 --agent codex / cursor
```

## `vcf-basics` {#vcf-basics}

View, query, and interpret VCF/BCF variant files with bcftools and cyvcf2. Use when inspecting variants, extracting fields with query format strings, converting VCF/BCF, or correctly reading a field -- QUAL (site) vs GQ (genotype) vs PL/GL likelihoods, AD vs DP and allele balance, GT phasing/ploidy/PS and missing-vs-hom-ref, INFO/FORMAT Number A/R/G semantics, symbolic alleles (&lt;DEL&gt;, &lt;NON_REF&gt;, spanning *) and END, or telling a raw gVCF apart from a filtered callset.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/vcf-basics)
- **许可证：**`MIT`
- **运行依赖：**命令 `bcftools`
- **别名：**`bio-vcf-basics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info vcf-basics
dnakit skills install vcf-basics --agent claude-code   # 或 --agent codex / cursor
```

## `vcf-manipulation` {#vcf-manipulation}

Combine, split, sort, intersect, and subset VCF/BCF files with bcftools merge, concat, isec, sort, view, and reheader. Use when merging different samples into a cohort VCF, concatenating per-chromosome or per-region call sets for the same samples, intersecting or complementing call sets from different callers, subsetting samples/regions, harmonizing sample names and

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/vcf-manipulation)
- **许可证：**`MIT`
- **运行依赖：**命令 `bcftools`
- **别名：**`bcftools-variant-manipulation`, `bio-vcf-manipulation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info vcf-manipulation
dnakit skills install vcf-manipulation --agent claude-code   # 或 --agent codex / cursor
```

## `vcf-statistics` {#vcf-statistics}

Compute and interpret VCF quality-control metrics (Ti/Tv, het/hom, novel/known, missingness, HWE, contamination, relatedness) with bcftools stats, vcftools, plot-vcfstats, and identity tools (somalier, peddy, KING). Use when judging whether a callset is trustworthy, diagnosing a low Ti/Tv or outlier het/hom sample, deciding whether an HWE deviation is error or biology, screening a cohort for sample swaps/contamination/wrong-sex before analysis, or comparing call sets before and after filtering. Not for applying filters (see variant-calling/filtering-best-practices) or normalizing representation (see variant-calling/variant-normalization).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/variant-calling/vcf-statistics)
- **许可证：**`MIT`
- **运行依赖：**命令 `bcftools`
- **别名：**`bio-vcf-statistics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info vcf-statistics
dnakit skills install vcf-statistics --agent claude-code   # 或 --agent codex / cursor
```
