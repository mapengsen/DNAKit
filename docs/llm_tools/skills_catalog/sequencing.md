# 测序与质控

本类共 **20** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`amplicon-qc`](#amplicon-qc) | clawbio | 16S/18S rRNA amplicon preprocessing — from raw paired-end FASTQ files through N removal and primer trimming, producing outputs ready for DADA2 quality filterin… |
| [`busco-status-interpretation`](#busco-status-interpretation) | sciagent | Guide to interpreting BUSCO completeness statuses: why Duplicated BUSCOs count as complete, parsing output files, computing/comparing completeness across prote… |
| [`instrument-data-to-allotrope`](#instrument-data-to-allotrope) | user335-anthropic | Convert laboratory instrument output files (PDF, CSV, Excel, TXT) to Allotrope Simple Model (ASM) JSON format or flattened 2D CSV. Use this skill when scientis… |
| [`long-read-sequencing-basecalling`](#long-read-sequencing-basecalling) | bioskills | Basecalls raw Oxford Nanopore signal (POD5/FAST5) into reads with Dorado, choosing the chemistry-matched model and accuracy tier (fast/hac/sup), requesting mod… |
| [`long-read-sequencing-clair3-variants`](#long-read-sequencing-clair3-variants) | bioskills | Calls germline small variants (SNPs and indels) from Oxford Nanopore and PacBio HiFi long reads with Clair3, a two-stage (pileup + full-alignment) deep-learnin… |
| [`long-read-sequencing-haplotype-phasing`](#long-read-sequencing-haplotype-phasing) | bioskills | Phases small variants, SVs, and methylation from Oxford Nanopore and PacBio long reads (read-backed/physical phasing) with WhatsHap, LongPhase, or HiPhase, and… |
| [`long-read-sequencing-long-read-alignment`](#long-read-sequencing-long-read-alignment) | bioskills | Aligns Oxford Nanopore and PacBio long reads (and assemblies) to a reference with minimap2 using the error-rate-matched preset (map-ont, lr:hq, map-hifi, map-p… |
| [`long-read-sequencing-long-read-qc`](#long-read-sequencing-long-read-qc) | bioskills | Assesses Oxford Nanopore and PacBio long-read quality with NanoPlot, cramino, NanoComp, pycoQC/toulligQC, and seqkit, and filters reads with chopper/Filtlong f… |
| [`long-read-sequencing-medaka-polishing`](#long-read-sequencing-medaka-polishing) | bioskills | Polishes Oxford Nanopore draft assemblies to higher consensus accuracy with medaka, a basecaller-model-specific neural consensus net, produces haploid variant… |
| [`long-read-sequencing-nanopore-methylation`](#long-read-sequencing-nanopore-methylation) | bioskills | Calls DNA base modifications (5mC, 5hmC, 6mA, 4mC) directly from Oxford Nanopore and PacBio HiFi long reads encoded as MM/ML SAM tags, piles them into per-site… |
| [`long-read-sequencing-structural-variants`](#long-read-sequencing-structural-variants) | bioskills | Detects structural variants (deletions, insertions, inversions, duplications, translocations) from Oxford Nanopore and PacBio long-read alignments with Sniffle… |
| [`multiqc-reporter`](#multiqc-reporter) | clawbio | Aggregates QC reports from any bioinformatics tool outputs (FastQC, fastp, STAR, Picard, samtools, etc.) into a single MultiQC HTML report plus a ClawBio markd… |
| [`pod5`](#pod5) | damlab | Inspect, merge, filter, subset, repack, and convert nanopore POD5 files using the pod5 CLI. Use when working with .pod5 files, converting fast5 to pod5 or pod5… |
| [`read-qc-adapter-trimming`](#read-qc-adapter-trimming) | bioskills | Removes sequencing adapters from FASTQ reads with Cutadapt and Trimmomatic, including paired-end read-through, small-RNA 3' adapters, amplicon primers, and anc… |
| [`read-qc-contamination-screening`](#read-qc-contamination-screening) | bioskills | Detects contamination in sequencing reads - cross-species (FastQ Screen, Kraken2), vector/PhiX/adapter, rRNA, and same-species cross-sample/index-hopping and s… |
| [`read-qc-fastp-workflow`](#read-qc-fastp-workflow) | bioskills | Runs all-in-one FASTQ preprocessing with fastp in a single pass - adapter trimming via paired-end overlap analysis, quality/length filtering, 2-color poly-G re… |
| [`read-qc-quality-filtering`](#read-qc-quality-filtering) | bioskills | Filters reads by quality, length, N content, and complexity with Trimmomatic, fastp, and Cutadapt, including sliding-window trimming, per-read unqualified-base… |
| [`read-qc-quality-reports`](#read-qc-quality-reports) | bioskills | Generates and interprets per-file and cross-sample QC reports from FASTQ data with FastQC, falco, and MultiQC, covering Phred quality, per-base composition, GC… |
| [`read-qc-umi-processing`](#read-qc-umi-processing) | bioskills | Extracts UMIs and collapses reads to original molecules with umi_tools (directional dedup) or builds error-corrected single-strand/duplex consensus reads with… |
| [`sample-qc-triage`](#sample-qc-triage) | clawbio | Deterministic multi-sample QC triage for identity, sex, contamination, and batch-shift outliers |

## `amplicon-qc` {#amplicon-qc}

16S/18S rRNA amplicon preprocessing — from raw paired-end FASTQ files through N removal and primer trimming, producing outputs ready for DADA2 quality filtering and denoising. Deliberately stops before quality-filtering decisions that require researcher judgment.

- **来源：**[Zabiulla; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/claw-amplicon-qc)
- **许可证：**`MIT`
- **运行依赖：**命令 `Rscript`, `cutadapt`, `seqkit`
- **别名：**`claw-amplicon-qc`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info amplicon-qc
dnakit skills install amplicon-qc --agent claude-code   # 或 --agent codex / cursor
```

## `busco-status-interpretation` {#busco-status-interpretation}

Guide to interpreting BUSCO completeness statuses: why Duplicated BUSCOs count as complete, parsing output files, computing/comparing completeness across proteomes/genomes, common counting mistakes. Use when running BUSCO QC, comparing assemblies, or reporting completeness. See also: prokka-genome-annotation for annotation workflows feeding BUSCO.

- **来源：**[jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/genomics-bioinformatics/qc/busco-status-interpretation)
- **许可证：**`CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info busco-status-interpretation
dnakit skills install busco-status-interpretation --agent claude-code   # 或 --agent codex / cursor
```

## `instrument-data-to-allotrope` {#instrument-data-to-allotrope}

Convert laboratory instrument output files (PDF, CSV, Excel, TXT) to Allotrope Simple Model (ASM) JSON format or flattened 2D CSV. Use this skill when scientists need to standardize instrument data for LIMS systems, data lakes, or downstream analysis. Supports auto-detection of instrument types. Outputs include full ASM JSON, flattened CSV for easy import, and exportable Python code for data engineers. Common triggers include converting instrument files, standardizing lab data, preparing data for upload to LIMS/ELN systems, or generating parser code for production pipelines.

- **来源：**[Anthropic and life-sciences skill contributors; leaf LICENSE.txt retained @ `e96556b`](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c/instrument-data-to-allotrope)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `allotropy`; 命令 `python3`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info instrument-data-to-allotrope
dnakit skills install instrument-data-to-allotrope --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-basecalling` {#long-read-sequencing-basecalling}

Basecalls raw Oxford Nanopore signal (POD5/FAST5) into reads with Dorado, choosing the chemistry-matched model and accuracy tier (fast/hac/sup), requesting modified bases (5mCG_5hmCG, 6mA, m6A) at basecall time, and handling duplex, demultiplexing, trimming, and HERRO read correction. Covers why the model+version is an irreversible analysis decision, why methylation cannot be recovered later, and why downstream polish/variant models must match the basecaller. Use when converting POD5/FAST5 to reads, picking a Dorado model for R9/R10 or RNA004, enabling methylation calling, basecalling duplex, demultiplexing barcoded runs, or correcting reads for assembly.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/basecalling)
- **许可证：**`MIT`
- **运行依赖：**命令 `dorado`
- **别名：**`bio-long-read-sequencing-basecalling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-basecalling
dnakit skills install long-read-sequencing-basecalling --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-clair3-variants` {#long-read-sequencing-clair3-variants}

Calls germline small variants (SNPs and indels) from Oxford Nanopore and PacBio HiFi long reads with Clair3, a two-stage (pileup + full-alignment) deep-learning caller, selecting the chemistry- and basecaller-version-matched model, enabling read-based phasing, and benchmarking against GIAB with stratification. Covers why the model string is the experiment (no auto-detection, silent degradation on mismatch), why ONT homopolymer/STR indels are the residual error whole-genome F1 hides, and the somatic/trio/RNA boundary to the ClairS/Clair3-Trio family. Use when calling germline SNVs/indels from ONT or HiFi BAMs, choosing a Clair3 model, phasing variants, or benchmarking long-read calls.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/clair3-variants)
- **许可证：**`MIT`
- **运行依赖：**命令 `Clair3`
- **别名：**`bio-long-read-sequencing-clair3-variants`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-clair3-variants
dnakit skills install long-read-sequencing-clair3-variants --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-haplotype-phasing` {#long-read-sequencing-haplotype-phasing}

Phases small variants, SVs, and methylation from Oxford Nanopore and PacBio long reads (read-backed/physical phasing) with WhatsHap, LongPhase, or HiPhase, and haplotags the BAM (HP/PS tags) for allele-resolved downstream analysis. Covers why phase blocks break at het-sparse gaps (read length x heterozygosity), why phasing the VCF is useless until the BAM is haplotagged, the GT-pipe/PS and read HP/PS tag spec, reporting block N50 with switch error, the diploid-assumption/CNV/haploid-region traps, trio phasing as the gold standard, and the boundary to statistical panel phasing. Use when phasing long-read variants, haplotagging reads for allele-specific methylation/expression or phased SVs, choosing WhatsHap vs LongPhase vs HiPhase, trio phasing, or assessing phasing quality.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/haplotype-phasing)
- **许可证：**`MIT`
- **运行依赖：**命令 `whatshap`
- **别名：**`bio-long-read-sequencing-haplotype-phasing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-haplotype-phasing
dnakit skills install long-read-sequencing-haplotype-phasing --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-long-read-alignment` {#long-read-sequencing-long-read-alignment}

Aligns Oxford Nanopore and PacBio long reads (and assemblies) to a reference with minimap2 using the error-rate-matched preset (map-ont, lr:hq, map-hifi, map-pb, splice/splice:hq, asm5/10/20, ava), producing a sorted/indexed BAM for variant, SV, methylation, or isoform analysis. Covers why the preset rewrites the scoring/chaining model, why SV calling rides on supplementary not secondary alignments, carrying MM/ML methylation tags through with -y, the multi-part-index MAPQ trap, and when to swap in Winnowmap/VACmap/lra/pbmm2. Use when mapping ONT or PacBio reads, choosing a minimap2 preset by platform/chemistry, preparing input for Clair3/medaka/Sniffles/modkit, aligning into repeats/centromeres, or spliced-aligning cDNA/Iso-Seq.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/long-read-alignment)
- **许可证：**`MIT`
- **运行依赖：**命令 `minimap2`
- **别名：**`bio-long-read-sequencing-long-read-alignment`, `minimap2`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-long-read-alignment
dnakit skills install long-read-sequencing-long-read-alignment --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-long-read-qc` {#long-read-sequencing-long-read-qc}

Assesses Oxford Nanopore and PacBio long-read quality with NanoPlot, cramino, NanoComp, pycoQC/toulligQC, and seqkit, and filters reads with chopper/Filtlong for the downstream goal. Covers why read-only Qscore is an uncalibrated posterior (real accuracy needs a reference BAM), why the sequencing_summary.txt is required for run-health metrics, intent-conditioned filtering (preserve long reads and small replicons for assembly, filter almost nothing for variant calling), the chimera/internal-adapter trap that fabricates SVs, and PacBio rq-based HiFi QC. Use when judging a long-read run, computing read N50 or percent identity, filtering reads before assembly or variant calling, comparing barcodes/runs, or reading run-health red flags.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/long-read-qc)
- **许可证：**`MIT`
- **运行依赖：**命令 `nanoplot`
- **别名：**`bio-long-read-sequencing-long-read-qc`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-long-read-qc
dnakit skills install long-read-sequencing-long-read-qc --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-medaka-polishing` {#long-read-sequencing-medaka-polishing}

Polishes Oxford Nanopore draft assemblies to higher consensus accuracy with medaka, a basecaller-model-specific neural consensus net, produces haploid variant calls (VCF) for microbial, mitochondrial, or viral samples, and generates amplicon/viral consensus sequences. Covers the model-matching footgun that silently degrades output, why Racon-first is obsolete and medaka runs directly on Flye output as a single pass, why HiFi must never be fed to medaka, the v1-&gt;v2 subcommand renames, and the precise medaka_variant deprecation. Use when polishing an ONT-only assembly, generating an amplicon/viral consensus, calling a haploid ONT consensus, or deciding whether medaka, dorado polish, or Clair3 is the right tool.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/medaka-polishing)
- **许可证：**`MIT`
- **运行依赖：**命令 `medaka`
- **别名：**`bio-long-read-sequencing-medaka-polishing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-medaka-polishing
dnakit skills install long-read-sequencing-medaka-polishing --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-nanopore-methylation` {#long-read-sequencing-nanopore-methylation}

Calls DNA base modifications (5mC, 5hmC, 6mA, 4mC) directly from Oxford Nanopore and PacBio HiFi long reads encoded as MM/ML SAM tags, piles them into per-site bedMethyl with modkit (or pb-CpG-tools for PacBio), and produces phased allele-specific methylation. Covers why methylation is a basecalling decision that cannot be recovered later, the MM/ML tag-drop failure that silently zeroes methylation through alignment, the MM ? vs . no-call semantics, 5mC/5hmC resolution vs bisulfite, modkit's 10th-percentile auto-threshold, and the haplotagged ASM workflow. Use when calling 5mC/5hmC/6mA from a modBAM, generating bedMethyl, preserving methylation tags through alignment, doing allele-specific or differential methylation, or QC-ing a modification BAM.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/nanopore-methylation)
- **许可证：**`MIT`
- **运行依赖：**命令 `modkit`
- **别名：**`bio-long-read-sequencing-nanopore-methylation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-nanopore-methylation
dnakit skills install long-read-sequencing-nanopore-methylation --agent claude-code   # 或 --agent codex / cursor
```

## `long-read-sequencing-structural-variants` {#long-read-sequencing-structural-variants}

Detects structural variants (deletions, insertions, inversions, duplications, translocations) from Oxford Nanopore and PacBio long-read alignments with Sniffles2, cuteSV, SVIM, and assembly-based callers, joint-genotypes cohorts via the Sniffles2 .snf workflow, and benchmarks with Truvari against GIAB. Covers why an SV call is a representation artifact (the tandem-repeat BED, aligner, and Truvari params set precision/recall as much as the caller), the cuteSV per-platform parameter trap, soft-clipped supplementary alignments as the SV substrate, and the somatic/mosaic boundary to Severus/nanomonsv. Use when calling germline or somatic SVs from ONT/HiFi reads, joint-genotyping a cohort, choosing or tuning an SV caller, or benchmarking SV calls.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/long-read-sequencing/structural-variants)
- **许可证：**`MIT`
- **运行依赖：**命令 `sniffles`
- **别名：**`bio-long-read-sequencing-structural-variants`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info long-read-sequencing-structural-variants
dnakit skills install long-read-sequencing-structural-variants --agent claude-code   # 或 --agent codex / cursor
```

## `multiqc-reporter` {#multiqc-reporter}

Aggregates QC reports from any bioinformatics tool outputs (FastQC, fastp, STAR, Picard, samtools, etc.) into a single MultiQC HTML report plus a ClawBio markdown summary with per-sample QC metrics.

- **来源：**[Cameron Lloyd; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/multiqc-reporter)
- **许可证：**`MIT`
- **运行依赖：**命令 `multiqc`, `python3`
- **别名：**`multiqc`, `multiqc-qc-reports`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info multiqc-reporter
dnakit skills install multiqc-reporter --agent claude-code   # 或 --agent codex / cursor
```

## `pod5` {#pod5}

Inspect, merge, filter, subset, repack, and convert nanopore POD5 files using the pod5 CLI. Use when working with .pod5 files, converting fast5 to pod5 or pod5 to fast5, extracting read summaries, filtering reads by read_id, subsetting by barcode or other metadata columns, merging multiple pod5 files, or updating pod5 files to the latest schema. Triggers on tasks involving .pod5, .fast5 conversion, nanopore raw signal files, read_id filtering, basecall input preparation, or Oxford Nanopore sequencing data.

- **来源：**[DamLabResources/damlab-skills contributors @ `4ad096e`](https://github.com/DamLabResources/damlab-skills/tree/4ad096e3b9b0116271d10bdb746641b60a1b2211/skills/pod5)
- **许可证：**`MIT`
- **运行依赖：**命令 `pod5`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info pod5
dnakit skills install pod5 --agent claude-code   # 或 --agent codex / cursor
```

## `read-qc-adapter-trimming` {#read-qc-adapter-trimming}

Removes sequencing adapters from FASTQ reads with Cutadapt and Trimmomatic, including paired-end read-through, small-RNA 3' adapters, amplicon primers, and anchored/linked adapters. Use when FastQC shows adapter content climbing toward the 3' end, when inserts are shorter than the read length (small-RNA, cfDNA, FFPE), or before assembly/k-mer analysis. For all-in-one trimming use fastp-workflow; for quality/length filtering use quality-filtering.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-qc/adapter-trimming)
- **许可证：**`MIT`
- **运行依赖：**命令 `cutadapt`
- **别名：**`bio-read-qc-adapter-trimming`, `cutadapt`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-qc-adapter-trimming
dnakit skills install read-qc-adapter-trimming --agent claude-code   # 或 --agent codex / cursor
```

## `read-qc-contamination-screening` {#read-qc-contamination-screening}

Detects contamination in sequencing reads - cross-species (FastQ Screen, Kraken2), vector/PhiX/adapter, rRNA, and same-species cross-sample/index-hopping and sample swaps (SNP fingerprints via verifyBamID2/NGSCheckMate/somalier). Use when suspecting cross-contamination, PDX host reads, microbial carry-over, or sample swaps, and to decide whether to report, filter, or align to a combined reference. For deep taxonomic profiling use metagenomics/kraken-classification.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-qc/contamination-screening)
- **许可证：**`MIT`
- **运行依赖：**命令 `fastq_screen`
- **别名：**`bio-read-qc-contamination-screening`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-qc-contamination-screening
dnakit skills install read-qc-contamination-screening --agent claude-code   # 或 --agent codex / cursor
```

## `read-qc-fastp-workflow` {#read-qc-fastp-workflow}

Runs all-in-one FASTQ preprocessing with fastp in a single pass - adapter trimming via paired-end overlap analysis, quality/length filtering, 2-color poly-G removal, base correction, optional dedup/UMI/merge, and HTML/JSON reports. Use when preprocessing bulk Illumina data and wanting one fast tool instead of separate Cutadapt, Trimmomatic, and FastQC steps. For precise small-RNA/amplicon adapters use adapter-trimming; for molecule-accurate UMI dedup use umi-processing.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-qc/fastp-workflow)
- **许可证：**`MIT`
- **运行依赖：**命令 `fastp`
- **别名：**`bio-read-qc-fastp-workflow`, `fastp`, `fastp-fastq-preprocessing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-qc-fastp-workflow
dnakit skills install read-qc-fastp-workflow --agent claude-code   # 或 --agent codex / cursor
```

## `read-qc-quality-filtering` {#read-qc-quality-filtering}

Filters reads by quality, length, N content, and complexity with Trimmomatic, fastp, and Cutadapt, including sliding-window trimming, per-read unqualified-base filtering, and 2-color poly-G removal. Use when reads have poor-quality tails, when an assembly or k-mer workflow needs clean input, or when a junk read subpopulation must be dropped. For adapter removal use adapter-trimming; for all-in-one preprocessing use fastp-workflow.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-qc/quality-filtering)
- **许可证：**`MIT`
- **运行依赖：**命令 `trimmomatic`
- **别名：**`bio-read-qc-quality-filtering`, `trimmomatic`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-qc-quality-filtering
dnakit skills install read-qc-quality-filtering --agent claude-code   # 或 --agent codex / cursor
```

## `read-qc-quality-reports` {#read-qc-quality-reports}

Generates and interprets per-file and cross-sample QC reports from FASTQ data with FastQC, falco, and MultiQC, covering Phred quality, per-base composition, GC, duplication, overrepresented sequences, and adapter content. Use when performing initial QC on raw sequencing reads, validating preprocessing, or judging a multi-sample cohort for outliers and batch effects. For long reads use NanoPlot; for adapter/quality remediation route to adapter-trimming, quality-filtering, or fastp-workflow.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-qc/quality-reports)
- **许可证：**`MIT`
- **运行依赖：**命令 `fastqc`
- **别名：**`bio-read-qc-quality-reports`, `bioqc-mcp`, `fastqc`, `ngs-quality-control`, `tooluniverse-fastq-qc`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-qc-quality-reports
dnakit skills install read-qc-quality-reports --agent claude-code   # 或 --agent codex / cursor
```

## `read-qc-umi-processing` {#read-qc-umi-processing}

Extracts UMIs and collapses reads to original molecules with umi_tools (directional dedup) or builds error-corrected single-strand/duplex consensus reads with fgbio. Use when the library has UMIs and accurate molecule counting or below-sequencer-floor error correction is needed - single-cell, low-input RNA-seq, targeted panels, and ctDNA/liquid-biopsy rare-variant detection. For UMI extraction during QC use fastp-workflow; do not dedup non-UMI bulk RNA-seq.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/read-qc/umi-processing)
- **许可证：**`MIT`
- **运行依赖：**命令 `umi_tools`
- **别名：**`bio-read-qc-umi-processing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info read-qc-umi-processing
dnakit skills install read-qc-umi-processing --agent claude-code   # 或 --agent codex / cursor
```

## `sample-qc-triage` {#sample-qc-triage}

Deterministic multi-sample QC triage for identity, sex, contamination, and batch-shift outliers

- **来源：**[ClawBio; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/sample-qc-triage)
- **许可证：**`MIT`
- **运行依赖：**命令 `python3`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info sample-qc-triage
dnakit skills install sample-qc-triage --agent claude-code   # 或 --agent codex / cursor
```
