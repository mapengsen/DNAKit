# 调控与表观基因组

本类共 **66** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`accessibility-aggregation`](#accessibility-aggregation) | encode-toolkit | Build comprehensive chromatin accessibility maps by aggregating ATAC-seq and DNase-seq narrowPeak data across multiple ENCODE experiments, donors, and labs. Us… |
| [`atac-seq-allele-specific-accessibility`](#atac-seq-allele-specific-accessibility) | bioskills | Detect allele-specific chromatin accessibility from ATAC-seq using WASP, GATK ASEReadCounter, or RASQUAL. Use when mapping cis-regulatory genetic variants from… |
| [`atac-seq-atac-peak-calling`](#atac-seq-atac-peak-calling) | bioskills | Call accessible chromatin regions from ATAC-seq BAM files using MACS3, MACS2, Genrich, or HMMRATAC. Use when identifying open chromatin from aligned ATAC-seq,… |
| [`atac-seq-atac-qc`](#atac-seq-atac-qc) | bioskills | ATAC-seq library quality control -- TSS enrichment, FRiP, fragment-size periodicity, library complexity (NRF/PBC1/PBC2), mitochondrial fraction, and ENCODE 4 t… |
| [`atac-seq-co-accessibility`](#atac-seq-co-accessibility) | bioskills | Infer cis-regulatory connections (peak-to-peak co-accessibility) from scATAC-seq using Cicero, ArchR getCoAccessibility, or SCENIC+. Use when linking enhancer… |
| [`atac-seq-consensus-peakset`](#atac-seq-consensus-peakset) | bioskills | Build a differential-ready consensus peakset from per-replicate ATAC-seq peaks using iterative overlap removal, fixed-width re-centering, and majority-rule ove… |
| [`atac-seq-deep-learning-atac`](#atac-seq-deep-learning-atac) | bioskills | Sequence-based deep learning for ATAC-seq using chromBPNet, BPNet, scBasset, or Enformer. Use when correcting Tn5 bias with neural networks beyond k-mer models… |
| [`atac-seq-differential-accessibility`](#atac-seq-differential-accessibility) | bioskills | Identify differentially accessible chromatin regions across conditions using DiffBind, csaw, DESeq2, or edgeR. Use when comparing ATAC-seq accessibility betwee… |
| [`atac-seq-enhancer-gene-linking`](#atac-seq-enhancer-gene-linking) | bioskills | Predict enhancer-gene regulatory connections from ATAC-seq using ABC, ENCODE-rE2G, HiChIP, or Cicero. Use when linking distal enhancers to target genes, choosi… |
| [`atac-seq-footprinting`](#atac-seq-footprinting) | bioskills | Detect transcription factor binding footprints in ATAC-seq using TOBIAS, HINT-ATAC, Wellington, or scprinter. Use when identifying bound TF sites within access… |
| [`atac-seq-motif-deviation`](#atac-seq-motif-deviation) | bioskills | Analyze TF motif accessibility variability across samples or single cells using chromVAR. Use when identifying TF motifs whose accessibility correlates with co… |
| [`atac-seq-nucleosome-positioning`](#atac-seq-nucleosome-positioning) | bioskills | Map nucleosome center positions, occupancy, and fuzziness from ATAC-seq fragment-size patterns using NucleoATAC, ATACseqQC, DANPOS3, or scprinter. Use when cha… |
| [`atac-seq-single-cell-atac`](#atac-seq-single-cell-atac) | bioskills | Process and analyze single-cell ATAC-seq data with Signac, ArchR, SnapATAC2, or Cell Ranger ATAC. Use when handling 10X scATAC or 10X Multiome (paired RNA+ATAC… |
| [`chip-atlas-diff-analysis`](#chip-atlas-diff-analysis) | agentic-science | Submit two groups of public ChIP-seq, ATAC-seq, DNase-seq, or Bisulfite-seq experiments to the ChIP-Atlas Diff Analysis API and recover differential peak regio… |
| [`chip-atlas-peak-enrichment`](#chip-atlas-peak-enrichment) | agentic-science | Submit a gene list to the ChIP-Atlas Enrichment Analysis API and recover the transcription factors, chromatin modifiers, and histone marks whose public ChIP-se… |
| [`chip-atlas-target-genes`](#chip-atlas-target-genes) | agentic-science | For a single transcription factor or chromatin regulator, download ChIP-Atlas's pre-computed per-gene MACS2 binding scores across every public ChIP-seq experim… |
| [`chipseq-allele-specific-binding`](#chipseq-allele-specific-binding) | bioskills | Detects allele-specific transcription factor or histone modification binding from heterozygous-variant ChIP-seq using WASP (reference-bias filter; mandatory up… |
| [`chipseq-chip-deep-learning`](#chipseq-chip-deep-learning) | bioskills | Trains and applies base-resolution deep learning models on ChIP-seq / ChIP-nexus / CUT&amp;RUN data. Uses BPNet (Avsec 2021 Nat Genet 53:354; soft motif syntax fro… |
| [`chipseq-chromatin-state-segmentation`](#chipseq-chromatin-state-segmentation) | bioskills | Segments the genome into chromatin states from combinatorial histone modification and chromatin factor ChIP-seq data. Uses ChromHMM (multivariate HMM on binari… |
| [`chipseq-cut-and-run-tag`](#chipseq-cut-and-run-tag) | bioskills | Analyzes CUT&amp;RUN (Skene Henikoff 2017) and CUT&amp;Tag (Kaya-Okur 2019) chromatin profiling data. Handles SEACR vs MACS2 peak calling (with the btaf375 2025 benchm… |
| [`chipseq-differential-binding`](#chipseq-differential-binding) | bioskills | Identifies differentially bound ChIP-seq regions between conditions using DiffBind, csaw (sliding windows), DESeq2/edgeR/PyDESeq2 on count matrices, NormR (con… |
| [`chipseq-motif-analysis`](#chipseq-motif-analysis) | bioskills | Discovers de novo motifs and tests known motif enrichment in ChIP-seq, ATAC-seq, or other peak sequences using HOMER, MEME-ChIP (STREME, CentriMo, TOMTOM, FIMO… |
| [`chipseq-peak-annotation`](#chipseq-peak-annotation) | bioskills | Annotates ChIP-seq peaks to genomic features, nearest genes, ENCODE candidate cis-regulatory elements (cCREs), and regulatory domains. Uses ChIPseeker (R), HOM… |
| [`chipseq-peak-calling`](#chipseq-peak-calling) | bioskills | Calls ChIP-seq peaks with MACS3, MACS2, HOMER, or SPP across narrow (TF) and broad (histone) modes. Handles input control matching, fragment-size modeling vs -… |
| [`chipseq-qc`](#chipseq-qc) | bioskills | Assesses ChIP-seq quality across antibody specificity, fragmentation, enrichment, replicate concordance, and library complexity. Computes FRiP, NSC/RSC (phanto… |
| [`chipseq-spike-in-normalization`](#chipseq-spike-in-normalization) | bioskills | Normalizes ChIP-seq data using exogenous spike-in (ChIP-Rx with Drosophila chromatin per Orlando 2014 / Egan 2016; E. coli carryover for CUT&amp;RUN/CUT&amp;Tag). Dist… |
| [`chipseq-super-enhancers`](#chipseq-super-enhancers) | bioskills | Identifies super-enhancers from H3K27ac, MED1, or BRD4 ChIP-seq using ROSE, ROSE2, LILY, HOMER -style super, and ENCODE dELS cross-referencing. Handles peak st… |
| [`chipseq-visualization`](#chipseq-visualization) | bioskills | Visualizes ChIP-seq data using deepTools (computeMatrix, plotHeatmap, plotProfile, bamCoverage, bamCompare), pyGenomeTracks (modern INI-driven track plots), Gv… |
| [`compare-biosamples`](#compare-biosamples) | encode-toolkit | Compare ENCODE experiments across different biosamples, tissues, or cell lines to identify tissue-specific regulatory patterns. Use when the user wants cross-t… |
| [`data-visualization-genome-tracks`](#data-visualization-genome-tracks) | bioskills | Build genome-browser-style multi-track figures with pyGenomeTracks (config-driven), Gviz (R), and IGV batch screenshotting. Covers BigWig coverage tracks, BED/… |
| [`dnakit-methylation-summary`](#dnakit-methylation-summary) | dnakit | Read modkit extended bedMethyl and compute coverage-weighted methylation summaries with DNAKit. Use for nanopore modification counts, including separate 5mC an… |
| [`epigenomics`](#epigenomics) | tooluniverse | Genomics and epigenomics analysis: DNA methylation (CpG, 5mC, 5hmC, bisulfite, RRBS), m6A RNA modification (MeRIP-seq), ChIP-seq peaks, ATAC-seq accessibility,… |
| [`epigenomics-chromatin`](#epigenomics-chromatin) | tooluniverse | Histone-modification ChIP-seq, ATAC-seq accessibility, chromatin state, and TF binding analysis from ENCODE, Roadmap Epigenomics, ChIP-Atlas. Use for chromatin… |
| [`gene-regulatory-networks`](#gene-regulatory-networks) | tooluniverse | Gene regulatory network analysis — TF-target inference (JASPAR motifs, ChIP-seq), motif scanning, eQTL integration, perturbation evidence (knockout/overexpress… |
| [`gene-regulatory-networks-multiomics-grn`](#gene-regulatory-networks-multiomics-grn) | bioskills | Build enhancer-driven gene regulatory networks (eGRNs) by integrating single-cell RNA-seq and ATAC-seq using SCENIC+, CellOracle base GRNs, Pando, FigR, DIRECT… |
| [`gene-regulatory-networks-scenic-regulons`](#gene-regulatory-networks-scenic-regulons) | bioskills | Infer transcription factor regulons from single-cell RNA-seq with pySCENIC by combining GRNBoost2 co-expression, cisTarget motif-enrichment pruning, and AUCell… |
| [`hi-c-analysis-compartment-analysis`](#hi-c-analysis-compartment-analysis) | bioskills | Detects A/B chromatin compartments from balanced Hi-C contact matrices via eigenvector decomposition of the distance-normalized, Pearson-correlated cis matrix… |
| [`hi-c-analysis-contact-pairs`](#hi-c-analysis-contact-pairs) | bioskills | Turns Hi-C/Micro-C FASTQ into a deduplicated, filtered .pairs file with pairtools and decides whether the library worked. Covers the bwa mem -SP5M / bwa-mem2 /… |
| [`hi-c-analysis-hic-data-io`](#hi-c-analysis-hic-data-io) | bioskills | Loads, converts, and manipulates Hi-C contact matrices in cooler format (.cool/.mcool/.scool) and Juicer .hic, using cooler (Python + CLI), hic2cool, and hictk… |
| [`hi-c-analysis-hic-differential`](#hi-c-analysis-hic-differential) | bioskills | Compares Hi-C contact maps between conditions across the right scale -- differential bin-pair contacts (multiHiCcompare, diffHic), differential A/B compartment… |
| [`hi-c-analysis-hic-visualization`](#hi-c-analysis-hic-visualization) | bioskills | Renders Hi-C contact matrices honestly and reproducibly with matplotlib, cooltools, HiCExplorer, pyGenomeTracks, FAN-C, CoolBox, and plotgardener. Covers the r… |
| [`hi-c-analysis-hichip-plac-loops`](#hi-c-analysis-hichip-plac-loops) | bioskills | Calls significant loops from protein-directed and targeted 3C assays (HiChIP, PLAC-seq, Capture Hi-C/PCHi-C, ChIA-PET) where the contact background is peak-anc… |
| [`hi-c-analysis-loop-calling`](#hi-c-analysis-loop-calling) | bioskills | Detects focal chromatin loops (point interactions / corner-dots) in balanced Hi-C and Micro-C contact maps and aggregates/validates a loop set. Covers de-novo… |
| [`hi-c-analysis-matrix-operations`](#hi-c-analysis-matrix-operations) | bioskills | Balances Hi-C contact matrices (ICE via cooler.balance_cooler, KR/SCALE/VC context), computes distance-decay expected with cooltools (expected_cis per-diagonal… |
| [`hi-c-analysis-tad-detection`](#hi-c-analysis-tad-detection) | bioskills | Detects TAD boundaries from balanced Hi-C contact matrices via the diamond-window insulation score (cooltools insulation) and HiCExplorer hicFindTADs, returnin… |
| [`hic-aggregation`](#hic-aggregation) | encode-toolkit | Build comprehensive chromatin contact maps by aggregating Hi-C loop calls (BEDPE) across multiple ENCODE experiments, donors, and labs. Use when the user wants… |
| [`histone-aggregation`](#histone-aggregation) | encode-toolkit | Build comprehensive histone mark maps by aggregating narrowPeak data across multiple ENCODE experiments, donors, and labs. Use when the user wants to answer "w… |
| [`methylation-aggregation`](#methylation-aggregation) | encode-toolkit | Build comprehensive DNA methylation maps by aggregating WGBS (Whole Genome Bisulfite Sequencing) data across multiple ENCODE experiments, donors, and labs. Use… |
| [`methylation-array-preprocessing`](#methylation-array-preprocessing) | bioskills | Turns raw Illumina Infinium methylation BeadChip IDATs (450K, EPIC, EPICv2) into a defensible beta/M matrix with sesame (openSesame/SigDF) or minfi (RGChannelS… |
| [`methylation-array-qc-filtering`](#methylation-array-qc-filtering) | bioskills | Performs probe filtering and sample-level QC on Illumina Infinium methylation arrays (450K / EPIC / EPICv2) to decide which probes and samples to trust. Drops… |
| [`methylation-bismark-alignment`](#methylation-bismark-alignment) | bioskills | Aligns bisulfite-converted (WGBS, RRBS, PBAT) and enzymatic (EM-seq) short reads to an in-silico C-&gt;T/G-&gt;A-converted reference with Bismark (Bowtie2 or HISAT2)… |
| [`methylation-calling`](#methylation-calling) | bioskills | Extracts per-cytosine methylation calls from aligned bisulfite/EM-seq reads with bismark_methylation_extractor (Bismark BAM) or the aligner-agnostic MethylDack… |
| [`methylation-cell-type-deconvolution`](#methylation-cell-type-deconvolution) | bioskills | Estimates cell-type composition from bulk DNA methylation and uses it to defuse the single biggest EWAS confounder. Covers reference-based deconvolution (House… |
| [`methylation-differential-cpg`](#methylation-differential-cpg) | bioskills | Tests individual CpG sites for differential methylation (DMC/DMP) from bisulfite sequencing counts or array/continuous beta-value matrices. Covers the count-vs… |
| [`methylation-dmr-detection`](#methylation-dmr-detection) | bioskills | Detects differentially methylated regions (DMRs) from short-read bisulfite (WGBS/RRBS), array, and long-read methylation count tables using dmrseq (permutation… |
| [`methylation-epigenetic-clocks`](#methylation-epigenetic-clocks) | bioskills | Computes DNA methylation age (DNAm age) and pace of aging by applying frozen elastic-net epigenetic clocks to a clean beta matrix with methylclock, dnaMethyAge… |
| [`methylation-ewas-design`](#methylation-ewas-design) | bioskills | Designs and defends an epigenome-wide association study (EWAS) on 450K/EPIC array or bisulfite methylation - the layer deciding whether a hit is credible. Cove… |
| [`methylation-methylkit`](#methylation-methylkit) | bioskills | Imports Bismark coverage or cytosine-report files into the methylKit object model, then runs the import-to-results spine - filterByCoverage, normalizeCoverage,… |
| [`pipeline-dnaseseq`](#pipeline-dnaseseq) | encode-toolkit | Execute ENCODE DNase-seq pipeline from FASTQ to hotspots and footprints. Child of pipeline-guide. Provides Nextflow execution with Docker and cloud deployment.… |
| [`quality-assessment`](#quality-assessment) | encode-toolkit | Evaluate ENCODE experiment quality using standard metrics and audit flags. Use when the user asks about data quality, wants to filter for high-quality experime… |
| [`regulatory-elements`](#regulatory-elements) | encode-toolkit | Discover and characterize regulatory elements (enhancers, promoters, silencers, insulators, super-enhancers) using ENCODE data and the cCRE catalog. Use when t… |
| [`regulatory-genomics`](#regulatory-genomics) | tooluniverse | Transcription factor binding, cis-regulatory elements (cCREs), chromatin accessibility, and regulatory annotation using JASPAR (motifs), ENCODE (cCREs, ChIP-se… |
| [`regulatory-variant-analysis`](#regulatory-variant-analysis) | tooluniverse | Non-coding/regulatory variant interpretation — GWAS association lookup, eQTL evidence (GTEx), chromatin state (ENCODE), regulatory variant scoring (RegulomeDB,… |
| [`scvi-tools`](#scvi-tools) | user335-anthropic | Deep learning for single-cell analysis using scvi-tools. This skill should be used when users need (1) data integration and batch correction with scVI/scANVI,… |
| [`ucsc-conservation-and-tfbs`](#ucsc-conservation-and-tfbs) | deepmind | Fetch Evolutionary Conservation scores (phyloP, phastCons) and Transcription Factor Binding Sites (TFBS) from the UCSC Genome Browser. Use when analyzing wheth… |
| [`variant-to-mechanism`](#variant-to-mechanism) | tooluniverse | End-to-end variant-to-mechanism analysis — trace a variant (rsID/coordinates) through regulatory context, target gene(s), molecular pathway(s), and phenotypic… |

## `accessibility-aggregation` {#accessibility-aggregation}

Build comprehensive chromatin accessibility maps by aggregating ATAC-seq and DNase-seq narrowPeak data across multiple ENCODE experiments, donors, and labs. Use when the user wants to answer "where is chromatin accessible in my tissue?" by combining peak calls into a union peak set. Handles cross-lab variation, ATAC vs DNase platform differences, and ENCODE blocklist filtering.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/accessibility-aggregation)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info accessibility-aggregation
dnakit skills install accessibility-aggregation --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-allele-specific-accessibility` {#atac-seq-allele-specific-accessibility}

Detect allele-specific chromatin accessibility from ATAC-seq using WASP, GATK ASEReadCounter, or RASQUAL. Use when mapping cis-regulatory genetic variants from heterozygous SNPs, separating cis from trans regulation, building chromatin QTL (caQTL) maps, validating GWAS variant function with allelic imbalance, or detecting reference allele mapping bias before downstream analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/allele-specific-accessibility)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-allele-specific-accessibility`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-allele-specific-accessibility
dnakit skills install atac-seq-allele-specific-accessibility --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-atac-peak-calling` {#atac-seq-atac-peak-calling}

Call accessible chromatin regions from ATAC-seq BAM files using MACS3, MACS2, Genrich, or HMMRATAC. Use when identifying open chromatin from aligned ATAC-seq, choosing between point-source vs HMM peak callers, applying ENCODE-style pseudoreplicate IDR, removing blacklist regions, or fixing 501bp consensus peaks for downstream differential analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/atac-peak-calling)
- **许可证：**`MIT`
- **运行依赖：**命令 `macs3`
- **别名：**`bio-atac-seq-atac-peak-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-atac-peak-calling
dnakit skills install atac-seq-atac-peak-calling --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-atac-qc` {#atac-seq-atac-qc}

ATAC-seq library quality control -- TSS enrichment, FRiP, fragment-size periodicity, library complexity (NRF/PBC1/PBC2), mitochondrial fraction, and ENCODE 4 thresholds. Use when assessing whether an ATAC-seq library passes ENCODE acceptance criteria, diagnosing transposition artefacts, comparing Omni-ATAC vs standard prep quality, or selecting which replicates to drop before peak calling.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/atac-qc)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-atac-qc`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-atac-qc
dnakit skills install atac-seq-atac-qc --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-co-accessibility` {#atac-seq-co-accessibility}

Infer cis-regulatory connections (peak-to-peak co-accessibility) from scATAC-seq using Cicero, ArchR getCoAccessibility, or SCENIC+. Use when linking enhancer accessibility to promoter accessibility, identifying enhancer-gene pairs from chromatin alone (without paired RNA), running gene-regulatory inference combining ATAC + RNA, or comparing predicted regulatory contacts against Hi-C/Micro-C ground truth.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/co-accessibility)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-co-accessibility`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-co-accessibility
dnakit skills install atac-seq-co-accessibility --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-consensus-peakset` {#atac-seq-consensus-peakset}

Build a differential-ready consensus peakset from per-replicate ATAC-seq peaks using iterative overlap removal, fixed-width re-centering, and majority-rule overlap. Use when generating a stable peak coordinate system for downstream differential accessibility, ML feature engineering, cross-sample comparison, or fixed-width peak counts; covers Corces 2018 iterative overlap (501 bp), DiffBind summit re-centering, and ENCODE consistency rules.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/consensus-peakset)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-consensus-peakset`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-consensus-peakset
dnakit skills install atac-seq-consensus-peakset --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-deep-learning-atac` {#atac-seq-deep-learning-atac}

Sequence-based deep learning for ATAC-seq using chromBPNet, BPNet, scBasset, or Enformer. Use when correcting Tn5 bias with neural networks beyond k-mer models, predicting per-base accessibility profiles, scoring in silico variant effects at GWAS or rare-variant SNPs, discovering motifs via DeepLIFT/TF-MoDISco from a trained model, or generating cell-type-specific accessibility predictions for unobserved cell states.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/deep-learning-atac)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `chrombpnet`, `dnakit`
- **别名：**`bio-atac-seq-deep-learning-atac`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-deep-learning-atac
dnakit skills install atac-seq-deep-learning-atac --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-differential-accessibility` {#atac-seq-differential-accessibility}

Identify differentially accessible chromatin regions across conditions using DiffBind, csaw, DESeq2, or edgeR. Use when comparing ATAC-seq accessibility between treatment groups, choosing between consensus-peak vs sliding-window approaches, picking the correct normalization (full library vs reads-in-peaks), correcting batch with SVA/RUVseq, or interpreting log2FC and FDR thresholds in a chromatin context.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/differential-accessibility)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-differential-accessibility`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-differential-accessibility
dnakit skills install atac-seq-differential-accessibility --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-enhancer-gene-linking` {#atac-seq-enhancer-gene-linking}

Predict enhancer-gene regulatory connections from ATAC-seq using ABC, ENCODE-rE2G, HiChIP, or Cicero. Use when linking distal enhancers to target genes, choosing between contact-aware (ABC, ENCODE-rE2G), accessibility-only (Cicero), and orthogonal (HiChIP H3K27ac, EpiMap) approaches, validating predictions against CRISPRi-FlowFISH gold-standard, or building cell-type-specific regulatory maps for fine-mapping or therapeutic target discovery.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/enhancer-gene-linking)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-enhancer-gene-linking`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-enhancer-gene-linking
dnakit skills install atac-seq-enhancer-gene-linking --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-footprinting` {#atac-seq-footprinting}

Detect transcription factor binding footprints in ATAC-seq using TOBIAS, HINT-ATAC, Wellington, or scprinter. Use when identifying bound TF sites within accessible regions, correcting Tn5 insertion bias before footprinting, choosing between cleavage-based and aggregate-based footprinters, or comparing differential TF activity between conditions.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/footprinting)
- **许可证：**`MIT`
- **运行依赖：**命令 `tobias`
- **别名：**`bio-atac-seq-footprinting`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-footprinting
dnakit skills install atac-seq-footprinting --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-motif-deviation` {#atac-seq-motif-deviation}

Analyze TF motif accessibility variability across samples or single cells using chromVAR. Use when identifying TF motifs whose accessibility correlates with conditions, computing per-sample motif z-scores after matched background correction, comparing to ArchR / Signac equivalents, or distinguishing motif-accessibility signal from per-site footprinting.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/motif-deviation)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-motif-deviation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-motif-deviation
dnakit skills install atac-seq-motif-deviation --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-nucleosome-positioning` {#atac-seq-nucleosome-positioning}

Map nucleosome center positions, occupancy, and fuzziness from ATAC-seq fragment-size patterns using NucleoATAC, ATACseqQC, DANPOS3, or scprinter. Use when characterizing nucleosome organization at promoters and enhancers, calling +1/-1 nucleosomes flanking NFRs, generating V-plots for chromatin structure visualization, or comparing nucleosome positioning between conditions.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/nucleosome-positioning)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-nucleosome-positioning`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-nucleosome-positioning
dnakit skills install atac-seq-nucleosome-positioning --agent claude-code   # 或 --agent codex / cursor
```

## `atac-seq-single-cell-atac` {#atac-seq-single-cell-atac}

Process and analyze single-cell ATAC-seq data with Signac, ArchR, SnapATAC2, or Cell Ranger ATAC. Use when handling 10X scATAC or 10X Multiome (paired RNA+ATAC) data, performing per-cell QC, choosing between ArchR/Signac/SnapATAC2 ecosystems, building per-cluster consensus peaksets, integrating with paired scRNA-seq, doublet detection (AMULET vs ArchR vs scDblFinder), or running pseudobulk differential accessibility per cluster.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/atac-seq/single-cell-atac)
- **许可证：**`MIT`
- **别名：**`bio-atac-seq-single-cell-atac`, `bio-single-cell-scatac-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info atac-seq-single-cell-atac
dnakit skills install atac-seq-single-cell-atac --agent claude-code   # 或 --agent codex / cursor
```

## `chip-atlas-diff-analysis` {#chip-atlas-diff-analysis}

Submit two groups of public ChIP-seq, ATAC-seq, DNase-seq, or Bisulfite-seq experiments to the ChIP-Atlas Diff Analysis API and recover differential peak regions (DPR via edgeR) or differentially methylated regions (DMR via metilene). Use when you have public experiment accessions and want a server-side comparison without downloading raw data.

- **来源：**[hubayirp/agentic-science contributors @ `213588a`](https://github.com/hubayirp/agentic-science/tree/213588a92bb74e705d795ba106e0119975877e8c/skills/chip-atlas-diff-analysis)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chip-atlas-diff-analysis
dnakit skills install chip-atlas-diff-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `chip-atlas-peak-enrichment` {#chip-atlas-peak-enrichment}

Submit a gene list to the ChIP-Atlas Enrichment Analysis API and recover the transcription factors, chromatin modifiers, and histone marks whose public ChIP-seq peaks cluster near those gene promoters. Use when you have a candidate gene set and want a ranked, FDR-controlled list of regulators across 430,000+ public ChIP-seq experiments.

- **来源：**[hubayirp/agentic-science contributors @ `213588a`](https://github.com/hubayirp/agentic-science/tree/213588a92bb74e705d795ba106e0119975877e8c/skills/chip-atlas-peak-enrichment)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chip-atlas-peak-enrichment
dnakit skills install chip-atlas-peak-enrichment --agent claude-code   # 或 --agent codex / cursor
```

## `chip-atlas-target-genes` {#chip-atlas-target-genes}

For a single transcription factor or chromatin regulator, download ChIP-Atlas's pre-computed per-gene MACS2 binding scores across every public ChIP-seq experiment, rank the targets, optionally subset to specific cell types, and integrate the pre-embedded STRING interaction scores. Use when you have one TF in mind and want a ranked list of its candidate target genes.

- **来源：**[hubayirp/agentic-science contributors @ `213588a`](https://github.com/hubayirp/agentic-science/tree/213588a92bb74e705d795ba106e0119975877e8c/skills/chip-atlas-target-genes)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chip-atlas-target-genes
dnakit skills install chip-atlas-target-genes --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-allele-specific-binding` {#chipseq-allele-specific-binding}

Detects allele-specific transcription factor or histone modification binding from heterozygous-variant ChIP-seq using WASP (reference-bias filter; mandatory upstream), RASQUAL (joint QTL + bias-corrected testing), BaalChIP (Bayesian beta-binomial with copy-number-aware overdispersion), and AlleleSeq (personalized diploid genome). Handles imprinted-locus awareness, X-inactivation artifacts, cancer copy-number imbalance, and integration with downstream caQTL / bQTL mapping. Use when identifying variants with allelic effects on TF binding, fine-mapping causal regulatory variants, validating deep-learning variant predictions, or characterizing cis-acting regulatory effects.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/allele-specific-binding)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-allele-specific-binding`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-allele-specific-binding
dnakit skills install chipseq-allele-specific-binding --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-chip-deep-learning` {#chipseq-chip-deep-learning}

Trains and applies base-resolution deep learning models on ChIP-seq / ChIP-nexus / CUT&amp;RUN data. Uses BPNet (Avsec 2021 Nat Genet 53:354; soft motif syntax from ChIP-nexus), chromBPNet (Pampari A et al 2024 bioRxiv; bias-factorized base-resolution profiles), EnFormer (Avsec 2021 Nat Methods 18:1196; 196 kb input, ~100 kb effective receptive field), DeepSEA (Zhou 2015; multi-task CNN), and JASPAR 2026 deep-learning collection (1259 BPNet ChIP models). Performs in silico mutagenesis for variant-effect prediction, DeepLIFT/Grad attribution, and TF-MoDISco motif discovery from attribution scores. Use when predicting variant effects on TF binding, discovering soft motif syntax / cooperativity, integrating ChIP-seq with sequence-only predictions, or applying precomputed JASPAR Deep Learning models to new variants.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/chip-deep-learning)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `chrombpnet`, `dnakit`
- **别名：**`bio-chipseq-chip-deep-learning`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-chip-deep-learning
dnakit skills install chipseq-chip-deep-learning --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-chromatin-state-segmentation` {#chipseq-chromatin-state-segmentation}

Segments the genome into chromatin states from combinatorial histone modification and chromatin factor ChIP-seq data. Uses ChromHMM (multivariate HMM on binarized signal, v1.27), Segway (Dynamic Bayesian Network on continuous signal), EpiSegMix (flexible-distribution HMM with duration modeling, 2024), EpiLogos (multi-biosample visualization), IDEAS (cell-type-aware joint), and full-stack ChromHMM (Vu Ernst 2022) for cross-cell-type segmentations. Handles state-count selection (15 vs 18 vs 25 states), binarization choice, OverlapEnrichment / NeighborhoodEnrichment downstream analysis, and cross-biosample integration. Use when learning chromatin states from a histone mark panel, characterizing learned states by genomic feature enrichment, or comparing chromatin landscapes across cell types.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/chromatin-state-segmentation)
- **许可证：**`MIT`
- **运行依赖：**命令 `ChromHMM`
- **别名：**`bio-chipseq-chromatin-state-segmentation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-chromatin-state-segmentation
dnakit skills install chipseq-chromatin-state-segmentation --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-cut-and-run-tag` {#chipseq-cut-and-run-tag}

Analyzes CUT&amp;RUN (Skene Henikoff 2017) and CUT&amp;Tag (Kaya-Okur 2019) chromatin profiling data. Handles SEACR vs MACS2 peak calling (with the btaf375 2025 benchmark guidance), pA-MNase vs pA-Tn5 vs pAG-Tn5 chimera differences, E. coli spike-in carryover normalization, IgG-only control logic (no input), characteristic fragment-size signatures (25-75 bp for CUT&amp;Tag), and lower depth requirements (5M reads typical vs 25M for ChIP). Use when calling peaks from CUT&amp;RUN/CUT&amp;Tag, scaling by E. coli spike-in carryover, choosing SEACR norm mode, or comparing CUT&amp;RUN/Tag results to traditional ChIP.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/cut-and-run-tag)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-cut-and-run-tag`, `pipeline-cutandrun`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-cut-and-run-tag
dnakit skills install chipseq-cut-and-run-tag --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-differential-binding` {#chipseq-differential-binding}

Identifies differentially bound ChIP-seq regions between conditions using DiffBind, csaw (sliding windows), DESeq2/edgeR/PyDESeq2 on count matrices, NormR (control-aware), or MAnorm2. Distinguishes three distinct normalization problems (composition bias, trended bias, global shifts) and matches each to its appropriate fix including spike-in scaling. Use when comparing ChIP-seq binding between experimental conditions, choosing normalization for global vs local changes, integrating spike-in data, or reconciling DiffBind/DESeq2 disagreement.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/differential-binding)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-differential-binding`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-differential-binding
dnakit skills install chipseq-differential-binding --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-motif-analysis` {#chipseq-motif-analysis}

Discovers de novo motifs and tests known motif enrichment in ChIP-seq, ATAC-seq, or other peak sequences using HOMER, MEME-ChIP (STREME, CentriMo, TOMTOM, FIMO), monaLisa, and AME. Handles background selection (GC-matched, dinucleotide-shuffled, Markov order-2, peak-flanks), motif databases (JASPAR 2024 CORE PWMs, JASPAR 2026 deep-learning collection, HOCOMOCO v12, HOMER built-in), centrally-enriched motif testing, and differential motif analysis. Use when identifying TF binding motifs in peaks, testing for known TF enrichment, scanning for motif instances, comparing motif content between conditions, or interpreting motifs from deep learning models.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/motif-analysis)
- **许可证：**`MIT`
- **运行依赖：**命令 `HOMER`
- **别名：**`bio-chipseq-motif-analysis`, `homer-motif-analysis`, `motif-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-motif-analysis
dnakit skills install chipseq-motif-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-peak-annotation` {#chipseq-peak-annotation}

Annotates ChIP-seq peaks to genomic features, nearest genes, ENCODE candidate cis-regulatory elements (cCREs), and regulatory domains. Uses ChIPseeker (R), HOMER annotatePeaks.pl (CLI), pyranges (Python), GREAT/rGREAT (regulatory domain gene-set enrichment), ChIP-Enrich (locus-length-adjusted), ENCODE SCREEN cCRE classification (PLS/pELS/dELS/CA-CTCF/CA-H3K4me3), and ENCODE-rE2G for cell-type-specific enhancer-gene linking. Handles nearest-TSS vs host-gene ambiguity, promoter window definition, and feature priority. Use when assigning genomic context to peaks, linking enhancer peaks to target genes, classifying peaks against ENCODE cCRE registry, or running gene-set enrichment on peak-associated genes.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/peak-annotation)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-peak-annotation`, `peak-annotation`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-peak-annotation
dnakit skills install chipseq-peak-annotation --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-peak-calling` {#chipseq-peak-calling}

Calls ChIP-seq peaks with MACS3, MACS2, HOMER, or SPP across narrow (TF) and broad (histone) modes. Handles input control matching, fragment-size modeling vs --nomodel, effective genome size, ENCODE-style IDR vs naive overlap, hyper-ChIPable artifacts, and aligner-specific shifts. Use when calling peaks from ChIP-seq alignments, choosing between narrow vs broad mode for a histone mark, deciding model vs nomodel for low-depth data, applying ENCODE pseudoreplicate IDR, or reconciling MACS vs HOMER vs SPP results.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/peak-calling)
- **许可证：**`MIT`
- **运行依赖：**命令 `macs3`
- **别名：**`bio-chipseq-peak-calling`, `macs3-peak-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-peak-calling
dnakit skills install chipseq-peak-calling --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-qc` {#chipseq-qc}

Assesses ChIP-seq quality across antibody specificity, fragmentation, enrichment, replicate concordance, and library complexity. Computes FRiP, NSC/RSC (phantompeakqualtools), library complexity (NRF/PBC1/PBC2), deepTools plotFingerprint (JS distance, AUC, synthetic JS), ChIPQC, IDR with ENCODE Nself/Nt rules, and detects hyper-ChIPable artifacts. Use when validating an antibody, diagnosing failed peak calls, deciding whether to proceed with downstream analysis, grading against ENCODE thresholds, or auditing replicate concordance.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/chipseq-qc)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-qc`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-qc
dnakit skills install chipseq-qc --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-spike-in-normalization` {#chipseq-spike-in-normalization}

Normalizes ChIP-seq data using exogenous spike-in (ChIP-Rx with Drosophila chromatin per Orlando 2014 / Egan 2016; E. coli carryover for CUT&amp;RUN/CUT&amp;Tag). Distinguishes RRPM from Rx-Input scaling, integrates with DiffBind / DESeq2 / edgeR / csaw via sizeFactors and DiffBind library-size vectors, applies the Patel et al 2024 *Nat Biotechnol* failure-mode framework, and validates that normalization is applied at the read level (not peak counts). Use when global signal shifts are expected (HDACi, BETi, EZH2i, dosage, target knockdown), when ChIPseqSpikeInFree detects post-hoc shifts, or when validating internal-control regions before publication.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/spike-in-normalization)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-spike-in-normalization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-spike-in-normalization
dnakit skills install chipseq-spike-in-normalization --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-super-enhancers` {#chipseq-super-enhancers}

Identifies super-enhancers from H3K27ac, MED1, or BRD4 ChIP-seq using ROSE, ROSE2, LILY, HOMER -style super, and ENCODE dELS cross-referencing. Handles peak stitching parameters, ranking choices, hockey-stick inflection, marker choice (H3K27ac vs MED1/BRD4), and cross-condition comparison with spike-in normalization. Constructs core regulatory circuitry (Saint-Andre 2016) from SE-encoded TFs. Use when identifying cell-identity / cancer-associated regulatory domains, comparing super-enhancers between conditions, identifying master transcription factor networks, or predicting BET-inhibitor responsiveness.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/super-enhancers)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-super-enhancers`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-super-enhancers
dnakit skills install chipseq-super-enhancers --agent claude-code   # 或 --agent codex / cursor
```

## `chipseq-visualization` {#chipseq-visualization}

Visualizes ChIP-seq data using deepTools (computeMatrix, plotHeatmap, plotProfile, bamCoverage, bamCompare), pyGenomeTracks (modern INI-driven track plots), Gviz (R browser-style), EnrichedHeatmap (ComplexHeatmap-based), ChIPseeker tag heatmaps, and IGV batch screenshots. Handles bigWig normalization choices (CPM, BPM, RPGC, spike-in scaled), bamCompare operations (log2 ratio, subtract) with SES scaling, k-means clustering of heatmaps for biological subgrouping, and spike-in-scaled tracks for global-shift experiments. Use when generating publication-quality ChIP-seq signal heatmaps, profile plots, genome-browser tracks, or comparing samples visually.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/chip-seq/chipseq-visualization)
- **许可证：**`MIT`
- **别名：**`bio-chipseq-visualization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info chipseq-visualization
dnakit skills install chipseq-visualization --agent claude-code   # 或 --agent codex / cursor
```

## `compare-biosamples` {#compare-biosamples}

Compare ENCODE experiments across different biosamples, tissues, or cell lines to identify tissue-specific regulatory patterns. Use when the user wants cross-tissue comparison, cell-type comparison, tissue-specific elements, differential chromatin, biosample matching, disease vs normal comparison, developmental time course, constitutive vs variable regulation, or multi-tissue data availability mapping. Handles batch effect detection, biosample hierarchy, and comparison design.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/compare-biosamples)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info compare-biosamples
dnakit skills install compare-biosamples --agent claude-code   # 或 --agent codex / cursor
```

## `data-visualization-genome-tracks` {#data-visualization-genome-tracks}

Build genome-browser-style multi-track figures with pyGenomeTracks (config-driven), Gviz (R), and IGV batch screenshotting. Covers BigWig coverage tracks, BED/peak overlays, gene-model rendering, Hi-C matrix tracks, BedPE link arcs, spike-in-aware normalization, and the bamCoverage --normalizeUsing trap. Use when producing publication figures of genomic loci with stacked aligned tracks (coverage, peaks, genes, interactions) for ChIP-seq, ATAC-seq, RNA-seq, Hi-C, or generic locus visualization.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/data-visualization/genome-tracks)
- **许可证：**`MIT`
- **别名：**`bio-data-visualization-genome-tracks`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info data-visualization-genome-tracks
dnakit skills install data-visualization-genome-tracks --agent claude-code   # 或 --agent codex / cursor
```

## `dnakit-methylation-summary` {#dnakit-methylation-summary}

Read modkit extended bedMethyl and compute coverage-weighted methylation summaries with DNAKit. Use for nanopore modification counts, including separate 5mC and 5hmC groups.

- **来源：**[DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-methylation-summary)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnakit-methylation-summary
dnakit skills install dnakit-methylation-summary --agent claude-code   # 或 --agent codex / cursor
```

## `epigenomics` {#epigenomics}

Genomics and epigenomics analysis: DNA methylation (CpG, 5mC, 5hmC, bisulfite, RRBS), m6A RNA modification (MeRIP-seq), ChIP-seq peaks, ATAC-seq accessibility, histone modifications, chromatin state, multi-omics integration. Combines pandas/scipy/pysam computation with ToolUniverse annotation tools. Use for genome-wide epigenomic statistics, methylation analysis, and chromatin-genome integration.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-epigenomics)
- **许可证：**`Apache-2.0`
- **别名：**`epigenome-profiling`, `tooluniverse-epigenomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info epigenomics
dnakit skills install epigenomics --agent claude-code   # 或 --agent codex / cursor
```

## `epigenomics-chromatin` {#epigenomics-chromatin}

Histone-modification ChIP-seq, ATAC-seq accessibility, chromatin state, and TF binding analysis from ENCODE, Roadmap Epigenomics, ChIP-Atlas. Use for chromatin-state-by-tissue queries, TF-binding-by-region, regulatory landscape mapping, and ENCODE-cCRE annotations. For DNA methylation use epigenomics; for RNA-seq use tooluniverse-rnaseq-deseq2.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-epigenomics-chromatin)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-epigenomics-chromatin`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info epigenomics-chromatin
dnakit skills install epigenomics-chromatin --agent claude-code   # 或 --agent codex / cursor
```

## `gene-regulatory-networks` {#gene-regulatory-networks}

Gene regulatory network analysis — TF-target inference (JASPAR motifs, ChIP-seq), motif scanning, eQTL integration, perturbation evidence (knockout/overexpression). Use for 'which TF regulates gene X', 'which genes does TF Y target', regulatory pathway reconstruction. Distinguishes direct (binding) vs indirect (co-expression) regulatory evidence.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-gene-regulatory-networks)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `dnakit`, `networkx`
- **别名：**`tooluniverse-gene-regulatory-networks`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info gene-regulatory-networks
dnakit skills install gene-regulatory-networks --agent claude-code   # 或 --agent codex / cursor
```

## `gene-regulatory-networks-multiomics-grn` {#gene-regulatory-networks-multiomics-grn}

Build enhancer-driven gene regulatory networks (eGRNs) by integrating single-cell RNA-seq and ATAC-seq using SCENIC+, CellOracle base GRNs, Pando, FigR, DIRECT-NET, TRIPOD, and scMEGA. Covers the accessibility-defines-enhancers principle, peak-to-gene linking and its cell-composition confound, the paired-vs-unpaired decision, and TF-region-gene eRegulon triplets. Use when analyzing 10x multiome or paired/unpaired scRNA+scATAC to infer cis-regulatory GRNs. For RNA-only regulons see scenic-regulons; for in silico TF perturbation see perturbation-simulation.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/gene-regulatory-networks/multiomics-grn)
- **许可证：**`MIT`
- **别名：**`bio-gene-regulatory-networks-multiomics-grn`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info gene-regulatory-networks-multiomics-grn
dnakit skills install gene-regulatory-networks-multiomics-grn --agent claude-code   # 或 --agent codex / cursor
```

## `gene-regulatory-networks-scenic-regulons` {#gene-regulatory-networks-scenic-regulons}

Infer transcription factor regulons from single-cell RNA-seq with pySCENIC by combining GRNBoost2 co-expression, cisTarget motif-enrichment pruning, and AUCell per-cell activity scoring. Covers the motif-pruning-as-directionality principle, regulon specificity scoring, run-to-run stability, and database/species matching. Use when identifying TF regulons, scoring TF activity per cell, finding master regulators of cell identity, or comparing regulon activity across conditions. For enhancer-driven multiomic GRNs see multiomics-grn; for bulk inference and VIPER protein-activity see grn-inference.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/gene-regulatory-networks/scenic-regulons)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `pySCENIC`
- **别名：**`bio-gene-regulatory-networks-scenic-regulons`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info gene-regulatory-networks-scenic-regulons
dnakit skills install gene-regulatory-networks-scenic-regulons --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-compartment-analysis` {#hi-c-analysis-compartment-analysis}

Detects A/B chromatin compartments from balanced Hi-C contact matrices via eigenvector decomposition of the distance-normalized, Pearson-correlated cis matrix with cooltools (eigs_cis), then orients (phases) the compartment eigenvector against a GC or gene-density track so the active (A) sign is not arbitrary. Covers the eigenvector-is-a-choice problem (per-arm view_df to remove the centromere gradient; picking the eigenvector by max correlation with activity, not by eigenvalue), GC phasing with bioframe.frac_gc, resolution choice (100kb-1Mb), saddle plots and saddle_strength for compartmentalization strength, the cohesin-loss-strengthens-compartments result, subcompartments (SNIPER/Calder/dcHiC), and cross-condition compartment switching. Use when calling A/B compartments, computing E1/eigenvectors, phasing the eigenvector, building saddle plots, choosing a compartment resolution, quantifying compartment strength, or comparing compartmentalization across conditions.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/compartment-analysis)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `cooltools`, `dnakit`
- **别名：**`bio-hi-c-analysis-compartment-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-compartment-analysis
dnakit skills install hi-c-analysis-compartment-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-contact-pairs` {#hi-c-analysis-contact-pairs}

Turns Hi-C/Micro-C FASTQ into a deduplicated, filtered .pairs file with pairtools and decides whether the library worked. Covers the bwa mem -SP5M / bwa-mem2 / chromap --preset hic alignment idiom (mates mapped as independent single-end reads), pairtools parse vs parse2 and the walks-policy choice (5unique pairwise vs all for Pore-C/Micro-C concatemers), pair-type classification (keep UU and rescued UC), dedup (PCR vs optical/by-tile), select by pair_type/MAPQ/distance, restriction-fragment handling (restrict, Arima dual-enzyme, Micro-C/DNase fragment-free), and allele-specific phasing (pairtools phase to two coolers). The library-QC decision uses % long-range cis as the one-number quality metric, trans as the noise floor, orientation balance as fragment-map-free dangling-end/self-circle QC, and % duplicates as a complexity proxy. Use when processing Hi-C/Micro-C/Omni-C reads into pairs, judging library quality, handling multi-enzyme or restriction-agnostic protocols, or generating allele-specific contacts.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/contact-pairs)
- **许可证：**`MIT`
- **运行依赖：**命令 `pairtools`
- **别名：**`bio-hi-c-analysis-contact-pairs`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-contact-pairs
dnakit skills install hi-c-analysis-contact-pairs --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-hic-data-io` {#hi-c-analysis-hic-data-io}

Loads, converts, and manipulates Hi-C contact matrices in cooler format (.cool/.mcool/.scool) and Juicer .hic, using cooler (Python + CLI), hic2cool, and hictk. Covers the single-resolution mcool URI (file.mcool::/resolutions/&lt;bp&gt;), the load-bearing divisive-vs-multiplicative weight-naming rule (KR/VC/VC_SQRT auto-divisive vs cooler's multiplicative weight), what survives .hic&lt;-&gt;.cool conversion (FRAG matrices and norm vectors do not), raw-vs-balanced coarsening, the .pairs upper-triangle/chromsize-order contract, and chrom-naming/bin-table provenance. Use when loading a cooler, converting .hic to .mcool, selecting a resolution, building a cooler from pairs or a matrix, coarsening/zoomifying, importing Juicer norm vectors, or debugging all-NaN balanced matrices and chr1-vs-1 empty fetches.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/hic-data-io)
- **许可证：**`MIT`
- **别名：**`bio-hi-c-analysis-hic-data-io`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-hic-data-io
dnakit skills install hi-c-analysis-hic-data-io --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-hic-differential` {#hi-c-analysis-hic-differential}

Compares Hi-C contact maps between conditions across the right scale -- differential bin-pair contacts (multiHiCcompare, diffHic), differential A/B compartments (dcHiC), differential TAD boundaries (delta insulation), and differential loops (diffloop, DiffHiChIP) -- with distance-stratified between-sample normalization, replicate-aware NB-GLM FDR, HiCRep SCC reproducibility gating, and CNV correction for cancer/aneuploid samples. Use when comparing Hi-C between treatment and control, finding differential contacts/compartments/boundaries/loops, normalizing two maps of unequal depth, choosing a replicate-aware test, gating replicates with SCC, or correcting copy-number artifacts before a tumor-vs-normal comparison.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/hic-differential)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `cooltools`, `dnakit`
- **别名：**`bio-hi-c-analysis-hic-differential`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-hic-differential
dnakit skills install hi-c-analysis-hic-differential --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-hic-visualization` {#hi-c-analysis-hic-visualization}

Renders Hi-C contact matrices honestly and reproducibly with matplotlib, cooltools, HiCExplorer, pyGenomeTracks, FAN-C, CoolBox, and plotgardener. Covers the raw/ICE-balanced/observed-over-expected transform choice, LogNorm vs symmetric-diverging colormaps with vmax/percentile clipping, resolution-to-feature matching (compartments 100-500kb, TADs 10-40kb, loops 5-10kb), square vs rotated-triangle track-stacking, NaN/white-stripe handling, virtual 4C, APA/saddle/on-diagonal pileups, two-condition side-by-side and log2-ratio maps, and interactive (HiGlass) vs scripted-static publication figures. Use when plotting a contact matrix, choosing a normalization or color scale, building a multi-track Hi-C figure, making a virtual 4C profile, piling up loops/boundaries, or comparing two conditions.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/hic-visualization)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `matplotlib`
- **别名：**`bio-hi-c-analysis-hic-visualization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-hic-visualization
dnakit skills install hi-c-analysis-hic-visualization --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-hichip-plac-loops` {#hi-c-analysis-hichip-plac-loops}

Calls significant loops from protein-directed and targeted 3C assays (HiChIP, PLAC-seq, Capture Hi-C/PCHi-C, ChIA-PET) where the contact background is peak-anchored and coverage-biased, so generic Hi-C loop callers (cooltools dots, Juicer HiCCUPS) use the wrong null. Covers FitHiChIP (config-driven coverage+distance-decay spline regression, peak-to-peak vs peak-to-all foreground, loose vs stringent background, coverage vs ICE bias), MAPS (positive Poisson regression on bias factors for PLAC-seq/HiChIP), hichipper (restriction-site-distance bias model + library QC), CHiCAGO (Delaporte two-component Brownian+technical background for asymmetric bait x other-end Capture Hi-C), the with/without separate-ChIP anchor decision, and differential loops via diffloop. Use when calling loops from HiChIP/PLAC-seq/Capture Hi-C, choosing FitHiChIP/MAPS/CHiCAGO, picking peak-to-all vs peak-to-peak, setting the loop FDR, supplying ChIP peaks as anchors, QCing a HiChIP library, or comparing loops between conditions.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/hichip-plac-loops)
- **许可证：**`MIT`
- **别名：**`bio-hi-c-analysis-hichip-plac-loops`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-hichip-plac-loops
dnakit skills install hi-c-analysis-hichip-plac-loops --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-loop-calling` {#hi-c-analysis-loop-calling}

Detects focal chromatin loops (point interactions / corner-dots) in balanced Hi-C and Micro-C contact maps and aggregates/validates a loop set. Covers de-novo calling with cooltools dots (HiCCUPS-style 4-background local enrichment with lambda-chunked FDR), chromosight (template-correlation), and Mustache (scale-space blob detection); aggregate peak analysis (APA) via cooltools pileup for confirmation; the depth/resolution prerequisite (de-novo needs ~5-10kb resolution = hundreds of millions to billions of valid pairs); consensus across callers and convergent-CTCF support as validation; and differential loops via union anchors plus chromosight quantify. Use when calling chromatin loops or dots from a cooler, deciding whether a map is deep enough to call de-novo vs running APA on known CTCF/cohesin anchors, building an aggregate peak pileup, comparing loops across conditions, or validating loop calls. For HiChIP/PLAC-seq/PCHi-C protein-anchored data use FitHiChIP/MAPS, not dots.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/loop-calling)
- **许可证：**`MIT`
- **别名：**`bio-hi-c-analysis-loop-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-loop-calling
dnakit skills install hi-c-analysis-loop-calling --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-matrix-operations` {#hi-c-analysis-matrix-operations}

Balances Hi-C contact matrices (ICE via cooler.balance_cooler, KR/SCALE/VC context), computes distance-decay expected with cooltools (expected_cis per-diagonal P(s), expected_trans scalar), builds observed/expected (O/E) matrices, and diagnoses polymer state from the P(s) log-derivative. Covers the within-matrix-vs-cross-sample distinction (balancing is NOT a normalizer), the equal-visibility assumption that CNV/aneuploidy violates (use raw counts for copy-number), cis-only balancing, mad_max/blacklist masking before balancing, multiplicative cooler weights vs divisive juicer weights, and the resolution-vs-depth budget. Use when balancing a .cool/.mcool, computing expected or P(s), making O/E matrices for compartments/loops, deciding ICE vs KR vs SCALE, choosing a resolution for a given depth, or troubleshooting NaN/all-NaN balanced matrices; route cross-sample comparison to hic-differential.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/matrix-operations)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `cooler`, `dnakit`
- **别名：**`bio-hi-c-analysis-matrix-operations`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-matrix-operations
dnakit skills install hi-c-analysis-matrix-operations --agent claude-code   # 或 --agent codex / cursor
```

## `hi-c-analysis-tad-detection` {#hi-c-analysis-tad-detection}

Detects TAD boundaries from balanced Hi-C contact matrices via the diamond-window insulation score (cooltools insulation) and HiCExplorer hicFindTADs, returning a continuous log2 insulation track, valley-prominence boundary_strength, and Li/Otsu-thresholded is_boundary flags across a list of window sizes. Covers the multi-scale window sweep (sub-TAD to compartment-domain), why the boundary is reproducible but the domain partition is not, cross-condition comparison via differential SCORE not differential partition, and the insulation-vs-compartment orthogonality. Use when calling TADs or domain boundaries, computing insulation scores, choosing a window size, ranking boundary strength, comparing boundaries across conditions, or annotating CTCF-backed boundaries; route domain rendering to hic-visualization and boundary-feature overlap to genome-intervals.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/hi-c-analysis/tad-detection)
- **许可证：**`MIT`
- **别名：**`bio-hi-c-analysis-tad-detection`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hi-c-analysis-tad-detection
dnakit skills install hi-c-analysis-tad-detection --agent claude-code   # 或 --agent codex / cursor
```

## `hic-aggregation` {#hic-aggregation}

Build comprehensive chromatin contact maps by aggregating Hi-C loop calls (BEDPE) across multiple ENCODE experiments, donors, and labs. Use when the user wants to answer "what regions are in 3D contact in my tissue?" by creating a union catalog of chromatin loops. Handles resolution-aware anchor matching, cross-lab variation, and Hi-C-specific quality metrics.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/hic-aggregation)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info hic-aggregation
dnakit skills install hic-aggregation --agent claude-code   # 或 --agent codex / cursor
```

## `histone-aggregation` {#histone-aggregation}

Build comprehensive histone mark maps by aggregating narrowPeak data across multiple ENCODE experiments, donors, and labs. Use when the user wants to answer "where is this histone mark present in my tissue?" by combining peak calls from multiple studies into a union peak set with confidence annotations. Handles cross-lab batch effects, broad vs narrow marks, and ENCODE blocklist filtering.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/histone-aggregation)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info histone-aggregation
dnakit skills install histone-aggregation --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-aggregation` {#methylation-aggregation}

Build comprehensive DNA methylation maps by aggregating WGBS (Whole Genome Bisulfite Sequencing) data across multiple ENCODE experiments, donors, and labs. Use when the user wants to answer "where is DNA methylated/unmethylated in my tissue?" by combining per-CpG methylation data into tissue-level methylation profiles. Handles coverage filtering, identifies hypomethylated regions (HMRs) and partially methylated domains (PMDs), and manages cross-lab variation.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/methylation-aggregation)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-aggregation
dnakit skills install methylation-aggregation --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-array-preprocessing` {#methylation-array-preprocessing}

Turns raw Illumina Infinium methylation BeadChip IDATs (450K, EPIC, EPICv2) into a defensible beta/M matrix with sesame (openSesame/SigDF) or minfi (RGChannelSet -&gt; MethylSet -&gt; GenomicRatioSet). Covers Type I vs Type II probe chemistry and why raw Type II beta is compressed, the signal-to-beta math (beta = M/(M+U+100)) and M-value logit, detection-p / pOOBAH masking including the out-of-band deletion-artifact catch, dye-bias correction, and the normalization decision (noob, funnorm, quantile, SWAN, BMIQ, dasen, sesame QCDPB). Use when reading IDATs, choosing a normalization for a 450K/EPIC/EPICv2 cohort, deciding beta vs M, masking failed probes, or producing the corrected matrix before testing. For probe/sample filtering, EPICv2 replicate collapse, and sample-identity QC see array-qc-filtering; for native long-read 5mC see long-read-sequencing/nanopore-methylation (a different platform).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/array-preprocessing)
- **许可证：**`MIT`
- **别名：**`bio-methylation-array-preprocessing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-array-preprocessing
dnakit skills install methylation-array-preprocessing --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-array-qc-filtering` {#methylation-array-qc-filtering}

Performs probe filtering and sample-level QC on Illumina Infinium methylation arrays (450K / EPIC / EPICv2) to decide which probes and samples to trust. Drops detection-p-failed and low-bead-count probes, removes cross-reactive/non-specific probes (Chen 2013 / Pidsley 2016 lists via maxprobes), excludes SNP-overlapping probes with dropLociWithSnps, and handles sex-chromosome probes. Collapses EPICv2 replicate probes with betasCollapseToPfx and harmonizes across array versions (EPICv2 hg38 vs 450K/EPIC hg19, intersect plus mLiftOver). Runs sample-identity QC: getSex sex prediction vs sample sheet for swap detection, rs-SNP fingerprint clustering for duplicates/swaps, and Sentrix chip/array-position batch diagnosis. Use when filtering methylation array probes, detecting sample swaps or mislabels, collapsing EPICv2 replicates, or merging 450K/EPIC/EPICv2 cohorts. For IDAT-to-corrected-beta normalization see array-preprocessing; for batch correction and study design see ewas-design.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/array-qc-filtering)
- **许可证：**`MIT`
- **别名：**`bio-methylation-array-qc-filtering`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-array-qc-filtering
dnakit skills install methylation-array-qc-filtering --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-bismark-alignment` {#methylation-bismark-alignment}

Aligns bisulfite-converted (WGBS, RRBS, PBAT) and enzymatic (EM-seq) short reads to an in-silico C-&gt;T/G-&gt;A-converted reference with Bismark (Bowtie2 or HISAT2), preparing the genome index, choosing the directional vs non-directional vs PBAT strand flag, deduplicating WGBS/EM-seq (never RRBS), and bounding bisulfite conversion efficiency with unmethylated lambda and methylated pUC19 spike-ins. Covers why the library protocol (not the aligner) decides whether calls are meaningful, why incomplete conversion masquerades as methylation, the 3-letter reduced-complexity mapping bias (50-70% efficiency is normal), and M-bias end-clipping. Use when aligning bisulfite or EM-seq reads, preparing a bisulfite genome, choosing the strand flag, or diagnosing low mapping efficiency. For methylation extraction see methylation-calling; for long-read MM/ML modification calling see long-read-sequencing/nanopore-methylation.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/bismark-alignment)
- **许可证：**`MIT`
- **运行依赖：**命令 `Bismark`
- **别名：**`bio-methylation-bismark-alignment`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-bismark-alignment
dnakit skills install methylation-bismark-alignment --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-calling` {#methylation-calling}

Extracts per-cytosine methylation calls from aligned bisulfite/EM-seq reads with bismark_methylation_extractor (Bismark BAM) or the aligner-agnostic MethylDackel/BISCUIT (bwa-meth BAM), producing the beta value M/(M+U) as a coverage file, bedGraph, or genome-wide cytosine report across CpG/CHG/CHH context. Covers conversion-rate QC as the first gate, the 5mC vs 5hmC summed caveat, variant-aware calling so a C/T SNP does not masquerade as unmethylation, paired-end --no_overlap double-counting, symmetric CpG dyad collapse, and the 0-based vs 1-based coordinate trap. Use when extracting methylation levels from a bisulfite/EM-seq alignment, choosing an extractor for a non-Bismark BAM, QC-ing conversion efficiency, or producing coverage/cytosine-report input for testing. For long-read MM/ML modification calling see long-read-sequencing/nanopore-methylation; for the upstream BAM see bismark-alignment; for per-CpG statistics see differential-cpg-testing.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/methylation-calling)
- **许可证：**`MIT`
- **运行依赖：**命令 `Bismark`
- **别名：**`bio-methylation-calling`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-calling
dnakit skills install methylation-calling --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-cell-type-deconvolution` {#methylation-cell-type-deconvolution}

Estimates cell-type composition from bulk DNA methylation and uses it to defuse the single biggest EWAS confounder. Covers reference-based deconvolution (Houseman constrained-projection, minfi estimateCellCounts2 with FlowSorted.Blood.EPIC + IDOL-optimized libraries, EpiDISH RPC/CBS/CP, 12-cell extended, cord-blood nRBC references, EpiSCORE/hepidish for solid tissue), reference-free correction (ReFACTor, RefFreeEWAS, SVA), using fractions as covariates vs the compositionality/collinearity trap, and cell-type-resolved EWAS (CellDMC, TCA, TOAST, omicwas, HIRE). Use when estimating blood/tissue cell fractions, adjusting an EWAS for composition, choosing a deconvolution reference, or attributing a methylation signal to a cell type. For the EWAS confounder-vs-mediator decision see ewas-design; for the IEAA cell-count adjustment of DNAm age see epigenetic-clocks; for clean beta input see array-preprocessing.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/cell-type-deconvolution)
- **许可证：**`MIT`
- **别名：**`bio-methylation-cell-type-deconvolution`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-cell-type-deconvolution
dnakit skills install methylation-cell-type-deconvolution --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-differential-cpg` {#methylation-differential-cpg}

Tests individual CpG sites for differential methylation (DMC/DMP) from bisulfite sequencing counts or array/continuous beta-value matrices. Covers the count-vs-continuous fork that dictates the model, beta-value vs M-value logit (Du 2010), beta-binomial overdispersion count models (DSS, methylKit, MOABS, RADMeth) for sequencing, limma moderated-t on M-values (eBayes trend/robust) for arrays, the bare-beta Welch t-test caveat, coverage-as-precision coupling, delta-beta effect size, BH-FDR with the neighboring-CpG dependence problem, EWAS genome-wide thresholds, and differential variability (DiffVar/iEVORA). Use when comparing per-CpG methylation between groups from WGBS/RRBS/targeted bisulfite or 450K/EPIC arrays, choosing a per-site test, or scanning for variance (not just mean) differences. For region-level aggregation see dmr-detection; for covariate/cell-fraction strategy and genomic inflation see ewas-design.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/differential-cpg-testing)
- **许可证：**`MIT`
- **别名：**`bio-methylation-differential-cpg`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-differential-cpg
dnakit skills install methylation-differential-cpg --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-dmr-detection` {#methylation-dmr-detection}

Detects differentially methylated regions (DMRs) from short-read bisulfite (WGBS/RRBS), array, and long-read methylation count tables using dmrseq (permutation region-FDR over the region selection), DSS callDMR (beta-binomial), methylKit tiles, bsseq BSmooth, DMRcate Gaussian-kernel smoothing, metilene, and comb-p. Covers why a DMR is DEFINED by arbitrary thresholds (min-CpGs, max-gap, delta-beta, q) and a smoothing bandwidth, why selecting extreme runs of CpGs then testing them on the same data is post-selection inference, why region q-values are not comparable across tools, and a single-sample domain-segmentation section (PMD, UMR/LMR, MethylSeekR, solo-WCGW) that must run before focal calling on cancer/aging genomes. Use when calling region-level methylation differences, choosing a DMR caller, controlling region-level FDR, or segmenting megabase methylation domains. For per-site testing see differential-cpg-testing; for the methylKit object model see methylkit-analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/dmr-detection)
- **许可证：**`MIT`
- **别名：**`bio-methylation-dmr-detection`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-dmr-detection
dnakit skills install methylation-dmr-detection --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-epigenetic-clocks` {#methylation-epigenetic-clocks}

Computes DNA methylation age (DNAm age) and pace of aging by applying frozen elastic-net epigenetic clocks to a clean beta matrix with methylclock, dnaMethyAge, or methylCIPHER. Covers the clock menu by question (chronological Horvath/Hannum/skin&amp;blood; health-mortality PhenoAge/GrimAge; DunedinPACE pace; pediatric/gestational; mitotic epiTOC), age acceleration (EAA/IEAA/EEAA) as the real endpoint, the principal-component (PC) clock fix for the per-CpG reliability crisis, and EPICv2 clock-CpG dropout with missing-CpG imputation bias. Use when estimating epigenetic age, computing age acceleration, choosing a clock for an outcome, assessing clock reliability, or porting a clock to EPICv2. A clock is a frozen predictor: do not GO-enrich its CpGs and do not train it here. For cell-count adjustment (IEAA) see cell-type-deconvolution; for predictor training/validation/leakage see machine-learning/model-validation; for survival modeling of age acceleration see clinical-biostatistics/survival-analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/epigenetic-clocks)
- **许可证：**`MIT`
- **别名：**`bio-methylation-epigenetic-clocks`, `methylation-clock`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-epigenetic-clocks
dnakit skills install methylation-epigenetic-clocks --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-ewas-design` {#methylation-ewas-design}

Designs and defends an epigenome-wide association study (EWAS) on 450K/EPIC array or bisulfite methylation - the layer deciding whether a hit is credible. Covers the confounding hierarchy (cell composition covariates as the dominant confounder, batch/Sentrix chip/array position, age/sex, smoking AHRR cg05575921, ancestry/mQTL, reverse causation), chip randomization (no-rescue theorem), surrogate variable analysis sva/SmartSVA, ComBat, RUVm, over-correction, genomic inflation lambda vs GWAS genomic control, BACON bias/inflation, genome-wide significance threshold 450K/EPIC, FWER vs FDR, pwrEWAS power, meta-analysis, EWAS Catalog/Atlas, methylation risk scores. Use when designing an EWAS, choosing a covariate set, randomizing a plate layout, interpreting lambda, applying BACON, setting a threshold, powering a study, or using an MRS. For the per-site test see differential-cpg-testing; for cell fractions see cell-type-deconvolution; for causal mQTL orientation see causal-genomics/mendelian-randomization.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/ewas-design)
- **许可证：**`MIT`
- **别名：**`bio-methylation-ewas-design`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-ewas-design
dnakit skills install methylation-ewas-design --agent claude-code   # 或 --agent codex / cursor
```

## `methylation-methylkit` {#methylation-methylkit}

Imports Bismark coverage or cytosine-report files into the methylKit object model, then runs the import-to-results spine - filterByCoverage, normalizeCoverage, unite/destrand, calculateDiffMeth, getMethylDiff - for both per-CpG (DMC) and fixed-tile (DMR) differential methylation, plus tileMethylCounts, PCA/correlation/clustering QC, and assocComp/removeComp batch handling. Covers the silent default traps that shape the false-positive rate: overdispersion='none' does no correction while 'MN' forces the F-test (ignoring test='Chisq'), adjust defaults to SLIM not BH, getMethylDiff defaults difference=25/qvalue=0.01, cov.bases=0 admits single-CpG tiles, and pool destroys biological replication. Use when importing bisulfite count tables, filtering/normalizing/uniting methylation samples, running methylKit differential testing, or QC-ing methylomes. For per-site test-choice (count vs continuous) see differential-cpg-testing; for selection-aware region FDR (dmrseq/DSS) see dmr-detection.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/methylation-analysis/methylkit-analysis)
- **许可证：**`MIT`
- **别名：**`bio-methylation-methylkit`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info methylation-methylkit
dnakit skills install methylation-methylkit --agent claude-code   # 或 --agent codex / cursor
```

## `pipeline-dnaseseq` {#pipeline-dnaseseq}

Execute ENCODE DNase-seq pipeline from FASTQ to hotspots and footprints. Child of pipeline-guide. Provides Nextflow execution with Docker and cloud deployment. Use when processing DNase-seq data, calling DNase hypersensitive sites, performing footprinting analysis. Trigger on: DNase-seq pipeline, DNase hypersensitive, DHS, Hotspot2, footprinting, DNase I, chromatin accessibility DNase.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/pipeline-dnaseseq)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info pipeline-dnaseseq
dnakit skills install pipeline-dnaseseq --agent claude-code   # 或 --agent codex / cursor
```

## `quality-assessment` {#quality-assessment}

Evaluate ENCODE experiment quality using standard metrics and audit flags. Use when the user asks about data quality, wants to filter for high-quality experiments, needs to interpret quality metrics (FRiP, NSC, RSC, NRF, IDR, TSS enrichment, fragment size), wants to understand ENCODE audit warnings, needs to compare quality across experiments, or is deciding whether data is usable for their analysis. Also use when the user mentions QC, quality control, or data filtering.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/quality-assessment)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info quality-assessment
dnakit skills install quality-assessment --agent claude-code   # 或 --agent codex / cursor
```

## `regulatory-elements` {#regulatory-elements}

Discover and characterize regulatory elements (enhancers, promoters, silencers, insulators, super-enhancers) using ENCODE data and the cCRE catalog. Use when the user wants to find candidate regulatory elements, identify active enhancers in a tissue, map promoter states, classify chromatin states with ChromHMM, identify super-enhancers with ROSE, understand the functional validation hierarchy (CRISPR &gt; MPRA &gt; reporter), or characterize non-coding genomic regions. Covers the full workflow from ENCODE cCRE lookup through chromatin state segmentation to functional validation and enhancer-gene linkage. Use this skill for ANY regulatory element discovery, classification, or characterization task.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/regulatory-elements)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info regulatory-elements
dnakit skills install regulatory-elements --agent claude-code   # 或 --agent codex / cursor
```

## `regulatory-genomics` {#regulatory-genomics}

Transcription factor binding, cis-regulatory elements (cCREs), chromatin accessibility, and regulatory annotation using JASPAR (motifs), ENCODE (cCREs, ChIP-seq), RegulomeDB (regulatory variant scoring), UCSC — plus sequence-based deep-learning prediction of regulatory activity and non-coding variant effects (AlphaGenome, Enformer, Borzoi, ChromBPNet, Evo 2). Use for regulatory element annotation, TF-binding-site prediction, regulatory-region functional impact assessment, and predicting how a non-coding variant or a raw DNA sequence affects expression/chromatin/accessibility. Use this whenever a user asks what regulates a gene, whether a SNP hits a regulatory element, or to predict a non-coding variant's functional effect from sequence.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-regulatory-genomics)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-regulatory-genomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info regulatory-genomics
dnakit skills install regulatory-genomics --agent claude-code   # 或 --agent codex / cursor
```

## `regulatory-variant-analysis` {#regulatory-variant-analysis}

Non-coding/regulatory variant interpretation — GWAS association lookup, eQTL evidence (GTEx), chromatin state (ENCODE), regulatory variant scoring (RegulomeDB, CADD), and TF-binding disruption. Use for non-coding GWAS hit interpretation, eQTL-based gene assignment, and regulatory mechanism reasoning. Distinct from coding-variant tools.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-regulatory-variant-analysis)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-regulatory-variant-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info regulatory-variant-analysis
dnakit skills install regulatory-variant-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `scvi-tools` {#scvi-tools}

Deep learning for single-cell analysis using scvi-tools. This skill should be used when users need (1) data integration and batch correction with scVI/scANVI, (2) ATAC-seq analysis with PeakVI, (3) CITE-seq multi-modal analysis with totalVI, (4) multiome RNA+ATAC analysis with MultiVI, (5) spatial transcriptomics deconvolution with DestVI, (6) label transfer and reference mapping with scANVI/scArches, (7) RNA velocity with veloVI, or (8) any deep learning-based single-cell method. Triggers include mentions of scVI, scANVI, totalVI, PeakVI, MultiVI, DestVI, veloVI, sysVI, scArches, variational autoencoder, VAE, batch correction, data integration, multi-modal, CITE-seq, multiome, reference mapping, latent space.

- **来源：**[Anthropic and life-sciences skill contributors; leaf LICENSE.txt retained @ `e96556b`](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c/scvi-tools)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `scvi`, `anndata`, `scanpy`; 命令 `python3`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info scvi-tools
dnakit skills install scvi-tools --agent claude-code   # 或 --agent codex / cursor
```

## `ucsc-conservation-and-tfbs` {#ucsc-conservation-and-tfbs}

Fetch Evolutionary Conservation scores (phyloP, phastCons) and Transcription Factor Binding Sites (TFBS) from the UCSC Genome Browser. Use when analyzing whether genomic variants or regions are evolutionarily conserved, functionally important, or bounded by TF regulators across major projects (ENCODE, JASPAR, ReMap).

- **来源：**[Copyright 2026 Google LLC @ `28b8482`](https://github.com/google-deepmind/science-skills/tree/28b8482603a420708c8896f6fe5e06c276d9933d/skills/ucsc_conservation_and_tfbs)
- **许可证：**`Apache-2.0 AND CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info ucsc-conservation-and-tfbs
dnakit skills install ucsc-conservation-and-tfbs --agent claude-code   # 或 --agent codex / cursor
```

## `variant-to-mechanism` {#variant-to-mechanism}

End-to-end variant-to-mechanism analysis — trace a variant (rsID/coordinates) through regulatory context, target gene(s), molecular pathway(s), and phenotypic consequences. Integrates 7+ databases across 3 evidence layers (regulatory, molecular, disease) for a mechanistic model. Use for GWAS-hit-to-mechanism, eQTL-causal-gene tracing, and full causal-chain reports.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-variant-to-mechanism)
- **许可证：**`Apache-2.0`
- **别名：**`tooluniverse-variant-to-mechanism`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info variant-to-mechanism
dnakit skills install variant-to-mechanism --agent claude-code   # 或 --agent codex / cursor
```
