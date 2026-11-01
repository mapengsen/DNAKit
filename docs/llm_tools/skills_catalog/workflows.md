# 组合流程与技能入口

本类共 **77** 个 Skills。每个 Skill 是一份让 Agent 按步骤完成分析的说明，分析所需的软件、参考数据和 API 需要另外准备。搜索、安装和使用方法见 [Skills](../skills.md)；网页版的 Skills 标签可以直接阅读每个 Skill 的完整内容。

| Skill | 来源 | 简介 |
| --- | --- | --- |
| [`batch-analysis`](#batch-analysis) | encode-toolkit | Guide for multi-experiment batch operations: QC screening, batch download, comparison, and report generation across many ENCODE experiments simultaneously. Use… |
| [`benchling-integration`](#benchling-integration) | sciagent | Benchling R&amp;D Python SDK: CRUD on registry entities (DNA, RNA, proteins, custom), inventory, ELN, workflow automation. Needs Benchling account and API key. Use… |
| [`bioconductor-bridge`](#bioconductor-bridge) | user335-clawbio | Bioconductor package discovery, workflow recommendation, setup inspection, and starter code generation grounded in official Bioconductor containers and BiocMan… |
| [`bioinfo-pitfalls`](#bioinfo-pitfalls) | plant-skills | Bioinformatics troubleshooting: HWE, multiple testing, aDNA damage, VCF merge traps, Slurm, array QC. Use when a pipeline looks significant but wrong, a tool e… |
| [`bioinformatics-methods-results-writer`](#bioinformatics-methods-results-writer) | damlab | Write Methods and Results sections for bioinformatics manuscripts from code, notebooks, logs, figures, and tables. Use when drafting or revising computational… |
| [`bioinformatics-project`](#bioinformatics-project) | omics | Structure reproducible bioinformatics projects with canonical layouts, restartable drivers, pinned environments, provenance, and lab notebooks. Use when starti… |
| [`cite-encode`](#cite-encode) | encode-toolkit | Generate proper ENCODE citations for publications, grants, and presentations. Use when the user needs to cite ENCODE data, create bibliography entries, write a… |
| [`data-provenance`](#data-provenance) | encode-toolkit | Track exact provenance for every operation on ENCODE data — tool versions, reference files, scripts, parameters, and timestamps — to enable publication-ready m… |
| [`dnakit`](#dnakit) | dnakit | Find and use curated DNA and genomics workflows from the DNAKit skill catalog. Use for DNA sequence processing, sequencing QC, alignment, assembly, variant ana… |
| [`dnakit-external-genomics`](#dnakit-external-genomics) | dnakit | Discover reviewed DNA projects, run registered external genomics analyses and connect to registered scientific MCP servers through DNAKit. Use when a workflow… |
| [`dnakit-plan`](#dnakit-plan) | dnakit | 在执行前规划多步骤 DNAKit 基因组学任务。适用于用户希望 Codex 检查 DNAKit 能力、确认关键选择，并保存经确认的 Markdown 计划供后续会话执行的情况；不用于立即执行任务或回答简单事实问题。 |
| [`dnakit-tool-generation`](#dnakit-tool-generation) | dnakit | 为 DNAKit 设计、实现和验证可复用的新工具。先检查现有工具和工作流能否满足需求，再实现一个或多个候选工具；适用于支持 DNAKit MCP 或等效工具接口的 Agent。 |
| [`dnanexus-integration`](#dnanexus-integration) | kdense | Build and operate reproducible genomics workloads on DNAnexus with the dx CLI, dxpy, apps/applets, native workflows, dxCompiler, and Nextflow. Use for DNAnexus… |
| [`encode-toolkit-setup`](#encode-toolkit-setup) | encode-toolkit | Set up ENCODE access through DNAKit, including MCP discovery, public queries, credentials and connection troubleshooting. |
| [`experimental-design-batch-design`](#experimental-design-batch-design) | bioskills | Designs genomics experiments so technical nuisance variation (batch, lane, plate, flow cell, operator, reagent lot, processing day) is balanced against the bio… |
| [`experimental-design-multiple-testing`](#experimental-design-multiple-testing) | bioskills | Controls error rates across thousands of simultaneous tests in genomics discovery using false-discovery-rate methods (Benjamini-Hochberg 1995; Benjamini-Yekuti… |
| [`experimental-design-power-analysis`](#experimental-design-power-analysis) | bioskills | Calculates statistical power for high-dimensional genomics experiments (bulk RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, proteomics) under negative-bi… |
| [`experimental-design-randomization-blocking`](#experimental-design-randomization-blocking) | bioskills | Structures biological experiments so inference is valid by construction, covering Fisher's principles (randomization, replication, local control), the experime… |
| [`experimental-design-sample-size`](#experimental-design-sample-size) | bioskills | Estimates the minimum biological replicates (or cells/events) for a target power at a target FDR in genomics experiments using ssizeRNA, PROPER, powsimR for sc… |
| [`foundation-housekeeping`](#foundation-housekeeping) | omics | Add schema-backed metadata validation, normalized Parquet tables, and a DuckDB catalog to a bioinformatics project. Use when an analysis needs LinkML/Pydantic… |
| [`galaxy-bioblend`](#galaxy-bioblend) | user335-galaxy | Use Galaxy REST or the BioBlend Python SDK to manage genomic analysis histories, upload or inspect datasets, discover tools and invoke existing workflows on a… |
| [`galaxy-bridge`](#galaxy-bridge) | user335-clawbio | Galaxy tool discovery, intelligent recommendation, and execution — 8,000+ bioinformatics tools from usegalaxy.org with multi-signal scoring and workflow sugges… |
| [`genomics-study-design`](#genomics-study-design) | research-suite | Design and audit non-model, polyploid, duplicated, or non-human genomic studies. Use when a task needs reference assemblies, truth sets, mapping, variant disco… |
| [`illumina-bridge`](#illumina-bridge) | clawbio | Import DRAGEN-exported Illumina result bundles into ClawBio for local tertiary analysis and downstream routing. |
| [`isabl-debug-analysis`](#isabl-debug-analysis) | isabl | Systematically debug a failed Isabl analysis. Use when an analysis has FAILED status or unexpected behavior. |
| [`isabl-merge-results`](#isabl-merge-results) | isabl | Aggregate results from multiple Isabl analyses into combined datasets. Use when merging VCFs, TSVs, or other outputs across samples or cohorts. |
| [`isabl-monitor-analyses`](#isabl-monitor-analyses) | isabl | Monitor and track the status of Isabl analyses. Use when checking job progress, finding failed analyses, or tracking pipeline execution. |
| [`isabl-project-report`](#isabl-project-report) | isabl | Generate a status report for an Isabl project. Use when summarizing project progress, identifying issues, or preparing updates. |
| [`isabl-query-data`](#isabl-query-data) | isabl | Help construct queries to retrieve data from Isabl using the Python SDK. Use when searching for experiments, analyses, samples, or other data. |
| [`isabl-run-pipeline`](#isabl-run-pipeline) | isabl | Submit and run multiple Isabl applications as a pipeline. Use when chaining apps together or running a sequence of analyses on samples. |
| [`isabl-submit-data`](#isabl-submit-data) | isabl | Submit sequencing data to Isabl by creating individuals, samples, and experiments. Use when registering new data files or importing from external sources. |
| [`isabl-write-app`](#isabl-write-app) | isabl | Guide through creating a new Isabl bioinformatics application. Use when building pipelines that integrate with the Isabl platform. |
| [`kiro-healthomics-genomics`](#kiro-healthomics-genomics) | user335-aws-kiro | Prepare DNA sequencing and variant workflows for AWS HealthOmics using selected upstream Kiro guidance. Use for WDL/Nextflow import, genomic file-format checks… |
| [`multi-omics-data-harmonization`](#multi-omics-data-harmonization) | bioskills | Harmonizes already-normalized per-omic matrices onto a common footing before joint integration - assembling a MultiAssayExperiment, choosing the per-omic varia… |
| [`multi-omics-integration`](#multi-omics-integration) | tooluniverse | Multi-omics integration — orchestrate per-layer analysis (transcriptomics, proteomics, epigenomics, genomics, metabolomics) then perform cross-omics correlatio… |
| [`multi-omics-mixomics-analysis`](#multi-omics-mixomics-analysis) | bioskills | Builds supervised and unsupervised multivariate integration across bulk omics blocks with mixOmics - sPLS for sparse pairwise correlation, DIABLO (block.splsda… |
| [`multi-omics-mofa-integration`](#multi-omics-mofa-integration) | bioskills | Discovers shared and view-specific latent factors across bulk multi-omics blocks (RNA-seq, proteomics, methylation) on a common sample axis with MOFA2's unsupe… |
| [`nextflow-development`](#nextflow-development) | user335-anthropic | Run nf-core bioinformatics pipelines (rnaseq, sarek, atacseq) on sequencing data. Use when analyzing RNA-seq, WGS/WES, or ATAC-seq data—either local FASTQs or… |
| [`nfcore-sarek-wrapper`](#nfcore-sarek-wrapper) | clawbio | ClawBio wrapper around nf-core/sarek 3.8.1 covering mapping through annotation for germline, tumor-only, and somatic paired analyses. |
| [`pipeline-guide`](#pipeline-guide) | encode-toolkit | Access ENCODE uniform analysis pipelines, generate user-specific Nextflow/WDL pipelines, manage compute resources, and integrate with cloud platforms. Use when… |
| [`plant-omics-db`](#plant-omics-db) | plant-skills | Build a journal-ready plant multi-omics database (genome, germplasm, variome, transcriptome, GWAS, tools). Use when the user wants a NAR/MP/PC class plant geno… |
| [`prefect-dask-nextflow`](#prefect-dask-nextflow) | omics | Design reproducible bioinformatics pipelines with Prefect plus Dask or Nextflow. Use when scaffolding local, distributed, or scheduler-backed workflows. |
| [`seq-wrangler`](#seq-wrangler) | clawbio | NGS read QC, alignment, and BAM processing pipeline. Wraps FastQC, BWA/Bowtie2/Minimap2, SAMtools, and MultiQC for automated read-to-BAM workflows. |
| [`single-cell-clustering`](#single-cell-clustering) | bioskills | Dimensionality reduction and graph-based clustering for single-cell RNA-seq with Scanpy (Python) and Seurat (R). Resolves which algorithm to use (Leiden vs Lou… |
| [`single-cell-data-io`](#single-cell-data-io) | bioskills | Read, write, create, and convert single-cell objects across AnnData (Python), Seurat (R), and SingleCellExperiment (R). Use when loading 10X Cell Ranger output… |
| [`single-cell-multimodal-integration`](#single-cell-multimodal-integration) | bioskills | Integrate multimodal single-cell data (CITE-seq RNA+protein, 10x Multiome RNA+ATAC, unpaired/diagonal RNA+ATAC) and choose the right joint method. Use when cla… |
| [`single-cell-preprocessing`](#single-cell-preprocessing) | bioskills | Quality control, ambient-RNA handling, normalization, and feature selection for single-cell RNA-seq using Scanpy (Python) and Seurat (R). Use when filtering lo… |
| [`track-experiments`](#track-experiments) | encode-toolkit | Track ENCODE experiments locally with publications, citations, and provenance. Use when the user wants to build a collection of experiments, manage citations,… |
| [`wgs-prs`](#wgs-prs) | clawbio | End-to-end WGS to polygenic risk score pipeline. Takes paired-end FASTQ files (or a pre-existing VCF) through nf-core/sarek for variant calling, applies VCF QC… |
| [`workflow-management-cwl-workflows`](#workflow-management-cwl-workflows) | bioskills | Authors portable, strongly-typed bioinformatics pipelines in the Common Workflow Language (CWL v1.2) as CommandLineTool/Workflow/ExpressionTool documents, vali… |
| [`workflow-management-nextflow-pipelines`](#workflow-management-nextflow-pipelines) | bioskills | Authors reproducible Nextflow DSL2 pipelines built on reactive dataflow, where processes communicate only through channels and execution order is not guarantee… |
| [`workflow-management-nf-core-pipelines`](#workflow-management-nf-core-pipelines) | bioskills | Runs and configures curated nf-core community Nextflow pipelines (rnaseq, sarek, atacseq, methylseq, ampliseq, taxprofiler, fetchngs) reproducibly, pinning the… |
| [`workflow-management-snakemake-workflows`](#workflow-management-snakemake-workflows) | bioskills | Authors reproducible bioinformatics pipelines with Snakemake - rules wired by output-file pattern, wildcards and expand() for sample fan-out, checkpoints for r… |
| [`workflow-management-wdl-workflows`](#workflow-management-wdl-workflows) | bioskills | Authors bioinformatics pipelines in WDL (Workflow Description Language) run by Cromwell or miniwdl, targeting the GATK/Broad and Terra/AnVIL/BioData Catalyst c… |
| [`workflow-methods-docwriter`](#workflow-methods-docwriter) | omics | Generate reproducible Methods from Nextflow, Snakemake, or CWL run artifacts. Use when documenting exact commands, versions, parameters, QC gates, provenance,… |
| [`workflows-atacseq-pipeline`](#workflows-atacseq-pipeline) | bioskills | Orchestrates the end-to-end bulk ATAC-seq pipeline from FASTQ to differential accessibility and TF footprints, chaining Nextera-aware fastp QC, Bowtie2 alignme… |
| [`workflows-causal-genomics-pipeline`](#workflows-causal-genomics-pipeline) | bioskills | End-to-end post-GWAS causal inference pipeline orchestrating heritability partitioning, genetic correlation, Mendelian randomization with CHP-aware sensitivity… |
| [`workflows-chipseq-pipeline`](#workflows-chipseq-pipeline) | bioskills | Orchestrates the end-to-end ChIP-seq pipeline from FASTQ to blacklist-filtered, annotated peaks, chaining fastp QC, Bowtie2 alignment, pre-dedup library-comple… |
| [`workflows-cnv-pipeline`](#workflows-cnv-pipeline) | bioskills | Orchestrates the copy-number pipeline from BAM to segmented, integer-called, annotated CNVs, forking on germline-vs-somatic - CNVkit (somatic exome/panel: cove… |
| [`workflows-crispr-editing-pipeline`](#workflows-crispr-editing-pipeline) | bioskills | Orchestrates an end-to-end CRISPR editing experiment design from target gene to delivery-ready, validatable constructs. Sequences guide design, off-target asse… |
| [`workflows-crispr-screen-pipeline`](#workflows-crispr-screen-pipeline) | bioskills | End-to-end pooled and single-cell CRISPR screen analysis from FASTQ to hit genes. Orchestrates library design QC, guide counting, six-stage screen QC (plasmid… |
| [`workflows-edna-pipeline`](#workflows-edna-pipeline) | bioskills | End-to-end eDNA metabarcoding from raw amplicons to community ecology. Covers QC, primer removal (mandatory before DADA2 filterAndTrim), denoising with OBITool… |
| [`workflows-fastq-to-variants`](#workflows-fastq-to-variants) | bioskills | Orchestrates the end-to-end germline short-variant pipeline from FASTQ to a filtered, normalized, benchmarked VCF, chaining QC/trim, BWA-MEM2 alignment, duplic… |
| [`workflows-genome-annotation-pipeline`](#workflows-genome-annotation-pipeline) | bioskills | Orchestrates genome annotation from assembled contigs to functional annotation, forking prokaryotic (Bakta one-step, genetic-code table from GTDB-Tk) vs eukary… |
| [`workflows-genome-assembly-pipeline`](#workflows-genome-assembly-pipeline) | bioskills | Orchestrates an end-to-end de novo genome assembly project, routing each step to the right genome-assembly skill rather than restating it. Profiles the genome… |
| [`workflows-gwas-pipeline`](#workflows-gwas-pipeline) | bioskills | Orchestrates the GWAS pipeline from genotypes to association results, chaining PLINK2 QC (variant-then-sample missingness, controls-only HWE, KING relatedness)… |
| [`workflows-hic-pipeline`](#workflows-hic-pipeline) | bioskills | End-to-end Hi-C analysis workflow from FASTQ to compartments, TADs, and loops, with the decision of WHICH features the sequencing depth can support. Covers pai… |
| [`workflows-liquid-biopsy-pipeline`](#workflows-liquid-biopsy-pipeline) | bioskills | Orchestrates the cell-free DNA / liquid-biopsy pipeline from plasma sequencing to tumor monitoring, forking tumor-naive (screening) vs tumor-informed (MRD), an… |
| [`workflows-longread-sv-pipeline`](#workflows-longread-sv-pipeline) | bioskills | Orchestrates an end-to-end long-read structural-variant pipeline - basecalling to minimap2 alignment (platform-matched preset) to Sniffles2/cuteSV/pbsv calling… |
| [`workflows-metagenomics-pipeline`](#workflows-metagenomics-pipeline) | bioskills | End-to-end shotgun metagenomics workflow from FASTQ to taxonomic and functional profiles, orchestrating controls/host depletion, Kraken2+Bracken classification… |
| [`workflows-methylation-pipeline`](#workflows-methylation-pipeline) | bioskills | Orchestrates the end-to-end bisulfite/EM-seq methylation pipeline from FASTQ to differentially methylated regions, chaining Trim Galore/fastp QC, Bismark align… |
| [`workflows-microbiome-pipeline`](#workflows-microbiome-pipeline) | bioskills | End-to-end 16S/ITS amplicon workflow from demultiplexed FASTQ to a consensus differential-abundance result, orchestrating cutadapt primer removal, per-run DADA… |
| [`workflows-multiome-pipeline`](#workflows-multiome-pipeline) | bioskills | Orchestrates the end-to-end 10x Multiome (paired scRNA + scATAC) pipeline from Cell Ranger ARC output to a jointly-embedded, annotated object, chaining per-mod… |
| [`workflows-neoantigen-pipeline`](#workflows-neoantigen-pipeline) | bioskills | Orchestrates neoantigen discovery from somatic variants to ranked vaccine candidates, chaining HLA typing (OptiType/arcasHLA + LOHHLA), VEP annotation (Wildtyp… |
| [`workflows-outbreak-pipeline`](#workflows-outbreak-pipeline) | bioskills | Orchestrates genomic-epidemiology outbreak investigation from pathogen isolates to transmission networks, forking bacterial (snippy -&gt; Gubbins recombination-ma… |
| [`workflows-somatic-variant-pipeline`](#workflows-somatic-variant-pipeline) | bioskills | Chains a somatic (tumor-normal) SNV/indel and structural-variant pipeline end to end with GATK Mutect2 (or Strelka2), wiring the somatic-specific machinery - p… |
| [`workflows-tcr-pipeline`](#workflows-tcr-pipeline) | bioskills | Orchestrates an end-to-end immune-repertoire pipeline from FASTQ to clonotypes, diversity, overlap, somatic hypermutation and lineages, routing on two forks. U… |

## `batch-analysis` {#batch-analysis}

Guide for multi-experiment batch operations: QC screening, batch download, comparison, and report generation across many ENCODE experiments simultaneously. Use when users need to process 5+ experiments together, create experiment comparison tables, perform batch quality checks, or generate summary reports. Trigger on: batch analysis, multiple experiments, bulk processing, experiment comparison, batch QC, multi-sample, batch download, experiment table, summary report, collection analysis.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/batch-analysis)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info batch-analysis
dnakit skills install batch-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `benchling-integration` {#benchling-integration}

Benchling R&amp;D Python SDK: CRUD on registry entities (DNA, RNA, proteins, custom), inventory, ELN, workflow automation. Needs Benchling account and API key. Use biopython for local sequence analysis; pubchem for chemical DBs.

- **来源：**[jaechang-hits; credited K-Dense-AI and Biomni contributors @ `fe505ca`](https://github.com/jaechang-hits/SciAgent-Skills/tree/fe505cae14d20b6c33be2e49666425be98f005bb/skills/lab-automation/benchling-integration)
- **许可证：**`CC-BY-4.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info benchling-integration
dnakit skills install benchling-integration --agent claude-code   # 或 --agent codex / cursor
```

## `bioconductor-bridge` {#bioconductor-bridge}

Bioconductor package discovery, workflow recommendation, setup inspection, and starter code generation grounded in official Bioconductor containers and BiocManager.

- **来源：**[Hiranyamaya Dash; Copyright 2026 Manuel Corpas (ClawBio license) @ `4848cc2`](https://github.com/ClawBio/ClawBio/tree/4848cc280c446093b1364a2445db3137d455cc17/skills/bioconductor-bridge)
- **许可证：**`MIT`
- **运行依赖：**命令 `python3`, `Rscript`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info bioconductor-bridge
dnakit skills install bioconductor-bridge --agent claude-code   # 或 --agent codex / cursor
```

## `bioinfo-pitfalls` {#bioinfo-pitfalls}

Bioinformatics troubleshooting: HWE, multiple testing, aDNA damage, VCF merge traps, Slurm, array QC. Use when a pipeline looks significant but wrong, a tool errors on macOS/HPC, or the user asks 踩坑 / 排错.

- **来源：**[Xuzhen-Li/bioinfo-agent-skills contributors @ `6dea4be`](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83/skills/bioinfo-pitfalls)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info bioinfo-pitfalls
dnakit skills install bioinfo-pitfalls --agent claude-code   # 或 --agent codex / cursor
```

## `bioinformatics-methods-results-writer` {#bioinformatics-methods-results-writer}

Write Methods and Results sections for bioinformatics manuscripts from code, notebooks, logs, figures, and tables. Use when drafting or revising computational biology prose, especially to preserve reproducibility details, figure/table references, and concise results reporting. Triggers on: "methods section", "results section", "manuscript", "write methods", "write results", "draft methods", "draft results", "publication", "paper", "scientific prose", "reproducibility", "methods and results".

- **来源：**[DamLabResources/damlab-skills contributors @ `4ad096e`](https://github.com/DamLabResources/damlab-skills/tree/4ad096e3b9b0116271d10bdb746641b60a1b2211/skills/bioinformatics-methods-results-writer)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info bioinformatics-methods-results-writer
dnakit skills install bioinformatics-methods-results-writer --agent claude-code   # 或 --agent codex / cursor
```

## `bioinformatics-project` {#bioinformatics-project}

Structure reproducible bioinformatics projects with canonical layouts, restartable drivers, pinned environments, provenance, and lab notebooks. Use when starting or reorganizing a genomics project or making a sequencing analysis rerunnable.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bioinformatics-project)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info bioinformatics-project
dnakit skills install bioinformatics-project --agent claude-code   # 或 --agent codex / cursor
```

## `cite-encode` {#cite-encode}

Generate proper ENCODE citations for publications, grants, and presentations. Use when the user needs to cite ENCODE data, create bibliography entries, write acknowledgment sections, or ensure compliance with ENCODE data use policy.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/cite-encode)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info cite-encode
dnakit skills install cite-encode --agent claude-code   # 或 --agent codex / cursor
```

## `data-provenance` {#data-provenance}

Track exact provenance for every operation on ENCODE data — tool versions, reference files, scripts, parameters, and timestamps — to enable publication-ready methods writing. Use when the user processes ENCODE files, runs any bioinformatics tool, creates filtered/merged datasets, runs pipelines, performs liftover, uses R/Python/Bash for analysis, or needs to document their analysis chain for reproducibility and publication. Also use when the user says "write me methods" to auto-generate methods sections from the provenance log. This skill implements comprehensive provenance documentation: every tool, every version, every reference file, every parameter, every accession — no shortcuts. Use this skill for ANY processing step, ANY file transformation, ANY analysis operation on ENCODE data.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/data-provenance)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info data-provenance
dnakit skills install data-provenance --agent claude-code   # 或 --agent codex / cursor
```

## `dnakit` {#dnakit}

Find and use curated DNA and genomics workflows from the DNAKit skill catalog. Use for DNA sequence processing, sequencing QC, alignment, assembly, variant analysis, GWAS, population genetics, epigenomics, CRISPR, and microbial genomics when the user needs an analysis workflow or wants to find an appropriate skill.

- **来源：**[DNAKit contributors @ `catalog`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit)
- **许可证：**`MIT`
- **别名：**`dna-skills`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnakit
dnakit skills install dnakit --agent claude-code   # 或 --agent codex / cursor
```

## `dnakit-external-genomics` {#dnakit-external-genomics}

Discover reviewed DNA projects, run registered external genomics analyses and connect to registered scientific MCP servers through DNAKit. Use when a workflow needs an optional executable or external service.

- **来源：**[DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-external-genomics)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnakit-external-genomics
dnakit skills install dnakit-external-genomics --agent claude-code   # 或 --agent codex / cursor
```

## `dnakit-plan` {#dnakit-plan}

在执行前规划多步骤 DNAKit 基因组学任务。适用于用户希望 Codex 检查 DNAKit 能力、确认关键选择，并保存经确认的 Markdown 计划供后续会话执行的情况；不用于立即执行任务或回答简单事实问题。

- **来源：**[DNAKit contributors @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-plan)
- **许可证：**`MIT`
- **运行依赖：**命令 `codex`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnakit-plan
dnakit skills install dnakit-plan --agent claude-code   # 或 --agent codex / cursor
```

## `dnakit-tool-generation` {#dnakit-tool-generation}

为 DNAKit 设计、实现和验证可复用的新工具。先检查现有工具和工作流能否满足需求，再实现一个或多个候选工具；适用于支持 DNAKit MCP 或等效工具接口的 Agent。

- **来源：**[DNAKit contributors @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-tool-generation)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnakit-tool-generation
dnakit skills install dnakit-tool-generation --agent claude-code   # 或 --agent codex / cursor
```

## `dnanexus-integration` {#dnanexus-integration}

Build and operate reproducible genomics workloads on DNAnexus with the dx CLI, dxpy, apps/applets, native workflows, dxCompiler, and Nextflow. Use for DNAnexus data transfers, dxapp.json development, execution monitoring, workflow import, and project automation.

- **来源：**[K-Dense Inc.; K-Dense Inc. and contributors @ `c1ed16d`](https://github.com/K-Dense-AI/scientific-agent-skills/tree/c1ed16d97dd61ff50a3bd46dd353e4a55fd77f34/skills/dnanexus-integration)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info dnanexus-integration
dnakit skills install dnanexus-integration --agent claude-code   # 或 --agent codex / cursor
```

## `encode-toolkit-setup` {#encode-toolkit-setup}

Set up ENCODE access through DNAKit, including MCP discovery, public queries, credentials and connection troubleshooting.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/setup)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info encode-toolkit-setup
dnakit skills install encode-toolkit-setup --agent claude-code   # 或 --agent codex / cursor
```

## `experimental-design-batch-design` {#experimental-design-batch-design}

Designs genomics experiments so technical nuisance variation (batch, lane, plate, flow cell, operator, reagent lot, processing day) is balanced against the biological variable of interest and therefore estimable rather than confounded, using constrained sample-to-batch assignment (designit, OSAT), the confounder/mediator/collider distinction, and the principle that no post-hoc correction recovers a fully confounded design. Covers detecting hidden batches with surrogate variable analysis, a decision table for downstream correction (ComBat-seq, RUVSeq, SVA) whose execution is deferred to differential-expression/batch-correction, and reproducibility metadata. Use when assigning samples to sequencing batches/lanes/plates, avoiding batch-condition confounding, deciding whether a design is salvageable by correction, choosing a correction method, or estimating the number of hidden batches. For the experimental unit, randomization, and blocking concepts see experimental-design/randomization-blocking.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/experimental-design/batch-design)
- **许可证：**`MIT`
- **别名：**`bio-experimental-design-batch-design`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info experimental-design-batch-design
dnakit skills install experimental-design-batch-design --agent claude-code   # 或 --agent codex / cursor
```

## `experimental-design-multiple-testing` {#experimental-design-multiple-testing}

Controls error rates across thousands of simultaneous tests in genomics discovery using false-discovery-rate methods (Benjamini-Hochberg 1995; Benjamini-Yekutieli 2001 for arbitrary dependence; Storey q-value with pi0 estimation; local FDR; independent filtering Bourgon 2010; covariate-weighted FDR via IHW Ignatiadis 2016), plus family-wise error control (Bonferroni, Holm) and the GWAS genome-wide threshold. Covers the FDR-versus-FWER choice as the discovery-versus-confirmatory distinction, the dependence assumptions behind BH (PRDS) versus BY, pi0 estimation, the independent-filtering and false-coverage-rate traps, and reproducibility ranking via IDR (Li 2011). Use when correcting p-values from genome-wide tests, choosing between BH/BY/q-value/Bonferroni, setting an FDR threshold, applying IHW or independent filtering, or interpreting q-values. For confirmatory trials with few pre-specified endpoints (closed testing, graphical/gatekeeping), see clinical-biostatistics/multiplicity-graphical.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/experimental-design/multiple-testing)
- **许可证：**`MIT`
- **别名：**`bio-experimental-design-multiple-testing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info experimental-design-multiple-testing
dnakit skills install experimental-design-multiple-testing --agent claude-code   # 或 --agent codex / cursor
```

## `experimental-design-power-analysis` {#experimental-design-power-analysis}

Calculates statistical power for high-dimensional genomics experiments (bulk RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, proteomics) under negative-binomial count models using RNASeqPower, PROPER, and simulation via powsimR, distinguishing per-gene from marginal (transcriptome-wide) power, the role of mean expression and dispersion, and the sequencing-depth-versus-replicate tradeoff. Covers simulation as the honest default for overdispersed counts, FDR-aware average power versus single-test power, observed/post-hoc power as an anti-pattern, and the winner's-curse / Type-S / Type-M consequences of underpowering. Use when planning replicate number for a sequencing experiment, deciding whether to add depth or samples, choosing closed-form versus simulation power, estimating power from pilot dispersions, or justifying replication in a grant. For clinical-trial power see clinical-biostatistics/power-and-sample-size; for the inverse sample-size question see experimental-design/sample-size.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/experimental-design/power-analysis)
- **许可证：**`MIT`
- **别名：**`bio-experimental-design-power-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info experimental-design-power-analysis
dnakit skills install experimental-design-power-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `experimental-design-randomization-blocking` {#experimental-design-randomization-blocking}

Structures biological experiments so inference is valid by construction, covering Fisher's principles (randomization, replication, local control), the experimental-vs-observational unit distinction and pseudoreplication (Hurlbert 1984; Lazic 2018), randomization mechanics (complete, restricted, stratified, rerandomization, run-order), blocking layouts (randomized complete block, Latin square, incomplete block), factorial designs and interactions, and the split-plot/nested error strata hidden inside multi-batch genomics. Use when deciding the experimental unit and what counts as a replicate, planning randomization and run order, choosing a blocked/factorial/split-plot/nested layout, avoiding pseudoreplication in cell-culture or animal studies, or specifying the random-effects structure of the analysis model. For assigning samples to sequencing batches/lanes/plates and batch-effect correction see experimental-design/batch-design; for regulated clinical-trial randomization see clinical-biostatistics.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/experimental-design/randomization-blocking)
- **许可证：**`MIT`
- **别名：**`bio-experimental-design-randomization-blocking`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info experimental-design-randomization-blocking
dnakit skills install experimental-design-randomization-blocking --agent claude-code   # 或 --agent codex / cursor
```

## `experimental-design-sample-size` {#experimental-design-sample-size}

Estimates the minimum biological replicates (or cells/events) for a target power at a target FDR in genomics experiments using ssizeRNA, PROPER, powsimR for scRNA-seq, and pilot-data dispersion estimation from DESeq2/edgeR. Covers the biological-versus-technical replication distinction (technical replicates do not add degrees of freedom for biological inference), replicate-number-versus-sequencing-depth budgeting, scRNA-seq sample-versus-cell allocation under a pseudobulk model, and the critique that "n=3" is a publication convention rather than a power calculation. Use when budgeting a sequencing experiment, writing the sample-size justification in a grant, estimating replicates from pilot data, allocating a fixed budget between samples and depth, or planning scRNA-seq cohort size. For clinical-trial sample size see clinical-biostatistics/power-and-sample-size; for the power-given-n direction see experimental-design/power-analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/experimental-design/sample-size)
- **许可证：**`MIT`
- **别名：**`bio-experimental-design-sample-size`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info experimental-design-sample-size
dnakit skills install experimental-design-sample-size --agent claude-code   # 或 --agent codex / cursor
```

## `foundation-housekeeping` {#foundation-housekeeping}

Add schema-backed metadata validation, normalized Parquet tables, and a DuckDB catalog to a bioinformatics project. Use when an analysis needs LinkML/Pydantic records or a queryable data catalog.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-foundation-housekeeping)
- **许可证：**`MIT`
- **别名：**`bio-foundation-housekeeping`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info foundation-housekeeping
dnakit skills install foundation-housekeeping --agent claude-code   # 或 --agent codex / cursor
```

## `galaxy-bioblend` {#galaxy-bioblend}

Use Galaxy REST or the BioBlend Python SDK to manage genomic analysis histories, upload or inspect datasets, discover tools and invoke existing workflows on a user-selected Galaxy server. Includes the original REST and SDK guides; server tools, account access and scientific workflows are external.

- **来源：**[Copyright 2025 Galaxy MCP Contributors @ `0f833c2`](https://github.com/nekrut/galaxy-claude-skills/tree/0f833c2f11fa624d9fd87c6d139dac404d6dd4ff)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `bioblend`; 命令 `python3`; 环境变量 `GALAXY_URL`, `GALAXY_API_KEY`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info galaxy-bioblend
dnakit skills install galaxy-bioblend --agent claude-code   # 或 --agent codex / cursor
```

## `galaxy-bridge` {#galaxy-bridge}

Galaxy tool discovery, intelligent recommendation, and execution — 8,000+ bioinformatics tools from usegalaxy.org with multi-signal scoring and workflow suggestions

- **来源：**[Copyright 2026 Manuel Corpas; Hiranyamaya Dash; Aviv Madar (individual skill notices retained) @ `4848cc2`](https://github.com/ClawBio/ClawBio/tree/4848cc280c446093b1364a2445db3137d455cc17/skills/galaxy-bridge)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `bioblend`; 命令 `python3`; 环境变量 `GALAXY_URL`, `GALAXY_API_KEY`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info galaxy-bridge
dnakit skills install galaxy-bridge --agent claude-code   # 或 --agent codex / cursor
```

## `genomics-study-design` {#genomics-study-design}

Design and audit non-model, polyploid, duplicated, or non-human genomic studies. Use when a task needs reference assemblies, truth sets, mapping, variant discovery, species-aware QC, paralogy or WGD, leakage, population or GWAS validation, trait simulation, or diagnosis of unexpected signals. Compose with evidence and statistics skills. Do not use for one command or generic biology.

- **来源：**[Ynzuses/research-skills-suite contributors @ `a5cb708`](https://github.com/Ynzuses/research-skills-suite/tree/a5cb708252a58bcfb771d46d39f5e38d64b8b2f6/skills/genomics-study-design)
- **许可证：**`Apache-2.0`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info genomics-study-design
dnakit skills install genomics-study-design --agent claude-code   # 或 --agent codex / cursor
```

## `illumina-bridge` {#illumina-bridge}

Import DRAGEN-exported Illumina result bundles into ClawBio for local tertiary analysis and downstream routing.

- **来源：**[ClawBio; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/illumina-bridge)
- **许可证：**`MIT`
- **运行依赖：**命令 `python3`; 环境变量 `ILLUMINA_ICA_API_KEY`, `ILLUMINA_ICA_BASE_URL`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info illumina-bridge
dnakit skills install illumina-bridge --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-debug-analysis` {#isabl-debug-analysis}

Systematically debug a failed Isabl analysis. Use when an analysis has FAILED status or unexpected behavior.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-debug-analysis)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-debug-analysis
dnakit skills install isabl-debug-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-merge-results` {#isabl-merge-results}

Aggregate results from multiple Isabl analyses into combined datasets. Use when merging VCFs, TSVs, or other outputs across samples or cohorts.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-merge-results)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-merge-results
dnakit skills install isabl-merge-results --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-monitor-analyses` {#isabl-monitor-analyses}

Monitor and track the status of Isabl analyses. Use when checking job progress, finding failed analyses, or tracking pipeline execution.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-monitor-analyses)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-monitor-analyses
dnakit skills install isabl-monitor-analyses --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-project-report` {#isabl-project-report}

Generate a status report for an Isabl project. Use when summarizing project progress, identifying issues, or preparing updates.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-project-report)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-project-report
dnakit skills install isabl-project-report --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-query-data` {#isabl-query-data}

Help construct queries to retrieve data from Isabl using the Python SDK. Use when searching for experiments, analyses, samples, or other data.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-query-data)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-query-data
dnakit skills install isabl-query-data --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-run-pipeline` {#isabl-run-pipeline}

Submit and run multiple Isabl applications as a pipeline. Use when chaining apps together or running a sequence of analyses on samples.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-run-pipeline)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-run-pipeline
dnakit skills install isabl-run-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-submit-data` {#isabl-submit-data}

Submit sequencing data to Isabl by creating individuals, samples, and experiments. Use when registering new data files or importing from external sources.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-submit-data)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-submit-data
dnakit skills install isabl-submit-data --agent claude-code   # 或 --agent codex / cursor
```

## `isabl-write-app` {#isabl-write-app}

Guide through creating a new Isabl bioinformatics application. Use when building pipelines that integrate with the Isabl platform.

- **来源：**[juanesarango/isabl-mcp contributors @ `85b3e29`](https://github.com/juanesarango/isabl-mcp/tree/85b3e298d9cfdf02cb82ebfdc780649a6811fd2e/skills/isabl-write-app)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `dnakit`, `isabl_cli`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info isabl-write-app
dnakit skills install isabl-write-app --agent claude-code   # 或 --agent codex / cursor
```

## `kiro-healthomics-genomics` {#kiro-healthomics-genomics}

Prepare DNA sequencing and variant workflows for AWS HealthOmics using selected upstream Kiro guidance. Use for WDL/Nextflow import, genomic file-format checks, container references and workflow/run inspection; requires a separately configured official AWS HealthOmics MCP server or AWS CLI.

- **来源：**[Amazon.com, Inc. or its affiliates @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc/kiro-life-sciences)
- **许可证：**`MIT-0`
- **运行依赖：**命令 `uvx`, `aws`; 环境变量 `AWS_PROFILE`, `AWS_REGION`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info kiro-healthomics-genomics
dnakit skills install kiro-healthomics-genomics --agent claude-code   # 或 --agent codex / cursor
```

## `multi-omics-data-harmonization` {#multi-omics-data-harmonization}

Harmonizes already-normalized per-omic matrices onto a common footing before joint integration - assembling a MultiAssayExperiment, choosing the per-omic variance-stabilizing transform, deciding per-view versus per-feature scaling, picking a cross-omic batch strategy, and triaging missing data (feature, value, or whole sample; MAR versus MNAR). Covers why a shared-latent integrator is blind to what an omic is so scaling silently decides which block dominates, why batch confounded with biology is irrecoverable and should be modeled as a covariate not scrubbed, and why stacking blocks and running one ComBat erases cross-omic signal. Use when preparing two or more omics for MOFA2, mixOmics, or SNF, deciding a transform or scaling, correcting batch across modalities, or handling missing omics per sample. For deep per-omic normalization see differential-expression, methylation-analysis, proteomics, metabolomics; for the method decision see integration-design; for fusion see mofa-integration, mixomics-analysis.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/multi-omics-integration/data-harmonization)
- **许可证：**`MIT`
- **别名：**`bio-multi-omics-data-harmonization`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info multi-omics-data-harmonization
dnakit skills install multi-omics-data-harmonization --agent claude-code   # 或 --agent codex / cursor
```

## `multi-omics-integration` {#multi-omics-integration}

Multi-omics integration — orchestrate per-layer analysis (transcriptomics, proteomics, epigenomics, genomics, metabolomics) then perform cross-omics correlation, multi-omics clustering, and pathway-level integration. Use for integrative systems-biology analysis, multi-modal disease characterization, and cross-omics biomarker discovery.

- **来源：**[ToolUniverse contributors, mims-harvard @ `7a0ceb2`](https://github.com/mims-harvard/ToolUniverse/tree/7a0ceb2f3a2ae28594e55ec27344e839df8bc9f7/skills/tooluniverse-multi-omics-integration)
- **许可证：**`Apache-2.0`
- **别名：**`integrative-analysis`, `multi-omics-pipeline`, `tooluniverse-multi-omics-integration`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info multi-omics-integration
dnakit skills install multi-omics-integration --agent claude-code   # 或 --agent codex / cursor
```

## `multi-omics-mixomics-analysis` {#multi-omics-mixomics-analysis}

Builds supervised and unsupervised multivariate integration across bulk omics blocks with mixOmics - sPLS for sparse pairwise correlation, DIABLO (block.splsda) for a multi-block discriminant signature, rCCA for regularized canonical correlation, and MINT for multi-study integration. Covers why these projection methods maximize covariance or correlation and not truth, why DIABLO's design matrix is the central correlation-versus-discrimination decision, why cross-validation must wrap keepX selection or the reported error is leaked, why balanced error rate is required under class imbalance, and why DIABLO needs matched samples while MINT handles multiple cohorts. Use when finding a cross-omic discriminant signature for a known outcome, selecting correlated features between two omics, tuning keepX, or integrating one omic across studies. For unsupervised factors see mofa-integration; for the method decision see integration-design; for cross-validation theory see machine-learning/model-validation.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/multi-omics-integration/mixomics-analysis)
- **许可证：**`MIT`
- **别名：**`bio-multi-omics-mixomics-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info multi-omics-mixomics-analysis
dnakit skills install multi-omics-mixomics-analysis --agent claude-code   # 或 --agent codex / cursor
```

## `multi-omics-mofa-integration` {#multi-omics-mofa-integration}

Discovers shared and view-specific latent factors across bulk multi-omics blocks (RNA-seq, proteomics, methylation) on a common sample axis with MOFA2's unsupervised Bayesian group factor model, then attributes per-view variance explained and interprets signed factor weights. Covers why a factor is an unsupervised axis of variance and not a pathway, why a factor that correlates with batch is a batch factor, why the per-view variance-explained table is the primary read-out rather than p-values, why raw counts in a Gaussian view make factor 1 the library-size factor, and why MOFA2 handles missing omics-per-sample natively. Use when integrating two or more bulk omics to find joint axes of variation, choosing factor count, labeling factors against metadata, or running enrichment on factor weights. For supervised discriminant integration see mixomics-analysis; for the method decision see integration-design; for single-cell see single-cell/multimodal-integration; for enrichment see pathway-analysis/gsea.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/multi-omics-integration/mofa-integration)
- **许可证：**`MIT`
- **别名：**`bio-multi-omics-mofa-integration`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info multi-omics-mofa-integration
dnakit skills install multi-omics-mofa-integration --agent claude-code   # 或 --agent codex / cursor
```

## `nextflow-development` {#nextflow-development}

Run nf-core bioinformatics pipelines (rnaseq, sarek, atacseq) on sequencing data. Use when analyzing RNA-seq, WGS/WES, or ATAC-seq data—either local FASTQs or public datasets from GEO/SRA. Triggers on nf-core, Nextflow, FASTQ analysis, variant calling, gene expression, differential expression, GEO reanalysis, GSE/GSM/SRR accessions, or samplesheet creation.

- **来源：**[Anthropic and life-sciences skill contributors; leaf LICENSE.txt retained @ `e96556b`](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c/nextflow-development)
- **许可证：**`Apache-2.0`
- **运行依赖：**Python 模块 `yaml`; 命令 `python3`, `nextflow`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info nextflow-development
dnakit skills install nextflow-development --agent claude-code   # 或 --agent codex / cursor
```

## `nfcore-sarek-wrapper` {#nfcore-sarek-wrapper}

ClawBio wrapper around nf-core/sarek 3.8.1 covering mapping through annotation for germline, tumor-only, and somatic paired analyses.

- **来源：**[ClawBio; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/nfcore-sarek-wrapper)
- **许可证：**`MIT`
- **运行依赖：**命令 `java`, `nextflow`, `python3`
- **别名：**`alterlab-nf-core-sarek`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info nfcore-sarek-wrapper
dnakit skills install nfcore-sarek-wrapper --agent claude-code   # 或 --agent codex / cursor
```

## `pipeline-guide` {#pipeline-guide}

Access ENCODE uniform analysis pipelines, generate user-specific Nextflow/WDL pipelines, manage compute resources, and integrate with cloud platforms. Use when the user wants to understand ENCODE pipelines, run pipelines on their own data, generate custom Nextflow workflows from ENCODE pipeline code, check compute requirements (CPU/GPU/memory), run pipelines in background, or integrate with Google Cloud, AWS, or other cloud platforms. Also use when the user asks about ENCODE pipeline outputs, processing standards, software versions, or wants to replicate ENCODE processing. Covers local execution, HPC, and cloud deployment with resource-aware scheduling. Use this skill for ANY pipeline execution, workflow generation, or compute resource management task involving ENCODE data.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/pipeline-guide)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info pipeline-guide
dnakit skills install pipeline-guide --agent claude-code   # 或 --agent codex / cursor
```

## `plant-omics-db` {#plant-omics-db}

Build a journal-ready plant multi-omics database (genome, germplasm, variome, transcriptome, GWAS, tools). Use when the user wants a NAR/MP/PC class plant genomics web database, module list, schema, or Flask/FastAPI prototype.

- **来源：**[Xuzhen-Li/bioinfo-agent-skills contributors @ `6dea4be`](https://github.com/Xuzhen-Li/bioinfo-agent-skills/tree/6dea4bee7aade611d22f2f1958612328bc61ec83/skills/plant-omics-db)
- **许可证：**`MIT`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info plant-omics-db
dnakit skills install plant-omics-db --agent claude-code   # 或 --agent codex / cursor
```

## `prefect-dask-nextflow` {#prefect-dask-nextflow}

Design reproducible bioinformatics pipelines with Prefect plus Dask or Nextflow. Use when scaffolding local, distributed, or scheduler-backed workflows.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-prefect-dask-nextflow)
- **许可证：**`MIT`
- **别名：**`bio-prefect-dask-nextflow`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info prefect-dask-nextflow
dnakit skills install prefect-dask-nextflow --agent claude-code   # 或 --agent codex / cursor
```

## `seq-wrangler` {#seq-wrangler}

NGS read QC, alignment, and BAM processing pipeline. Wraps FastQC, BWA/Bowtie2/Minimap2, SAMtools, and MultiQC for automated read-to-BAM workflows.

- **来源：**[Daniel Garbozo; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/seq-wrangler)
- **许可证：**`MIT`
- **别名：**`seqkit`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info seq-wrangler
dnakit skills install seq-wrangler --agent claude-code   # 或 --agent codex / cursor
```

## `single-cell-clustering` {#single-cell-clustering}

Dimensionality reduction and graph-based clustering for single-cell RNA-seq with Scanpy (Python) and Seurat (R). Resolves which algorithm to use (Leiden vs Louvain), how many PCs and neighbors to set, how to sweep and validate resolution, when a split is over-clustering, and why post-clustering marker p-values are not valid inference. Use when clustering cells, choosing a clustering resolution, deciding whether two clusters are one population, building a UMAP/tSNE, or judging whether clusters are real.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/single-cell/clustering)
- **许可证：**`MIT`
- **别名：**`bio-single-cell-clustering`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info single-cell-clustering
dnakit skills install single-cell-clustering --agent claude-code   # 或 --agent codex / cursor
```

## `single-cell-data-io` {#single-cell-data-io}

Read, write, create, and convert single-cell objects across AnnData (Python), Seurat (R), and SingleCellExperiment (R). Use when loading 10X Cell Ranger output (raw vs filtered), importing or exporting h5ad/RDS/h5mu/zarr, building AnnData or Seurat objects from matrices, moving objects between Python and R, or debugging lost layers, transposed matrices, or mangled gene names during conversion.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/single-cell/data-io)
- **许可证：**`MIT`
- **别名：**`bio-single-cell-data-io`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info single-cell-data-io
dnakit skills install single-cell-data-io --agent claude-code   # 或 --agent codex / cursor
```

## `single-cell-multimodal-integration` {#single-cell-multimodal-integration}

Integrate multimodal single-cell data (CITE-seq RNA+protein, 10x Multiome RNA+ATAC, unpaired/diagonal RNA+ATAC) and choose the right joint method. Use when classifying an integration task by anchor structure (paired vs unpaired), denoising CITE-seq ADT background before joint embedding, picking between WNN, totalVI, MultiVI, MOFA+, GLUE, or Seurat v5 bridge integration, or diagnosing why a modality dominates a joint clustering.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/single-cell/multimodal-integration)
- **许可证：**`MIT`
- **别名：**`bio-single-cell-multimodal-integration`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info single-cell-multimodal-integration
dnakit skills install single-cell-multimodal-integration --agent claude-code   # 或 --agent codex / cursor
```

## `single-cell-preprocessing` {#single-cell-preprocessing}

Quality control, ambient-RNA handling, normalization, and feature selection for single-cell RNA-seq using Scanpy (Python) and Seurat (R). Use when filtering low-quality cells with MAD-adaptive thresholds, setting tissue-aware mito cutoffs, removing ambient RNA (SoupX/CellBender/DecontX), choosing a normalization (shifted-log vs scran vs sctransform vs Pearson residuals), selecting highly variable genes, or deciding whether to scale and regress out covariates.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/single-cell/preprocessing)
- **许可证：**`MIT`
- **别名：**`bio-single-cell-preprocessing`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info single-cell-preprocessing
dnakit skills install single-cell-preprocessing --agent claude-code   # 或 --agent codex / cursor
```

## `track-experiments` {#track-experiments}

Track ENCODE experiments locally with publications, citations, and provenance. Use when the user wants to build a collection of experiments, manage citations, compare experiments, or track data provenance.

- **来源：**[ammawla/encode-toolkit contributors @ `b1d90bc`](https://github.com/ammawla/encode-toolkit/tree/b1d90bc9e7ead40d49b8c89f0b9b9e7b8e592a5b/skills/track-experiments)
- **许可证：**`AGPL-3.0-only`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info track-experiments
dnakit skills install track-experiments --agent claude-code   # 或 --agent codex / cursor
```

## `wgs-prs` {#wgs-prs}

End-to-end WGS to polygenic risk score pipeline. Takes paired-end FASTQ files (or a pre-existing VCF) through nf-core/sarek for variant calling, applies VCF QC (normalisation, hard filtering, Ti/Tv and Het/Hom checks), then computes polygenic risk scores via the PGS Catalog. Fills the FASTQ to VCF gap upstream of the gwas-prs skill.

- **来源：**[David de Lorenzo; ClawBio contributors @ `1fcecb7`](https://github.com/ClawBio/ClawBio/tree/1fcecb7e9531210baabe6d5cdbdd08de1ff3ab8d/skills/wgs-prs)
- **许可证：**`MIT`
- **运行依赖：**命令 `nextflow`, `python3`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info wgs-prs
dnakit skills install wgs-prs --agent claude-code   # 或 --agent codex / cursor
```

## `workflow-management-cwl-workflows` {#workflow-management-cwl-workflows}

Authors portable, strongly-typed bioinformatics pipelines in the Common Workflow Language (CWL v1.2) as CommandLineTool/Workflow/ExpressionTool documents, validated with cwltool and run at scale on Toil/Arvados/Calrissian. Use when deciding CWL (portability/provenance/regulated) vs Nextflow/WDL/Snakemake; declaring secondaryFiles for indexed companions (.bai/.fai/.dict/.tbi and the caret rule); putting resources/containers under requirements (must-hold) vs hints (advisory) to avoid silent OOM; choosing scatterMethod (dotproduct vs flat_/nested_crossproduct); preferring $(...) parameter refs over $\{...\} JavaScript for portability; pinning DockerRequirement images; or emitting a CWLProv provenance object for audited/clinical settings.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflow-management/cwl-workflows)
- **许可证：**`MIT`
- **运行依赖：**命令 `cwltool`
- **别名：**`bio-workflow-management-cwl-workflows`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflow-management-cwl-workflows
dnakit skills install workflow-management-cwl-workflows --agent claude-code   # 或 --agent codex / cursor
```

## `workflow-management-nextflow-pipelines` {#workflow-management-nextflow-pipelines}

Authors reproducible Nextflow DSL2 pipelines built on reactive dataflow, where processes communicate only through channels and execution order is not guaranteed. Use when deciding channel/dataflow (Nextflow) vs rule-based (Snakemake) authoring; wiring queue vs value channels and fixing shared-reference exhaustion with .first(); composing DSL2 modules and subworkflows with take/main/emit; selecting container/conda profiles and pinning images by digest for portability across local/SLURM/LSF/AWS Batch/Google Batch/Kubernetes executors; diagnosing why -resume misses the cache (nondeterministic input order, mtime on network filesystems, mutable :latest tags) with cache 'lenient' and -dump-hashes; managing work/ vs publishDir and dynamic retry escalation; and choosing whether to adopt an nf-core community pipeline or author from scratch.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflow-management/nextflow-pipelines)
- **许可证：**`MIT`
- **运行依赖：**命令 `Nextflow`
- **别名：**`bio-workflow-management-nextflow-pipelines`, `nextflow`, `nextflow-dsl2-generator`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflow-management-nextflow-pipelines
dnakit skills install workflow-management-nextflow-pipelines --agent claude-code   # 或 --agent codex / cursor
```

## `workflow-management-nf-core-pipelines` {#workflow-management-nf-core-pipelines}

Runs and configures curated nf-core community Nextflow pipelines (rnaseq, sarek, atacseq, methylseq, ampliseq, taxprofiler, fetchngs) reproducibly, pinning the pipeline revision with -r and selecting a container engine and institutional config via -profile. Use when deciding to adopt a community pipeline versus author one from scratch; picking a pipeline and pinning its -r revision; selecting -profile test/docker/singularity/conda plus an institutional config from nf-core/configs; building and validating a samplesheet CSV against the pipeline schema (nf-schema); choosing --genome/iGenomes versus custom references; configuring resources and max_memory for SLURM/AWS Batch; using -resume and -stub; and reading MultiQC outputs.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflow-management/nf-core-pipelines)
- **许可证：**`MIT`
- **运行依赖：**命令 `nf-core`
- **别名：**`bio-workflow-management-nf-core-pipelines`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflow-management-nf-core-pipelines
dnakit skills install workflow-management-nf-core-pipelines --agent claude-code   # 或 --agent codex / cursor
```

## `workflow-management-snakemake-workflows` {#workflow-management-snakemake-workflows}

Authors reproducible bioinformatics pipelines with Snakemake - rules wired by output-file pattern, wildcards and expand() for sample fan-out, checkpoints for runtime-unknown outputs, resource/retry escalation, and conda/container software deployment on HPC and cloud. Use when deciding rule-based (Snakemake) vs channel/dataflow (Nextflow) authoring; wiring rules by OUTPUT-file pattern rather than imperative order; using wildcards + expand() for sample fan-out and constraining them to stop silent mis-routing; adding checkpoints when the set of outputs is unknown until a step runs (dynamic DAG); diagnosing why a job reran (or did not) under the mtime-plus-provenance trigger set; escalating memory on retry for OOM-killed jobs; and porting a Snakemake 7 `--cluster`/remote-provider command to the Snakemake 8+ executor-plugin and storage-plugin model (snakemake-executor-plugin-slurm) with `--software-deployment-method`.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflow-management/snakemake-workflows)
- **许可证：**`MIT`
- **运行依赖：**Python 模块 `Snakemake`, `dnakit`
- **别名：**`bio-workflow-management-snakemake-workflows`, `snakemake`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflow-management-snakemake-workflows
dnakit skills install workflow-management-snakemake-workflows --agent claude-code   # 或 --agent codex / cursor
```

## `workflow-management-wdl-workflows` {#workflow-management-wdl-workflows}

Authors bioinformatics pipelines in WDL (Workflow Description Language) run by Cromwell or miniwdl, targeting the GATK/Broad and Terra/AnVIL/BioData Catalyst cloud ecosystem, with tasks, workflows, scatter-gather parallelism, structs, and a runtime block that sizes the cloud VM. Use when deciding to target Terra/AnVIL/GATK/WARP (chosen for the ecosystem, not the language); sizing runtime disks dynamically for a fresh-per-task cloud VM (ceil(size(f)*factor)+buffer); choosing preemptible vs on-demand VMs by task length and idempotency; picking Cromwell (production, cloud, call-caching) vs miniwdl (local dev, miniwdl check linting, readable errors); enabling and debugging call-caching silent-miss modes; pinning Docker by digest for reproducibility and cache stability; or scattering an array for parallel fan-out.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflow-management/wdl-workflows)
- **许可证：**`MIT`
- **运行依赖：**命令 `cromwell`
- **别名：**`bio-workflow-management-wdl-workflows`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflow-management-wdl-workflows
dnakit skills install workflow-management-wdl-workflows --agent claude-code   # 或 --agent codex / cursor
```

## `workflow-methods-docwriter` {#workflow-methods-docwriter}

Generate reproducible Methods from Nextflow, Snakemake, or CWL run artifacts. Use when documenting exact commands, versions, parameters, QC gates, provenance, and outputs.

- **来源：**[Copyright (c) 2026 Frederik Schulz @ `7a4eb5b`](https://github.com/fmschulz/omics-skills/tree/7a4eb5bace3f874c3be2ec697e9d2e695c374637/skills/bio-workflow-methods-docwriter)
- **许可证：**`MIT`
- **别名：**`bio-workflow-methods-docwriter`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflow-methods-docwriter
dnakit skills install workflow-methods-docwriter --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-atacseq-pipeline` {#workflows-atacseq-pipeline}

Orchestrates the end-to-end bulk ATAC-seq pipeline from FASTQ to differential accessibility and TF footprints, chaining Nextera-aware fastp QC, Bowtie2 alignment, chrM removal, dedup, a single Tn5 +4/-5 shift, MACS3 peak calling, Corces fixed-width consensus, DiffBind/csaw differential accessibility, and TOBIAS footprinting. Use when committing the reference build + blacklist once, recognizing ATAC has NO input control (the shift-extend model IS the background), applying the Tn5 shift exactly once (never combining -f BAMPE with --shift/--extsize), removing chrM before calling, building a fixed-width consensus so per-sample counts are comparable, or choosing MACS3 vs Genrich vs HMMRATAC. Hands mechanism to the atac-seq component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/atacseq-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'read-qc/fastp-workflow', 'id': 'read-qc-fastp-workflow', 'status': 'bundled'}`, `{'upstream_path': 'read-alignment/bowtie2-alignment', 'id': 'read-alignment-bowtie2-alignment', 'status': 'bundled'}`, `{'upstream_path': 'alignment-files/duplicate-handling', 'id': 'duplicate-handling', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/atac-peak-calling', 'id': 'atac-seq-atac-peak-calling', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/atac-qc', 'id': 'atac-seq-atac-qc', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/consensus-peakset', 'id': 'atac-seq-consensus-peakset', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/differential-accessibility', 'id': 'atac-seq-differential-accessibility', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/footprinting', 'id': 'atac-seq-footprinting', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/motif-deviation', 'id': 'atac-seq-motif-deviation', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/nucleosome-positioning', 'id': 'atac-seq-nucleosome-positioning', 'status': 'bundled'}`
- **别名：**`atac-seq`, `bio-workflows-atacseq-pipeline`, `pipeline-atacseq`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-atacseq-pipeline
dnakit skills install workflows-atacseq-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-causal-genomics-pipeline` {#workflows-causal-genomics-pipeline}

End-to-end post-GWAS causal inference pipeline orchestrating heritability partitioning, genetic correlation, Mendelian randomization with CHP-aware sensitivity (CAUSE / LHC-MR), colocalization, fine-mapping with SuSiE / FOCUS, mediation, TWAS triangulation, cis-pQTL drug-target MR, effector-gene prioritization (L2G / PoPS / cS2G), and GenomicSEM common-factor GWAS. Use when triangulating causal inference across multiple complementary methods, prioritizing tissues via stratified LDSC, nominating or de-risking drug targets, mapping a lead SNP to a candidate effector gene, modeling shared genetic architecture across correlated traits, or producing a STROBE-MR-compliant publication-grade evidence battery from GWAS summary statistics.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/causal-genomics-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'causal-genomics/mendelian-randomization', 'id': 'causal-genomics-mendelian-randomization', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/colocalization-analysis', 'id': 'causal-genomics-colocalization-analysis', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/fine-mapping', 'id': 'causal-genomics-fine-mapping', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/pleiotropy-detection', 'id': 'causal-genomics-pleiotropy-detection', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/mediation-analysis', 'id': 'causal-genomics-mediation-analysis', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/transcriptome-wide-association', 'id': 'causal-genomics-transcriptome-wide-association', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/heritability-partitioning', 'id': 'causal-genomics-heritability-partitioning', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/proteome-mr-drug-target', 'id': 'causal-genomics-proteome-mr-drug-target', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/effector-gene-prioritization', 'id': 'causal-genomics-effector-gene-prioritization', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/genetic-correlation', 'id': 'causal-genomics-genetic-correlation', 'status': 'bundled'}`, `{'upstream_path': 'causal-genomics/genomic-sem', 'id': 'causal-genomics-genomic-sem', 'status': 'bundled'}`
- **别名：**`bio-workflows-causal-genomics-pipeline`, `causal-genomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-causal-genomics-pipeline
dnakit skills install workflows-causal-genomics-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-chipseq-pipeline` {#workflows-chipseq-pipeline}

Orchestrates the end-to-end ChIP-seq pipeline from FASTQ to blacklist-filtered, annotated peaks, chaining fastp QC, Bowtie2 alignment, pre-dedup library-complexity QC (NRF/PBC), duplicate removal, chrM + ENCODE-blacklist filtering, MACS3 peak calling against a matched input, IDR/consensus reproducibility, deepTools signal tracks, and ChIPseeker annotation. Use when committing the reference build + blacklist version + effective genome size once, pairing each IP with its matched control, computing complexity metrics BEFORE dedup, choosing narrow vs broad and MACS3 vs SEACR/Genrich, keeping per-replicate peaks for IDR, or avoiding depth-normalization that erases a spike-in global shift. Hands mechanism to the chip-seq component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/chipseq-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'read-qc/fastp-workflow', 'id': 'read-qc-fastp-workflow', 'status': 'bundled'}`, `{'upstream_path': 'read-alignment/bowtie2-alignment', 'id': 'read-alignment-bowtie2-alignment', 'status': 'bundled'}`, `{'upstream_path': 'alignment-files/duplicate-handling', 'id': 'duplicate-handling', 'status': 'bundled'}`, `{'upstream_path': 'chip-seq/chipseq-qc', 'id': 'chipseq-qc', 'status': 'bundled'}`, `{'upstream_path': 'chip-seq/peak-calling', 'id': 'chipseq-peak-calling', 'status': 'bundled'}`, `{'upstream_path': 'chip-seq/peak-annotation', 'id': 'chipseq-peak-annotation', 'status': 'bundled'}`, `{'upstream_path': 'chip-seq/differential-binding', 'id': 'chipseq-differential-binding', 'status': 'bundled'}`, `{'upstream_path': 'chip-seq/chipseq-visualization', 'id': 'chipseq-visualization', 'status': 'bundled'}`, `{'upstream_path': 'chip-seq/motif-analysis', 'id': 'chipseq-motif-analysis', 'status': 'bundled'}`
- **别名：**`bio-workflows-chipseq-pipeline`, `chip-seq`, `pipeline-chipseq`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-chipseq-pipeline
dnakit skills install workflows-chipseq-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-cnv-pipeline` {#workflows-cnv-pipeline}

Orchestrates the copy-number pipeline from BAM to segmented, integer-called, annotated CNVs, forking on germline-vs-somatic - CNVkit (somatic exome/panel: coverage -&gt; assay-matched reference/PoN -&gt; fix -&gt; segment -&gt; purity/ploidy-aware call), GATK gCNV (germline rare-CNV cohort), and allele-specific callers (ASCAT/FACETS/PURPLE) for purity/ploidy. Use when committing the build + target/access BED + PoN once (assay-matched), building the reference from normals BEFORE segmenting, fitting purity/ploidy BEFORE integer calls in tumors, centering on the true (non-diploid) mode before GISTIC2 recurrence, or routing cfDNA to ichorCNA. Hands mechanism to the copy-number component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/cnv-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'copy-number/cnvkit-analysis', 'id': 'copy-number-cnvkit-analysis', 'status': 'bundled'}`, `{'upstream_path': 'copy-number/gatk-cnv', 'id': 'copy-number-gatk-cnv', 'status': 'bundled'}`, `{'upstream_path': 'copy-number/copy-ratio-segmentation', 'id': 'copy-number-copy-ratio-segmentation', 'status': 'bundled'}`, `{'upstream_path': 'copy-number/allele-specific-copy-number', 'id': 'copy-number-allele-specific-copy-number', 'status': 'bundled'}`, `{'upstream_path': 'copy-number/cnv-visualization', 'id': 'copy-number-cnv-visualization', 'status': 'bundled'}`, `{'upstream_path': 'copy-number/cnv-annotation', 'id': 'copy-number-cnv-annotation', 'status': 'bundled'}`, `{'upstream_path': 'copy-number/recurrent-cnv', 'id': 'copy-number-recurrent-cnv', 'status': 'bundled'}`
- **别名：**`bio-workflows-cnv-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-cnv-pipeline
dnakit skills install workflows-cnv-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-crispr-editing-pipeline` {#workflows-crispr-editing-pipeline}

Orchestrates an end-to-end CRISPR editing experiment design from target gene to delivery-ready, validatable constructs. Sequences guide design, off-target assessment, edit-modality selection (knockout, base editing, prime editing, HDR knock-in), and template/donor design, with a QC checkpoint at each handoff. Use when designing a complete CRISPR experiment for knockout, point correction, or tagging and the order of operations, the modality decision, and the cross-cutting traps are needed rather than a single step. Defers each step's mechanics to the genome-engineering skills.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/crispr-editing-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'genome-engineering/grna-design', 'id': 'genome-engineering-grna-design', 'status': 'bundled'}`, `{'upstream_path': 'genome-engineering/off-target-prediction', 'id': 'genome-engineering-off-target-prediction', 'status': 'bundled'}`, `{'upstream_path': 'genome-engineering/base-editing-design', 'id': 'genome-engineering-base-editing-design', 'status': 'bundled'}`, `{'upstream_path': 'genome-engineering/prime-editing-design', 'id': 'genome-engineering-prime-editing-design', 'status': 'bundled'}`, `{'upstream_path': 'genome-engineering/hdr-template-design', 'id': 'genome-engineering-hdr-template-design', 'status': 'bundled'}`
- **别名：**`bio-workflows-crispr-editing-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-crispr-editing-pipeline
dnakit skills install workflows-crispr-editing-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-crispr-screen-pipeline` {#workflows-crispr-screen-pipeline}

End-to-end pooled and single-cell CRISPR screen analysis from FASTQ to hit genes. Orchestrates library design QC, guide counting, six-stage screen QC (plasmid Gini, replicate Pearson, CEGv2 PR-AUC, copy-number artifact), method-appropriate hit calling across MAGeCK RRA/MLE, BAGEL2, drugZ, JACKS, and Chronos, cancer-cell-line copy-number correction (CRISPRcleanR / Chronos), batch correction for multi-batch screens, and the specialized branches for combinatorial paralog screens, single-cell Perturb-seq, base-editor variant-function screens, prime-editor screens, and in vivo bottleneck-aware screens. Use when analyzing any pooled CRISPR screen end-to-end, matching the hit-calling method to the experimental design, integrating copy-number correction into the pipeline, or branching the workflow for single-cell, combinatorial, base-editor, prime-editor, or in vivo variants.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/crispr-screen-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'crispr-screens/library-design', 'id': 'crispr-screens-library-design', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/screen-qc', 'id': 'crispr-screens-screen-qc', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/mageck-analysis', 'id': 'crispr-screens-mageck-analysis', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/bagel-essentiality', 'id': 'crispr-screens-bagel-essentiality', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/drugz-chemogenomic', 'id': 'crispr-screens-drugz-chemogenomic', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/jacks-analysis', 'id': 'crispr-screens-jacks-analysis', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/hit-calling', 'id': 'crispr-screens-hit-calling', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/copy-number-correction', 'id': 'crispr-screens-copy-number-correction', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/batch-correction', 'id': 'crispr-screens-batch-correction', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/crispresso-editing', 'id': 'crispr-screens-crispresso-editing', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/base-editing-analysis', 'id': 'crispr-screens-base-editing-analysis', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/prime-editing-screens', 'id': 'crispr-screens-prime-editing-screens', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/perturb-seq-analysis', 'id': 'crispr-screens-perturb-seq-analysis', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/combinatorial-screens', 'id': 'crispr-screens-combinatorial-screens', 'status': 'bundled'}`, `{'upstream_path': 'crispr-screens/in-vivo-screens', 'id': 'crispr-screens-in-vivo-screens', 'status': 'bundled'}`
- **别名：**`bio-workflows-crispr-screen-pipeline`, `pooled-crispr-screens`, `tooluniverse-crispr-screen-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-crispr-screen-pipeline
dnakit skills install workflows-crispr-screen-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-edna-pipeline` {#workflows-edna-pipeline}

End-to-end eDNA metabarcoding from raw amplicons to community ecology. Covers QC, primer removal (mandatory before DADA2 filterAndTrim), denoising with OBITools3 v3 (obi stats plural; DMS-based) or DADA2 ASVs (Callahan 2017), decontam combined method as screening-not-classifier (Davis 2018), tag-jumping (Schnell 2015) with a platform-dependent baseline (NovaSeq patterned flow cells ~10x MiSeq), Hill-number effective species counts with coverage-based rarefaction (Jost 2006; Chao &amp; Jost 2012; doubling rule), beta-diversity decomposition with MANDATORY PERMANOVA + PERMDISP pair (Anderson &amp; Walsh 2013), constrained ordination, and the read-counts-not-abundance critique (Lamb 2019). Use when processing eDNA samples for biodiversity assessment, deciding ASV vs OTU, configuring OBITools3 v3, interpreting decontam screening, or reporting community comparisons with the dispersion confound check.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/edna-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'ecological-genomics/edna-metabarcoding', 'id': 'ecological-genomics-edna-metabarcoding', 'status': 'bundled'}`, `{'upstream_path': 'ecological-genomics/biodiversity-metrics', 'id': 'ecological-genomics-biodiversity-metrics', 'status': 'bundled'}`, `{'upstream_path': 'ecological-genomics/community-ecology', 'id': 'ecological-genomics-community-ecology', 'status': 'bundled'}`, `{'upstream_path': 'read-qc/quality-reports', 'id': 'read-qc-quality-reports', 'status': 'bundled'}`
- **别名：**`bio-workflows-edna-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-edna-pipeline
dnakit skills install workflows-edna-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-fastq-to-variants` {#workflows-fastq-to-variants}

Orchestrates the end-to-end germline short-variant pipeline from FASTQ to a filtered, normalized, benchmarked VCF, chaining QC/trim, BWA-MEM2 alignment, duplicate marking, optional BQSR, calling (bcftools/GATK HaplotypeCaller/DeepVariant/DRAGEN), normalization, site+genotype filtering, annotation, and hap.py/vcfeval benchmarking. Use when deciding the pipeline-wide reference-genome commitment (GRCh38 analysis set vs T2T, ALT/decoy handling), sequencing the steps in the defensible order (normalize BEFORE annotate, filter site- then genotype-level), choosing the calling engine and single-sample vs cohort joint-calling, picking a filtering strategy by cohort size, or benchmarking stratified within GIAB confident regions. Hands off mechanism to the variant-calling and read-alignment component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/fastq-to-variants)
- **许可证：**`MIT`
- **运行依赖：**命令 `bcftools`
- **依赖的 Skills：**`{'upstream_path': 'read-qc/fastp-workflow', 'id': 'read-qc-fastp-workflow', 'status': 'bundled'}`, `{'upstream_path': 'read-alignment/bwa-alignment', 'id': 'read-alignment-bwa-alignment', 'status': 'bundled'}`, `{'upstream_path': 'alignment-files/alignment-sorting', 'id': 'alignment-sorting', 'status': 'bundled'}`, `{'upstream_path': 'alignment-files/duplicate-handling', 'id': 'duplicate-handling', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/variant-calling', 'id': 'variant-calling', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/joint-calling', 'id': 'variant-calling-joint-calling', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/variant-normalization', 'id': 'variant-normalization', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/filtering-best-practices', 'id': 'variant-calling-filtering-best-practices', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/variant-annotation', 'id': 'variant-annotation', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/vcf-statistics', 'id': 'vcf-statistics', 'status': 'bundled'}`
- **别名：**`bio-workflows-fastq-to-variants`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-fastq-to-variants
dnakit skills install workflows-fastq-to-variants --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-genome-annotation-pipeline` {#workflows-genome-annotation-pipeline}

Orchestrates genome annotation from assembled contigs to functional annotation, forking prokaryotic (Bakta one-step, genetic-code table from GTDB-Tk) vs eukaryotic (RepeatMask -&gt; BRAKER3 -&gt; functional -&gt; ncRNA), then eggNOG/InterProScan functional assignment and Infernal/tRNAscan ncRNA. Use when committing the pro-vs-eukaryotic path and the genetic-code table from taxonomy (never guessing), annotating ONLY a decontaminated QC-passed assembly (CheckM2 before prokaryotic annotation is non-negotiable), committing the evidence set (RNA-seq + protein drives BRAKER3 training), soft-masking with a curated repeat library before gene prediction, or pinning the tool + DB version for any pangenome comparison. Hands mechanism to the genome-annotation component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/genome-annotation-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'genome-annotation/prokaryotic-annotation', 'id': 'genome-annotation-prokaryotic-annotation', 'status': 'bundled'}`, `{'upstream_path': 'genome-annotation/eukaryotic-gene-prediction', 'id': 'genome-annotation-eukaryotic-gene-prediction', 'status': 'bundled'}`, `{'upstream_path': 'genome-annotation/repeat-annotation', 'id': 'genome-annotation-repeat-annotation', 'status': 'bundled'}`, `{'upstream_path': 'genome-annotation/functional-annotation', 'id': 'genome-annotation-functional-annotation', 'status': 'bundled'}`, `{'upstream_path': 'genome-annotation/ncrna-annotation', 'id': 'genome-annotation-ncrna-annotation', 'status': 'bundled'}`, `{'upstream_path': 'genome-annotation/annotation-qc', 'id': 'genome-annotation-annotation-qc', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/assembly-qc', 'id': 'genome-assembly-assembly-qc', 'status': 'bundled'}`
- **别名：**`bio-workflows-genome-annotation-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-genome-annotation-pipeline
dnakit skills install workflows-genome-annotation-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-genome-assembly-pipeline` {#workflows-genome-assembly-pipeline}

Orchestrates an end-to-end de novo genome assembly project, routing each step to the right genome-assembly skill rather than restating it. Profiles the genome first (k-mer spectrum -&gt; size, heterozygosity, ploidy), QCs reads, chooses an assembly path by data type (SPAdes for Illumina, Flye for noisy long reads, hifiasm for HiFi, metaFlye for communities), polishes only when needed, decontaminates, scaffolds with Hi-C, and finishes with three-axis QC (contiguity + completeness + correctness). Use when assembling a genome from raw reads and deciding which assembler, whether to polish, and how to prove the result is good.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/genome-assembly-pipeline)
- **许可证：**`MIT`
- **运行依赖：**命令 `Flye`
- **依赖的 Skills：**`{'upstream_path': 'genome-assembly/genome-profiling', 'id': 'genome-assembly-genome-profiling', 'status': 'bundled'}`, `{'upstream_path': 'read-qc/fastp-workflow', 'id': 'read-qc-fastp-workflow', 'status': 'bundled'}`, `{'upstream_path': 'long-read-sequencing/long-read-qc', 'id': 'long-read-sequencing-long-read-qc', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/short-read-assembly', 'id': 'genome-assembly-short-read-assembly', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/long-read-assembly', 'id': 'genome-assembly-long-read-assembly', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/hifi-assembly', 'id': 'genome-assembly-hifi-assembly', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/metagenome-assembly', 'id': 'genome-assembly-metagenome-assembly', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/assembly-polishing', 'id': 'genome-assembly-assembly-polishing', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/contamination-detection', 'id': 'genome-assembly-contamination-detection', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/scaffolding', 'id': 'genome-assembly-scaffolding', 'status': 'bundled'}`, `{'upstream_path': 'genome-assembly/assembly-qc', 'id': 'genome-assembly-assembly-qc', 'status': 'bundled'}`
- **别名：**`bio-assembly-qc`, `bio-workflows-genome-assembly-pipeline`, `genome-assembly`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-genome-assembly-pipeline
dnakit skills install workflows-genome-assembly-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-gwas-pipeline` {#workflows-gwas-pipeline}

Orchestrates the GWAS pipeline from genotypes to association results, chaining PLINK2 QC (variant-then-sample missingness, controls-only HWE, KING relatedness), panel harmonization + joint phasing/imputation to dosages, long-range-LD-excluded PCA, and an engine chosen by sample structure (PLINK2-GLM / regenie / SAIGE / BOLT-LMM), with LDSC-intercept diagnostics. Use when committing the genome build + ancestry-matched imputation panel once (ancestry match &gt; panel size), running the strand/allele harmonization gate (drop intermediate-frequency palindromes), imputing cases+controls TOGETHER on dosages, excluding long-range-LD regions before PCA, choosing an LMM when relatedness/structure is present (PCs cannot remove a covariance), or separating polygenicity from confounding via the LDSC intercept. Hands mechanism to the population-genetics and phasing-imputation component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/gwas-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'population-genetics/plink-basics', 'id': 'population-genetics-plink-basics', 'status': 'bundled'}`, `{'upstream_path': 'phasing-imputation/reference-panels', 'id': 'phasing-imputation-reference-panels', 'status': 'bundled'}`, `{'upstream_path': 'phasing-imputation/haplotype-phasing', 'id': 'phasing-imputation-haplotype-phasing', 'status': 'bundled'}`, `{'upstream_path': 'phasing-imputation/genotype-imputation', 'id': 'phasing-imputation-genotype-imputation', 'status': 'bundled'}`, `{'upstream_path': 'phasing-imputation/imputation-qc', 'id': 'phasing-imputation-imputation-qc', 'status': 'bundled'}`, `{'upstream_path': 'population-genetics/population-structure', 'id': 'population-genetics-population-structure', 'status': 'bundled'}`, `{'upstream_path': 'population-genetics/association-testing', 'id': 'population-genetics-association-testing', 'status': 'bundled'}`, `{'upstream_path': 'population-genetics/rare-variant-association', 'id': 'population-genetics-rare-variant-association', 'status': 'bundled'}`, `{'upstream_path': 'population-genetics/linkage-disequilibrium', 'id': 'population-genetics-linkage-disequilibrium', 'status': 'bundled'}`
- **别名：**`bio-workflows-gwas-pipeline`, `gwas-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-gwas-pipeline
dnakit skills install workflows-gwas-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-hic-pipeline` {#workflows-hic-pipeline}

End-to-end Hi-C analysis workflow from FASTQ to compartments, TADs, and loops, with the decision of WHICH features the sequencing depth can support. Covers pairtools read-pair processing and library QC, cooler matrices, ICE balancing and distance-decay expected, A/B compartments, TAD boundaries, loop calling, and the routing of HiChIP/PLAC-seq/Capture Hi-C to protein-directed loop callers. Use when processing Hi-C data end to end, deciding a resolution for a given depth, or choosing between bulk-Hi-C and protein-directed loop calling.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/hic-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'hi-c-analysis/contact-pairs', 'id': 'hi-c-analysis-contact-pairs', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/hic-data-io', 'id': 'hi-c-analysis-hic-data-io', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/matrix-operations', 'id': 'hi-c-analysis-matrix-operations', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/compartment-analysis', 'id': 'hi-c-analysis-compartment-analysis', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/tad-detection', 'id': 'hi-c-analysis-tad-detection', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/loop-calling', 'id': 'hi-c-analysis-loop-calling', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/hic-visualization', 'id': 'hi-c-analysis-hic-visualization', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/hic-differential', 'id': 'hi-c-analysis-hic-differential', 'status': 'bundled'}`, `{'upstream_path': 'hi-c-analysis/hichip-plac-loops', 'id': 'hi-c-analysis-hichip-plac-loops', 'status': 'bundled'}`
- **别名：**`bio-workflows-hic-pipeline`, `hi-c-3d-genomics`, `pipeline-hic`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-hic-pipeline
dnakit skills install workflows-hic-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-liquid-biopsy-pipeline` {#workflows-liquid-biopsy-pipeline}

Orchestrates the cell-free DNA / liquid-biopsy pipeline from plasma sequencing to tumor monitoring, forking tumor-naive (screening) vs tumor-informed (MRD), and chaining pre-analytic QC, UMI/duplex error-suppression (fgbio), fragment QC, ichorCNA tumor fraction (sWGS) or VarDict low-VAF calling (panel), CHIP subtraction against matched WBC, optional fragmentomics/methylation, and longitudinal tracking. Use when treating pre-analytics as the irreversible sensitivity ceiling (tube/time-to-plasma/hemolysis), running error-suppression BEFORE calling (single-strand consensus does not remove deamination; only duplex does), reporting a VAF only with input genome-equivalents (TF ~ 2x VAF only for clonal-het-diploid), subtracting CHIP before reporting somatic, or keeping tube/panel/pipeline identical across a longitudinal MRD series. Hands mechanism to the liquid-biopsy component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/liquid-biopsy-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'liquid-biopsy/cfdna-preprocessing', 'id': 'cfdna-preprocessing', 'status': 'bundled'}`, `{'upstream_path': 'liquid-biopsy/analytical-validation', 'id': 'analytical-validation', 'status': 'bundled'}`, `{'upstream_path': 'liquid-biopsy/ctdna-mutation-detection', 'id': 'ctdna-mutation-detection', 'status': 'bundled'}`, `{'upstream_path': 'liquid-biopsy/tumor-fraction-estimation', 'id': 'tumor-fraction-estimation', 'status': 'bundled'}`, `{'upstream_path': 'liquid-biopsy/fragment-analysis', 'id': 'fragment-analysis', 'status': 'bundled'}`, `{'upstream_path': 'liquid-biopsy/methylation-based-detection', 'id': 'methylation-based-detection', 'status': 'bundled'}`, `{'upstream_path': 'liquid-biopsy/longitudinal-monitoring', 'id': 'longitudinal-monitoring', 'status': 'bundled'}`
- **别名：**`bio-workflows-liquid-biopsy-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-liquid-biopsy-pipeline
dnakit skills install workflows-liquid-biopsy-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-longread-sv-pipeline` {#workflows-longread-sv-pipeline}

Orchestrates an end-to-end long-read structural-variant pipeline - basecalling to minimap2 alignment (platform-matched preset) to Sniffles2/cuteSV/pbsv calling to optional assembly-based calling (dipcall/PAV) to two-step .snf cohort merging to Truvari benchmarking - chaining ONT and PacBio HiFi runs while handing the SV signal mechanism off to the component skills. Use when running a long-read SV workflow from reads to a benchmarked callset, choosing the minimap2 preset and SV caller by platform and goal, deciding when long reads are worth it for the insertions and repeat-mediated SVs short reads physically miss, building a joint-genotyped cohort with the two-step .snf design, or parameterizing a Truvari benchmark against GIAB HG002 Tier 1 plus CMRG. Not for the SV signal mechanism itself (see variant-calling/structural-variant-calling) or short-read SV.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/longread-sv-pipeline)
- **许可证：**`MIT`
- **运行依赖：**命令 `Sniffles`
- **依赖的 Skills：**`{'upstream_path': 'long-read-sequencing/basecalling', 'id': 'long-read-sequencing-basecalling', 'status': 'bundled'}`, `{'upstream_path': 'long-read-sequencing/long-read-alignment', 'id': 'long-read-sequencing-long-read-alignment', 'status': 'bundled'}`, `{'upstream_path': 'long-read-sequencing/long-read-qc', 'id': 'long-read-sequencing-long-read-qc', 'status': 'bundled'}`, `{'upstream_path': 'long-read-sequencing/structural-variants', 'id': 'long-read-sequencing-structural-variants', 'status': 'bundled'}`
- **别名：**`bio-workflows-longread-sv-pipeline`, `long-read-genomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-longread-sv-pipeline
dnakit skills install workflows-longread-sv-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-metagenomics-pipeline` {#workflows-metagenomics-pipeline}

End-to-end shotgun metagenomics workflow from FASTQ to taxonomic and functional profiles, orchestrating controls/host depletion, Kraken2+Bracken classification, MetaPhlAn marker profiling, and HUMAnN functional profiling. Covers the controls-first ordering, why Kraken2 read counts are not abundances and MetaPhlAn cell fractions do not equal Bracken read fractions, and the consistent-pipeline framing. Use when profiling shotgun metagenomic samples end to end, or chaining classification, abundance, and function. For resistome see metagenomics/amr-detection; for strains see metagenomics/strain-tracking; for assembly see genome-assembly/metagenome-assembly.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/metagenomics-pipeline)
- **许可证：**`MIT`
- **运行依赖：**命令 `Kraken2`
- **依赖的 Skills：**`{'upstream_path': 'read-qc/fastp-workflow', 'id': 'read-qc-fastp-workflow', 'status': 'bundled'}`, `{'upstream_path': 'metagenomics/contamination-controls', 'id': 'metagenomics-contamination-controls', 'status': 'bundled'}`, `{'upstream_path': 'metagenomics/kraken-classification', 'id': 'metagenomics-kraken', 'status': 'bundled'}`, `{'upstream_path': 'metagenomics/metaphlan-profiling', 'id': 'metagenomics-metaphlan', 'status': 'bundled'}`, `{'upstream_path': 'metagenomics/abundance-estimation', 'id': 'metagenomics-abundance', 'status': 'bundled'}`, `{'upstream_path': 'metagenomics/functional-profiling', 'id': 'metagenomics-functional-profiling', 'status': 'bundled'}`, `{'upstream_path': 'metagenomics/metagenome-visualization', 'id': 'metagenomics-visualization', 'status': 'bundled'}`
- **别名：**`bio-workflows-metagenomics-pipeline`, `claw-metagenomics`, `metagenomics`, `metagenomics-workflow`, `tooluniverse-metagenomics-analysis`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-metagenomics-pipeline
dnakit skills install workflows-metagenomics-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-methylation-pipeline` {#workflows-methylation-pipeline}

Orchestrates the end-to-end bisulfite/EM-seq methylation pipeline from FASTQ to differentially methylated regions, chaining Trim Galore/fastp QC, Bismark alignment + deduplication, methylation calling, methylKit coverage-filtering/normalization, and selection-aware DMR detection (dmrseq/DSS). Use when gating the run on bisulfite conversion (lambda + pUC19 controls) BEFORE any beta value, committing the genome build + library directionality once, keeping mate-overlap deduplicated (--no_overlap), M-bias-trimming from the plot, filtering coverage before testing, choosing a count model (beta-binomial/DSS) over a bare-beta t-test, or using a region-selection-aware FDR (dmrseq/DSS) rather than raw methylKit tiles. Hands mechanism to the methylation-analysis component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/methylation-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'read-qc/fastp-workflow', 'id': 'read-qc-fastp-workflow', 'status': 'bundled'}`, `{'upstream_path': 'methylation-analysis/bismark-alignment', 'id': 'methylation-bismark-alignment', 'status': 'bundled'}`, `{'upstream_path': 'methylation-analysis/methylation-calling', 'id': 'methylation-calling', 'status': 'bundled'}`, `{'upstream_path': 'methylation-analysis/methylkit-analysis', 'id': 'methylation-methylkit', 'status': 'bundled'}`, `{'upstream_path': 'methylation-analysis/differential-cpg-testing', 'id': 'methylation-differential-cpg', 'status': 'bundled'}`, `{'upstream_path': 'methylation-analysis/dmr-detection', 'id': 'methylation-dmr-detection', 'status': 'bundled'}`
- **别名：**`bio-workflows-methylation-pipeline`, `methylation-analysis`, `pipeline-wgbs`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-methylation-pipeline
dnakit skills install workflows-methylation-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-microbiome-pipeline` {#workflows-microbiome-pipeline}

End-to-end 16S/ITS amplicon workflow from demultiplexed FASTQ to a consensus differential-abundance result, orchestrating cutadapt primer removal, per-run DADA2 ASV inference (learnErrors/mergeSequenceTables/removeBimeraDenovo), region-matched taxonomy assignment, a SEPP/Greengenes2 tree, alpha/beta diversity at a declared sampling depth (phyloseq/vegan, adonis2 paired with betadisper), compositional DA as a consensus of &gt;=2 tools (ALDEx2/ANCOM-BC2) on unrarefied counts, and optional PICRUSt2 functional prediction gated on NSTI. Covers the stage-ordering decisions (primers before truncation, per-run error model, rarefy for diversity not DA, predicted potential not activity) and defers each per-step choice to the six microbiome skills. Use when staging an amplicon study end to end or chaining ASV inference, taxonomy, diversity, and differential abundance. For shotgun reads see workflows/metagenomics-pipeline.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/microbiome-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'read-qc/adapter-trimming', 'id': 'read-qc-adapter-trimming', 'status': 'bundled'}`, `{'upstream_path': 'microbiome/amplicon-processing', 'id': 'microbiome-amplicon-processing', 'status': 'bundled'}`, `{'upstream_path': 'microbiome/taxonomy-assignment', 'id': 'microbiome-taxonomy-assignment', 'status': 'bundled'}`, `{'upstream_path': 'microbiome/diversity-analysis', 'id': 'microbiome-diversity-analysis', 'status': 'bundled'}`, `{'upstream_path': 'microbiome/differential-abundance', 'id': 'microbiome-differential-abundance', 'status': 'bundled'}`, `{'upstream_path': 'microbiome/functional-prediction', 'id': 'microbiome-functional-prediction', 'status': 'bundled'}`
- **别名：**`bio-workflows-microbiome-pipeline`, `microbiome-amplicon`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-microbiome-pipeline
dnakit skills install workflows-microbiome-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-multiome-pipeline` {#workflows-multiome-pipeline}

Orchestrates the end-to-end 10x Multiome (paired scRNA + scATAC) pipeline from Cell Ranger ARC output to a jointly-embedded, annotated object, chaining per-modality QC, AMULET fragment-based ATAC doublet detection, per-modality normalization (RNA SCT/PCA; ATAC TF-IDF/LSI), WNN (or MultiVI) integration, joint clustering, RNA-based annotation, and LinkPeaks peak-to-gene linking. Use when enforcing the shared cell-barcode intersection between modalities (cellranger-ARC not -atac), keeping per-modality QC/doublets before the joint embedding, dropping the depth-correlated LSI component, annotating identity from RNA (ATAC is regulatory state), treating peak-to-gene links as correlational hypotheses, or aggregating to pseudobulk for cross-condition DE. Hands mechanism to the single-cell and atac-seq component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/multiome-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'single-cell/data-io', 'id': 'single-cell-data-io', 'status': 'bundled'}`, `{'upstream_path': 'single-cell/preprocessing', 'id': 'single-cell-preprocessing', 'status': 'bundled'}`, `{'upstream_path': 'single-cell/clustering', 'id': 'single-cell-clustering', 'status': 'bundled'}`, `{'upstream_path': 'single-cell/multimodal-integration', 'id': 'single-cell-multimodal-integration', 'status': 'bundled'}`, `{'upstream_path': 'single-cell/scatac-analysis', 'id': 'atac-seq-single-cell-atac', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/single-cell-atac', 'id': 'atac-seq-single-cell-atac', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/co-accessibility', 'id': 'atac-seq-co-accessibility', 'status': 'bundled'}`, `{'upstream_path': 'atac-seq/motif-deviation', 'id': 'atac-seq-motif-deviation', 'status': 'bundled'}`
- **别名：**`bio-workflows-multiome-pipeline`, `multiome-scatac`, `single-cell-encode`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-multiome-pipeline
dnakit skills install workflows-multiome-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-neoantigen-pipeline` {#workflows-neoantigen-pipeline}

Orchestrates neoantigen discovery from somatic variants to ranked vaccine candidates, chaining HLA typing (OptiType/arcasHLA + LOHHLA), VEP annotation (Wildtype+Frameshift plugins) + expression/readcount annotation, proximal-variant phasing, pVACseq MHC-I/II binding, CCF/clonality, and immunogenicity/quality ranking. Use when recognizing that binding is single-digit PPV and the critical steps are downstream (full-resolution HLA + LOH gating, proximal-variant phasing, clonality from purity+CN not raw VAF, expression), sequencing normalize+annotate -&gt; phase -&gt; HLA -&gt; binding -&gt; quality in the defensible order, dropping candidates on LOH-lost alleles, supplying --phased-proximal-variants-vcf so the mutant peptide is real, or ranking WITHIN patient rather than a fixed IC50 threshold. Hands mechanism to the immunoinformatics component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/neoantigen-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'clinical-databases/hla-typing', 'id': 'clinical-databases-hla-typing', 'status': 'bundled'}`, `{'upstream_path': 'immunoinformatics/mhc-binding-prediction', 'id': 'immunoinformatics-mhc-binding-prediction', 'status': 'bundled'}`, `{'upstream_path': 'immunoinformatics/mhc-class-ii-prediction', 'id': 'immunoinformatics-mhc-class-ii-prediction', 'status': 'bundled'}`, `{'upstream_path': 'immunoinformatics/neoantigen-prediction', 'id': 'immunoinformatics-neoantigen-prediction', 'status': 'bundled'}`, `{'upstream_path': 'immunoinformatics/immunogenicity-scoring', 'id': 'immunoinformatics-immunogenicity-scoring', 'status': 'bundled'}`, `{'upstream_path': 'immunoinformatics/epitope-prediction', 'id': 'immunoinformatics-epitope-prediction', 'status': 'bundled'}`
- **别名：**`bio-workflows-neoantigen-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-neoantigen-pipeline
dnakit skills install workflows-neoantigen-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-outbreak-pipeline` {#workflows-outbreak-pipeline}

Orchestrates genomic-epidemiology outbreak investigation from pathogen isolates to transmission networks, forking bacterial (snippy -&gt; Gubbins recombination-masking -&gt; IQ-TREE -&gt; TreeTime -&gt; TransPhylo) vs viral (Nextstrain/augur), with parallel MLST typing (cgMLST delegated to epidemiological-genomics/pathogen-typing) and AMR surveillance. Use when committing ONE reference genome for SNP calling (every isolate and distance inherits its coordinates), applying MANDATORY Gubbins recombination-masking on core.full.aln before the tree for recombining bacteria (skipping it inflates the clock 2-5x), gating time-scaling on a temporal-signal test (TempEst R2 &gt;= 0.3), using a pathogen- AND population-specific cluster threshold rather than a universal SNP cutoff, or pinning pangolin-data/Nextclade/Freyja versions for the viral route. Hands mechanism to the epidemiological-genomics component skills; not a re-teach of any single step.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/outbreak-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'epidemiological-genomics/pathogen-typing', 'id': 'epidemiological-genomics-pathogen-typing', 'status': 'bundled'}`, `{'upstream_path': 'epidemiological-genomics/amr-surveillance', 'id': 'epidemiological-genomics-amr-surveillance', 'status': 'bundled'}`, `{'upstream_path': 'epidemiological-genomics/phylodynamics', 'id': 'epidemiological-genomics-phylodynamics', 'status': 'bundled'}`, `{'upstream_path': 'epidemiological-genomics/transmission-inference', 'id': 'epidemiological-genomics-transmission-inference', 'status': 'bundled'}`, `{'upstream_path': 'epidemiological-genomics/variant-surveillance', 'id': 'epidemiological-genomics-variant-surveillance', 'status': 'bundled'}`
- **别名：**`bio-workflows-outbreak-pipeline`, `pathogen-epi-genomics`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-outbreak-pipeline
dnakit skills install workflows-outbreak-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-somatic-variant-pipeline` {#workflows-somatic-variant-pipeline}

Chains a somatic (tumor-normal) SNV/indel and structural-variant pipeline end to end with GATK Mutect2 (or Strelka2), wiring the somatic-specific machinery - panel-of-normals and gnomAD germline-resource priors, GetPileupSummaries/CalculateContamination, and LearnReadOrientationModel FFPE/oxoG orientation-bias filtering fed into FilterMutectCalls. Use when calling somatic mutations from a tumor-normal pair (or tumor-only with PoN caveats), deciding which artifact filter removes which class of false positive, reasoning about VAF/purity/ploidy and clonal-vs-subclonal detection, adding somatic SV/CNV or TMB/MSI/signatures, or routing variants to AMP/ASCO/CAP tier and oncogenicity interpretation (never germline ACMG).

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/somatic-variant-pipeline)
- **许可证：**`MIT`
- **依赖的 Skills：**`{'upstream_path': 'read-alignment/bwa-alignment', 'id': 'read-alignment-bwa-alignment', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/gatk-variant-calling', 'id': 'gatk-variant-calling', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/filtering-best-practices', 'id': 'variant-calling-filtering-best-practices', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/structural-variant-calling', 'id': 'variant-calling-structural-variant-calling', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/variant-annotation', 'id': 'variant-annotation', 'status': 'bundled'}`, `{'upstream_path': 'variant-calling/clinical-interpretation', 'id': 'clinical-databases-acmg-classification', 'status': 'bundled'}`, `{'upstream_path': 'copy-number/cnvkit-analysis', 'id': 'copy-number-cnvkit-analysis', 'status': 'bundled'}`
- **别名：**`bio-workflows-somatic-variant-pipeline`, `gatk4-somatic-variant-caller`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-somatic-variant-pipeline
dnakit skills install workflows-somatic-variant-pipeline --agent claude-code   # 或 --agent codex / cursor
```

## `workflows-tcr-pipeline` {#workflows-tcr-pipeline}

Orchestrates an end-to-end immune-repertoire pipeline from FASTQ to clonotypes, diversity, overlap, somatic hypermutation and lineages, routing on two forks. Use when deciding bulk vs single-cell (bulk amplicon/RNA-seq -&gt; MiXCR analyze preset -&gt; VDJtools/immunarch depth-normalized diversity and overlap -&gt; figures; 10x paired VDJ -&gt; MiXCR 10x preset or Cell Ranger -&gt; scirpy gene-expression integration, chain QC, clonotype clusters); and TCR vs BCR (TCR -&gt; exact CDR3-nt+V/J clonotypes, VDJtools diversity is fine; BCR -&gt; somatic hypermutation makes exact clonotypes wrong -&gt; Immcantation distToNearest/findThreshold clonal clustering, germline reconstruction, SHM, Dowser lineages); selecting the MiXCR 4.x preset by chemistry and activating its license; downsampling to equal depth before comparing diversity or overlap; and optionally annotating antigen specificity.

- **来源：**[GPTomics bioSkills contributors @ `d91ed3d`](https://github.com/GPTomics/bioSkills/tree/d91ed3d563019e649dc854c56ccd62551359488a/workflows/tcr-pipeline)
- **许可证：**`MIT`
- **运行依赖：**命令 `MiXCR`
- **依赖的 Skills：**`{'upstream_path': 'tcr-bcr-analysis/mixcr-analysis', 'id': 'tcr-bcr-analysis-mixcr-analysis', 'status': 'bundled'}`, `{'upstream_path': 'tcr-bcr-analysis/vdjtools-analysis', 'id': 'tcr-bcr-analysis-vdjtools-analysis', 'status': 'bundled'}`, `{'upstream_path': 'tcr-bcr-analysis/immcantation-analysis', 'id': 'tcr-bcr-analysis-immcantation-analysis', 'status': 'bundled'}`, `{'upstream_path': 'tcr-bcr-analysis/scirpy-analysis', 'id': 'tcr-bcr-analysis-scirpy-analysis', 'status': 'bundled'}`, `{'upstream_path': 'tcr-bcr-analysis/repertoire-visualization', 'id': 'tcr-bcr-analysis-repertoire-visualization', 'status': 'bundled'}`, `{'upstream_path': 'tcr-bcr-analysis/specificity-annotation', 'id': 'tcr-bcr-analysis-specificity-annotation', 'status': 'bundled'}`
- **别名：**`bio-workflows-tcr-pipeline`
- **核对：**元数据 通过，资源 已校验，执行 未运行

```bash
dnakit skills info workflows-tcr-pipeline
dnakit skills install workflows-tcr-pipeline --agent claude-code   # 或 --agent codex / cursor
```
