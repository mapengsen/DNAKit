# 安装

## 从 PyPI 安装

安装 `0.2.0`：

```bash
pip install dnakit==0.2.0
```

地址为：[pypi.org/project/dnakit](https://pypi.org/project/dnakit/)

## 其他平台

项目仓库已提供以下平台适配文件：

- [Bioconda 配方](https://github.com/mapengsen/DNAKit/tree/main/packaging/bioconda)；
- [GNU Guix 包定义](https://github.com/mapengsen/DNAKit/tree/main/packaging/guix)；
- [Galaxy Tool Shed 包装器](https://github.com/mapengsen/DNAKit/tree/main/galaxy/dnakit)。

这些适配文件固定对应 `0.1.1`，正式安装命令需要等相应平台审核并发布后才能使用。

## 开发环境

在项目根目录创建并激活完整 Conda 开发环境：

```bash
conda env create -f environment-dev.yml
conda activate dnakit-dev
```

开发环境包含测试、类型检查、文档和构建工具。可选功能的额外依赖见
[`pyproject.toml`](https://github.com/mapengsen/DNAKit/blob/main/pyproject.toml)。

## 可选生物信息后端

安装 VEP、ClinVar、dbSNP、gnomAD、dN/dS 和
Golden Gate 功能所需的可选后端：

```bash
python -m pip install "dnakit[external-tools]"
```

默认的 `pip install dnakit` 不会安装或加载该后端。

## 第三方 pip 适配器

覆盖矩阵中的 scikit-bio、seqme、DNA Chisel、pySBOL3、scadnano、polars-bio、DnaCauldron、DNA Features Viewer、gtars、GeneBlocks、sourmash 和 Cutadapt 等工具通过显式可选适配器接入统一接口；fastp、FastQC 和 MultiQC 通过用户提供的外部可执行程序接入：

```bash
python -m pip install "dnakit[bio-skibio,bio-seqme,bio-dnachisel,bio-biotite,bio-edlib,bio-geneblocks,bio-sourmash,bio-cutadapt,bio-pyranges,bio-pyfastx,bio-dnaio,bio-mappy,bio-pybigwig,bio-pybedtools,bio-gffutils,bio-pysbol3,bio-scadnano,bio-polarsbio,bio-dnacauldron,bio-dna-features-viewer,bio-gtars,bio-pyfaidx,bio-vcfpy,bio-cyvcf2]"
```

随后可使用 `dnakit.integrations.beta_diversity_scikit_bio()`、
`dnakit.integrations.permanova_scikit_bio()`、`dnakit.integrations.evaluate_seqme()` 和
`dnakit.integrations.optimize_with_dnachisel()`、`dnakit.integrations.summarize_sbol3_document()` 和
`dnakit.integrations.summarize_scadnano_design()`、`dnakit.integrations.overlap_polars_bio_bed()` 和
`dnakit.integrations.simulate_dnacauldron_assembly()` 和
`dnakit.integrations.summarize_dna_features_viewer()` 和
`dnakit.integrations.summarize_gtars_regions()` 以及
`dnakit.integrations.diff_geneblocks_sequences()` 以及
`dnakit.integrations.compare_sourmash_sequences()`；Cutadapt、fastp、FastQC 和
MultiQC 的 FASTQ 预处理与 QC 接口是 `dnakit.sequencing_qc.run_cutadapt()`、
`run_fastp()`、`run_fastqc()` 和 `run_multiqc()`，需要用户单独提供相应可执行程序和输出路径。
seqme 的模型 extras、模型权重和
DNA Chisel 的外部数据仍需按上游条款单独处理。

## 深度学习性质预测

安装 Enformer 轨道以及 NT Revised/Genomic Benchmarks 任务 checkpoint 推理依赖：

```bash
python -m pip install "dnakit[neural,neural-enformer]"
```

27 个任务 checkpoint 不随安装包分发，预测时会自动下载当前任务缺失的权重，已有文件直接复用。默认缓存目录为 `./ckpt/enformer-benchmarks/`，也可在调用时使用 `checkpoint_dir` 或
`checkpoint_path` 指定。缓存规则及各模型依赖见 [FAQ：性质预测安装与权重缓存](faq/tutorial_integrations.md#property-prediction)。调用示例按任务分为[调控轨道预测](api/features/23_prediction_regulatory_tracks.md)、[表观遗传判别](api/features/23_prediction_chromatin_classification.md)和[启动子与增强子判别](api/features/23_prediction_promoter_enhancer.md)。

## Agent 与 MCP

如需让 Codex、Claude、Cursor 等支持 MCP 的 Agent 调用 DNAKit：

```bash
python -m pip install "dnakit[agent]"
```

安装后使用 `dnakit-mcp` 启动本地 MCP 服务，具体配置见 [大模型工具 → MCP](llm_tools/mcp.md)。

# 仓库地址

[github.com/mapengsen/DNAKit](https://github.com/mapengsen/DNAKit)

## 生物信息学 Python API 可选依赖

安装当前源码的常用基因组计算能力（Python ≥3.10）：

```bash
python -m pip install -e ".[genomics]"
```

`genomics` 一并安装 NumPy、Pyrodigal、pysam 和 PyFastANI，用于基因型 PCA、基因组关系矩阵与育种预测、原核 CDS 预测、BAM 统计和基因组 ANI 计算。这些接口不要求另装 PLINK、Prodigal、samtools 或 FastANI 命令行程序；接口范围和示例见[基因组教程](api/features/38_genomics_extensions.md#pip-genomics)。

本次新增能力尚未发布到 PyPI；发布后可使用 `python -m pip install "dnakit[genomics]"`。二进制 wheel 的平台覆盖由各依赖决定，缺少匹配 wheel 时可能需要编译器；这里不包含全部 43 项外部分析，也不下载数据库或模型权重。

当前源码的新增基因组功能可独立安装：`.[genomics-genes]` 用于 Pyrodigal，`.[dna-storage]` 用于 Chamaeleo，`.[dna-nanotechnology]` 用于 scadnano；育种矩阵使用 `.[bio-analysis]`，MCP 客户端使用 `.[agent]`。

例如安装本地源码和基因识别依赖：

```bash
python -m pip install -e ".[genomics-genes]"
```

外部命令、数据库和模型权重仍需按具体分析配置，见[外部基因组分析与 MCP](api/features/39_external_genomics.md)。

当前工作区可使用 `python -m pip install -e ".[bio]"` 安装通用生物信息学 Python 分析依赖，或仅安装 `.[bio-analysis]`。新增基因组计算依赖使用上面的 `.[genomics]`。MAFFT、BLAST、HMMER、BWA、samtools、bcftools 等外部命令接口所需的可执行程序需单独安装。详见[依赖与实现方式](faq/bioinformatics.md#dependencies)。

## 安装与使用入口 {#source-python-usage}

这些修改尚未发布到 PyPI。在 DNAKit 源码根目录安装：

```bash
python -m pip install -e "."
```

Python 直接使用 `import dnakit`，不需要初始化 ToolUniverse，也不要求先安装 Skill。科学软件、模型权重、API 凭据和参考数据仍按具体任务准备。来源和许可见[第三方声明](acknowledgements.md#third-party-notices)。
