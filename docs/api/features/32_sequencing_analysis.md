# 测序数据分析

## 1) 统计 BAM 并生成区域共识 {#bio-14}

- **作用：** 查看比对记录统计、各位点碱基深度与多数共识。
- **计算方法：** 解析 BAM flags；筛选质量合格的 A/C/G/T 比对，平票或深度不足返回 N。
- **API：** `dnakit.bioinformatics.reads.alignment_summary(path)`; `dnakit.bioinformatics.reads.coverage_consensus(path, region, min_base_quality=20, min_mapping_quality=20, min_depth=1)`
- **输入：** 已排序并索引的 `aligned.bam`；需要 `bio-reads`。共识保持参考正链坐标，不加入插入/删除。
- **输出：** `alignment_summary` 返回比对数量字典；`coverage_consensus` 返回共识序列和覆盖信息。
- **示例代码：**

```python
from dnakit.bioinformatics.reads import Region, alignment_summary, coverage_consensus

summary = alignment_summary("aligned.bam")
result = coverage_consensus("aligned.bam", Region("chr1", 100, 180))
print(summary)
print(result["consensus"])
```

- **示例结果：**

`summary` 为计数字典；`result` 含区域、A/C/G/T 计数、深度与共识序列。需要真实 BAM，此处未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#reads)。

## 2) 一键运行单样本 DNA-seq 流程 {#bio-20}

- **作用：** 从 reads 依次得到质控报告、BAM 和 VCF。
- **计算方法：** 串联 FastQC → BWA → samtools → bcftools，并保留各阶段结果。
- **API：** `dnakit.bioinformatics.ngs.dna_seq_workflow(reference, first, output_directory, second=None, ploidy=2, qc=True, config=None)`
- **输入：** 带 `.fai` 的参考 FASTA 和 FASTQ；对应四个程序在 PATH 中；新输出目录。
- **输出：** 返回按阶段名索引的 `ToolResult` 字典，包含 VCF 等产物路径。
- **示例代码：**

```python
from dnakit.bioinformatics import ExecutionConfig
from dnakit.bioinformatics.ngs import dna_seq_workflow

stages = dna_seq_workflow("reference.fa", "R1.fq.gz", "new_dna_analysis",
                          second="R2.fq.gz", ploidy=2, config=ExecutionConfig(threads=4))
print(stages["variants"].artifacts["variants.vcf.gz"])
```

- **示例结果：**

返回阶段名到 `ToolResult` 的字典：默认包含 `qc`、`mapping`、`bam`、`variants`。外部程序示例未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#reads)。

## 3) 运行单样本 RNA-seq 流程 {#bio-21}

- **作用：** 从 RNA-seq reads 得到 BAM、转录本和基因计数。
- **计算方法：** 依次运行 HISAT2、samtools、StringTie 和 featureCounts。
- **API：** `dnakit.bioinformatics.ngs.rna_seq_workflow(reference, annotation, first, output_directory, second=None, strandedness=0)`
- **输入：** 参考 FASTA、GTF 和 FASTQ；上述程序在 PATH 中；`strandedness` 为 0/1/2；新输出目录。
- **输出：** 返回按阶段名索引的 `ToolResult` 字典，包含表达计数表等产物路径。
- **示例代码：**

```python
from dnakit.bioinformatics.ngs import rna_seq_workflow

stages = rna_seq_workflow("reference.fa", "genes.gtf", "R1.fq.gz",
                          "new_rna_analysis", second="R2.fq.gz", strandedness=0)
print(stages["rna"].artifacts["counts.tsv"])
```

- **示例结果：**

返回含 `bam`、`rna` 的字典；RNA 阶段含 `transcripts.gtf` 和 `counts.tsv`。外部程序示例未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#reads)。

## 4) 从 ChIP-seq BAM 检测富集峰 {#bio-22}

- **作用：** 比较 treatment/control BAM，得到富集峰与覆盖轨迹。
- **计算方法：** 使用 MACS3 callpeak 进行峰检测。
- **API：** `dnakit.bioinformatics.ngs.chip_seq_workflow(treatment, control, output_directory, genome_size=..., paired_end=False, qvalue=0.05)`
- **输入：** 已处理的 treatment/control BAM、有效基因组大小、文库单双端信息；MACS3 在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
- **示例代码：**

```python
from dnakit.bioinformatics.ngs import chip_seq_workflow

result = chip_seq_workflow("treatment.bam", "control.bam", "new_chip_analysis",
                           genome_size=2700000000, paired_end=True)
print(sorted(result.artifacts))
```

- **示例结果：**

返回 `ToolResult`；`artifacts` 包含 MACS3 生成的峰文件和 coverage 文件。外部程序示例未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#reads)。

## 5) 生成 Sanger 参考坐标共识 {#bio-56}

- **作用：** 合并一个或多个 Sanger reads，检查构建序列。
- **计算方法：** 对裁剪后的 read 做双向局部比对，筛选质量、identity 和覆盖度后逐位投票。
- **API：** `dnakit.bioinformatics.cloning.sanger_consensus(traces, reference, min_quality=20, min_identity=0.8, min_read_coverage=0.8, min_depth=1)`
- **输入：** 含 `id/sequence/qualities` 的 trace 字典和参考序列；需要 `bio-analysis`。N 表示平票或不足覆盖；共识不加入 indel。
- **输出：** 返回字典，包含共识序列、深度数组和变异信息。
- **示例代码：**

```python
from dnakit.bioinformatics.cloning import sanger_consensus

trace = {"id": "synthetic", "sequence": "AACCGGTT", "qualities": (30,) * 8}
result = sanger_consensus([trace], "AACCGGTT")
print(result["consensus"])
print(result["depth"].tolist())
```

- **示例结果：**

```text
AACCGGTT
[1, 1, 1, 1, 1, 1, 1, 1]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#optimization)。
