# 序列距离与比对

比较两条 DNA 序列：“差异有多大”的序列距离，“碱基如何对应”的成对序列比对。

## 使用 SPOA 生成共识序列

- **作用：** 调用用户安装的 SPOA，对 FASTA/FASTQ reads 生成共识 FASTA。
- **API：** `dnakit.alignment.consensus_spoa(source, output_directory, algorithm=2, strand_ambiguous=False, ...)`。
- **输出：** 新目录中的 `consensus.fasta`、标准错误和结构化执行信封；只固定暴露共识结果模式。
- **引用与许可：** Christopher Lee (2003), [10.1093/bioinformatics/btg109](https://doi.org/10.1093/bioinformatics/btg109)；SPOA 官方仓库 MIT。完整图、比对和打分参数没有全部重导出。

## Picard 标记重复

- **作用：** 调用调用者安装的 Picard `MarkDuplicates`，生成标记重复后的 BAM 和重复统计文件。
- **API：** `dnakit.ngs.mark_duplicates_picard(source, output_directory, remove_duplicates=False, assume_sorted=True, ...)`。
- **输出：** 新目录中的 `marked_duplicates.bam`、`duplicate_metrics.txt`、诊断日志和 JSON 可序列化执行信封。
- **引用与许可：** Broad Institute, *Picard Toolkit* (2019)，[官方仓库](https://github.com/broadinstitute/picard)；MIT。DNAKit 只暴露固定的 MarkDuplicates 参数，不捆绑 Picard 或 Java。

<span id="1"></span>

## 1) Hamming distance {#11-sim-006-hamming-distance}

- **作用：** 逐个位置比较两条等长序列，返回不同碱基的数量和具体位置，用于统计点突变或快速判断等长序列差异。
- **计算方法：** 逐位置计数不同字符；输入必须等长。

- **API：** `dnakit.similarity.hamming_distance(left[必须], right[必须], max_distance[可选])`。

- **输入：** 必填两条等长序列；可选 `max_distance` 审计阈值。
- **输出：** 返回 `DistanceResult`，示例展示 `distance`、`mismatches`。

- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import hamming_distance

result = hamming_distance(DNASequence("ACGT"), DNASequence("ACCT"))
print(result.distance, [item.position for item in result.mismatches])
```

- **示例结果：**

```text
1.0 [2]
```

## 2) Edit distance {#12-sim-007-edit-distance}

- **作用：** 计算把一条序列变成另一条序列所需的最小替换、插入和删除代价，并可返回最优编辑路径，用于量化不同长度序列的差异。
- **计算方法：** 通过动态规划求替换、插入和删除的最小累计代价；需要时回溯编辑路径。

- **API：** `dnakit.similarity.edit_distance(left[必须], right[必须], substitution_cost[可选], insertion_cost[可选], deletion_cost[可选], max_distance[可选], return_path[可选], max_cells[可选])`。

- **输入：** 必填两条序列；可选操作代价、最大距离、编辑路径和 DP 上限。
- **输出：** 返回 `DistanceResult`，示例展示 `distance`、`edit_path`。

- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import edit_distance

result = edit_distance(DNASequence("ACGT"), DNASequence("AGT"), return_path=True)
print(result.distance, [step.operation for step in result.edit_path or ()])
```

- **示例结果：**

```text
1.0 ['match', 'delete', 'match', 'match']
```

<span id="2-sim-008"></span>

## 3) 全局比对 {#21}

- **作用：** 从头到尾对齐两条完整序列，返回加入 Gap 后的对应关系、得分、Identity 和 Coverage，适合比较长度接近且整体相关的序列。
- **计算方法：** 动态规划优化覆盖两条完整序列的累计匹配、错配和 Gap 得分。

- **API：** `dnakit.alignment.align_pairwise(query[必须], target[必须], config[可选])`。
- **输入：** 两条 DNA 序列；可选 `AlignmentConfig(mode="global")` 设置匹配/错配分数、线性或 affine Gap 分数及 `max_cells`。返回 `AlignmentResult`，含加入 Gap 后的对应序列、得分、一致性和覆盖度。
- **输出：** 返回 `AlignmentResult`，包含对齐后的序列、分数和坐标。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

result = align_pairwise(
    DNASequence("ACGT"),
    DNASequence("AGT"),
    config=AlignmentConfig(mode="global", gap_score=-1),
)
print(result.aligned_query, result.aligned_target, result.score, result.identity)
```

- **示例结果：**

```text
ACGT A-GT 2.0 0.75
```

## 4) 局部比对 {#22}

- **作用：** 寻找两条序列中得分最高的局部对齐区域，返回该片段在两条原序列中的坐标、Identity 和 Coverage，适合发现共同片段。
- **计算方法：** 动态规划允许得分从零重新开始，并回溯最高分的局部片段。

- **API：** `dnakit.alignment.align_pairwise(query[必须], target[必须], config[可选])`。
- **输入：** 两条 DNA 序列；可选 `AlignmentConfig(mode="local")` 设置匹配/错配分数、线性或 affine Gap 分数及 `max_cells`。返回 `AlignmentResult`，含加入 Gap 后的对应序列、得分、一致性和覆盖度。
- **输出：** 返回 `AlignmentResult`，包含对齐后的序列、分数和坐标。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

result = align_pairwise(
    DNASequence("TTACGTAA"),
    DNASequence("GGACGTCC"),
    config=AlignmentConfig(mode="local", mismatch_score=-2, gap_score=-2),
)
print(
    result.aligned_query,
    result.aligned_target,
    result.identity,
    result.query_coverage,
    result.target_coverage,
)
```

- **示例结果：**

```text
ACGT ACGT 1.0 0.5 0.5
```

## 5) 半全局比对 {#23}

- **作用：** 在两端未对齐区域不扣 Gap 分，返回核心对齐及其原序列坐标，适合把引物、扩增子或短序列与长序列中的目标区域进行比较。
- **计算方法：** 动态规划免除两端未对齐区域的 Gap 惩罚，回溯核心匹配区域。

- **API：** `dnakit.alignment.align_pairwise(query[必须], target[必须], config[可选])`。
- **输入：** 两条 DNA 序列；可选 `AlignmentConfig(mode="semi_global")` 设置匹配/错配分数、线性或 affine Gap 分数及 `max_cells`。返回 `AlignmentResult`，含加入 Gap 后的对应序列、得分、一致性和覆盖度。
- **输出：** 返回 `AlignmentResult`，包含对齐后的序列、分数和坐标。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.alignment import AlignmentConfig, align_pairwise

result = align_pairwise(
    DNASequence("ACGT"),
    DNASequence("TTACGTAA"),
    config=AlignmentConfig(mode="semi_global", mismatch_score=-2, gap_score=-2),
)
print(
    result.aligned_query,
    result.aligned_target,
    result.query_coverage,
    result.target_coverage,
    (result.target_start, result.target_end),
)
```

- **示例结果：**

```text
ACGT ACGT 1.0 0.5 (2, 6)
```

## 6) 计算多序列比对的距离 {#bio-1}

- **作用：** 把已经对齐的序列转换为样本间距离矩阵。
- **计算方法：** 逐对统计不同碱基比例；默认忽略含缺口的位置。
- **API：** `dnakit.bioinformatics.alignment.MultipleAlignment.from_sequences(sequences)`; `MultipleAlignment.distances(gap_policy='pairwise_delete')`
- **输入：** 唯一 ID 对应的等宽序列；允许缺口 `-`。
- **输出：** 返回多序列比对的成对距离矩阵。
- **示例代码：**

```python
from dnakit.bioinformatics.alignment import MultipleAlignment

msa = MultipleAlignment.from_sequences({"a": "ACGT", "b": "ACTT", "c": "TCGT"})
print(msa.ids)
print(msa.distances().tolist())
```

- **示例结果：**

```text
('a', 'b', 'c')
[[0.0, 0.25, 0.25], [0.25, 0.0, 0.5], [0.25, 0.5, 0.0]]
```

依赖与说明：[FAQ](../../faq/bioinformatics.md#alignment)。

## 7) 调用 MAFFT、MUSCLE 或 Clustal Omega 进行比对 {#bio-3}

- **作用：** 把未对齐的 DNA 序列转为多序列比对。
- **计算方法：** 执行选定的外部程序，再验证返回序列的 ID 和内容。
- **API：** `dnakit.bioinformatics.alignment.align_sequences(sequences, output_directory, method='mafft', config=None)`
- **Agent/MCP API：** `dnakit.alignment.align_sequences(sequences, output_directory, method='mafft', ...)`；返回 JSON 可序列化的阶段信封，固定支持 `mafft`、`muscle` 和 `clustalo`。
- **输入：** `DNA([...])` 或 ID→序列字典、新输出目录；PATH 中需有 MAFFT、MUSCLE 5 或 Clustal Omega，并安装 `bio-formats`。
- **输出：** Python API 返回外部工具生成的 `MultipleAlignment`；Agent/MCP 门面返回包含对齐结果、产物路径和日志的 JSON 信封。
- **示例代码：**

```python
from dnakit import DNA
from dnakit.bioinformatics.alignment import align_sequences

msa = align_sequences(DNA(["ACGTACGT", "ACGACGT"]), "new_msa", method="mafft")
print(msa.ids, msa.sequences)
```

- **示例结果：**

返回 `MultipleAlignment`；`ids` 保留输入 ID，`sequences` 为等宽比对序列。具体缺口位置由外部程序决定，此处未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#alignment)。

## 8) 调用 NCBI BLAST+、minimap2 和 BEDTools 的统一门面 {#bio-16}

- **作用：** 通过统一的 JSON/MCP 安全接口执行本地核酸搜索、长读长比对和 BED 区间交集。
- **API：** `dnakit.alignment.build_blast_database(fasta, output_directory, ...)`、`dnakit.alignment.local_blast(query, database, output_directory, ...)`、`dnakit.alignment.map_long_reads(reference, reads, output_directory, ...)`、`dnakit.alignment.intersect_bed(first, second, output_directory, ...)`。
- **输入：** 本地 FASTA/FASTQ/BED 路径和新输出目录；BLAST+、minimap2 或 BEDTools 由调用者安装并置于 PATH，或通过 `program` 指定可执行文件。
- **输出：** JSON 可序列化的阶段信封，包含后端、产物路径、日志、固定操作和来源信息；外部程序不会自动安装或下载。
- **引用与许可：** NCBI BLAST+ 使用 Camacho et al. (2009), DOI `10.1186/1471-2105-10-421`；minimap2 使用 Li (2018), DOI `10.1093/bioinformatics/bty191`（MIT）；BEDTools 使用 Quinlan & Hall (2010), DOI `10.1093/bioinformatics/btq033`（MIT）。
- **边界：** 只开放固定操作和有界参数子集，不宣称重导出上游完整 CLI。

```python
from dnakit.alignment import intersect_bed, map_long_reads

mapping = map_long_reads("reference.fa", "reads.fastq", "new_minimap", preset="map-ont")
overlap = intersect_bed("a.bed", "b.bed", "new_intersection")
print(mapping["result"]["artifacts"], overlap["result"]["artifacts"])
```

## 9) 将 reads 比对到参考基因组 {#bio-15}

- **作用：** 将单端或双端 reads 转为 SAM 比对文件。
- **计算方法：** 在新目录建立私有索引，然后运行 BWA MEM、BWA-MEM2 或 Bowtie2。
- **API：** `dnakit.bioinformatics.ngs.map_reads(reference, first, output_directory, second=None, method='bwa', config=None)`
- **输入：** 参考 FASTA、R1/R2 FASTQ；所选比对程序在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
- **Agent/MCP API：** `dnakit.ngs.map_reads(reference, read1, output_directory, read2=None, method='bwa', ...)`；`method` 可选 `bwa`、`bwa-mem2`、`bowtie2`；返回 JSON 可序列化的阶段信封，内部包含产物路径与日志。
- **引用与许可：** BWA 使用 Li & Durbin (2009), DOI `10.1093/bioinformatics/btp324`（GPL-3.0）；BWA-MEM2 使用 Vasimuddin et al. (2019), DOI `10.1109/IPDPS.2019.00041`（MIT）；Bowtie2 使用 Langmead & Salzberg (2012), DOI `10.1038/nmeth.1923`（GPL-3.0）。后端由调用者安装。
- **示例代码：**

```python
from dnakit.bioinformatics import ExecutionConfig
from dnakit.bioinformatics.ngs import map_reads

result = map_reads("reference.fa", "R1.fq.gz", "new_mapping",
                   second="R2.fq.gz", config=ExecutionConfig(threads=4))
print(result.artifacts["aligned.sam"])
```

- **示例结果：**

返回 `ToolResult`；`artifacts['aligned.sam']` 是 SAM 文件路径。外部程序示例未执行。

依赖与说明：[FAQ](../../faq/bioinformatics.md#reads)。
