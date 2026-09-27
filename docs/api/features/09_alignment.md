# 序列距离与比对

比较两条 DNA 序列：“差异有多大”的序列距离，“碱基如何对应”的成对序列比对。

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
- **输入：** `DNA([...])` 或 ID→序列字典、新输出目录；PATH 中需有 MAFFT、MUSCLE 5 或 Clustal Omega，并安装 `bio-formats`。
- **输出：** 返回外部工具生成的 `MultipleAlignment`，包含对齐后的序列。
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

## 8) 将 reads 比对到参考基因组 {#bio-15}

- **作用：** 将单端或双端 reads 转为 SAM 比对文件。
- **计算方法：** 在新目录建立私有索引，然后运行 BWA MEM 或 Bowtie2。
- **API：** `dnakit.bioinformatics.ngs.map_reads(reference, first, output_directory, second=None, method='bwa', config=None)`
- **输入：** 参考 FASTA、R1/R2 FASTQ；所选比对程序在 PATH 中；新输出目录。
- **输出：** 返回 `ToolResult`，包含后端 `backend`、产物路径 `artifacts` 和运行日志 `log`。
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
