# 序列距离与比对

本页用于比较两条 DNA 序列，包括“差异有多大”的序列距离，以及“碱基如何对应”的成对序列比对。

除特别说明外，IUPAC 字符按字面符号比较，输入中的 Gap 不会被静默删除。

## 1. 序列距离

序列距离返回一个差异数值，不生成加入 Gap 后的对齐序列。

### 1.1 SIM-006 · Hamming distance

**作用：** 逐个位置比较两条等长序列，返回不同碱基的数量和具体位置，用于统计点突变或快速判断等长序列差异。

**API：** `dnakit.similarity.hamming_distance(left[必须], right[必须], max_distance[可选])`。

**输入：** 必填两条等长序列；可选 `max_distance` 审计阈值。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import hamming_distance

result = hamming_distance(DNASequence("ACGT"), DNASequence("ACCT"))
print(result.distance, [item.position for item in result.mismatches])
```

**示例结果：**

```text
1.0 [2]
```

**限制：** 序列必须等长且无 Gap；不计算插入和删除，IUPAC 不按集合兼容关系解释。

### 1.2 SIM-007 · Edit distance

**作用：** 计算把一条序列变成另一条序列所需的最小替换、插入和删除代价，并可返回最优编辑路径，用于量化不同长度序列的差异。

**API：** `dnakit.similarity.edit_distance(left[必须], right[必须], substitution_cost[可选], insertion_cost[可选], deletion_cost[可选], max_distance[可选], return_path[可选], max_cells[可选])`。

**输入：** 必填两条序列；可选操作代价、最大距离、编辑路径和 DP 上限。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import edit_distance

result = edit_distance(DNASequence("ACGT"), DNASequence("AGT"), return_path=True)
print(result.distance, [step.operation for step in result.edit_path or ()])
```

**示例结果：**

```text
1.0 ['match', 'delete', 'match', 'match']
```

**限制：** 不接受带 Gap 的输入；`max_distance` 只标记是否超过阈值，不截断精确结果。

## 2. 成对序列比对（SIM-008）

成对序列比对会在适当位置加入 Gap，生成两条相互对应的序列，并计算得分、一致性和覆盖度。

**API：** `dnakit.alignment.align_pairwise(query[必须], target[必须], config[可选])`。

**配置：** 使用 `dnakit.alignment.AlignmentConfig` 设置模式、匹配/错配分数、线性或 affine Gap 分数及 `max_cells`，返回 `dnakit.alignment.AlignmentResult`。

### 2.1 全局比对

**作用：** 从头到尾对齐两条完整序列，返回加入 Gap 后的对应关系、得分、Identity 和 Coverage，适合比较长度接近且整体相关的序列。

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

```text
ACGT A-GT 2.0 0.75
```

### 2.2 局部比对

**作用：** 寻找两条序列中得分最高的局部对齐区域，返回该片段在两条原序列中的坐标、Identity 和 Coverage，适合发现共同片段。

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

```text
ACGT ACGT 1.0 0.5 0.5
```

### 2.3 半全局比对

**作用：** 在两端未对齐区域不扣 Gap 分，返回核心对齐及其原序列坐标，适合把引物、扩增子或短序列与长序列中的目标区域进行比较。

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

```text
ACGT ACGT 1.0 0.5 (2, 6)
```

**公共结果指标：**

- `aligned_query`、`aligned_target`：加入比对 Gap 后的序列。
- `score`：按匹配、错配和 Gap 参数得到的比对分数。
- `matches`、`mismatches`、`insertions`、`deletions`：各类比对列数量。
- `identity`：完全匹配列数占全部比对列的比例。
- `query_coverage`、`target_coverage`：参与比对的碱基分别占原序列的比例。
- `query_start/end`、`target_start/end`：比对区域在原序列中的范围。

**限制：** 当前只接受线性、无 Gap 输入；IUPAC 按字面字符计分；计算规模受 DP 单元上限约束。Identity 和 Coverage 会随比对模式及 Gap 参数改变。
