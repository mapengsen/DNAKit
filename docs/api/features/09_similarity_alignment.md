# 相似度计算

计算 DNA 序列之间的 k-mer、指纹和 sketch 相似度，并生成批量相似度矩阵及参考库相似度结果。

## 1) k-mer相似度 {#1-sim-010-k-mer}

- **作用：** 把序列转换为 k-mer 集合或计数向量，计算 Jaccard、Containment 或 Cosine 相似度，用于量化局部组成的重合程度；也可比较两组序列的平均相似度。
- **计算方法：** 单条对单条时，对 k-mer 集合计算交并比/包含率，或对计数向量计算余弦相似度。输入多条序列时，逐条比较左侧每条序列与右侧所有序列，`value` 为全部配对分数的算术平均：`sum(sim(query, reference)) / (查询数 × 参考数)`。这不是每条查询的最近邻分数。

- **API：** `dnakit.similarity.kmer_similarity(left[必须], right[必须], k[必须], metric[可选], mode[可选], canonical[可选], overlapping[可选], zero_vector_policy[可选], max_pairwise_comparisons[可选])`、`dnakit.similarity.kmer_vector_similarity(left[必须], right[必须], metric[可选], mode[可选], zero_vector_policy[可选])`。

- **输入：** `left` 为查询序列，`right` 为参考序列；两侧均可用单条 `DNASequence`/`DNARecord`、`DNA([...])`、`DNASet` 或由 `DNASequence`/`DNARecord` 组成的列表。集合不能为空。必填 `k`；可选 metric、set/count、canonical、重叠及零向量策略。默认最多计算 1,000,000 对，可用 `max_pairwise_comparisons` 调整。
- **输出：** 返回 `SimilarityResult`。单条对单条时保留原有 `components["shared_weight"]` 等字段；多序列时 `value` 为整体平均，`components["query_means"]` 按查询顺序给出各自与全部参考序列的平均分，并包含 `query_count`、`reference_count`、`pair_count`。

- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import kmer_similarity

result = kmer_similarity(DNASequence("AAA"), DNASequence("AAC"), k=1)
print(result.value, result.components["shared_weight"])
```

- **示例结果：**

```text
0.5 1.0
```

- **多序列示例：**

```python
from dnakit import DNA
from dnakit.similarity import kmer_similarity

queries = DNA(["AAAA", "GGGG"])
references = DNA(["AAAA", "AAAT"])
result = kmer_similarity(queries, references, k=2)
print(result.value, result.components["query_means"])
```

```text
0.375 (0.75, 0.0)
```

## 2) 指纹相似度 {#2-sim-011}

- **作用：** 比较两个 schema 一致的 DNA 指纹或数值向量，计算 Tanimoto、Jaccard、Cosine 或数值距离，返回可用于排序、聚类和阈值判断的分数。
- **计算方法：** 先校验两侧 schema 相同，再按所选位集合或数值向量公式比较。

- **API：** `dnakit.similarity.fingerprint_similarity(left[必须], right[必须], metric[可选], weights[可选], zero_vector_policy[可选])`。

- **输入：** 必填两个同构指纹/向量；可选 metric、特征权重和零向量策略。
- **输出：** 返回 `SimilarityResult`，包含相似度数值 `value` 和比较方法 `method`。

- **示例代码：**

```python
from dnakit.similarity import fingerprint_similarity

result = fingerprint_similarity((1.0, 1.0), (1.0, 0.0), metric="tanimoto")
print(result.value)
```

- **示例结果：**

```text
0.5
```

## 3) Sketch相似度 {#3-sim-012-sketch}

- **作用：** 比较两个 MinHash/FracMinHash Sketch 中保留的哈希值，近似估计 Jaccard 或 Containment，并报告兼容性参数，适合大序列快速比较。
- **计算方法：** 先检查 k、seed 和抽样参数兼容，再用共享哈希估计集合交并比或包含率。

- **API：** `dnakit.fingerprints.minhash(value[必须], k[必须], num_hashes[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`、`dnakit.fingerprints.fracminhash(value[必须], k[必须], scaled[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`，以及 `dnakit.similarity.sketch_similarity(left[必须], right[必须], metric[可选], min_shared_hashes[可选])`。

- **输入：** 必填两个兼容 sketch；可选 metric 和最少共享 hash 数。
- **输出：** `minhash`/`fracminhash` 返回 `SketchResult`；`sketch_similarity` 返回 `SketchSimilarityResult`，包含相似度 `value` 和共享哈希数。

- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import minhash
from dnakit.similarity import sketch_similarity

left = minhash(DNASequence("ACGTAC"), k=2, num_hashes=100)
right = minhash(DNASequence("ACGTTC"), k=2, num_hashes=100)
result = sketch_similarity(left, right, min_shared_hashes=1)
print(result.value, result.shared_hash_count)
```

- **示例结果：**

```text
0.4 2
```

## 4) Dashing sketch相似度 {#4-sim-013-dashing-sketch}

- **作用：** 调用外部 Dashing 为多条序列构建高性能 Sketch，计算近似 Jaccard 矩阵或 Top-k 邻居，适合原生穷举计算难以承受的数据规模。
- **计算方法：** 调用显式配置的 Dashing 可执行程序生成 sketch，并解析矩阵或 Top-k 结果。

- **API：** `dnakit.similarity.DashingAdapter(executable_path[必须])`、`dnakit.similarity.DashingAdapter.matrix(inputs[必须], k[可选], mode[可选], sketch_size_log2[可选], canonical[可选], threads[可选], temp_dir[可选], output_path[可选], overwrite[可选], timeout_seconds[可选], max_items[可选], max_input_bytes[可选], max_output_bytes[可选], max_capture_bytes[可选], max_sketch_memory_bytes[可选])`、`dnakit.similarity.DashingAdapter.top_k(inputs[必须], top_k[必须], k[可选], mode[可选], sketch_size_log2[可选], canonical[可选], threads[可选], temp_dir[可选], output_path[可选], overwrite[可选], timeout_seconds[可选], max_items[可选], max_input_bytes[可选], max_output_bytes[可选], max_capture_bytes[可选], max_sketch_memory_bytes[可选])`。

- **输入：** 必填调用方显式配置的可执行文件和序列/FASTA/FASTQ；可选 k、模式、sketch 大小、canonical、线程、临时目录和输出路径。
- **输出：** 返回 Dashing sketch 比较的相似度结果。

- **示例代码：**

```python
import os
from pathlib import Path

from dnakit import DNASequence
from dnakit.similarity import DashingAdapter

configured = os.environ.get("DNAKIT_DASHING_EXECUTABLE")
if configured is None:
    print("跳过：请先显式配置 DNAKIT_DASHING_EXECUTABLE")
else:
    result = DashingAdapter(Path(configured)).matrix(
        (DNASequence("AACCGG"), DNASequence("AACCTT")), k=2, mode="exact"
    )
    print(result.values)
```

- **示例结果：**

```text
跳过：请先显式配置 DNAKIT_DASHING_EXECUTABLE
```

## 5) Reference similarity {#6-eval-011-reference-similarity}

- **作用：** 为每条查询序列在版本化参考库中找到最相似记录，返回参考 ID、相似度、覆盖度和排名，用于参考归属及新颖性分析。
- **计算方法：** 按指定相似度方法逐条检索参考库，并按分数与覆盖度筛选和排序。
- **API：** `dnakit.evaluation.evaluate_reference_similarity(queries[必须], reference[必须], config[可选])`；`config` 使用 `dnakit.evaluation.ReferenceSearchConfig`。
- **输入：** query 和 `ReferenceLibrary`；可选相似度与 coverage 配置。
- **输出：** 返回 `EvaluationReport`，整体指标在 `metrics`，逐条结果在 `entries`。
- **示例代码：** 先按[参考库示例准备](12_evaluation.md#eval-008-novelty)构造 `queries` 和 `reference`，再运行下面代码。

```python
from dnakit.evaluation import evaluate_reference_similarity

report = evaluate_reference_similarity(queries, reference)
print(report.metrics["mean_nearest_similarity"])
```

- **示例结果：**

```text
0.5
```
