# 相似度计算

计算 DNA 序列之间的 k-mer、指纹和 sketch 相似度，并生成批量相似度矩阵及参考库相似度结果。

## 1) SIM-010 · k-mer相似度

**作用：** 把两条序列转换为 k-mer 集合或计数向量，计算 Jaccard、Containment 或 Cosine 相似度，用于量化局部组成的重合程度。

**API：** `dnakit.similarity.kmer_similarity(left[必须], right[必须], k[必须], metric[可选], mode[可选], canonical[可选], overlapping[可选], zero_vector_policy[可选])`、`dnakit.similarity.kmer_vector_similarity(left[必须], right[必须], metric[可选], mode[可选], zero_vector_policy[可选])`。

**输入：** 必填两条序列和 `k`；可选 metric、set/count、canonical、重叠及零向量策略。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.similarity import kmer_similarity

result = kmer_similarity(DNASequence("AAA"), DNASequence("AAC"), k=1)
print(result.value, result.components["shared_weight"])
```

**示例结果：**

```text
0.5 1.0
```

## 2) SIM-011 · 指纹相似度

**作用：** 比较两个 schema 一致的 DNA 指纹或数值向量，计算 Tanimoto、Jaccard、Cosine 或数值距离，返回可用于排序、聚类和阈值判断的分数。

**API：** `dnakit.similarity.fingerprint_similarity(left[必须], right[必须], metric[可选], weights[可选], zero_vector_policy[可选])`。

**输入：** 必填两个同构指纹/向量；可选 metric、特征权重和零向量策略。

**示例代码：**

```python
from dnakit.similarity import fingerprint_similarity

result = fingerprint_similarity((1.0, 1.0), (1.0, 0.0), metric="tanimoto")
print(result.value)
```

**示例结果：**

```text
0.5
```

## 3) SIM-012 · Sketch相似度

**作用：** 比较两个 MinHash/FracMinHash Sketch 中保留的哈希值，近似估计 Jaccard 或 Containment，并报告兼容性参数，适合大序列快速比较。

**API：** `dnakit.fingerprints.minhash(value[必须], k[必须], num_hashes[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`、`dnakit.fingerprints.fracminhash(value[必须], k[必须], scaled[可选], canonical[可选], seed[可选], max_hashes[可选], max_unique_hashes[可选])`，以及 `dnakit.similarity.sketch_similarity(left[必须], right[必须], metric[可选], min_shared_hashes[可选])`。

**输入：** 必填两个兼容 sketch；可选 metric 和最少共享 hash 数。

**示例代码：**

```python
from dnakit import DNASequence
from dnakit.fingerprints import minhash
from dnakit.similarity import sketch_similarity

left = minhash(DNASequence("ACGTAC"), k=2, num_hashes=100)
right = minhash(DNASequence("ACGTTC"), k=2, num_hashes=100)
result = sketch_similarity(left, right, min_shared_hashes=1)
print(result.value, result.shared_hash_count)
```

**示例结果：**

```text
0.4 2
```

## 4) SIM-013 · Dashing sketch相似度

**作用：** 调用外部 Dashing 为多条序列构建高性能 Sketch，计算近似 Jaccard 矩阵或 Top-k 邻居，适合原生穷举计算难以承受的数据规模。

**API：** `dnakit.similarity.DashingAdapter(executable_path[必须])`、`dnakit.similarity.DashingAdapter.matrix(inputs[必须], k[可选], mode[可选], sketch_size_log2[可选], canonical[可选], threads[可选], temp_dir[可选], output_path[可选], overwrite[可选], timeout_seconds[可选], max_items[可选], max_input_bytes[可选], max_output_bytes[可选], max_capture_bytes[可选], max_sketch_memory_bytes[可选])`、`dnakit.similarity.DashingAdapter.top_k(inputs[必须], top_k[必须], k[可选], mode[可选], sketch_size_log2[可选], canonical[可选], threads[可选], temp_dir[可选], output_path[可选], overwrite[可选], timeout_seconds[可选], max_items[可选], max_input_bytes[可选], max_output_bytes[可选], max_capture_bytes[可选], max_sketch_memory_bytes[可选])`。

**输入：** 必填调用方显式配置的可执行文件和序列/FASTA/FASTQ；可选 k、模式、sketch 大小、canonical、线程、临时目录和输出路径。

**示例代码：**

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

**示例结果：**

```text
跳过：请先显式配置 DNAKIT_DASHING_EXECUTABLE
```

## 6) `EVAL-011` Reference similarity

- **作用：** 为每条查询序列在版本化参考库中找到最相似记录，返回参考 ID、相似度、覆盖度和排名，用于参考归属及新颖性分析。
- **API**：`dnakit.evaluation.evaluate_reference_similarity(queries[必须], reference[必须], config[可选])`；`config` 使用 `dnakit.evaluation.ReferenceSearchConfig`。
- **输入**：query 和 `ReferenceLibrary`；可选相似度与 coverage 配置。
- **示例代码**：先按[参考库示例准备](12_evaluation.md#eval-008-novelty)构造 `queries` 和 `reference`，再运行下面代码。

```python
from dnakit.evaluation import evaluate_reference_similarity

report = evaluate_reference_similarity(queries, reference)
print(report.metrics["mean_nearest_similarity"])
```

- **示例结果：**

```text
0.5
```
