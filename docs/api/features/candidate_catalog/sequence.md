# 通用序列操作、k-mer 与编码

本类收录 **52** 项有文档的候选功能，其中 52 项有可在本地调用的正式入口（0 项因依赖二进制格式仅部分实现；0 项为受限的外部后端；0 项仅为依赖适配层）。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) apply_reverse_complement_augmentation

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `apply_reverse_complement_augmentation` 功能。候选编号 CAP-0583，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `apply_reverse_complement_augmentation` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.apply_reverse_complement_augmentation(sequence[required], include_original[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`include_original`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `apply_reverse_complement_augmentation` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import apply_reverse_complement_augmentation

result = apply_reverse_complement_augmentation(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0583",
  "operation": "apply_reverse_complement_augmentation",
  "result_type": "list",
  "result_preview": "['ATGGCC', 'GGCCAT']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) back_transcribe_rna_to_dna

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `back_transcribe_rna_to_dna` 功能。候选编号 CAP-0584，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `back_transcribe_rna_to_dna` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.back_transcribe_rna_to_dna(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `back_transcribe_rna_to_dna` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import back_transcribe_rna_to_dna

result = back_transcribe_rna_to_dna(
    sequence="AUGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0584",
  "operation": "back_transcribe_rna_to_dna",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) build_minimizer_index

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `build_minimizer_index` 功能。候选编号 CAP-0585，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `build_minimizer_index` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.build_minimizer_index(sequence[required], k[optional], window[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`k`、`window`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `build_minimizer_index` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import build_minimizer_index

result = build_minimizer_index(
    sequence="ATGGCC",
    k=3,
    window=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0585",
  "operation": "build_minimizer_index",
  "result_type": "dict",
  "result_preview": "{'ATG': [0], 'GCC': [3]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) calculate_sequence_metrics

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `calculate_sequence_metrics` 功能。候选编号 CAP-0586，来自 27 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_sequence_metrics` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `metric` 用于选择一种支持的变体（取值见源目录记录）。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.calculate_sequence_metrics(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `calculate_sequence_metrics` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import calculate_sequence_metrics

result = calculate_sequence_metrics(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0586",
  "operation": "calculate_sequence_metrics",
  "result_type": "dict",
  "result_preview": "{'length': 6, 'gc_fraction': 0.6666666666666666, 'gc_percent': 66.66666666666666, 'n_count': 0, 'ambiguous_count': 0, 'entropy_bits': 1.9182958340544896, 'composition': {'A': 1, 'C': 2, 'G': 2, 'T': 1}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) cluster_sequences_by_identity

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `cluster_sequences_by_identity` 功能。候选编号 CAP-0587，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `cluster_sequences_by_identity` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.cluster_sequences_by_identity(sequences[required], threshold[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequences`；可选参数：`threshold`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `cluster_sequences_by_identity` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import cluster_sequences_by_identity

result = cluster_sequences_by_identity(
    sequences=["ACGT", "ACGA"],
    threshold=0.5
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0587",
  "operation": "cluster_sequences_by_identity",
  "result_type": "list",
  "result_preview": "[[0, 1]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 6) compare_motif_alleles

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `compare_motif_alleles` 功能。候选编号 CAP-0588，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `compare_motif_alleles` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.compare_motif_alleles(left[required], right[required], motif[required], start[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`left`、`right`、`motif`；可选参数：`start`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `compare_motif_alleles` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import compare_motif_alleles

result = compare_motif_alleles(
    left="ACGT",
    right="AGGT",
    motif="GCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0588",
  "operation": "compare_motif_alleles",
  "result_type": "dict",
  "result_preview": "{'start': 0, 'end': 3, 'motif': 'GCC', 'left': 'ACG', 'right': 'AGG', 'left_matches': False, 'right_matches': False, 'changed': True}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 7) compare_sequence_sketches

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `compare_sequence_sketches` 功能。候选编号 CAP-0589，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `compare_sequence_sketches` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.compare_sequence_sketches(left[required], right[required], k[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`left`、`right`；可选参数：`k`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `compare_sequence_sketches` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import compare_sequence_sketches

result = compare_sequence_sketches(
    left="ACGT",
    right="AGGT",
    k=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0589",
  "operation": "compare_sequence_sketches",
  "result_type": "dict",
  "result_preview": "{'k': 3, 'left_size': 2, 'right_size': 2, 'intersection_size': 0, 'jaccard': 0.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 8) complement_dna

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `complement_dna` 功能。候选编号 CAP-0590，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `complement_dna` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`complement`。

- **API：** `dnakit.sequence.complement_dna(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `complement_dna` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import complement_dna

result = complement_dna(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0590",
  "operation": "complement_dna",
  "result_type": "str",
  "result_preview": "'TACCGG'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) count_canonical_kmers

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `count_canonical_kmers` 功能。候选编号 CAP-0591，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `count_canonical_kmers` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。已有同名函数证据：`canonical_kmer`。

- **API：** `dnakit.sequence.count_canonical_kmers(sequence[required], k[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`k`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `count_canonical_kmers` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import count_canonical_kmers

result = count_canonical_kmers(
    sequence="ATGGCC",
    k=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0591",
  "operation": "count_canonical_kmers",
  "result_type": "dict",
  "result_preview": "{'ATG': 1, 'CCA': 1, 'GCC': 2}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 10) count_kmers

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `count_kmers` 功能。候选编号 CAP-0592，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `count_kmers` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.count_kmers(sequence[required], k[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`k`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `count_kmers` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import count_kmers

result = count_kmers(
    sequence="ATGGCC",
    k=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0592",
  "operation": "count_kmers",
  "result_type": "dict",
  "result_preview": "{'ATG': 1, 'GCC': 1, 'GGC': 1, 'TGG': 1}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 11) count_kmers_batch

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `count_kmers_batch` 功能。候选编号 CAP-0593，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `count_kmers_batch` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.count_kmers_batch(sequences[required], k[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequences`；可选参数：`k`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `count_kmers_batch` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import count_kmers_batch

result = count_kmers_batch(
    sequences=["ACGT", "ACGA"],
    k=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0593",
  "operation": "count_kmers_batch",
  "result_type": "list",
  "result_preview": "[{'ACG': 1, 'CGT': 1}, {'ACG': 1, 'CGA': 1}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 12) deduplicate_sequences

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `deduplicate_sequences` 功能。候选编号 CAP-0594，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `deduplicate_sequences` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.deduplicate_sequences(sequences[required], canonical[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequences`；可选参数：`canonical`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `deduplicate_sequences` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import deduplicate_sequences

result = deduplicate_sequences(
    sequences=["ACGT", "ACGA"]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0594",
  "operation": "deduplicate_sequences",
  "result_type": "list",
  "result_preview": "['ACGT', 'ACGA']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 13) digest_sequence_in_silico

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `digest_sequence_in_silico` 功能。候选编号 CAP-0595，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `digest_sequence_in_silico` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.digest_sequence_in_silico(sequence[required], enzymes[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`、`enzymes`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `digest_sequence_in_silico` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import digest_sequence_in_silico

result = digest_sequence_in_silico(
    sequence="ATGGCC",
    enzymes={"EcoRI": "GAATTC"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0595",
  "operation": "digest_sequence_in_silico",
  "result_type": "dict",
  "result_preview": "{'sites': {'EcoRI': []}, 'cut_positions': [0, 6], 'fragment_lengths': [6]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) extract_cds_from_transcript

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `extract_cds_from_transcript` 功能。候选编号 CAP-0596，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_cds_from_transcript` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.extract_cds_from_transcript(transcript[required], intervals[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`transcript`、`intervals`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `extract_cds_from_transcript` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import extract_cds_from_transcript

result = extract_cds_from_transcript(
    transcript="ATGGCC",
    intervals=[(0, 6)]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0596",
  "operation": "extract_cds_from_transcript",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) generate_kmer_table

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `generate_kmer_table` 功能。候选编号 CAP-0598，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `generate_kmer_table` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.generate_kmer_table(sequences[required], k[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequences`；可选参数：`k`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `generate_kmer_table` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import generate_kmer_table

result = generate_kmer_table(
    sequences=["ACGT", "ACGA"],
    k=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0598",
  "operation": "generate_kmer_table",
  "result_type": "list",
  "result_preview": "[{'sequence_index': 0, 'kmer': 'ACG', 'count': 1}, {'sequence_index': 0, 'kmer': 'CGT', 'count': 1}, {'sequence_index': 1, 'kmer': 'ACG', 'count': 1}, {'sequence_index': 1, 'kmer': 'CGA', 'count': 1}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) generate_markov_sequence

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `generate_markov_sequence` 功能。候选编号 CAP-0599，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `generate_markov_sequence` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.generate_markov_sequence(sequence[required], length[required], order[optional], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`、`length`；可选参数：`order`、`seed`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `generate_markov_sequence` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import generate_markov_sequence

result = generate_markov_sequence(
    sequence="ATGGCC",
    length=12,
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0599",
  "operation": "generate_markov_sequence",
  "result_type": "str",
  "result_preview": "'ATGCCCCCCCCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 17) generate_random_dna

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `generate_random_dna` 功能。候选编号 CAP-0600，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `generate_random_dna` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.generate_random_dna(length[required], seed[optional], alphabet[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`length`；可选参数：`seed`、`alphabet`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `generate_random_dna` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import generate_random_dna

result = generate_random_dna(
    length=12,
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0600",
  "operation": "generate_random_dna",
  "result_type": "str",
  "result_preview": "'GGGGGGGGGGGG'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 18) mask_low_complexity

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `mask_low_complexity` 功能。候选编号 CAP-0601，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `mask_low_complexity` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.mask_low_complexity(sequence[required], window[optional], entropy_threshold[optional], mask[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`window`、`entropy_threshold`、`mask`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `mask_low_complexity` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import mask_low_complexity

result = mask_low_complexity(
    sequence="ATGGCC",
    window=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0601",
  "operation": "mask_low_complexity",
  "result_type": "str",
  "result_preview": "'ANNNNN'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 19) merge_overlapping_windows

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `merge_overlapping_windows` 功能。候选编号 CAP-0602，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `merge_overlapping_windows` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.merge_overlapping_windows(windows[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`windows`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `merge_overlapping_windows` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import merge_overlapping_windows

result = merge_overlapping_windows(
    windows=[(0, 3), (2, 5)]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0602",
  "operation": "merge_overlapping_windows",
  "result_type": "list",
  "result_preview": "[[0, 5]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 20) mutate_sequence_at_positions

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `mutate_sequence_at_positions` 功能。候选编号 CAP-0603，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `mutate_sequence_at_positions` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.mutate_sequence_at_positions(sequence[required], mutations[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`、`mutations`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `mutate_sequence_at_positions` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import mutate_sequence_at_positions

result = mutate_sequence_at_positions(
    sequence="ATGGCC",
    mutations={1: "T"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0603",
  "operation": "mutate_sequence_at_positions",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 21) normalize_dna_sequence

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `normalize_dna_sequence` 功能。候选编号 CAP-0604，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `normalize_dna_sequence` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.normalize_dna_sequence(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `normalize_dna_sequence` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import normalize_dna_sequence

result = normalize_dna_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0604",
  "operation": "normalize_dna_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 22) normalize_protein_sequence

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `normalize_protein_sequence` 功能。候选编号 CAP-0605，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `normalize_protein_sequence` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.normalize_protein_sequence(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `normalize_protein_sequence` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import normalize_protein_sequence

result = normalize_protein_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0605",
  "operation": "normalize_protein_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) normalize_rna_sequence

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `normalize_rna_sequence` 功能。候选编号 CAP-0606，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `normalize_rna_sequence` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.normalize_rna_sequence(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `normalize_rna_sequence` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import normalize_rna_sequence

result = normalize_rna_sequence(
    sequence="AUGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0606",
  "operation": "normalize_rna_sequence",
  "result_type": "str",
  "result_preview": "'AUGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) one_hot_encode_dna

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `one_hot_encode_dna` 功能。候选编号 CAP-0607，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `one_hot_encode_dna` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`one_hot_encode`。

- **API：** `dnakit.sequence.one_hot_encode_dna(sequence[required], include_ambiguous[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`include_ambiguous`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `one_hot_encode_dna` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import one_hot_encode_dna

result = one_hot_encode_dna(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0607",
  "operation": "one_hot_encode_dna",
  "result_type": "list",
  "result_preview": "[[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 25) one_hot_encode_rna

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `one_hot_encode_rna` 功能。候选编号 CAP-0608，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `one_hot_encode_rna` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`one_hot_encode`。

- **API：** `dnakit.sequence.one_hot_encode_rna(sequence[required], include_ambiguous[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`include_ambiguous`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `one_hot_encode_rna` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import one_hot_encode_rna

result = one_hot_encode_rna(
    sequence="AUGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0608",
  "operation": "one_hot_encode_rna",
  "result_type": "list",
  "result_preview": "[[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 26) pad_sequence

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `pad_sequence` 功能。候选编号 CAP-0609，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `pad_sequence` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`ena_sequence|plot_sequence|sequence`。

- **API：** `dnakit.sequence.pad_sequence(sequence[required], length[required], symbol[optional], side[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`、`length`；可选参数：`symbol`、`side`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `pad_sequence` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import pad_sequence

result = pad_sequence(
    sequence="ATGGCC",
    length=12
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0609",
  "operation": "pad_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCCNNNNNN'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 27) predict_coding_regions

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `predict_coding_regions` 功能。候选编号 CAP-0610，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `predict_coding_regions` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.predict_coding_regions(sequence[required], min_aa[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`min_aa`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `predict_coding_regions` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import predict_coding_regions

result = predict_coding_regions(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0610",
  "operation": "predict_coding_regions",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 28) remove_sequence_gaps

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `remove_sequence_gaps` 功能。候选编号 CAP-0611，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `remove_sequence_gaps` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.remove_sequence_gaps(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `remove_sequence_gaps` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import remove_sequence_gaps

result = remove_sequence_gaps(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0611",
  "operation": "remove_sequence_gaps",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 29) remove_sequence_whitespace

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `remove_sequence_whitespace` 功能。候选编号 CAP-0612，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `remove_sequence_whitespace` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.remove_sequence_whitespace(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `remove_sequence_whitespace` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import remove_sequence_whitespace

result = remove_sequence_whitespace(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0612",
  "operation": "remove_sequence_whitespace",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 30) replace_ambiguous_bases

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `replace_ambiguous_bases` 功能。候选编号 CAP-0613，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `replace_ambiguous_bases` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.replace_ambiguous_bases(sequence[required], replacement[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`replacement`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `replace_ambiguous_bases` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import replace_ambiguous_bases

result = replace_ambiguous_bases(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0613",
  "operation": "replace_ambiguous_bases",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 31) reverse_complement_dna

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `reverse_complement_dna` 功能。候选编号 CAP-0614，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `reverse_complement_dna` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`reverse_complement|reverse_complement_record|reverse_complement_search`。

- **API：** `dnakit.sequence.reverse_complement_dna(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `reverse_complement_dna` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import reverse_complement_dna

result = reverse_complement_dna(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0614",
  "operation": "reverse_complement_dna",
  "result_type": "str",
  "result_preview": "'GGCCAT'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 32) reverse_sequence

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `reverse_sequence` 功能。候选编号 CAP-0615，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `reverse_sequence` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.reverse_sequence(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `reverse_sequence` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import reverse_sequence

result = reverse_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0615",
  "operation": "reverse_sequence",
  "result_type": "str",
  "result_preview": "'CCGGTA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 33) scan_pwm_batch

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `scan_pwm_batch` 功能。候选编号 CAP-0617，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `scan_pwm_batch` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.scan_pwm_batch(sequences[required], pwm[required], threshold[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequences`、`pwm`；可选参数：`threshold`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `scan_pwm_batch` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import scan_pwm_batch

result = scan_pwm_batch(
    sequences=["ACGT", "ACGA"],
    pwm=[{"A": 1, "C": 0, "G": 0, "T": 0}],
    threshold=0.5
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0617",
  "operation": "scan_pwm_batch",
  "result_type": "list",
  "result_preview": "[[{'start': 0, 'end': 1, 'sequence': 'A', 'score': 1.0}], [{'start': 0, 'end': 1, 'sequence': 'A', 'score': 1.0}, {'start': 3, 'end': 4, 'sequence': 'A', 'score': 1.0}]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) score_motif_hit

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `score_motif_hit` 功能。候选编号 CAP-0618，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `score_motif_hit` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.score_motif_hit(sequence[required], motif[required], position[optional], match_score[optional], mismatch_score[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`、`motif`；可选参数：`position`、`match_score`、`mismatch_score`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `score_motif_hit` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import score_motif_hit

result = score_motif_hit(
    sequence="ATGGCC",
    motif="GCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0618",
  "operation": "score_motif_hit",
  "result_type": "dict",
  "result_preview": "{'position': 0, 'motif': 'GCC', 'observed': 'ATG', 'score': 0.0, 'matches': False}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 35) select_representative_sequences

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `select_representative_sequences` 功能。候选编号 CAP-0619，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `select_representative_sequences` 操作，并返回附带来源与诊断信息的结构化结果。无外部依赖的实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。它直接处理显式给出的 DNA 字符串，不依赖参考数据库、训练好的模型或外部可执行程序。审核依据：`no_exact_public_match`。已有同名函数证据：`select_representatives`。

- **API：** `dnakit.sequence.select_representative_sequences(sequences[required], threshold[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequences`；可选参数：`threshold`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `select_representative_sequences` 根据显式给出的 DNA 字符串返回确定性的 Python 结果。它不查询参考数据库、不加载训练好的模型，也不调用外部可执行程序。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import select_representative_sequences

result = select_representative_sequences(
    sequences=["ACGT", "ACGA"],
    threshold=0.5
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0619",
  "operation": "select_representative_sequences",
  "result_type": "list",
  "result_preview": "[0]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 36) simulate_indels

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `simulate_indels` 功能。候选编号 CAP-0620，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `simulate_indels` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.simulate_indels(sequence[required], count[optional], seed[optional], max_length[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`count`、`seed`、`max_length`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `simulate_indels` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import simulate_indels

result = simulate_indels(
    sequence="ATGGCC",
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0620",
  "operation": "simulate_indels",
  "result_type": "str",
  "result_preview": "'ATTGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 37) simulate_point_mutations

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `simulate_point_mutations` 功能。候选编号 CAP-0621，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `simulate_point_mutations` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.simulate_point_mutations(sequence[required], rate[optional], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`rate`、`seed`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `simulate_point_mutations` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import simulate_point_mutations

result = simulate_point_mutations(
    sequence="ATGGCC",
    rate=0.5,
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0621",
  "operation": "simulate_point_mutations",
  "result_type": "str",
  "result_preview": "'CATTTA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 38) simulate_recombination

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `simulate_recombination` 功能。候选编号 CAP-0622，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `simulate_recombination` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.simulate_recombination(left[required], right[required], crossover[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`left`、`right`；可选参数：`crossover`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `simulate_recombination` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import simulate_recombination

result = simulate_recombination(
    left="ACGT",
    right="AGGT"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0622",
  "operation": "simulate_recombination",
  "result_type": "str",
  "result_preview": "'ACGT'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 39) simulate_sequence_evolution

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `simulate_sequence_evolution` 功能。候选编号 CAP-0623，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `simulate_sequence_evolution` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.simulate_sequence_evolution(sequence[required], generations[optional], mutation_rate[optional], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`generations`、`mutation_rate`、`seed`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `simulate_sequence_evolution` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import simulate_sequence_evolution

result = simulate_sequence_evolution(
    sequence="ATGGCC",
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0623",
  "operation": "simulate_sequence_evolution",
  "result_type": "list",
  "result_preview": "['ATGGCC', 'ATGGCC']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 40) split_sequence_windows

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `split_sequence_windows` 功能。候选编号 CAP-0624，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_sequence_windows` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.split_sequence_windows(sequence[required], window[required], step[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`、`window`；可选参数：`step`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `split_sequence_windows` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import split_sequence_windows

result = split_sequence_windows(
    sequence="ATGGCC",
    window=3,
    step=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0624",
  "operation": "split_sequence_windows",
  "result_type": "list",
  "result_preview": "[{'start': 0, 'end': 3, 'sequence': 'ATG'}, {'start': 3, 'end': 6, 'sequence': 'GCC'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 41) tokenize_dna_bpe

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `tokenize_dna_bpe` 功能。候选编号 CAP-0625，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `tokenize_dna_bpe` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.tokenize_dna_bpe(sequence[required], merges[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`merges`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `tokenize_dna_bpe` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import tokenize_dna_bpe

result = tokenize_dna_bpe(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0625",
  "operation": "tokenize_dna_bpe",
  "result_type": "list",
  "result_preview": "['A', 'T', 'G', 'G', 'C', 'C']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 42) tokenize_dna_kmer

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `tokenize_dna_kmer` 功能。候选编号 CAP-0626，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `tokenize_dna_kmer` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.tokenize_dna_kmer(sequence[required], k[optional], stride[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`k`、`stride`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `tokenize_dna_kmer` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import tokenize_dna_kmer

result = tokenize_dna_kmer(
    sequence="ATGGCC",
    k=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0626",
  "operation": "tokenize_dna_kmer",
  "result_type": "list",
  "result_preview": "['ATG', 'TGG', 'GGC', 'GCC']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 43) transcribe_dna_to_rna

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `transcribe_dna_to_rna` 功能。候选编号 CAP-0627，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `transcribe_dna_to_rna` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.transcribe_dna_to_rna(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `transcribe_dna_to_rna` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import transcribe_dna_to_rna

result = transcribe_dna_to_rna(
    sequence="AUGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0627",
  "operation": "transcribe_dna_to_rna",
  "result_type": "str",
  "result_preview": "'AUGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 44) translate_codon

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `translate_codon` 功能。候选编号 CAP-0628，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `translate_codon` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.translate_codon(codon[required], unknown[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`codon`；可选参数：`unknown`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `translate_codon` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import translate_codon

result = translate_codon(
    codon="ATG"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0628",
  "operation": "translate_codon",
  "result_type": "str",
  "result_preview": "'M'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 45) translate_dna_to_protein

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `translate_dna_to_protein` 功能。候选编号 CAP-0629，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `translate_dna_to_protein` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.translate_dna_to_protein(sequence[required], frame[optional], stop_at_stop[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`frame`、`stop_at_stop`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `translate_dna_to_protein` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import translate_dna_to_protein

result = translate_dna_to_protein(
    sequence="ATGGCC",
    frame=0
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0629",
  "operation": "translate_dna_to_protein",
  "result_type": "str",
  "result_preview": "'MA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 46) translate_orf

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `translate_orf` 功能。候选编号 CAP-0630，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `translate_orf` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。已有同名函数证据：`translate`。

- **API：** `dnakit.sequence.translate_orf(sequence[required], start_codon[optional], stop_codons[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`start_codon`、`stop_codons`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `translate_orf` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import translate_orf

result = translate_orf(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0630",
  "operation": "translate_orf",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 47) translate_rna_to_protein

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `translate_rna_to_protein` 功能。候选编号 CAP-0631，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `translate_rna_to_protein` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.translate_rna_to_protein(sequence[required], frame[optional], stop_at_stop[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`frame`、`stop_at_stop`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `translate_rna_to_protein` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import translate_rna_to_protein

result = translate_rna_to_protein(
    sequence="AUGGCC",
    frame=0
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0631",
  "operation": "translate_rna_to_protein",
  "result_type": "str",
  "result_preview": "'MA'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 48) translate_six_frames

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `translate_six_frames` 功能。候选编号 CAP-0632，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `translate_six_frames` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.translate_six_frames(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `translate_six_frames` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import translate_six_frames

result = translate_six_frames(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0632",
  "operation": "translate_six_frames",
  "result_type": "dict",
  "result_preview": "{'+1': 'MA', '+2': 'W', '+3': 'G', '-1': 'GH', '-2': 'A', '-3': 'P'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 49) trim_sequence

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `trim_sequence` 功能。候选编号 CAP-0633，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `trim_sequence` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.trim_sequence(sequence[required], left[optional], right[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`left`、`right`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `trim_sequence` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import trim_sequence

result = trim_sequence(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0633",
  "operation": "trim_sequence",
  "result_type": "str",
  "result_preview": "'ATGGCC'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 50) validate_dna_alphabet

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `validate_dna_alphabet` 功能。候选编号 CAP-0634，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `validate_dna_alphabet` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.validate_dna_alphabet(sequence[required], allow_ambiguous[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`allow_ambiguous`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `validate_dna_alphabet` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import validate_dna_alphabet

result = validate_dna_alphabet(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0634",
  "operation": "validate_dna_alphabet",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'length': 6, 'invalid_symbols': [], 'alphabet': 'iupac'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 51) validate_protein_alphabet

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `validate_protein_alphabet` 功能。候选编号 CAP-0635，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `validate_protein_alphabet` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.validate_protein_alphabet(sequence[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：无。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `validate_protein_alphabet` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import validate_protein_alphabet

result = validate_protein_alphabet(
    sequence="ATGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0635",
  "operation": "validate_protein_alphabet",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'length': 6, 'invalid_symbols': [], 'alphabet': 'protein'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 52) validate_rna_alphabet

- **作用：** 为“通用序列操作、k-mer 与编码”类分析提供可复用的 `validate_rna_alphabet` 功能。候选编号 CAP-0636，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `validate_rna_alphabet` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/sequence/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.sequence.validate_rna_alphabet(sequence[required], allow_ambiguous[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sequence`；可选参数：`allow_ambiguous`。序列参数按操作接受 DNA、RNA 或蛋白质文本。

- **输出：** 本地函数 `validate_rna_alphabet` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequence import validate_rna_alphabet

result = validate_rna_alphabet(
    sequence="AUGGCC"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0636",
  "operation": "validate_rna_alphabet",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'length': 6, 'invalid_symbols': [], 'alphabet': 'iupac'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```
