# FASTA 与 FASTQ 记录处理

本类收录 **66** 项有文档的候选功能，其中 66 项有可在本地调用的正式入口（0 项因依赖二进制格式仅部分实现；0 项为受限的外部后端；0 项仅为依赖适配层）。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) build_fasta_index

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `build_fasta_index` 功能。候选编号 CAP-0176，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `build_fasta_index` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`generic_operation_name`。已有同名函数证据：`build_fasta_index|build_fastq_index|build_sketch_index|dnakit.io.build_fasta_index|load_fasta_index`。

- **API：** `dnakit.fasta_fastq.build_fasta_index(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `build_fasta_index` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import build_fasta_index

result = build_fasta_index(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0176",
  "operation": "build_fasta_index",
  "result_type": "dict",
  "result_preview": "{'a': {'record_index': 0, 'length': 4}, 'b': {'record_index': 1, 'length': 4}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) build_fastq_index

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `build_fastq_index` 功能。候选编号 CAP-0177，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `build_fastq_index` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`generic_operation_name`。已有同名函数证据：`build_fasta_index|build_fastq_index|build_sketch_index|dnakit.io.build_fastq_index|load_fastq_index`。

- **API：** `dnakit.fasta_fastq.build_fastq_index(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `build_fastq_index` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import build_fastq_index

result = build_fastq_index(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0177",
  "operation": "build_fastq_index",
  "result_type": "dict",
  "result_preview": "{'a': {'record_index': 0, 'length': 4, 'quality_length': 4}, 'b': {'record_index': 1, 'length': 4, 'quality_length': 4}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) calculate_gc_profile

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `calculate_gc_profile` 功能。候选编号 CAP-0178，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_gc_profile` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.calculate_gc_profile(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `calculate_gc_profile` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_gc_profile

result = calculate_gc_profile(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0178",
  "operation": "calculate_gc_profile",
  "result_type": "list",
  "result_preview": "[{'id': 'a', 'length': 4, 'gc_fraction': 0.5, 'gc_percent': 50.0}, {'id': 'b', 'length': 4, 'gc_fraction': 0.5, 'gc_percent': 50.0}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) calculate_length_distribution

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `calculate_length_distribution` 功能。候选编号 CAP-0179，来自 2 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_length_distribution` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta、fastq）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.calculate_length_distribution(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `calculate_length_distribution` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_length_distribution

result = calculate_length_distribution(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0179",
  "operation": "calculate_length_distribution",
  "result_type": "dict",
  "result_preview": "{'count': 2, 'min': 4, 'max': 4, 'mean': 4.0, 'counts': {4: 2}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) calculate_n_content_profile

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `calculate_n_content_profile` 功能。候选编号 CAP-0180，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_n_content_profile` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.calculate_n_content_profile(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `calculate_n_content_profile` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_n_content_profile

result = calculate_n_content_profile(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0180",
  "operation": "calculate_n_content_profile",
  "result_type": "list",
  "result_preview": "[{'id': 'a', 'n_count': 0, 'n_fraction': 0.0}, {'id': 'b', 'n_count': 0, 'n_fraction': 0.0}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 6) calculate_quality_distribution

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `calculate_quality_distribution` 功能。候选编号 CAP-0181，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_quality_distribution` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.calculate_quality_distribution(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `calculate_quality_distribution` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_quality_distribution

result = calculate_quality_distribution(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0181",
  "operation": "calculate_quality_distribution",
  "result_type": "dict",
  "result_preview": "{'count': 8, 'min': 40, 'max': 40, 'mean': 40.0, 'counts': {40: 8}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 7) compare_fasta_files

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `compare_fasta_files` 功能。候选编号 CAP-0182，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `compare_fasta_files` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.compare_fasta_files(left[required], right[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`left`、`right`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `compare_fasta_files` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import compare_fasta_files

result = compare_fasta_files(
    left='>a\nACGT\n>b\nGGAA\n',
    right='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0182",
  "operation": "compare_fasta_files",
  "result_type": "dict",
  "result_preview": "{'left_only': [], 'right_only': [], 'common': ['a', 'b'], 'different_sequences': [], 'identical': True}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 8) convert_from_fasta

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `convert_from_fasta` 功能。候选编号 CAP-0183，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `convert_from_fasta` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.convert_from_fasta(source[required], quality[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`quality`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `convert_from_fasta` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_from_fasta

result = convert_from_fasta(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0183",
  "operation": "convert_from_fasta",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) convert_from_tabular

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `convert_from_tabular` 功能。候选编号 CAP-0184，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `convert_from_tabular` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.convert_from_tabular(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `convert_from_tabular` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_from_tabular

result = convert_from_tabular(
    source='a\tACGT\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0184",
  "operation": "convert_from_tabular",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 10) convert_to_fasta

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `convert_to_fasta` 功能。候选编号 CAP-0185，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `convert_to_fasta` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.convert_to_fasta(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `convert_to_fasta` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_to_fasta

result = convert_to_fasta(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0185",
  "operation": "convert_to_fasta",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 11) convert_to_tabular

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `convert_to_tabular` 功能。候选编号 CAP-0186，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `convert_to_tabular` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.convert_to_tabular(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `convert_to_tabular` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_to_tabular

result = convert_to_tabular(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0186",
  "operation": "convert_to_tabular",
  "result_type": "list",
  "result_preview": "[{'id': 'a', 'sequence': 'ACGT'}, {'id': 'b', 'sequence': 'GGAA'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 12) count_reads

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `count_reads` 功能。候选编号 CAP-0187，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `count_reads` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.count_reads(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `count_reads` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import count_reads

result = count_reads(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0187",
  "operation": "count_reads",
  "result_type": "int",
  "result_preview": "2",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 13) count_records

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `count_records` 功能。候选编号 CAP-0188，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `count_records` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.count_records(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `count_records` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import count_records

result = count_records(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0188",
  "operation": "count_records",
  "result_type": "int",
  "result_preview": "2",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) deinterleave_pairs

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `deinterleave_pairs` 功能。候选编号 CAP-0189，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `deinterleave_pairs` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.deinterleave_pairs(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `deinterleave_pairs` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import deinterleave_pairs

result = deinterleave_pairs(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0189",
  "operation": "deinterleave_pairs",
  "result_type": "tuple",
  "result_preview": "([FastqRecord(identifier='a', sequence='ACGT', quality='IIII')], [FastqRecord(identifier='b', sequence='GGAA', quality='IIII')])",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) estimate_duplication

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `estimate_duplication` 功能。候选编号 CAP-0190，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `estimate_duplication` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.estimate_duplication(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `estimate_duplication` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import estimate_duplication

result = estimate_duplication(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0190",
  "operation": "estimate_duplication",
  "result_type": "dict",
  "result_preview": "{'records': 2, 'unique_sequences': 2, 'duplicate_records': 0, 'duplicate_fraction': 0.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) estimate_error_rate

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `estimate_error_rate` 功能。候选编号 CAP-0191，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `estimate_error_rate` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.estimate_error_rate(source[required], expected_quality[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`expected_quality`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `estimate_error_rate` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import estimate_error_rate

result = estimate_error_rate(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0191",
  "operation": "estimate_error_rate",
  "result_type": "float",
  "result_preview": "0.0",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 17) extract_cds

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_cds` 功能。候选编号 CAP-0192，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_cds` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_cds(source[required], identifier[required], intervals[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`identifier`、`intervals`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_cds` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_cds

result = extract_cds(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    intervals=[(0, 3)]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0192",
  "operation": "extract_cds",
  "result_type": "tuple",
  "result_preview": "('a:cds', 'ACG')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 18) extract_exon

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_exon` 功能。候选编号 CAP-0193，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_exon` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_exon(source[required], identifier[required], start[required], end[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`identifier`、`start`、`end`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_exon` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_exon

result = extract_exon(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    start=0,
    end=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0193",
  "operation": "extract_exon",
  "result_type": "tuple",
  "result_preview": "('a:0-3', 'ACG')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 19) extract_flanking

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_flanking` 功能。候选编号 CAP-0194，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_flanking` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_flanking(source[required], identifier[required], start[required], end[required], flank[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`identifier`、`start`、`end`；可选参数：`flank`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_flanking` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_flanking

result = extract_flanking(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    start=0,
    end=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0194",
  "operation": "extract_flanking",
  "result_type": "tuple",
  "result_preview": "('a:flank', 'ACGT')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 20) extract_introns

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_introns` 功能。候选编号 CAP-0195，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_introns` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_introns(source[required], identifier[required], exons[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`identifier`、`exons`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_introns` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_introns

result = extract_introns(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    exons=[(0, 2), (4, 6)]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0195",
  "operation": "extract_introns",
  "result_type": "list",
  "result_preview": "[('a:intron:1', 'GT')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 21) extract_promoter

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_promoter` 功能。候选编号 CAP-0196，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_promoter` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_promoter(source[required], identifier[required], tss[required], upstream[optional], downstream[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`identifier`、`tss`；可选参数：`upstream`、`downstream`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_promoter` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_promoter

result = extract_promoter(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    tss=2
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0196",
  "operation": "extract_promoter",
  "result_type": "tuple",
  "result_preview": "('a:promoter', 'ACGT')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 22) extract_read_names

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_read_names` 功能。候选编号 CAP-0197，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_read_names` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_read_names(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_read_names` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_read_names

result = extract_read_names(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0197",
  "operation": "extract_read_names",
  "result_type": "list",
  "result_preview": "['a', 'b']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) extract_region

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_region` 功能。候选编号 CAP-0198，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_region` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_region(source[required], identifier[required], start[required], end[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`identifier`、`start`、`end`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_region` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_region

result = extract_region(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    start=0,
    end=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0198",
  "operation": "extract_region",
  "result_type": "tuple",
  "result_preview": "('a:0-3', 'ACG')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) extract_regions_from_bed

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `extract_regions_from_bed` 功能。候选编号 CAP-0199，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `extract_regions_from_bed` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.extract_regions_from_bed(source[required], bed[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`bed`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `extract_regions_from_bed` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_regions_from_bed

result = extract_regions_from_bed(
    source='>a\nACGT\n>b\nGGAA\n',
    bed='a\t0\t3\tregion1\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0199",
  "operation": "extract_regions_from_bed",
  "result_type": "list",
  "result_preview": "[('region1', 'ACG')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 25) filter

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `filter` 功能。候选编号 CAP-0200，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.filter(source[required], min_length[optional], max_length[optional], ids[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`min_length`、`max_length`、`ids`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `filter` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import filter

result = filter(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0200",
  "operation": "filter",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 26) filter_by_length

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `filter_by_length` 功能。候选编号 CAP-0201，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_by_length` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.filter_by_length(source[required], min_length[optional], max_length[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`min_length`、`max_length`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `filter_by_length` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import filter_by_length

result = filter_by_length(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0201",
  "operation": "filter_by_length",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 27) filter_by_quality

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `filter_by_quality` 功能。候选编号 CAP-0202，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_by_quality` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.filter_by_quality(source[required], min_quality[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`min_quality`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `filter_by_quality` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import filter_by_quality

result = filter_by_quality(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    min_quality=20
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0202",
  "operation": "filter_by_quality",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 28) get_ids

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `get_ids` 功能。候选编号 CAP-0203，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `get_ids` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.get_ids(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `get_ids` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import get_ids

result = get_ids(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0203",
  "operation": "get_ids",
  "result_type": "list",
  "result_preview": "['a', 'b']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 29) get_lengths

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `get_lengths` 功能。候选编号 CAP-0204，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `get_lengths` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.get_lengths(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `get_lengths` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import get_lengths

result = get_lengths(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0204",
  "operation": "get_lengths",
  "result_type": "dict",
  "result_preview": "{'a': 4, 'b': 4}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 30) get_sequence_by_id

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `get_sequence_by_id` 功能。候选编号 CAP-0205，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `get_sequence_by_id` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。已有同名函数证据：`sequence_by_id`。

- **API：** `dnakit.fasta_fastq.get_sequence_by_id(source[required], identifier[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`identifier`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `get_sequence_by_id` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import get_sequence_by_id

result = get_sequence_by_id(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0205",
  "operation": "get_sequence_by_id",
  "result_type": "str",
  "result_preview": "'ACGT'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 31) interleave_pairs

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `interleave_pairs` 功能。候选编号 CAP-0206，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `interleave_pairs` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.interleave_pairs(first[required], second[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`first`、`second`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `interleave_pairs` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import interleave_pairs

result = interleave_pairs(
    first='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    second='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0206",
  "operation": "interleave_pairs",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 32) iterate

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `iterate` 功能。候选编号 CAP-0207，来自 2 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `iterate` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta、fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。已有同名函数证据：`literature`。

- **API：** `dnakit.fasta_fastq.iterate(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `iterate` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import iterate

result = iterate(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0207",
  "operation": "iterate",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 33) mask_regions

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `mask_regions` 功能。候选编号 CAP-0208，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `mask_regions` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.mask_regions(source[required], regions[required], symbol[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`regions`；可选参数：`symbol`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `mask_regions` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import mask_regions

result = mask_regions(
    source='>a\nACGT\n>b\nGGAA\n',
    regions={"a": [(0, 2)]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0208",
  "operation": "mask_regions",
  "result_type": "list",
  "result_preview": "[('a', 'NNGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) merge

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `merge` 功能。候选编号 CAP-0209，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `merge` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.merge(sources[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sources`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `merge` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import merge

result = merge(
    sources=['>a\nACGT\n>b\nGGAA\n']
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0209",
  "operation": "merge",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 35) merge_lanes

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `merge_lanes` 功能。候选编号 CAP-0210，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `merge_lanes` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.merge_lanes(sources[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`sources`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `merge_lanes` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import merge_lanes

result = merge_lanes(
    sources=['@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n']
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0210",
  "operation": "merge_lanes",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 36) query_fasta_index

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `query_fasta_index` 功能。候选编号 CAP-0211，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `query_fasta_index` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.query_fasta_index(index[required], identifier[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`index`、`identifier`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `query_fasta_index` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import query_fasta_index

result = query_fasta_index(
    index={"a": {"record_index": 0, "length": 4}},
    identifier="a"
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0211",
  "operation": "query_fasta_index",
  "result_type": "dict",
  "result_preview": "{'record_index': 0, 'length': 4}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 37) remove_duplicate_ids

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `remove_duplicate_ids` 功能。候选编号 CAP-0212，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `remove_duplicate_ids` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.remove_duplicate_ids(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `remove_duplicate_ids` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_duplicate_ids

result = remove_duplicate_ids(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0212",
  "operation": "remove_duplicate_ids",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 38) remove_duplicate_sequences

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `remove_duplicate_sequences` 功能。候选编号 CAP-0213，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `remove_duplicate_sequences` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.remove_duplicate_sequences(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `remove_duplicate_sequences` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_duplicate_sequences

result = remove_duplicate_sequences(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0213",
  "operation": "remove_duplicate_sequences",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 39) remove_poly_a

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `remove_poly_a` 功能。候选编号 CAP-0214，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `remove_poly_a` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.remove_poly_a(source[required], minimum[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`minimum`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `remove_poly_a` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_poly_a

result = remove_poly_a(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0214",
  "operation": "remove_poly_a",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 40) remove_poly_g

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `remove_poly_g` 功能。候选编号 CAP-0215，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `remove_poly_g` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.remove_poly_g(source[required], minimum[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`minimum`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `remove_poly_g` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_poly_g

result = remove_poly_g(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0215",
  "operation": "remove_poly_g",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 41) rename_headers

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `rename_headers` 功能。候选编号 CAP-0216，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `rename_headers` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.rename_headers(source[required], prefix[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`prefix`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `rename_headers` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import rename_headers

result = rename_headers(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0216",
  "operation": "rename_headers",
  "result_type": "list",
  "result_preview": "[('record_1', 'ACGT'), ('record_2', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 42) repair_fasta_headers

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `repair_fasta_headers` 功能。候选编号 CAP-0217，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `repair_fasta_headers` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.repair_fasta_headers(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `repair_fasta_headers` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import repair_fasta_headers

result = repair_fasta_headers(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0217",
  "operation": "repair_fasta_headers",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 43) repair_pairs

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `repair_pairs` 功能。候选编号 CAP-0218，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `repair_pairs` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.repair_pairs(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `repair_pairs` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import repair_pairs

result = repair_pairs(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0218",
  "operation": "repair_pairs",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 44) reverse_complement_records

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `reverse_complement_records` 功能。候选编号 CAP-0219，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `reverse_complement_records` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。已有同名函数证据：`reverse_complement|reverse_complement_record|reverse_complement_search`。

- **API：** `dnakit.fasta_fastq.reverse_complement_records(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `reverse_complement_records` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import reverse_complement_records

result = reverse_complement_records(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0219",
  "operation": "reverse_complement_records",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'TTCC')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 45) sample_records

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `sample_records` 功能。候选编号 CAP-0221，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `sample_records` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.sample_records(source[required], size[required], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`size`；可选参数：`seed`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `sample_records` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import sample_records

result = sample_records(
    source='>a\nACGT\n>b\nGGAA\n',
    size=1,
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0221",
  "operation": "sample_records",
  "result_type": "list",
  "result_preview": "[('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 46) sanitize_headers

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `sanitize_headers` 功能。候选编号 CAP-0222，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `sanitize_headers` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.sanitize_headers(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `sanitize_headers` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import sanitize_headers

result = sanitize_headers(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0222",
  "operation": "sanitize_headers",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 47) shuffle_reads

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `shuffle_reads` 功能。候选编号 CAP-0223，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `shuffle_reads` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.shuffle_reads(source[required], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`seed`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `shuffle_reads` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import shuffle_reads

result = shuffle_reads(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0223",
  "operation": "shuffle_reads",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 48) shuffle_records

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `shuffle_records` 功能。候选编号 CAP-0224，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `shuffle_records` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.shuffle_records(source[required], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`seed`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `shuffle_records` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import shuffle_records

result = shuffle_records(
    source='>a\nACGT\n>b\nGGAA\n',
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0224",
  "operation": "shuffle_records",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 49) sort

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `sort` 功能。候选编号 CAP-0225，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `sort` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.sort(source[required], by[optional], reverse[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`by`、`reverse`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `sort` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import sort

result = sort(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0225",
  "operation": "sort",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 50) split_by_barcode

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `split_by_barcode` 功能。候选编号 CAP-0226，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_by_barcode` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.split_by_barcode(source[required], barcodes[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`barcodes`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `split_by_barcode` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_barcode

result = split_by_barcode(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    barcodes={"sample1": "AC"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0226",
  "operation": "split_by_barcode",
  "result_type": "dict",
  "result_preview": "{'sample1': [FastqRecord(identifier='a', sequence='ACGT', quality='IIII')], 'unmatched': [FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 51) split_by_chromosome

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `split_by_chromosome` 功能。候选编号 CAP-0227，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_by_chromosome` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.split_by_chromosome(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `split_by_chromosome` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_chromosome

result = split_by_chromosome(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0227",
  "operation": "split_by_chromosome",
  "result_type": "dict",
  "result_preview": "{'a': [('a', 'ACGT')], 'b': [('b', 'GGAA')]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 52) split_by_length

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `split_by_length` 功能。候选编号 CAP-0228，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_by_length` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.split_by_length(source[required], threshold[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`threshold`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `split_by_length` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_length

result = split_by_length(
    source='>a\nACGT\n>b\nGGAA\n',
    threshold=4
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0228",
  "operation": "split_by_length",
  "result_type": "dict",
  "result_preview": "{'short': [], 'long': [('a', 'ACGT'), ('b', 'GGAA')]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 53) split_by_record

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `split_by_record` 功能。候选编号 CAP-0229，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_by_record` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.split_by_record(source[required], records_per_part[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`records_per_part`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `split_by_record` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_record

result = split_by_record(
    source='>a\nACGT\n>b\nGGAA\n',
    records_per_part=1
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0229",
  "operation": "split_by_record",
  "result_type": "list",
  "result_preview": "[[('a', 'ACGT')], [('b', 'GGAA')]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 54) split_by_size

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `split_by_size` 功能。候选编号 CAP-0230，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_by_size` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.split_by_size(source[required], bases_per_part[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`bases_per_part`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `split_by_size` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_size

result = split_by_size(
    source='>a\nACGT\n>b\nGGAA\n',
    bases_per_part=4
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0230",
  "operation": "split_by_size",
  "result_type": "list",
  "result_preview": "[[('a', 'ACGT')], [('b', 'GGAA')]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 55) split_by_window

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `split_by_window` 功能。候选编号 CAP-0231，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_by_window` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.split_by_window(source[required], window[required], step[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`window`；可选参数：`step`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `split_by_window` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_window

result = split_by_window(
    source='>a\nACGT\n>b\nGGAA\n',
    window=3,
    step=2
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0231",
  "operation": "split_by_window",
  "result_type": "list",
  "result_preview": "[('a:0-3', 'ACG'), ('a:2-4', 'GT'), ('b:0-3', 'GGA'), ('b:2-4', 'AA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 56) stats

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `stats` 功能。候选编号 CAP-0232，来自 2 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `stats` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta、fastq）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.stats(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `stats` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import stats

result = stats(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0232",
  "operation": "stats",
  "result_type": "dict",
  "result_preview": "{'records': 2, 'bases': 8, 'lengths': {'count': 2, 'min': 4, 'max': 4, 'mean': 4.0, 'counts': {4: 2}}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 57) subsample

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `subsample` 功能。候选编号 CAP-0233，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `subsample` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。已有同名函数证据：`sample`。

- **API：** `dnakit.fasta_fastq.subsample(source[required], fraction[required], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`fraction`；可选参数：`seed`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `subsample` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import subsample

result = subsample(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    fraction=0.5,
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0233",
  "operation": "subsample",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 58) subsample_records

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `subsample_records` 功能。候选编号 CAP-0234，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `subsample_records` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.subsample_records(source[required], fraction[required], seed[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`fraction`；可选参数：`seed`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `subsample_records` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import subsample_records

result = subsample_records(
    source='>a\nACGT\n>b\nGGAA\n',
    fraction=0.5,
    seed=7
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0234",
  "operation": "subsample_records",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 59) synchronize_pairs

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `synchronize_pairs` 功能。候选编号 CAP-0235，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `synchronize_pairs` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.synchronize_pairs(first[required], second[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`first`、`second`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `synchronize_pairs` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import synchronize_pairs

result = synchronize_pairs(
    first='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    second='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0235",
  "operation": "synchronize_pairs",
  "result_type": "tuple",
  "result_preview": "([FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')], [FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')])",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 60) trim_adapters

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `trim_adapters` 功能。候选编号 CAP-0236，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `trim_adapters` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.trim_adapters(source[required], adapters[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`adapters`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `trim_adapters` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import trim_adapters

result = trim_adapters(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    adapters=["AAA"]
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0236",
  "operation": "trim_adapters",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 61) trim_low_quality_ends

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `trim_low_quality_ends` 功能。候选编号 CAP-0237，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `trim_low_quality_ends` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.trim_low_quality_ends(source[required], minimum_quality[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`minimum_quality`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `trim_low_quality_ends` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import trim_low_quality_ends

result = trim_low_quality_ends(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    minimum_quality=20
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0237",
  "operation": "trim_low_quality_ends",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 62) unmask_regions

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `unmask_regions` 功能。候选编号 CAP-0238，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `unmask_regions` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.unmask_regions(source[required], regions[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`、`regions`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `unmask_regions` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import unmask_regions

result = unmask_regions(
    source='>a\nACGT\n>b\nGGAA\n',
    regions={"a": [(0, 2)]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0238",
  "operation": "unmask_regions",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 63) unwrap_fasta_lines

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `unwrap_fasta_lines` 功能。候选编号 CAP-0239，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `unwrap_fasta_lines` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.unwrap_fasta_lines(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `unwrap_fasta_lines` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import unwrap_fasta_lines

result = unwrap_fasta_lines(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0239",
  "operation": "unwrap_fasta_lines",
  "result_type": "str",
  "result_preview": "'>a\\nACGT\\n>b\\nGGAA\\n'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 64) validate

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `validate` 功能。候选编号 CAP-0240，来自 2 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `validate` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta、fastq）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`generic_operation_name`。已有同名函数证据：`dnakit.standardize.validate|validate|validate_set`。

- **API：** `dnakit.fasta_fastq.validate(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `validate` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import validate

result = validate(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0240",
  "operation": "validate",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'format': 'fasta', 'records': 2, 'duplicate_ids': []}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 65) validate_sequence_ids

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `validate_sequence_ids` 功能。候选编号 CAP-0241，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `validate_sequence_ids` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。无外部依赖的实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。它处理显式给出的 FASTA/FASTQ 记录或文本，不调用比对工具、参考基因组或外部数据库。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.validate_sequence_ids(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `validate_sequence_ids` 根据显式给出的 FASTA/FASTQ 输入返回 Python 原生值（记录、映射或摘要列表）。它不推断注释，也不调用外部后端。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import validate_sequence_ids

result = validate_sequence_ids(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0241",
  "operation": "validate_sequence_ids",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'ids': ['a', 'b'], 'duplicates': []}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 66) wrap_fasta_lines

- **作用：** 为“FASTA 与 FASTQ 记录处理”类分析提供可复用的 `wrap_fasta_lines` 功能。候选编号 CAP-0242，来自 1 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `wrap_fasta_lines` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `format` 用于选择一种支持的变体（fasta）。实现位于 `src/dnakit/fasta_fastq/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.fasta_fastq.wrap_fasta_lines(source[required], width[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：`width`。FASTA/FASTQ 参数接受文本、路径或已读入的记录。

- **输出：** 本地函数 `wrap_fasta_lines` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。

- **示例代码：**

```python
from pathlib import Path

from dnakit.fasta_fastq import wrap_fasta_lines

result = wrap_fasta_lines(
    source='>a\nACGT\n>b\nGGAA\n',
    width=4
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0242",
  "operation": "wrap_fasta_lines",
  "result_type": "str",
  "result_preview": "'>a\\nACGT\\n>b\\nGGAA\\n'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```
