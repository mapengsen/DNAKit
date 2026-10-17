# 读段比对与 BAM 统计

本类收录 **34** 项有文档的功能，均已在本地实现（原先唯一的占位函数 `calculate_clipping_profile` 现在由 CIGAR 字符串计算结果）。其中 5 项为部分实现，因为二进制 BAM/CRAM 输入需要可选包 `pysam`，其余处理 SAM 文本或记录映射。示例输出均为真实运行结果；`input_summary` 简写为 `{...}`。

## 1) align_long_reads

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `align_long_reads` 功能。候选编号 CAP-0001，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `align_long_reads` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.align_long_reads(input[required], backend[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`backend`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `align_long_reads` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import align_long_reads

result = align_long_reads(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0001",
  "operation": "align_long_reads",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) align_paired_reads

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `align_paired_reads` 功能。候选编号 CAP-0002，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `align_paired_reads` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.align_paired_reads(input[required], backend[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`backend`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `align_paired_reads` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import align_paired_reads

result = align_paired_reads(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0002",
  "operation": "align_paired_reads",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) align_reads

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `align_reads` 功能。候选编号 CAP-0003，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `align_reads` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.align_reads(input[required], backend[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`backend`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `align_reads` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import align_reads

result = align_reads(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0003",
  "operation": "align_reads",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) align_to_graph

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `align_to_graph` 功能。候选编号 CAP-0004，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `align_to_graph` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.align_to_graph(input[required], backend[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`backend`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `align_to_graph` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import align_to_graph

result = align_to_graph(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0004",
  "operation": "align_to_graph",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) bam_to_cram

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `bam_to_cram` 功能。候选编号 CAP-0005，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `bam_to_cram` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/alignment/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。此条目仅部分实现：本地示例处理 SAM 文本；二进制 BAM/CRAM 需要经过验证的 pysam 或同等后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.bam_to_cram(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `bam_to_cram` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。二进制 BAM/CRAM 的解析或转换尚未验证，仍待完成。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import bam_to_cram

result = bam_to_cram(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0005",
  "operation": "bam_to_cram",
  "result_type": "bytes",
  "result_preview": "<bytes length=34>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 6) bam_to_sam

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `bam_to_sam` 功能。候选编号 CAP-0006，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `bam_to_sam` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/alignment/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。此条目仅部分实现：本地示例处理 SAM 文本；二进制 BAM/CRAM 需要经过验证的 pysam 或同等后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.bam_to_sam(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `bam_to_sam` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。二进制 BAM/CRAM 的解析或转换尚未验证，仍待完成。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import bam_to_sam

result = bam_to_sam(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0006",
  "operation": "bam_to_sam",
  "result_type": "str",
  "result_preview": "'r1\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 7) build_index

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `build_index` 功能。候选编号 CAP-0007，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `build_index` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`generic_operation_name`。已有同名函数证据：`build_index|dnakit.download.build_index`。

- **API：** `dnakit.alignment.build_index(input[required], k[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`k`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `build_index` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import build_index

result = build_index(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"},
    k=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0007",
  "operation": "build_index",
  "result_type": "dict",
  "result_preview": "{'ACG': [0, 4], 'CGT': [1, 5], 'GTA': [2], 'TAC': [3]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 8) calculate_base_depth

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_base_depth` 功能。候选编号 CAP-0008，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_base_depth` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_base_depth(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_base_depth` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_base_depth

result = calculate_base_depth(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0008",
  "operation": "calculate_base_depth",
  "result_type": "dict",
  "result_preview": "{'chr1': {0: 1, 1: 1, 2: 1, 3: 1}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) calculate_clipping_profile

- **作用：** 汇总 SAM 文本或类 SAM 记录中已比对 reads 左右两端的软剪切和硬剪切。

- **处理方法：** 剪切信息取自 CIGAR 字符串：左端和右端连续的 `S`（软）与 `H`（硬）操作，其中左、右指参考序列方向（而不是 read 的 5'/3' 端）。未比对的 reads（flag 为 4 或参考为 `*`）被排除。`clipped_reads` 统计含任何 `S` 或 `H` 操作的已比对 reads，`clipped_fraction` 为其占已比对 reads 的比例。`clip_length_histogram` 把一条 read 的剪切碱基总数（软加硬、两端合计）映射到对应的已比对 reads 数，包含 0。只使用 CIGAR，不检查 read 序列。

- **API：** `dnakit.alignment.calculate_clipping_profile(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：SAM 文本或 SAM 文件路径（忽略 `@` 头部行）、记录映射列表（`query_name`、`flag`、`reference`、`position`、`cigar` 等），或带 `sam` 或 `records` 的映射。无法读取的记录会抛 `InputFormatError`（`INVALID_ALIGNMENT_INPUT`）。`options`：`per_read`（默认 true：包含逐 read 表，最多 10,000 行）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `reads`、`mapped_reads`、`clipped_reads`、`clipped_fraction`、碱基总数 `left_soft`、`right_soft`、`left_hard`、`right_hard`、`soft_clipped_bases`、`hard_clipped_bases`、`clip_length_histogram`（JSON 对象的键是字符串形式的剪切碱基数），以及（除非 `per_read` 为 false）`per_read`（`query_name` 与每条 read 的四个剪切长度）。示例中 `r3` 未比对，`r4` 没有剪切。

- **示例代码：**

```python
import json

from dnakit.alignment import calculate_clipping_profile

sam = (
    "@HD\tVN:1.6\n"
    "r1\t0\tchr1\t10\t60\t3S10M2S\t*\t0\t0\tACGTACGTACGTACG\tIIIIIIIIIIIIIII\n"
    "r2\t16\tchr1\t20\t60\t5H10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
    "r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
    "r4\t0\tchr1\t30\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
)
result = calculate_clipping_profile(sam)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "alignment",
  "operation": "calculate_clipping_profile",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "mapped_reads": 3,
    "clipped_reads": 2,
    "clipped_fraction": 0.6666666666666666,
    "left_soft": 3,
    "right_soft": 2,
    "left_hard": 5,
    "right_hard": 0,
    "soft_clipped_bases": 5,
    "hard_clipped_bases": 5,
    "clip_length_histogram": {
      "0": 1,
      "5": 2
    },
    "per_read": [
      {
        "query_name": "r1",
        "left_soft": 3,
        "right_soft": 2,
        "left_hard": 0,
        "right_hard": 0
      },
      {
        "query_name": "r2",
        "left_soft": 0,
        "right_soft": 0,
        "left_hard": 5,
        "right_hard": 0
      },
      {
        "query_name": "r4",
        "left_soft": 0,
        "right_soft": 0,
        "left_hard": 0,
        "right_hard": 0
      }
    ]
  }
}
```

## 10) calculate_coverage

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_coverage` 功能。候选编号 CAP-0010，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_coverage` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_coverage(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_coverage` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_coverage

result = calculate_coverage(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0010",
  "operation": "calculate_coverage",
  "result_type": "dict",
  "result_preview": "{'covered_bases': 4, 'total_depth': 4, 'mean_depth': 1.0, 'by_reference': {'chr1': {0: 1, 1: 1, 2: 1, 3: 1}}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 11) calculate_depth

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_depth` 功能。候选编号 CAP-0011，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_depth` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_depth(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_depth` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_depth

result = calculate_depth(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0011",
  "operation": "calculate_depth",
  "result_type": "dict",
  "result_preview": "{'covered_bases': 4, 'total_depth': 4, 'mean_depth': 1.0, 'by_reference': {'chr1': {0: 1, 1: 1, 2: 1, 3: 1}}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 12) calculate_insert_size

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_insert_size` 功能。候选编号 CAP-0012，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_insert_size` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_insert_size(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_insert_size` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_insert_size

result = calculate_insert_size(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0012",
  "operation": "calculate_insert_size",
  "result_type": "dict",
  "result_preview": "{'count': 0, 'mean': 0.0, 'values': []}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 13) calculate_mapping_rate

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_mapping_rate` 功能。候选编号 CAP-0013，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_mapping_rate` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_mapping_rate(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_mapping_rate` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_mapping_rate

result = calculate_mapping_rate(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0013",
  "operation": "calculate_mapping_rate",
  "result_type": "float",
  "result_preview": "1.0",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) calculate_mismatch_rate

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_mismatch_rate` 功能。候选编号 CAP-0014，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_mismatch_rate` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_mismatch_rate(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_mismatch_rate` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_mismatch_rate

result = calculate_mismatch_rate(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0014",
  "operation": "calculate_mismatch_rate",
  "result_type": "float",
  "result_preview": "0.0",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) calculate_read_group_stats

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_read_group_stats` 功能。候选编号 CAP-0015，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_read_group_stats` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_read_group_stats(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_read_group_stats` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_read_group_stats

result = calculate_read_group_stats(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0015",
  "operation": "calculate_read_group_stats",
  "result_type": "dict",
  "result_preview": "{'unknown': {'records': 1, 'mapped': 1}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) calculate_soft_clipping

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_soft_clipping` 功能。候选编号 CAP-0016，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_soft_clipping` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_soft_clipping(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_soft_clipping` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_soft_clipping

result = calculate_soft_clipping(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0016",
  "operation": "calculate_soft_clipping",
  "result_type": "dict",
  "result_preview": "{'r1': 0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 17) calculate_strand_balance

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_strand_balance` 功能。候选编号 CAP-0017，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_strand_balance` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_strand_balance(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_strand_balance` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_strand_balance

result = calculate_strand_balance(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0017",
  "operation": "calculate_strand_balance",
  "result_type": "dict",
  "result_preview": "{'forward': 1, 'reverse': 0, 'forward_fraction': 1.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 18) calculate_strand_bias

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `calculate_strand_bias` 功能。候选编号 CAP-0018，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `calculate_strand_bias` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.calculate_strand_bias(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `calculate_strand_bias` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import calculate_strand_bias

result = calculate_strand_bias(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0018",
  "operation": "calculate_strand_bias",
  "result_type": "dict",
  "result_preview": "{'bias': 1.0, 'forward': 1, 'reverse': 0, 'forward_fraction': 1.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 19) call_consensus_from_pileup

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `call_consensus_from_pileup` 功能。候选编号 CAP-0019，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `call_consensus_from_pileup` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.call_consensus_from_pileup(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `call_consensus_from_pileup` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import call_consensus_from_pileup

result = call_consensus_from_pileup(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0019",
  "operation": "call_consensus_from_pileup",
  "result_type": "dict",
  "result_preview": "{'chr1': 'ACGT'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 20) filter_bam_by_flags

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `filter_bam_by_flags` 功能。候选编号 CAP-0020，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_bam_by_flags` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.filter_bam_by_flags(input[required], include_flags[optional], exclude_flags[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`include_flags`、`exclude_flags`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `filter_bam_by_flags` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import filter_bam_by_flags

result = filter_bam_by_flags(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0020",
  "operation": "filter_bam_by_flags",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 21) filter_bam_by_mapping_quality

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `filter_bam_by_mapping_quality` 功能。候选编号 CAP-0021，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `filter_bam_by_mapping_quality` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.filter_bam_by_mapping_quality(input[required], minimum[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`minimum`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `filter_bam_by_mapping_quality` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import filter_bam_by_mapping_quality

result = filter_bam_by_mapping_quality(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0021",
  "operation": "filter_bam_by_mapping_quality",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 22) generate_pileup

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `generate_pileup` 功能。候选编号 CAP-0022，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `generate_pileup` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.generate_pileup(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `generate_pileup` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import generate_pileup

result = generate_pileup(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0022",
  "operation": "generate_pileup",
  "result_type": "list",
  "result_preview": "[{'reference': 'chr1', 'position': 0, 'depth': 1}, {'reference': 'chr1', 'position': 1, 'depth': 1}, {'reference': 'chr1', 'position': 2, 'depth': 1}, {'reference': 'chr1', 'position': 3, 'depth': 1}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) index_alignments

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `index_alignments` 功能。候选编号 CAP-0023，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `index_alignments` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.index_alignments(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `index_alignments` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import index_alignments

result = index_alignments(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0023",
  "operation": "index_alignments",
  "result_type": "dict",
  "result_preview": "{'chr1': [0]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) mark_duplicates

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `mark_duplicates` 功能。候选编号 CAP-0024，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `mark_duplicates` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.mark_duplicates(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `mark_duplicates` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import mark_duplicates

result = mark_duplicates(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0024",
  "operation": "mark_duplicates",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 25) merge_bam_files

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `merge_bam_files` 功能。候选编号 CAP-0025，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `merge_bam_files` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.merge_bam_files(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `merge_bam_files` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import merge_bam_files

result = merge_bam_files(
    input=['@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n']
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0025",
  "operation": "merge_bam_files",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 26) parse_bam_records

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `parse_bam_records` 功能。候选编号 CAP-0026，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `parse_bam_records` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/alignment/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。此条目仅部分实现：本地示例处理 SAM 文本；二进制 BAM/CRAM 需要经过验证的 pysam 或同等后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.parse_bam_records(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `parse_bam_records` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。二进制 BAM/CRAM 的解析或转换尚未验证，仍待完成。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import parse_bam_records

result = parse_bam_records(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0026",
  "operation": "parse_bam_records",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 27) parse_cram_records

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `parse_cram_records` 功能。候选编号 CAP-0027，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `parse_cram_records` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/alignment/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。此条目仅部分实现：本地示例处理 SAM 文本；二进制 BAM/CRAM 需要经过验证的 pysam 或同等后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.parse_cram_records(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `parse_cram_records` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。二进制 BAM/CRAM 的解析或转换尚未验证，仍待完成。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import parse_cram_records

result = parse_cram_records(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0027",
  "operation": "parse_cram_records",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 28) parse_sam_records

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `parse_sam_records` 功能。候选编号 CAP-0028，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `parse_sam_records` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.parse_sam_records(source[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`source`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `parse_sam_records` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import parse_sam_records

result = parse_sam_records(
    source='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0028",
  "operation": "parse_sam_records",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 29) plot_alignment_view

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `plot_alignment_view` 功能。候选编号 CAP-0029，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `plot_alignment_view` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。已有同名函数证据：`plot_alignment`。

- **API：** `dnakit.alignment.plot_alignment_view(input[required], width[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`width`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `plot_alignment_view` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import plot_alignment_view

result = plot_alignment_view(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n',
    width=3
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0029",
  "operation": "plot_alignment_view",
  "result_type": "str",
  "result_preview": "'<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"24\" height=\"30\"><text x=\"0\" y=\"18\">r1 chr1:1 ACGT</text></svg>'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 30) plot_coverage_track

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `plot_coverage_track` 功能。候选编号 CAP-0030，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `plot_coverage_track` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.plot_coverage_track(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `plot_coverage_track` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import plot_coverage_track

result = plot_coverage_track(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0030",
  "operation": "plot_coverage_track",
  "result_type": "str",
  "result_preview": "'<svg xmlns=\"http://www.w3.org/2000/svg\"><text>mean_depth=1.0</text></svg>'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 31) remove_duplicates

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `remove_duplicates` 功能。候选编号 CAP-0031，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `remove_duplicates` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.remove_duplicates(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `remove_duplicates` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import remove_duplicates

result = remove_duplicates(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0031",
  "operation": "remove_duplicates",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 32) sam_to_bam

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `sam_to_bam` 功能。候选编号 CAP-0032，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `sam_to_bam` 操作，并返回附带来源与诊断信息的结构化结果。实现位于 `src/dnakit/alignment/__init__.py`，并通过正式的 DNAKit 命名空间导出。函数的 JSON 调用方案兼容时，也可以通过 Agent/MCP 调用。此条目仅部分实现：本地示例处理 SAM 文本；二进制 BAM/CRAM 需要经过验证的 pysam 或同等后端。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.sam_to_bam(input[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `sam_to_bam` 直接返回 Python 原生值（例如列表、映射、字符串、字节或记录对象），而不是正式的 DNAKit 结果封装。下方示例输出给出实际的返回类型和截断后的预览。二进制 BAM/CRAM 的解析或转换尚未验证，仍待完成。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import sam_to_bam

result = sam_to_bam(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0032",
  "operation": "sam_to_bam",
  "result_type": "bytes",
  "result_preview": "<bytes length=34>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 33) sort_alignments

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `sort_alignments` 功能。候选编号 CAP-0033，来自 10 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `sort_alignments` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `backend` 用于选择一种支持的变体（bowtie2、bwa、graphmap、hisat2、minimap2、ngmlr、pbmm2、star、vg、winnowmap）。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`merged_dimension`。已有同名函数证据：`plot_alignment`。

- **API：** `dnakit.alignment.sort_alignments(input[required], by[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`；可选参数：`by`。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `sort_alignments` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import sort_alignments

result = sort_alignments(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0033",
  "operation": "sort_alignments",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) split_bam_by_region

- **作用：** 为“读段比对与 BAM 统计”类分析提供可复用的 `split_bam_by_region` 功能。候选编号 CAP-0034，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `split_bam_by_region` 操作，并返回附带来源与诊断信息的结构化结果。确定性的本地实现位于 `src/dnakit/alignment/_candidates.py`，并通过正式的 DNAKit 命名空间导出。它处理类 SAM 文本或已读入的记录；二进制 BAM/CRAM 解析、基于参考索引的比对和外部比对工具不在此范围内。审核依据：`no_exact_public_match`。比对中未发现名称完全相同的公开函数。

- **API：** `dnakit.alignment.split_bam_by_region(input[required], regions[required])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`input`、`regions`；可选参数：无。比对输入接受 SAM 文本、记录映射或本地的 read/reference 封装。

- **输出：** 本地函数 `split_bam_by_region` 根据类 SAM 文本或已读入的记录返回确定性的 Python 结果。它不解析二进制 BAM/CRAM，也不调用基于参考索引的外部比对工具。

- **示例代码：**

```python
from pathlib import Path

from dnakit.alignment import split_bam_by_region

result = split_bam_by_region(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n',
    regions={"chr1": [(0, 4)]}
)
print(result)
```

- **示例输出：**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0034",
  "operation": "split_bam_by_region",
  "result_type": "dict",
  "result_preview": "{'chr1': [{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```
