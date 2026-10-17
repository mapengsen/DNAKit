# 测序质量控制

本类收录 **50** 项功能。其中 46 项为 DNAKit 原生的纯 Python 实现：42 项根据所给读段、SAM 文本、片段表或计数实际计算并返回 `status: "ok"`，4 项（`basecall_nanopore_reads`、`calculate_phasing_qc`、`calculate_phasing_qc_tenx`、`correct_long_reads`）始终返回 `status: "not_implemented"`，并给出 `required_dependency` 与 `message`，因为它们需要碱基识别器、Illumina InterOp 文件或长读段校正程序。其余 4 项（`run_cutadapt`、`run_fastp`、`run_fastqc`、`run_multiqc`）是调用调用者自行安装的外部程序的受限封装。原生指标仅基于所给读段（不含比对），遵循 FastQC、Picard、umi_tools 和 Cell Ranger 的定义，其中若干在各自小节中标明为启发式。这些函数通过正式的 DNAKit 命名空间导出，JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。原生函数的示例输出均为真实运行结果，较长的数组和字符串以 `...` 省略。

## 1) assign_sample_barcodes

- **作用：** 将每条读段的条形码与样本表比对（允许一定错配），把读段分配到样本。

- **处理方法：** 条形码取自读段起始端、末端或 Illumina 头注释（注释末尾的索引，双索引为 `i7+i5`），与每个样本条形码计算汉明距离，`N` 计为错配。若最近样本的距离不超过 `max_mismatches`，且最优样本比次优样本至少近 `min_distance_gap` 个错配，则分配给该样本；没有任何条形码落在容忍范围内的读段为 `unassigned`，不满足间隔规则的为 `ambiguous`。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.assign_sample_barcodes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列）。不接受读段对映射（抛出 `INVALID_INPUT`）；双端数据请只传入带样本条形码的那条读段（通常是索引读段或 read 1）。选项：`samples`（必填；样本 -> 条形码的映射，双索引写作 `i7+i5`）、`barcode_location`（`read_start` 默认、`read_end` 或 `header`）、`barcode_length`（默认取样本表条形码长度）、`max_mismatches`（默认 1）、`min_distance_gap`（默认 1）、`max_listed`（返回的逐读段记录数，默认 100000）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`barcode_location`、`max_mismatches`、`min_distance_gap`、`reads_per_sample`、`assigned`、`unassigned`、`ambiguous`、`assigned_fraction`、`min_distance_between_sample_barcodes` 以及 `assignments`（每条读段的 `id`、`barcode`、`sample`、`distance`、`status`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import assign_sample_barcodes

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("a", "ACGTAAAA", "IIIIJJJJ"),
    ("b", "ACGAAAAA", "IIIIJJJJ"),
    ("c", "TTTTCCCC", "IIIIJJJJ"),
    ("d", "GGGGCCCC", "IIIIJJJJ"),
)
result = assign_sample_barcodes(reads, options={"samples": {"S1": "ACGT", "S2": "TTTT"}, "barcode_length": 4})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "assign_sample_barcodes",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 92,
    "sequence_bases": 32,
    "lengths": [8, 8, 8, 8],
    "sha256": "c5f05df413cd46a9d8fecf917a84d8586235800f704f7497bc3ea55d1dfcc2ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "barcode_location": "read_start",
    "max_mismatches": 1,
    "min_distance_gap": 1,
    "reads_per_sample": {
      "S1": 2,
      "S2": 1
    },
    "assigned": 3,
    "unassigned": 1,
    "ambiguous": 0,
    "assigned_fraction": 0.75,
    "min_distance_between_sample_barcodes": 3,
    "assignments": [
      {
        "id": "a",
        "barcode": "ACGT",
        "sample": "S1",
        "distance": 0,
        "status": "assigned"
      },
      {
        "id": "b",
        "barcode": "ACGA",
        "sample": "S1",
        "distance": 1,
        "status": "assigned"
      },
      "..."
    ]
  }
}
```

## 2) basecall_nanopore_reads

- **作用：** 原定用神经网络碱基识别器将 Oxford Nanopore 原始信号（FAST5/POD5）转换为碱基。

- **处理方法：** 未实现。碱基识别需要用神经网络模型把原始电流信号转换为碱基，而 DNAKit 未内置任何碱基识别器或模型。该函数始终返回 `status: "not_implemented"`，不做任何计算。

- **API：** `dnakit.sequencing_qc.basecall_nanopore_reads(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 可为路径或任意对象，不参与任何计算；`options` 可省略，会被忽略。

- **输出：** 始终返回封装，含 `status: "not_implemented"`、`module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary`、`warnings`、`required_dependency`（"Oxford Nanopore basecaller (dorado/guppy) and trained model weights"）和说明原因的 `message`，没有 `result` 键。

- **示例代码：**

```python
from dnakit.sequencing_qc import basecall_nanopore_reads

result = basecall_nanopore_reads("run1.pod5")
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "basecall_nanopore_reads",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 9,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "10557ad64d8a6c03ccd022c14b535ba4280940dd255f1a2f5cede803ce64fa27"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Oxford Nanopore basecaller (dorado/guppy) and trained model weights",
  "message": "Basecalling converts raw FAST5/POD5 current signal to bases with a neural-network model; no basecaller or model is bundled."
}
```

## 3) calculate_adapter_content

- **作用：** 按 FastQC 风格的 k-mer 规则计算沿读段位置累积的接头含量。

- **处理方法：** 当接头前 min(`kmer_size`, 接头长度) 个碱基出现在读段中（汉明距离不超过 `max_mismatches`）时，该读段命中该接头。位置 p 处的累积百分比为首次命中起点不晚于 p（从 1 起计）的读段所占比例，与 FastQC 一致。被读段末端截断的接头检测不到。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_adapter_content(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。选项：`adapters`（名称 -> 序列或序列列表的映射；默认为 FastQC 的 Illumina 通用接头、small RNA 3' 与 5' 接头、Nextera 和 SOLiD small RNA 接头）、`kmer_size`（默认 12）、`max_mismatches`（默认 0）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`kmer_size`、`max_mismatches`、`reads_with_any_adapter`、`fraction_reads_with_any_adapter`、`adapters`（每个接头的 `sequence`、`kmer_used`、`reads_with_adapter`、`fraction_reads_with_adapter`、`max_cumulative_percent`、`cumulative_percent_by_position`）以及 `definition` 说明。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_adapter_content

reads = ["AAAAAGATCGGAAGAGCTT", "A" * 19]
result = calculate_adapter_content(reads, options={"adapters": {"Illumina Universal Adapter": "AGATCGGAAGAG"}})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_adapter_content",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 46,
    "sequence_bases": 38,
    "lengths": [
      19,
      19
    ],
    "sha256": "fbd53e50eb19d3d8570f3ff850b3a816ffdeb3ed41eb8f31c591d1957cee9b9c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 2,
    "kmer_size": 12,
    "max_mismatches": 0,
    "reads_with_any_adapter": 1,
    "fraction_reads_with_any_adapter": 0.5,
    "adapters": {
      "Illumina Universal Adapter": {
        "sequence": "AGATCGGAAGAG",
        "kmer_used": "AGATCGGAAGAG",
        "reads_with_adapter": 1,
        "fraction_reads_with_adapter": 0.5,
        "max_cumulative_percent": 50.0,
        "cumulative_percent_by_position": [
          0.0,
          0.0,
          "..."
        ]
      }
    },
    "definition": "A read hits an adapter when the first min(kmer_size, len(adapter)) bases of the adapter occur in the read (Hamming <= max_mismatches). The cumulative percentage..."
  }
}
```

## 4) calculate_adapter_content_tenx

- **作用：** 计算 10x Chromium read 2（cDNA）中的接头、模板转换寡核苷酸（TSO）和 poly-A 含量。

- **处理方法：** 对 read 2 使用与 `calculate_adapter_content` 相同的 k-mer 规则；默认接头列表另外加入 10x 模板转换寡核苷酸（AAGCAGTGGTATCAACGCAGAGTACATGGG）和 12 聚体 poly-A。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_adapter_content_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有与 `calculate_adapter_content` 相同的选项：`adapters`、`kmer_size`（12）、`max_mismatches`（0）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`region`（`rna_read (read 2)`）以及 `calculate_adapter_content` 对 read 2 给出的字段（`reads`、`kmer_size`、`max_mismatches`、`reads_with_any_adapter`、`fraction_reads_with_any_adapter`、`adapters`、`definition`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_adapter_content_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_adapter_content_tenx(tenx, options={"adapters": {"Illumina Universal Adapter": "AGATCGGAAGAG"}})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_adapter_content_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "region": "rna_read (read 2)",
    "reads": 3,
    "kmer_size": 12,
    "max_mismatches": 0,
    "reads_with_any_adapter": 2,
    "fraction_reads_with_any_adapter": 0.666667,
    "adapters": {
      "Illumina Universal Adapter": {
        "sequence": "AGATCGGAAGAG",
        "kmer_used": "AGATCGGAAGAG",
        "reads_with_adapter": 2,
        "fraction_reads_with_adapter": 0.666667,
        "max_cumulative_percent": 66.666667,
        "cumulative_percent_by_position": [
          0.0,
          0.0,
          "..."
        ]
      }
    },
    "definition": "A read hits an adapter when the first min(kmer_size, len(adapter)) bases of the adapter occur in the read (Hamming <= max_mismatches). The cumulative percentage..."
  }
}
```

## 5) calculate_barcode_quality

- **作用：** 计算读段条形码（样本条形码或细胞条形码）的测序质量与有效性。

- **处理方法：** 条形码取自读段起始端、末端或头注释。给出不同条形码的数量、含 `N` 的条形码、出现最多的条形码、条形码计数的香农熵（比特），以及（有质量值时）平均质量、Q30 比例和逐位置平均质量。给出 `whitelist` 时，每个条形码被归为完全匹配、已校正（`max_mismatches` 内存在唯一的白名单条形码）、有歧义（多个距离相同）或无效。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_barcode_quality(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列）。不接受读段对映射（抛出 `INVALID_INPUT`）；双端数据请只传入带条形码的那条读段。选项：`barcode_location`（`read_start` 默认、`read_end` 或 `header`）、`barcode_length`（默认 8；`header` 时不使用）、`whitelist`（有效条形码的列表或换行分隔的文本）、`max_mismatches`（默认 1，配合白名单使用）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`barcode_source`、`distinct_barcodes`、`reads_with_N_in_barcode`、`top_barcodes`、`barcode_entropy_bits`、`mean_quality`、`q30_fraction`、`mean_quality_by_position`，以及（给出白名单时）`whitelist`（`size`、`max_mismatches`、`exact_match_fraction`、`corrected_fraction`、`ambiguous_fraction`、`invalid_fraction`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_barcode_quality

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("a", "ACGTTTTT", "IIII####"),
    ("b", "ACGATTTT", "IIII####"),
    ("c", "TTTTAAAA", "IIII####"),
    ("d", "ACGTCCCC", "IIII####"),
)
result = calculate_barcode_quality(reads, options={"barcode_length": 4, "whitelist": ["ACGT", "TTTT"]})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_barcode_quality",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 92,
    "sequence_bases": 32,
    "lengths": [8, 8, 8, 8],
    "sha256": "1ea787f978a89222ad3f525a54d47607e60402c92e5b3e1e32f8f3633b41e325"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "barcode_source": "read_start",
    "distinct_barcodes": 3,
    "reads_with_N_in_barcode": 0,
    "top_barcodes": [
      {
        "barcode": "ACGT",
        "reads": 2
      },
      {
        "barcode": "ACGA",
        "reads": 1
      },
      "..."
    ],
    "barcode_entropy_bits": 1.5,
    "mean_quality": 40.0,
    "q30_fraction": 1.0,
    "mean_quality_by_position": [
      40.0,
      40.0,
      "..."
    ],
    "whitelist": {
      "size": 2,
      "max_mismatches": 1,
      "exact_match_fraction": 0.75,
      "corrected_fraction": 0.25,
      "ambiguous_fraction": 0.0,
      "invalid_fraction": 0.0
    }
  }
}
```

## 6) calculate_barcode_quality_tenx

- **作用：** 计算 read 1 中 10x 细胞条形码和 UMI 的质量以及白名单有效性。

- **处理方法：** 按化学版本的布局把 read 1 拆为细胞条形码和 UMI。细胞条形码的报告同 `calculate_barcode_quality`（可带白名单）；UMI 给出不同 UMI 数、`N` 数、平均质量、Q30 比例以及同聚物 UMI 的读段数。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_barcode_quality_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有选项：`whitelist`（有效细胞条形码）、`max_mismatches`（默认 1）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`cell_barcode`（`calculate_barcode_quality` 的报告）以及 `umi`（`distinct_umis`、`reads_with_N_in_umi`、`mean_quality`、`q30_fraction`、`homopolymer_umi_reads`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_barcode_quality_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_barcode_quality_tenx(tenx, options={"whitelist": [barcode]})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_barcode_quality_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "cell_barcode": {
      "reads": 3,
      "barcode_source": "read_start",
      "distinct_barcodes": 2,
      "reads_with_N_in_barcode": 0,
      "top_barcodes": [
        {
          "barcode": "ACGTACGTACGTACGT",
          "reads": 2
        },
        {
          "barcode": "TTTTTTTTTTTTTTTT",
          "reads": 1
        }
      ],
      "barcode_entropy_bits": 0.918296,
      "mean_quality": 40.0,
      "q30_fraction": 1.0,
      "mean_quality_by_position": [
        40.0,
        40.0,
        "..."
      ],
      "whitelist": {
        "size": 1,
        "max_mismatches": 1,
        "exact_match_fraction": 0.666667,
        "corrected_fraction": 0.0,
        "ambiguous_fraction": 0.0,
        "invalid_fraction": 0.333333
      }
    },
    "umi": {
      "distinct_umis": 2,
      "reads_with_N_in_umi": 0,
      "mean_quality": 33.333333,
      "q30_fraction": 0.666667,
      "homopolymer_umi_reads": 3
    }
  }
}
```

## 7) calculate_basecalling_qc

- **作用：** 计算碱基识别质量指标：无法识别（N）比例、Q30 比例、每条读段的期望错误数以及 Illumina 通过过滤比例。

- **处理方法：** 无法识别比例为 `N` 碱基所占比例，Q30 比例为 Phred Q >= 30 的碱基所占比例，每条读段期望错误数为各读段所有碱基 10^(-Q/10) 之和的平均值。统计平均 Q 低于 `min_mean_quality` 的读段。通过过滤比例取自 CASAVA 头标志（`1:Y:0:...` 表示被过滤，`N` 表示通过），没有该标志时为 `null`。仪器端指标（phasing、信号强度）需要厂商文件，不予计算。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_basecalling_qc(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。必须带有 Phred 质量值。选项：`min_mean_quality`（默认 20）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`bases`、`no_call_fraction`、`reads_with_no_calls`、`mean_q`、`q30_fraction`、`expected_errors_per_read`、`min_mean_q`、`reads_below_min_mean_q`、`fraction_reads_below_min_mean_q`、`illumina_filter_flag_reads`、`pass_filter_fraction`（存在该标志时另有 `filtered_reads`）。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_basecalling_qc

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_basecalling_qc(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_basecalling_qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "bases": 40,
    "no_call_fraction": 0.05,
    "reads_with_no_calls": 1,
    "mean_q": 33.0,
    "q30_fraction": 0.7,
    "expected_errors_per_read": 0.5257,
    "min_mean_q": 20.0,
    "reads_below_min_mean_q": 0,
    "fraction_reads_below_min_mean_q": 0.0,
    "illumina_filter_flag_reads": 0,
    "pass_filter_fraction": null
  }
}
```

## 8) calculate_basecalling_qc_tenx

- **作用：** 计算 10x 读段中 Cell Ranger 风格的碱基识别指标（条形码、UMI 和 RNA 读段的 Q30）。

- **处理方法：** 把 read 1 拆为条形码和 UMI 区域，read 2 为 RNA 读段；每个区域给出 `calculate_basecalling_qc` 的各项指标（无法识别比例、Q30 比例、每条读段期望错误数、低于平均质量阈值的读段、通过过滤比例）。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_basecalling_qc_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。必须带有 Phred 质量值。另有选项：`min_mean_quality`（默认 20）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`，以及 `regions`（含 `barcode`、`umi`、`rna_read`，各自包含 `calculate_basecalling_qc` 的字段）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_basecalling_qc_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_basecalling_qc_tenx(tenx)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_basecalling_qc_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "reads": 3,
        "bases": 48,
        "no_call_fraction": 0.0,
        "reads_with_no_calls": 0,
        "mean_q": 40.0,
        "q30_fraction": 1.0,
        "expected_errors_per_read": 0.0016,
        "min_mean_q": 20.0,
        "reads_below_min_mean_q": 0,
        "fraction_reads_below_min_mean_q": 0.0,
        "illumina_filter_flag_reads": 0,
        "pass_filter_fraction": null
      },
      "umi": {
        "reads": 3,
        "bases": 36,
        "no_call_fraction": 0.0,
        "reads_with_no_calls": 0,
        "mean_q": 33.333333,
        "q30_fraction": 0.666667,
        "expected_errors_per_read": 0.0408,
        "min_mean_q": 20.0,
        "reads_below_min_mean_q": 0,
        "fraction_reads_below_min_mean_q": 0.0,
        "illumina_filter_flag_reads": 0,
        "pass_filter_fraction": null
      },
      "rna_read": {
        "reads": 3,
        "bases": 90,
        "no_call_fraction": 0.0,
        "reads_with_no_calls": 0,
        "mean_q": 33.333333,
        "q30_fraction": 0.666667,
        "expected_errors_per_read": 0.102,
        "min_mean_q": 20.0,
        "reads_below_min_mean_q": 0,
        "fraction_reads_below_min_mean_q": 0.0,
        "illumina_filter_flag_reads": 0,
        "pass_filter_fraction": null
      }
    }
  }
}
```

## 9) calculate_contamination

- **作用：** 用与用户提供参考序列的精确 k-mer 包含关系筛查读段污染。

- **处理方法：** 当读段至少 `min_kmer_fraction` 比例的规范化 k-mer（不区分链）出现在某参考序列中时，该读段命中该参考序列。结果给出各参考序列的命中数、无命中读段数以及命中多个参考序列的读段数。给出 `expected_reference` 时，污染比例为命中其他参考序列而未命中预期参考序列的读段占被筛查读段的比例。这是 k-mer 筛查而非比对，只能检测调用者提供了序列的污染物。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_contamination(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。选项：`references`（必填；名称 -> 序列的映射，或 FASTA 文本）、`kmer_size`（默认 21）、`min_kmer_fraction`（默认 0.5）、`expected_reference`（可选的名称）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads_screened`、`kmer_size`、`min_kmer_fraction`、`hits`（每个参考序列的 `reads` 与 `fraction`）、`reads_without_hit`、`reads_with_multiple_hits`、`definition`，以及（给出 `expected_reference` 时）`expected_reference`、`contamination_fraction` 和 `contamination_definition`。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_contamination

host = "ACGTTGCAAGGCTTAACCGGATCGATTCGAGGCTATTGCCAGTCAAGGTACCTGATGCCGTAAGCTTGGAACCTAGGTCCA"
other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

reads = [host[:30], host[20:50], reverse_complement(host[10:40]), other[:30]]
result = calculate_contamination(
    reads,
    options={"references": {"host": host, "bug": other}, "kmer_size": 12, "expected_reference": "host"},
)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_contamination",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 136,
    "sequence_bases": 120,
    "lengths": [30, 30, 30, 30],
    "sha256": "c2c4d13828683299b138663672c38c92457692686454ffa93ca8a805435c7942"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_screened": 4,
    "kmer_size": 12,
    "min_kmer_fraction": 0.5,
    "hits": {
      "host": {
        "reads": 3,
        "fraction": 0.75
      },
      "bug": {
        "reads": 1,
        "fraction": 0.25
      }
    },
    "reads_without_hit": 0,
    "reads_with_multiple_hits": 0,
    "definition": "A read hits a reference when at least min_kmer_fraction of its canonical k-mers occur in that reference (exact k-mer containment, both strands). It is a k-mer s...",
    "expected_reference": "host",
    "contamination_fraction": 0.25,
    "contamination_definition": "Fraction of screened reads that hit at least one non-expected reference and not the expected one."
  }
}
```

## 10) calculate_contamination_tenx

- **作用：** 用与所提供参考序列的 k-mer 包含关系筛查 10x 数据 cDNA 读段（read 2）的污染。

- **处理方法：** 对 read 2 的筛查方式与 `calculate_contamination` 完全相同（规范化 k-mer、`min_kmer_fraction` 规则、可选的 `expected_reference`）；read 1 仅提供布局和读段数。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_contamination_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有与 `calculate_contamination` 相同的选项：`references`（必填）、`kmer_size`（21）、`min_kmer_fraction`（0.5）、`expected_reference`。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`region`（`rna_read (read 2)`）以及 `calculate_contamination` 的字段（`reads_screened`、`hits`、`reads_without_hit`、`reads_with_multiple_hits` 等）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_contamination_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
host = "GATTACAGATCGGAAGAGCGGTTCAGCAGGTTACCAGTCAGGCCTAAGGTCCATTGGCC"
result = calculate_contamination_tenx(tenx, options={"references": {"host": host}, "kmer_size": 12})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_contamination_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "region": "rna_read (read 2)",
    "reads_screened": 3,
    "kmer_size": 12,
    "min_kmer_fraction": 0.5,
    "hits": {
      "host": {
        "reads": 2,
        "fraction": 0.666667
      }
    },
    "reads_without_hit": 1,
    "reads_with_multiple_hits": 0,
    "definition": "A read hits a reference when at least min_kmer_fraction of its canonical k-mers occur in that reference (exact k-mer containment, both strands). It is a k-mer s..."
  }
}
```

## 11) calculate_duplicate_rate

- **作用：** 计算读段的完全序列重复率：(读段数 - 不同序列数) / 读段数。

- **处理方法：** 按完全序列比较读段（可只比较前 `prefix_length` 个碱基，也可把读段与其反向互补序列视为相同）。对于读段对，比较键为两条读段的拼接。读段位置未知，因此这是序列重复，而不是 Picard 标记的比对重复。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_duplicate_rate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或读段对 `{"read1": ..., "read2": ...}`（每对一个键）。选项：`prefix_length`（只比较前 N 个碱基；默认 0 即整条读段）、`reverse_complement_aware`（默认 `False`）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `paired`、`reads`、`distinct_sequences`、`duplicate_reads`、`duplicate_rate`、`unique_fraction`、`most_common`（序列及计数）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_duplicate_rate

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_duplicate_rate(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_duplicate_rate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "paired": false,
    "reads": 4,
    "distinct_sequences": 3,
    "duplicate_reads": 1,
    "duplicate_rate": 0.25,
    "unique_fraction": 0.75,
    "most_common": [
      {
        "sequence": "ACGTACGTAC",
        "count": 2
      },
      {
        "sequence": "NNACGTACGT",
        "count": 1
      },
      "..."
    ],
    "definition": "duplicate_rate = (reads - distinct sequences) / reads; sequences are compared exactly (read pairs by the concatenation of both mates)."
  }
}
```

## 12) calculate_duplication

- **作用：** 计算读段的 FastQC 风格序列重复水平直方图。

- **处理方法：** 序列的重复水平指该完全相同序列的拷贝数，按 FastQC 的方式分箱：1-9、>10（10-49）、>50、>100、>500、>1k、>5k 和 >10k+。每个水平分别给出占总读段和占不同序列的百分比，以及去重后剩余读段的百分比。不做 FastQC 式的对未见读段的外推。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_duplication(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或读段对 `{"read1": ..., "read2": ...}`（每对一个键）。选项：`prefix_length`（默认 0）、`reverse_complement_aware`（默认 `False`）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `paired`、`reads`、`distinct_sequences`、`percent_remaining_after_deduplication`、`levels`（16 行，含 `level`、`distinct_sequences`、`percent_of_total_reads`、`percent_of_deduplicated`）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_duplication

reads = ["AAAA"] * 3 + ["CCCC", "GGGG"]
result = calculate_duplication(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_duplication",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 40,
    "sequence_bases": 20,
    "lengths": [4, 4, 4, 4, 4],
    "sha256": "9a2abe2e853265a2f1ad0d543707ac64ac86a065f17d998ffdb72db0f8f0dc35"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "paired": false,
    "reads": 5,
    "distinct_sequences": 3,
    "percent_remaining_after_deduplication": 60.0,
    "levels": [
      {
        "level": "1",
        "distinct_sequences": 2,
        "percent_of_total_reads": 40.0,
        "percent_of_deduplicated": 66.666667
      },
      {
        "level": "2",
        "distinct_sequences": 0,
        "percent_of_total_reads": 0.0,
        "percent_of_deduplicated": 0.0
      },
      "..."
    ],
    "definition": "Level = number of copies of an identical sequence; bins 1-9, >10 (10-49), >50, >100, >500, >1k, >5k, >10k+ as in FastQC. No extrapolation beyond the supplied re..."
  }
}
```

## 13) calculate_duplication_tenx

- **作用：** 以细胞条形码 + UMI + read 2 起始碱基为键，计算 10x 读段的重复情况。

- **处理方法：** 细胞条形码、UMI 和 read 2 前 `cdna_prefix_length` 个碱基都相同的两条读段视为重复。这是测序饱和度的启发式替代指标：Cell Ranger 使用比对得到的 (条形码, UMI, 基因)，此处无法获得。重复水平的分箱方式同 `calculate_duplication`。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_duplication_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有选项：`cdna_prefix_length`（默认 25）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`reads`、`distinct_sequences`、`percent_remaining_after_deduplication`、`levels`、`definition`、`duplicate_rate` 以及 `sequencing_saturation_proxy`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_duplication_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_duplication_tenx(tenx)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_duplication_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 3,
    "distinct_sequences": 2,
    "percent_remaining_after_deduplication": 66.666667,
    "levels": [
      {
        "level": "1",
        "distinct_sequences": 1,
        "percent_of_total_reads": 33.333333,
        "percent_of_deduplicated": 50.0
      },
      {
        "level": "2",
        "distinct_sequences": 1,
        "percent_of_total_reads": 66.666667,
        "percent_of_deduplicated": 50.0
      },
      "..."
    ],
    "definition": "Reads are duplicates when cell barcode, UMI and the first 25 bases of read 2 are identical. This is a proxy: Cell Ranger's sequencing saturation uses (barcode, ...",
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "duplicate_rate": 0.333333,
    "sequencing_saturation_proxy": 0.333333
  }
}
```

## 14) calculate_error_profile

- **作用：** 根据 Phred 质量值计算读段逐位置的期望错误率曲线。

- **处理方法：** 某位置的期望错误率为覆盖该位置的读段的 10^(-Q/10) 的平均值；结果还给出沿读段累积的期望错误数、每条读段期望错误数和逐位置的 `N` 比例。这是质量值所声称的错误率，而非观测到的错配率（观测错配率请用带参考序列的 `estimate_read_error_model`）。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_error_profile(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。必须带有 Phred 质量值。选项：`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `source`（`phred_quality`）、`overall_expected_error_rate`、`overall_q_equivalent`、`expected_errors_per_read`、`per_position`（`position`、`bases`、`expected_error_rate`、`cumulative_expected_errors`、`n_fraction`）以及 `note`。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_error_profile

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_error_profile(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_error_profile",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "phred_quality",
    "overall_expected_error_rate": 0.05256999999999999,
    "overall_q_equivalent": 12.79262,
    "expected_errors_per_read": 0.5257,
    "per_position": [
      {
        "position": 1,
        "bases": 4,
        "expected_error_rate": 0.25255,
        "cumulative_expected_errors": 0.25255,
        "n_fraction": 0.25
      },
      {
        "position": 2,
        "bases": 4,
        "expected_error_rate": 0.25255,
        "cumulative_expected_errors": 0.5051,
        "n_fraction": 0.25
      },
      "..."
    ],
    "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
  }
}
```

## 15) calculate_error_profile_tenx

- **作用：** 计算 10x 读段条形码、UMI 和 RNA 区域由 Phred 质量值推导的期望错误率曲线。

- **处理方法：** 把 read 1 拆为条形码和 UMI 区域，read 2 为 RNA 读段；每个区域给出 `calculate_error_profile` 的结果，即逐位置平均的期望错误率 10^(-Q/10)。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_error_profile_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。必须带有 Phred 质量值。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`，以及 `regions`（含 `barcode`、`umi`、`rna_read`，各自包含 `calculate_error_profile` 的字段）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_error_profile_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_error_profile_tenx(tenx)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_error_profile_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "source": "phred_quality",
        "overall_expected_error_rate": 0.0001,
        "overall_q_equivalent": 40.0,
        "expected_errors_per_read": 0.0016,
        "per_position": [
          {
            "...": "5 keys"
          },
          {
            "...": "5 keys"
          },
          "..."
        ],
        "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
      },
      "umi": {
        "source": "phred_quality",
        "overall_expected_error_rate": 0.0034000000000000002,
        "overall_q_equivalent": 24.685211,
        "expected_errors_per_read": 0.0408,
        "per_position": [
          {
            "...": "5 keys"
          },
          {
            "...": "5 keys"
          },
          "..."
        ],
        "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
      },
      "rna_read": {
        "source": "phred_quality",
        "overall_expected_error_rate": 0.0034000000000000007,
        "overall_q_equivalent": 24.685211,
        "expected_errors_per_read": 0.10200000000000002,
        "per_position": [
          {
            "...": "5 keys"
          },
          {
            "...": "5 keys"
          },
          "..."
        ],
        "note": "Expected error = 10^(-Q/10) from the reported Phred scores; not an observed mismatch rate."
      }
    }
  }
}
```

## 16) calculate_gc_bias

- **作用：** 计算读段 GC 含量的分布及其相对拟合正态分布的偏离（FastQC 约定）。

- **处理方法：** 每条读段的 GC 百分比（G+C 占 A+C+G+T 的比例）按 1% 分箱，并与具有观测均值和标准差的正态分布比较。`deviation_from_normal_percent` 为各分箱 |观测值 - 理论计数| 之和除以读段数再乘以 100（标准差为 0 时为 `null`）。同时给出按读段位置的 GC%。给出 `reference_gc` 时报告均值相对期望 GC 的偏移。这是读段组成谱，不是 Picard 的覆盖度-GC 偏好（后者需要比对）。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_gc_bias(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。选项：`reference_gc`（可选，期望 GC 比例，取值 [0, 1]）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads_with_gc`、`mean_gc_percent`、`stdev_gc_percent`、`mode_gc_percent`、`histogram`（GC% -> 读段数）、`theoretical_normal_counts`、`deviation_from_normal_percent`、`gc_percent_by_position`、`definition`，以及（给出 `reference_gc` 时）`reference_gc_percent` 和 `gc_shift_vs_reference_percent`。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_gc_bias

reads = ["GGGG", "AAAA", "GCAT", "GGAT"]
result = calculate_gc_bias(reads, options={"reference_gc": 0.4})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_gc_bias",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 32,
    "sequence_bases": 16,
    "lengths": [4, 4, 4, 4],
    "sha256": "2a486830d70f56db6c17aebf31a0bb6fc9ec7c1df90a1cb1dfbeced28e260d66"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_with_gc": 4,
    "mean_gc_percent": 50.0,
    "stdev_gc_percent": 40.824829,
    "mode_gc_percent": 50,
    "histogram": {
      "0": 1,
      "50": 2,
      "100": 1
    },
    "theoretical_normal_counts": [
      0.0185,
      0.019,
      "..."
    ],
    "deviation_from_normal_percent": 174.590226,
    "gc_percent_by_position": [
      75.0,
      75.0,
      "..."
    ],
    "definition": "Per-read GC% (G+C over A+C+G+T) histogrammed in 1% bins; the theoretical curve is a normal distribution with the observed mean and standard deviation (FastQC co...",
    "reference_gc_percent": 40.0,
    "gc_shift_vs_reference_percent": 10.0
  }
}
```

## 17) calculate_gc_bias_tenx

- **作用：** 计算 10x cDNA 读段（read 2）的 GC 含量分布。

- **处理方法：** 对 read 2 的分析方式与 `calculate_gc_bias` 完全相同（1% 分箱、拟合正态曲线、偏离度、逐位置 GC%）。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_gc_bias_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有选项：`reference_gc`（可选的 GC 比例，取值 [0, 1]）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`region`（`rna_read (read 2)`）以及 `calculate_gc_bias` 的字段。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_gc_bias_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_gc_bias_tenx(tenx, options={"reference_gc": 0.4})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_gc_bias_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "region": "rna_read (read 2)",
    "reads_with_gc": 3,
    "mean_gc_percent": 68.888889,
    "stdev_gc_percent": 26.943013,
    "mode_gc_percent": 53,
    "histogram": {
      "53": 2,
      "100": 1
    },
    "theoretical_normal_counts": [
      0.0017,
      0.0019,
      "..."
    ],
    "deviation_from_normal_percent": 183.455679,
    "gc_percent_by_position": [
      100.0,
      33.333333,
      "..."
    ],
    "definition": "Per-read GC% (G+C over A+C+G+T) histogrammed in 1% bins; the theoretical curve is a normal distribution with the observed mean and standard deviation (FastQC co...",
    "reference_gc_percent": 40.0,
    "gc_shift_vs_reference_percent": 28.888889
  }
}
```

## 18) calculate_index_hopping

- **作用：** 根据读段头中的样本索引，计算双索引读段的索引跳跃（index hopping）率。

- **处理方法：** 从 Illumina 头注释末尾读取 `i7+i5` 索引（如 `1:N:0:ACGTAC+TTGGCA`）。若读段的 i7 与 i5 分别出现在样本表中、但并非同属一个样本，则为跳跃读段；跳跃率 = 跳跃读段 / (跳跃读段 + 具有预期配对的读段)。索引必须完全匹配（不容错配）。索引不在样本表中的读段单独作为未知索引报告，不计入该比率。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_index_hopping(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），其头注释须以 `i7+i5` 结尾。不接受读段对映射（抛出 `INVALID_INPUT`）；双端数据请只传入其中一条读段（两条读段的索引相同）。选项：`samples`（必填；样本 -> `i7+i5` 的映射）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`expected_combination_reads`、`index_hopped_reads`、`unknown_index_reads`、`index_hopping_rate`、`unknown_index_fraction`、`hopped_combinations`（`index`、`reads`）、`reads_per_expected_sample` 以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_index_hopping

records = (
    [("s1", "A", "I", "1:N:0:AAAA+CCCC")] * 6
    + [("s2", "A", "I", "1:N:0:GGGG+TTTT")] * 2
    + [("hop", "A", "I", "1:N:0:AAAA+TTTT")]
    + [("unk", "A", "I", "1:N:0:CCCC+AAAA")]
)
reads = "".join(f"@{name} {comment}\n{seq}\n+\n{qual}\n" for name, seq, qual, comment in records)
result = calculate_index_hopping(reads, options={"samples": {"S1": "AAAA+CCCC", "S2": "GGGG+TTTT"}})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_index_hopping",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 10,
    "text_length": 262,
    "sequence_bases": 10,
    "lengths": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "sha256": "3be23d01c8b0b3436149a94c2fe9a3fb2a272880ae45bc5fdf7650fe52f6b2ee"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 10,
    "expected_combination_reads": 8,
    "index_hopped_reads": 1,
    "unknown_index_reads": 1,
    "index_hopping_rate": 0.111111,
    "unknown_index_fraction": 0.1,
    "hopped_combinations": [
      {
        "index": "AAAA+TTTT",
        "reads": 1
      }
    ],
    "reads_per_expected_sample": {
      "S1": 6,
      "S2": 2
    },
    "definition": "A read is index-hopped when its i7 and i5 are each used in the sample sheet but the pair is not. rate = hopped / (hopped + expected-pair reads). Reads whose ind..."
  }
}
```

## 19) calculate_index_hopping_tenx

- **作用：** 根据读段头中的样本索引，计算双索引 10x 文库的索引跳跃率。

- **处理方法：** 计算方式与 `calculate_index_hopping` 相同（`i7+i5` 精确匹配，跳跃率 = 跳跃读段 / (跳跃读段 + 预期配对读段)）。适用于 10x 双索引试剂盒；单索引试剂盒每个样本只有一个索引，没有跳跃信号，会报错。

- **API：** `dnakit.sequencing_qc.calculate_index_hopping_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），其头注释须以 `i7+i5` 结尾。不接受读段对映射（抛出 `INVALID_INPUT`）；双端数据请只传入其中一条读段。选项同 `calculate_index_hopping`：`samples`（必填）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 与 `calculate_index_hopping` 相同的字段（`reads`、`expected_combination_reads`、`index_hopped_reads`、`unknown_index_reads`、`index_hopping_rate`、`unknown_index_fraction`、`hopped_combinations`、`reads_per_expected_sample`、`definition`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_index_hopping_tenx

records = (
    [("s1", "A", "I", "1:N:0:AAAA+CCCC")] * 6
    + [("s2", "A", "I", "1:N:0:GGGG+TTTT")] * 2
    + [("hop", "A", "I", "1:N:0:AAAA+TTTT")]
    + [("unk", "A", "I", "1:N:0:CCCC+AAAA")]
)
reads = "".join(f"@{name} {comment}\n{seq}\n+\n{qual}\n" for name, seq, qual, comment in records)
result = calculate_index_hopping_tenx(reads, options={"samples": {"S1": "AAAA+CCCC", "S2": "GGGG+TTTT"}})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_index_hopping_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 10,
    "text_length": 262,
    "sequence_bases": 10,
    "lengths": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    "sha256": "3be23d01c8b0b3436149a94c2fe9a3fb2a272880ae45bc5fdf7650fe52f6b2ee"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 10,
    "expected_combination_reads": 8,
    "index_hopped_reads": 1,
    "unknown_index_reads": 1,
    "index_hopping_rate": 0.111111,
    "unknown_index_fraction": 0.1,
    "hopped_combinations": [
      {
        "index": "AAAA+TTTT",
        "reads": 1
      }
    ],
    "reads_per_expected_sample": {
      "S1": 6,
      "S2": 2
    },
    "definition": "A read is index-hopped when its i7 and i5 are each used in the sample sheet but the pair is not. rate = hopped / (hopped + expected-pair reads). Reads whose ind..."
  }
}
```

## 20) calculate_insert_size

- **作用：** 根据读段对重叠、SAM 文本、10x fragments 表或大小列表，计算插入片段（fragment）大小的统计量。

- **处理方法：** 根据输入类型，大小来源于 read 1 与 read 2 反向互补序列的重叠（无法检测到重叠的读段对，即插入片段长于两条读段之和者，被排除并计数）、已比对的第一条 SAM 读段的绝对模板长度（TLEN，默认只用 proper pair）、fragments 表（end - start，按 count 加权）或直接给出的数值。超过 `max_insert_size` 的大小被排除并计数。统计量包括数量、范围、均值、中位数、标准差、绝对中位差、众数和百分位数。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_insert_size(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为读段对 `{"read1": ..., "read2": ...}`（片段长度由重叠推得，即 read 1 与 read 2 长度之和减去重叠长度，其中 read 1 与 read 2 的反向互补序列重叠）、SAM 文本或路径（已比对的第一条读段的 |TLEN|）、片段/类 BED 的 TSV（`chrom`、`start`、`end`，另可有 `barcode`、`count`；大小 = end - start，按 count 加权）或数值列表。选项：`min_overlap`（默认 10）、`max_mismatch_rate`（默认 0.1）、`max_pairs`（按重叠定大小的读段对上限，默认 200000）、`proper_pairs_only`（SAM，默认 `True`）、`max_insert_size`（默认 100000）、`phred_offset`、`max_reads`。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `source`（`read-pair overlap`、`SAM TLEN`、`fragments (end - start)` 或 `numbers`）、`count`、`excluded_above_max`、`min`、`max`、`mean`、`median`、`stdev`、`median_absolute_deviation`、`mode`、`percentiles`（`p5`、`p25`、`p75`、`p95`）以及与来源相关的计数（`pairs`、`pairs_without_overlap`、`non_proper_pairs_skipped` 等）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_insert_size

fragment = "GATTACAGCTAGGCTTACGGATCCATGCAA"

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

pairs = {"read1": [fragment[:20]], "read2": [reverse_complement(fragment[10:])]}
result = calculate_insert_size(pairs, options={"min_overlap": 10, "max_mismatch_rate": 0.0})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_insert_size",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "52677fd537fe6e846e979fc7d0a79c8e4dfe046926d05023ec12df1cfd1de470"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "read-pair overlap",
    "pairs": 1,
    "pairs_without_overlap": 0,
    "min_overlap": 10,
    "max_mismatch_rate": 0.0,
    "count": 1,
    "excluded_above_max": 0,
    "min": 30,
    "max": 30,
    "mean": 30.0,
    "median": 30.0,
    "stdev": null,
    "median_absolute_deviation": 0.0,
    "mode": 30,
    "percentiles": {
      "p5": 30.0,
      "p25": 30.0,
      "p75": 30.0,
      "p95": 30.0
    }
  }
}
```

## 21) calculate_insert_size_distribution

- **作用：** 计算插入片段（fragment）大小的直方图及累积比例。

- **处理方法：** 大小的获取方式同 `calculate_insert_size`，随后按 `bin_width` 分箱（默认 1，范围很宽时自动加宽）。每个分箱给出计数、占全部大小的比例和累积比例；`end` 为该分箱内的最大尺寸。

- **API：** `dnakit.sequencing_qc.calculate_insert_size_distribution(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 与选项同 `calculate_insert_size`，另有 `bin_width`（默认 1）。`input` 为读段对 `{"read1": ..., "read2": ...}`（片段长度由重叠推得，即 read 1 与 read 2 长度之和减去重叠长度，其中 read 1 与 read 2 的反向互补序列重叠）、SAM 文本或路径（已比对的第一条读段的 |TLEN|）、片段/类 BED 的 TSV（`chrom`、`start`、`end`，另可有 `barcode`、`count`；大小 = end - start，按 count 加权）或数值列表。选项：`min_overlap`（默认 10）、`max_mismatch_rate`（默认 0.1）、`max_pairs`（按重叠定大小的读段对上限，默认 200000）、`proper_pairs_only`（SAM，默认 `True`）、`max_insert_size`（默认 100000）、`phred_offset`、`max_reads`。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `calculate_insert_size` 的汇总统计量，以及 `bin_width` 和 `histogram`（每个分箱的 `start`、`end`、`count`、`fraction`、`cumulative_fraction`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_insert_size_distribution

result = calculate_insert_size_distribution([100, 100, 200, 300], options={"bin_width": 100})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_insert_size_distribution",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 20,
    "sequence_bases": 700,
    "lengths": [100, 100, 200, 300],
    "sha256": "1cc3ddba72c516c92c1b4e49a8cce6e66fe4f916a8034bfc62bccecfb408d5e6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "numbers",
    "count": 4,
    "excluded_above_max": 0,
    "min": 100,
    "max": 300,
    "mean": 175.0,
    "median": 150.0,
    "stdev": 95.742711,
    "median_absolute_deviation": 50.0,
    "mode": 100,
    "percentiles": {
      "p5": 100.0,
      "p25": 100.0,
      "p75": 225.0,
      "p95": 285.0
    },
    "bin_width": 100,
    "histogram": [
      {
        "start": 100,
        "end": 199,
        "count": 2,
        "fraction": 0.5,
        "cumulative_fraction": 0.5
      },
      {
        "start": 200,
        "end": 299,
        "count": 1,
        "fraction": 0.25,
        "cumulative_fraction": 0.75
      },
      {
        "start": 300,
        "end": 399,
        "count": 1,
        "fraction": 0.25,
        "cumulative_fraction": 1.0
      }
    ]
  }
}
```

## 22) calculate_insert_size_tenx

- **作用：** 根据 fragments.tsv 表、SAM 文本或大小列表，计算 10x 数据的片段大小分布。

- **处理方法：** 计算方式与 `calculate_insert_size` 相同；对 10x scATAC 的 `fragments.tsv`（`chrom`、`start`、`end`、`barcode`、`count`），大小为 end - start 并按 count 加权。10x scRNA 的 read 1 含条形码和 UMI，因此 scRNA 文库没有读段对插入片段大小。

- **API：** `dnakit.sequencing_qc.calculate_insert_size_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 10x scATAC 的 `fragments.tsv` 文本或路径、SAM 文本、数值列表或读段对；选项同 `calculate_insert_size`（`min_overlap`、`max_mismatch_rate`、`max_pairs`、`proper_pairs_only`、`max_insert_size`、`phred_offset`、`max_reads`）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `source`（`read-pair overlap`、`SAM TLEN`、`fragments (end - start)` 或 `numbers`）、`count`、`excluded_above_max`、`min`、`max`、`mean`、`median`、`stdev`、`median_absolute_deviation`、`mode`、`percentiles`（`p5`、`p25`、`p75`、`p95`）以及与来源相关的计数（`pairs`、`pairs_without_overlap`、`non_proper_pairs_skipped` 等）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_insert_size_tenx

fragments = "chr1\t100\t250\tAAAC\t2\nchr1\t100\t300\tAAAC\t1\n"
result = calculate_insert_size_tenx(fragments)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_insert_size_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 40,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5a53a6078f92a73b442f3887fee115319bb05572c40c635d00d176af98c64a00"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "fragments (end - start)",
    "count": 3,
    "excluded_above_max": 0,
    "min": 150,
    "max": 200,
    "mean": 166.666667,
    "median": 150.0,
    "stdev": 28.867513,
    "median_absolute_deviation": 0.0,
    "mode": 150,
    "percentiles": {
      "p5": 150.0,
      "p25": 150.0,
      "p75": 175.0,
      "p95": 195.0
    }
  }
}
```

## 23) calculate_phasing_qc

- **作用：** 原定从测序仪的 InterOp 运行指标计算 Illumina 的 phasing 与 pre-phasing 百分比。

- **处理方法：** 未实现。phasing 与 pre-phasing 是存储在测序仪 InterOp 二进制文件中的仪器端指标，无法由 FASTQ 读段推出。该函数始终返回 `status: "not_implemented"`，不做任何计算。

- **API：** `dnakit.sequencing_qc.calculate_phasing_qc(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 可为路径或任意对象，不参与任何计算；`options` 可省略，会被忽略。

- **输出：** 始终返回封装，含 `status: "not_implemented"`、`module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary`、`warnings`、`required_dependency`（"Illumina InterOp run-folder files (phasing/prephasing metrics)"）和说明原因的 `message`，没有 `result` 键。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_phasing_qc

result = calculate_phasing_qc("run_folder/InterOp")
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_phasing_qc",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6dff2776f630ad4ba95c41e73226a139851dee3b5bd6678bac6016b7c1ca929c"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Illumina InterOp run-folder files (phasing/prephasing metrics)",
  "message": "Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads."
}
```

## 24) calculate_phasing_qc_tenx

- **作用：** 原定从测序仪的 InterOp 运行指标计算 10x 文库测序运行的 Illumina phasing 与 pre-phasing 百分比。

- **处理方法：** 未实现。phasing 与 pre-phasing 是存储在测序仪 InterOp 二进制文件中的仪器端指标，无法由 FASTQ 读段推出。该函数始终返回 `status: "not_implemented"`，不做任何计算。

- **API：** `dnakit.sequencing_qc.calculate_phasing_qc_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 可为路径或任意对象，不参与任何计算；`options` 可省略，会被忽略。

- **输出：** 始终返回封装，含 `status: "not_implemented"`、`module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary`、`warnings`、`required_dependency`（"Illumina InterOp run-folder files (phasing/prephasing metrics)"）和说明原因的 `message`，没有 `result` 键。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_phasing_qc_tenx

result = calculate_phasing_qc_tenx("run_folder/InterOp")
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_phasing_qc_tenx",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 18,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6dff2776f630ad4ba95c41e73226a139851dee3b5bd6678bac6016b7c1ca929c"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "Illumina InterOp run-folder files (phasing/prephasing metrics)",
  "message": "Phasing and pre-phasing are instrument-side metrics stored in the sequencer's InterOp binaries and cannot be derived from FASTQ reads."
}
```

## 25) calculate_qscore_summary

- **作用：** 计算全局 Q 值汇总：平均和中位 Q、Q10-Q40 的碱基比例以及按读段 Q 阈值统计的读段数。

- **处理方法：** 碱基层面的统计量为 Phred Q 的平均值、中位数、最小值、最大值，平均错误概率对应的 Q，以及 Q 不低于 10、20、30、40 的碱基比例。读段的 Q 默认（`probability`）为其平均错误概率的 -10 log10（与 MinKNOW、NanoPlot 的报告方式一致），也可取碱基 Q 的算术平均（`arithmetic`）；并统计 Q 不低于 7、10、15、20、30 的读段数。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_qscore_summary(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。必须带有 Phred 质量值。选项：`read_q_method`（`probability` 默认或 `arithmetic`）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `bases`、`mean_q`、`q_of_mean_error_probability`、`median_q`、`min_q`、`max_q`、`fraction_bases_at_least`（`Q10`、`Q20`、`Q30`、`Q40`）、`read_q_method`、`reads`、`mean_read_q`、`median_read_q` 以及 `reads_with_mean_q_at_least`（`Q7`、`Q10`、`Q15`、`Q20`、`Q30`）。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_qscore_summary

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_qscore_summary(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_qscore_summary",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "bases": 40,
    "mean_q": 33.0,
    "q_of_mean_error_probability": 12.79262,
    "median_q": 40.0,
    "min_q": 0,
    "max_q": 40,
    "fraction_bases_at_least": {
      "Q10": 0.95,
      "Q20": 0.95,
      "Q30": 0.7,
      "Q40": 0.7
    },
    "read_q_method": "probability",
    "reads": 4,
    "mean_read_q": 26.746991,
    "median_read_q": 30.0,
    "reads_with_mean_q_at_least": {
      "Q7": 3,
      "Q10": 3,
      "Q15": 3,
      "Q20": 3,
      "Q30": 2
    }
  }
}
```

## 26) calculate_quality

- **作用：** 计算逐位置的碱基质量统计量（均值、中位数、四分位数、第 10/90 百分位数）以及每条读段平均质量的分布。

- **处理方法：** 对每个读段位置，用覆盖该位置的读段的 Phred Q 值给出均值、中位数、下四分位数和上四分位数、第 10 与第 90 百分位数、最小值和最大值。百分位数在顺序统计量之间线性插值，并四舍五入到 6 位小数，因此不会出现 6.000000000000001 这类浮点误差。每条读段平均质量为 Q 的算术平均（FastQC），汇总为均值、中位数、众数和整数值直方图。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_quality(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。必须带有 Phred 质量值。选项：`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads_with_quality`、`per_position`（每个位置的 `position`、`bases`、`mean`、`median`、`q1`、`q3`、`p10`、`p90`、`min`、`max`）以及 `per_read_mean_quality`（`mean`、`median`、`mode`、`distribution`）。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_quality

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_quality(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_quality",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_with_quality": 4,
    "per_position": [
      {
        "position": 1,
        "bases": 4,
        "mean": 25.0,
        "median": 30.0,
        "q1": 15.0,
        "q3": 40.0,
        "p10": 6.0,
        "p90": 40.0,
        "min": 0,
        "max": 40
      },
      {
        "position": 2,
        "bases": 4,
        "mean": 25.0,
        "median": 30.0,
        "q1": 15.0,
        "q3": 40.0,
        "p10": 6.0,
        "p90": 40.0,
        "min": 0,
        "max": 40
      },
      "..."
    ],
    "per_read_mean_quality": {
      "mean": 33.0,
      "median": 36.0,
      "mode": 40,
      "distribution": {
        "20": 1,
        "32": 1,
        "40": 2
      }
    }
  }
}
```

## 27) calculate_quality_tenx

- **作用：** 计算 10x 读段条形码、UMI 和 RNA 区域的逐位置质量统计量。

- **处理方法：** 把 read 1 拆为条形码和 UMI 区域，read 2 为 RNA 读段；每个区域给出 `calculate_quality` 的结果（逐位置均值、中位数、四分位数、第 10/90 百分位数以及每条读段平均质量分布）。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_quality_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。必须带有 Phred 质量值。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`，以及 `regions`（含 `barcode`、`umi`、`rna_read`，各自包含 `reads_with_quality`、`per_position` 与 `per_read_mean_quality`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_quality_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_quality_tenx(tenx)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_quality_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "reads_with_quality": 3,
        "per_position": [
          {
            "...": "10 keys"
          },
          {
            "...": "10 keys"
          },
          "..."
        ],
        "per_read_mean_quality": {
          "mean": 40.0,
          "median": 40.0,
          "mode": 40,
          "distribution": {
            "...": "1 keys"
          }
        }
      },
      "umi": {
        "reads_with_quality": 3,
        "per_position": [
          {
            "...": "10 keys"
          },
          {
            "...": "10 keys"
          },
          "..."
        ],
        "per_read_mean_quality": {
          "mean": 33.333333,
          "median": 40.0,
          "mode": 40,
          "distribution": {
            "...": "2 keys"
          }
        }
      },
      "rna_read": {
        "reads_with_quality": 3,
        "per_position": [
          {
            "...": "10 keys"
          },
          {
            "...": "10 keys"
          },
          "..."
        ],
        "per_read_mean_quality": {
          "mean": 33.333333,
          "median": 40.0,
          "mode": 40,
          "distribution": {
            "...": "2 keys"
          }
        }
      }
    }
  }
}
```

## 28) calculate_read_length

- **作用：** 计算读段长度分布：最小、最大、均值、中位数、标准差、N50、N90、L50 及直方图。

- **处理方法：** N50 为满足“长度不小于 L 的读段包含至少一半总碱基”的长度 L（N90 同理，为 90%），L50 为从最长读段起累加达到一半碱基所需的读段数。长度按 `bin_width` 分箱（默认 1，范围很宽时自动确定）。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_read_length(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。选项：`bin_width`（默认 1，范围很宽时自动确定）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`bases`、`min`、`max`、`mean`、`median`、`stdev`、`n50`、`n90`、`l50`、`bin_width` 以及 `histogram`（`start`、`end`、`count`）。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_read_length

reads = ["A" * n for n in (100, 200, 300, 400, 500)]
result = calculate_read_length(reads, options={"bin_width": 100})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_read_length",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 1520,
    "sequence_bases": 1500,
    "lengths": [100, 200, 300, 400, 500],
    "sha256": "0212d55422b398acab34c76a9e36b0a6991a4cdefef16a25c8f9b01bea4e23e9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 5,
    "bases": 1500,
    "min": 100,
    "max": 500,
    "mean": 300.0,
    "median": 300.0,
    "stdev": 158.113883,
    "n50": 400,
    "n90": 200,
    "l50": 2,
    "bin_width": 100,
    "histogram": [
      {
        "start": 100,
        "end": 199,
        "count": 1
      },
      {
        "start": 200,
        "end": 299,
        "count": 1
      },
      "..."
    ]
  }
}
```

## 29) calculate_read_length_tenx

- **作用：** 计算 10x 数据 read 1 与 read 2 的读段长度分布。

- **处理方法：** 对每条读段（read 1、read 2）分别给出 `calculate_read_length` 的统计量（最小、最大、均值、中位数、标准差、N50、N90、L50、直方图）。

- **API：** `dnakit.sequencing_qc.calculate_read_length_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为读段对 `{"read1": ..., "read2": ...}`（FASTQ 文本、路径或读段列表）。选项：`bin_width`、`phred_offset`、`max_reads`。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `read_pairs`，以及 `read1` 与 `read2`，各自包含 `calculate_read_length` 的字段。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_read_length_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_read_length_tenx(tenx)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_read_length_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "read_pairs": 3,
    "read1": {
      "reads": 3,
      "bases": 84,
      "min": 28,
      "max": 28,
      "mean": 28.0,
      "median": 28.0,
      "stdev": 0.0,
      "n50": 28,
      "n90": 28,
      "l50": 2,
      "bin_width": 1,
      "histogram": [
        {
          "start": 28,
          "end": 28,
          "count": 3
        }
      ]
    },
    "read2": {
      "reads": 3,
      "bases": 90,
      "min": 30,
      "max": 30,
      "mean": 30.0,
      "median": 30.0,
      "stdev": 0.0,
      "n50": 30,
      "n90": 30,
      "l50": 2,
      "bin_width": 1,
      "histogram": [
        {
          "start": 30,
          "end": 30,
          "count": 3
        }
      ]
    }
  }
}
```

## 30) calculate_yield

- **作用：** 计算测序产量：读段数、碱基数、平均长度、N 碱基数、GC 比例以及 Q20/Q30 碱基数。

- **处理方法：** 统计读段数和碱基数、平均读长、`N` 碱基数、GC 比例（G+C 占 A+C+G+T 的比例），有质量值时还统计 Q20、Q30 及以上的碱基数及其比例。给出 `genome_size` 时，期望覆盖度为碱基数除以基因组大小。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.calculate_yield(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。选项：`genome_size`（可选；增加 `expected_coverage`）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`bases`、`mean_read_length`、`n_bases`、`gc_fraction`、`q20_bases`、`q30_bases`、`q20_fraction`、`q30_fraction`（无质量值时为 `null`），以及（给出 `genome_size` 时）`genome_size` 和 `expected_coverage`。输入为读段对时，结果含 `paired: true`、`read_pairs`、`total_bases`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_yield

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = calculate_yield(reads, options={"genome_size": 20})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_yield",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 4,
    "bases": 40,
    "mean_read_length": 10.0,
    "n_bases": 2,
    "gc_fraction": 0.5,
    "q20_bases": 38,
    "q30_bases": 28,
    "q20_fraction": 0.95,
    "q30_fraction": 0.7,
    "genome_size": 20.0,
    "expected_coverage": 2.0
  }
}
```

## 31) calculate_yield_tenx

- **作用：** 按区域（条形码、UMI、RNA 读段）计算 10x 读段的产量，并统计条形码数量。

- **处理方法：** 把 read 1 拆为条形码和 UMI 区域，read 2 为 RNA 读段；每个区域给出 `calculate_yield` 的产量字段。此外给出不同细胞条形码数和每个条形码的中位读段数；给出 `whitelist` 时还给出有效条形码比例（与白名单完全相同，或在 `max_mismatches` 内有唯一的白名单近邻）。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.calculate_yield_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有选项：`whitelist`（有效细胞条形码）、`max_mismatches`（默认 1）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`regions`（`barcode`、`umi`、`rna_read`，各自包含 `calculate_yield` 的字段）、`total_bases_read1_read2`、`distinct_barcodes`、`median_reads_per_barcode`，以及（给出白名单时）`valid_barcode_fraction` 和 `valid_barcode_definition`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import calculate_yield_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = calculate_yield_tenx(tenx, options={"whitelist": [barcode]})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "calculate_yield_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "regions": {
      "barcode": {
        "reads": 3,
        "bases": 48,
        "mean_read_length": 16.0,
        "n_bases": 0,
        "gc_fraction": 0.333333,
        "q20_bases": 48,
        "q30_bases": 48,
        "q20_fraction": 1.0,
        "q30_fraction": 1.0
      },
      "umi": {
        "reads": 3,
        "bases": 36,
        "mean_read_length": 12.0,
        "n_bases": 0,
        "gc_fraction": 0.333333,
        "q20_bases": 36,
        "q30_bases": 24,
        "q20_fraction": 1.0,
        "q30_fraction": 0.666667
      },
      "rna_read": {
        "reads": 3,
        "bases": 90,
        "mean_read_length": 30.0,
        "n_bases": 0,
        "gc_fraction": 0.688889,
        "q20_bases": 90,
        "q30_bases": 60,
        "q20_fraction": 1.0,
        "q30_fraction": 0.666667
      }
    },
    "total_bases_read1_read2": 174,
    "distinct_barcodes": 2,
    "median_reads_per_barcode": 1.5,
    "valid_barcode_fraction": 0.666667,
    "valid_barcode_definition": "Barcode equals a whitelist entry or has a unique whitelist neighbour within 1 mismatches."
  }
}
```

## 32) collapse_umi_reads

- **作用：** 将共享同一 UMI（及读段起始）的读段折叠为每个分子一条一致序列读段。

- **处理方法：** 由于没有比对位置，读段按 (UMI, 前 `start_length` 个碱基) 分组，读段起始只是比对位置的近似；`group_by` 为 `umi` 时仅按 UMI 分组。启用 `correct_errors` 时，先用 `correct_umi_errors` 的 directional 网络合并 UMI。一致序列在每个位置按 Phred 质量加权投票，只取被至少一半组内读段覆盖的位置。一致序列质量 = clip(一致碱基 Q 之和 - 不一致碱基 Q 之和, 0, 93)。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.collapse_umi_reads(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为读段名以 `<umi_separator><UMI>` 结尾的 FASTQ 读段（umi_tools 约定），例如 `read1_ACGTACGT`；格式可为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列）。不接受读段对映射（抛出 `INVALID_INPUT`）；双端数据请对带 UMI 读段名的那条读段做合并。选项：`umi_separator`（默认 `_`）、`group_by`（`umi_and_start` 默认或 `umi`）、`start_length`（默认 10）、`correct_errors`（默认 `False`）、`threshold`（校正用的汉明距离，默认 1）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads_in`、`molecules_out`、`collapse_ratio`、`group_by`、`start_length`、`error_corrected`、`groups`（每组的 `umi`、`start`、`reads`、`read_names`、`consensus_sequence`、`consensus_quality`）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import collapse_umi_reads

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1_AAAA", "ACGTACGT", "IIIIIIII"),
    ("r2_AAAA", "ACGTACGT", "IIIIIIII"),
    ("r3_AAAA", "ACGAACGT", "III5IIII"),
    ("r4_CCCC", "ACGTACGT", "IIIIIIII"),
)
result = collapse_umi_reads(reads, options={"group_by": "umi"})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "collapse_umi_reads",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 116,
    "sequence_bases": 32,
    "lengths": [8, 8, 8, 8],
    "sha256": "e7ea522ad8a504ec9f45d37aebd7c1e1e0964c3e20b71854138375e0984a97f9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads_in": 4,
    "molecules_out": 2,
    "collapse_ratio": 2.0,
    "group_by": "umi",
    "start_length": null,
    "error_corrected": false,
    "groups": [
      {
        "umi": "AAAA",
        "start": null,
        "reads": 3,
        "read_names": [
          "r1_AAAA",
          "r2_AAAA",
          "..."
        ],
        "consensus_sequence": "ACGTACGT",
        "consensus_quality": "~~~]~~~~"
      },
      {
        "umi": "CCCC",
        "start": null,
        "reads": 1,
        "read_names": [
          "r4_CCCC"
        ],
        "consensus_sequence": "ACGTACGT",
        "consensus_quality": "IIIIIIII"
      }
    ],
    "definition": "Reads are grouped by (UMI, first start_length bases) because no alignment position is available; the consensus is a per-position vote weighted by Phred quality ..."
  }
}
```

## 33) correct_long_reads

- **作用：** 原定通过基于重叠的一致序列（canu/ratatosk 风格）校正长读段中的测序错误。

- **处理方法：** 未实现。长读段自校正需要全对全重叠检测和一致序列生成，DNAKit 未原生实现这些功能。该函数始终返回 `status: "not_implemented"`，不做任何计算。

- **API：** `dnakit.sequencing_qc.correct_long_reads(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 可为路径或任意对象，不参与任何计算；`options` 可省略，会被忽略。

- **输出：** 始终返回封装，含 `status: "not_implemented"`、`module`、`operation`、`provenance`（`implementation: "not_implemented"`）、`input_summary`、`warnings`、`required_dependency`（"long-read error-correction program (canu, ratatosk or similar)"）和说明原因的 `message`，没有 `result` 键。

- **示例代码：**

```python
from dnakit.sequencing_qc import correct_long_reads

result = correct_long_reads("long_reads.fastq")
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "correct_long_reads",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 16,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b8567bebd311c7484fd0d8dbde61c467727df0b1d64337a7ef96d0dac8922bf8"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "long-read error-correction program (canu, ratatosk or similar)",
  "message": "Long-read self-correction needs all-vs-all overlap detection and consensus calling that DNAKit does not implement natively."
}
```

## 34) correct_umi_errors

- **作用：** 在汉明距离网络上对 UMI 聚类，校正 UMI 的测序错误（UMI-tools 方法）。

- **处理方法：** 汉明距离不超过 `threshold` 的 UMI 相互连接，并按 UMI-tools 的三种方法之一分组。`directional`（默认）：当 count(a) >= 2 x count(b) - 1 时存在边 a -> b，从丰度最高的 UMI 起做广度优先搜索扩展分组。`adjacency`：丰度最高的 UMI 加上其邻居，在每个连通分量内重复。`cluster`：连通分量。每组中丰度最高的 UMI 为代表，其他成员被校正为它。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.correct_umi_errors(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 UMI 字符串列表、UMI -> 读段数的映射，或读段名以 `_<UMI>` 结尾的读段（最多 5000 个不同 UMI）。选项：`method`（`directional` 默认、`adjacency` 或 `cluster`）、`threshold`（汉明距离，默认 1）、`umi_separator`（默认 `_`）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `method`、`threshold`、`distinct_umis_before`、`distinct_umis_after`、`reads`、`groups`（每组的 `representative`、`members`、`reads`、`member_counts`）、`corrections`（被校正的 UMI -> 代表 UMI）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import correct_umi_errors

counts = {"AAAAAA": 100, "AAAAAT": 1, "CCCCCC": 50, "CCCCCG": 30}
result = correct_umi_errors(counts)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "correct_umi_errors",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3fefe1abde49ead8296fe12856ea08439dda7deb9e9f8853bab0cb087b034aba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "directional",
    "threshold": 1,
    "distinct_umis_before": 4,
    "distinct_umis_after": 3,
    "reads": 181,
    "groups": [
      {
        "representative": "AAAAAA",
        "members": [
          "AAAAAT"
        ],
        "reads": 101,
        "member_counts": {
          "AAAAAA": 100,
          "AAAAAT": 1
        }
      },
      {
        "representative": "CCCCCC",
        "members": [],
        "reads": 50,
        "member_counts": {
          "CCCCCC": 50
        }
      },
      "..."
    ],
    "corrections": {
      "AAAAAT": "AAAAAA"
    },
    "definition": "UMI-tools networks over Hamming distance <= threshold. directional: edge a->b when count(a) >= 2*count(b)-1, groups grown by breadth-first search from the most ..."
  }
}
```

## 35) demultiplex_barcoded_reads

- **作用：** 按条形码把读段拆分为各样本的读段集合，条形码属于读段时将其剪除。

- **处理方法：** 条形码的分配方式同 `assign_sample_barcodes`（汉明距离，受 `max_mismatches` 和 `min_distance_gap` 约束）。已分配的读段按样本返回；对 `read_start` 与 `read_end` 条形码，`trim_barcode` 为真时剪除条形码碱基（及其质量值）。没有唯一匹配的读段归入 `unassigned`，不做剪除。给出 `output_directory` 时，每个样本写出一个 `<sample>.fastq`，且不会覆盖已有文件。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.demultiplex_barcoded_reads(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 与选项同 `assign_sample_barcodes`（`samples` 必填，另有 `barcode_location`、`barcode_length`、`max_mismatches`、`min_distance_gap`；不接受读段对映射，会抛出 `INVALID_INPUT`；双端数据请对带条形码的那条读段做拆分），另有 `trim_barcode`（`read_start` 与 `read_end` 时默认 `True`）、`output_directory`（可选，须为已存在的目录）和 `max_listed`（每个样本返回的读段数，默认 1000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `assign_sample_barcodes` 的分配汇总（`reads`、`reads_per_sample`、`assigned`、`unassigned`、`ambiguous`、`assigned_fraction` 等）、`barcode_trimmed`、`samples`（每个样本及 `unassigned`：`reads` 和含 `id`、`sequence`、`quality` 的 `records`）以及 `written_files`（未给出 `output_directory` 时为 `null`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import demultiplex_barcoded_reads

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("a", "ACGTAAAA", "IIIIJJJJ"),
    ("b", "TTTTCCCC", "IIIIJJJJ"),
    ("c", "GGGGCCCC", "IIIIJJJJ"),
)
result = demultiplex_barcoded_reads(reads, options={"samples": {"S1": "ACGT", "S2": "TTTT"}, "barcode_length": 4})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "demultiplex_barcoded_reads",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 69,
    "sequence_bases": 24,
    "lengths": [8, 8, 8],
    "sha256": "554e09db1fbcb9792fc57c3561da72f69cb5cf168c635f51759ab0ce0ac285a1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 3,
    "barcode_location": "read_start",
    "max_mismatches": 1,
    "min_distance_gap": 1,
    "reads_per_sample": {
      "S1": 1,
      "S2": 1
    },
    "assigned": 2,
    "unassigned": 1,
    "ambiguous": 0,
    "assigned_fraction": 0.666667,
    "min_distance_between_sample_barcodes": 3,
    "barcode_trimmed": true,
    "samples": {
      "S1": {
        "reads": 1,
        "records": [
          {
            "...": "3 keys"
          }
        ]
      },
      "S2": {
        "reads": 1,
        "records": [
          {
            "...": "3 keys"
          }
        ]
      },
      "unassigned": {
        "reads": 1,
        "records": [
          {
            "...": "3 keys"
          }
        ]
      }
    },
    "written_files": null
  }
}
```

## 36) detect_cross_sample_contamination

- **作用：** 标记共享 k-mer 比例异常高的样本对（可能的样本间交叉污染）。

- **处理方法：** 为每个样本建立规范化 k-mer 集合，对每个有序样本对计算 containment[a][b] = |K(a) 与 K(b) 的交集| / |K(a)|，containment 不低于 `threshold` 的样本对被标记。启发式：仅当样本预期互不相关（不同物种或不同个体）时才有意义；相关样本本来就会共享 k-mer。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.detect_cross_sample_contamination(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为样本名 -> 读段的映射（至少两个样本；每个样本的格式可为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列））。选项：`kmer_size`（默认 21）、`threshold`（containment，取值 [0, 1]，默认 0.1）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `samples`（样本 -> 不同 k-mer 数）、`kmer_size`、`threshold`、`containment`（嵌套映射 a -> b -> 比例）、`flagged_pairs`（`sample`、`found_in`、`containment`，按 containment 降序）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import detect_cross_sample_contamination

host = "ACGTTGCAAGGCTTAACCGGATCGATTCGAGGCTATTGCCAGTCAAGGTACCTGATGCCGTAAGCTTGGAACCTAGGTCCA"
other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
samples = {"A": [host[:40]], "B": [host[20:60]], "C": [other[:40]]}
result = detect_cross_sample_contamination(samples, options={"kmer_size": 12, "threshold": 0.2})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "detect_cross_sample_contamination",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 153,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8435281fcc502972cb1960fbe982a9f256340d91c17229f1e34dadc196a800c1"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": {
      "A": 29,
      "B": 29,
      "C": 29
    },
    "kmer_size": 12,
    "threshold": 0.2,
    "containment": {
      "A": {
        "B": 0.310345,
        "C": 0.0
      },
      "B": {
        "A": 0.310345,
        "C": 0.0
      },
      "C": {
        "A": 0.0,
        "B": 0.0
      }
    },
    "flagged_pairs": [
      {
        "sample": "A",
        "found_in": "B",
        "containment": 0.310345
      },
      {
        "sample": "B",
        "found_in": "A",
        "containment": 0.310345
      }
    ],
    "definition": "containment[a][b] = |K(a) intersect K(b)| / |K(a)| over canonical k-mers. Heuristic: samples that are expected to be biologically unrelated should share almost ..."
  }
}
```

## 37) detect_platform

- **作用：** 根据读段名格式（基于规则）判断测序平台。

- **处理方法：** 将每条读段名与已知格式匹配：Illumina（instrument:run:flowcell:lane:tile:x:y）、Oxford Nanopore（`runid=`、`ch=`、UUID）、PacBio（movie/zmw/ccs）、MGI/BGI（L#C#R#）和 Ion Torrent；经 SRA 重命名的读段不含平台信息。当多数格式覆盖至少一半读段时，取其为平台，否则为 `unknown`。不推断仪器型号。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.detect_platform(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。选项：`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `platform`（`illumina`、`oxford_nanopore`、`pacbio`、`mgi_dnbseq`、`ion_torrent`、`sra_renamed` 或 `unknown`）、`confidence_fraction`、`read_name_votes`、`read_length`（`min`、`median`、`max`、`uniform`）、`evidence`（Illumina 时为 `instrument_id`、`run_number`、`flowcell`、`lane`）以及 `definition` 说明。输入为读段对时，结果含 `paired: true`，上述字段分别放在 `read1` 和 `read2` 下。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import detect_platform

reads = "@A00123:45:HXYZ123DSXX:1:1101:1000:2000 1:N:0:ACGT\nACGT\n+\nIIII\n"
result = detect_platform(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "detect_platform",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 63,
    "sequence_bases": 4,
    "lengths": [4],
    "sha256": "104d5bfb8cfd56a46818297271296941ac13efeb1f9873b9e9ad65abc7e35c8e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "platform": "illumina",
    "confidence_fraction": 1.0,
    "read_name_votes": {
      "illumina": 1
    },
    "read_length": {
      "min": 4,
      "median": 4.0,
      "max": 4,
      "uniform": true
    },
    "evidence": {
      "instrument_id": "A00123",
      "run_number": "45",
      "flowcell": "HXYZ123DSXX",
      "lane": "1"
    },
    "definition": "Rule-based: the read-name (header) format of each read is matched against known layouts (Illumina instrument:run:flowcell:lane:tile:x:y, Oxford Nanopore runid=/..."
  }
}
```

## 38) detect_platform_tenx

- **作用：** 根据 read 2 的读段名判断 10x 读段的测序平台，并由 read 1 长度推断 Chromium 条形码+UMI 布局。

- **处理方法：** 平台判断方式同 `detect_platform`，依据 read 2 的读段名。只使用 read 1 的众数长度（read 1 序列本身不与化学版本列表或白名单比较）：24、26 或 28 分别对应 14+10、16+10 或 16+12 bp 的条形码+UMI 布局及兼容的 Chromium 化学版本。这只是长度启发式：不执行 Cell Ranger 基于白名单的化学版本自动判断，其他 read 1 长度则不给出布局。

- **API：** `dnakit.sequencing_qc.detect_platform_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为读段对 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错。选项：`phred_offset`（33）、`max_reads`（1000000）。其他 10x 函数的化学版本选项在此不使用：布局由 read 1 长度推断。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `detect_platform` 对 read 2 给出的字段（`platform`、`confidence_fraction`、`read_name_votes`、`read_length`、`evidence`、`definition`）以及 `tenx`（`read1_modal_length`、`layout`、`compatible_chemistries`、`note`；长度未知时 `layout` 与 `compatible_chemistries` 为 `null`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import detect_platform_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

name = "A00123:45:HXYZ123DSXX:1:1101:1000:2000"
tenx = {
    "read1": fastq((name, "ACGTACGTACGTACGT" + "A" * 12, "I" * 28)),
    "read2": fastq((name, "ACGT" * 5, "I" * 20)),
}
result = detect_platform_tenx(tenx)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "detect_platform_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 218,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1b27ae92e3ff96f5c0a9fe0740c51141c29144ec1cc927301189bc3b8dc3d586"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "platform": "illumina",
    "confidence_fraction": 1.0,
    "read_name_votes": {
      "illumina": 1
    },
    "read_length": {
      "min": 20,
      "median": 20.0,
      "max": 20,
      "uniform": true
    },
    "evidence": {
      "instrument_id": "A00123",
      "run_number": "45",
      "flowcell": "HXYZ123DSXX",
      "lane": "1"
    },
    "definition": "Rule-based: the read-name (header) format of each read is matched against known layouts (Illumina instrument:run:flowcell:lane:tile:x:y, Oxford Nanopore runid=/...",
    "tenx": {
      "read1_modal_length": 28,
      "layout": "16 bp barcode + 12 bp UMI",
      "compatible_chemistries": "Chromium 3' v3 / v3.1 / 5' HT",
      "note": "Chemistry is inferred only from the modal read-1 length (24/26/28). Cell Ranger's whitelist-based autodetection needs the 10x barcode whitelists and is not perf..."
    }
  }
}
```

## 39) estimate_library_complexity

- **作用：** 用 Lander-Waterman / Picard 模型估计文库复杂度（不同分子的数量）。

- **处理方法：** 设总读段数为 N、不同读段数为 C，文库大小 X 满足 C / X = 1 - exp(-N / X)（Picard EstimateLibraryComplexity），数值求解；并预测在当前测序深度的 0.5x、1x、2x、5x、10x 下的期望唯一读段数 X (1 - exp(-k N / X))。所有读段都不同时 X 无界（`null`）。给出读段时，用完全序列重复定义不同读段集合（读段对使用两条读段）。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.estimate_library_complexity(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为读段（FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或读段对）或映射 `{"total_reads": N, "unique_reads": C}`（也可放在 `options` 中）。选项：`prefix_length`、`reverse_complement_aware`、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `source`（`counts` 或读段）、`total_reads`、`unique_reads`、`duplicate_rate`、`estimated_library_size`（无界时为 `null`）、`predicted_unique_reads`（`0.5x`、`1x`、`2x`、`5x`、`10x`）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import estimate_library_complexity

result = estimate_library_complexity({"total_reads": 1_000_000, "unique_reads": 632_121})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "estimate_library_complexity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "34f4ef63ded2460579d54590e66569cec3c6acf812aa7ca34e25176a55cfe90a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "counts",
    "total_reads": 1000000,
    "unique_reads": 632121,
    "duplicate_rate": 0.367879,
    "estimated_library_size": 1000001.669581,
    "predicted_unique_reads": {
      "0.5x": 393469.49089,
      "1x": 632121.0,
      "2x": 864665.708484,
      "5x": 993263.655084,
      "10x": 999956.268817
    },
    "definition": "Lander-Waterman / Picard EstimateLibraryComplexity: find X with unique/X = 1 - exp(-total/X)."
  }
}
```

## 40) estimate_read_error_model

- **作用：** 估计读段错误模型：相对参考序列的观测错配率，或由 Phred 质量值推导的期望错误率曲线。

- **处理方法：** 给出 `reference` 时，读段通过 k-mer 种子投票（两条链，取前 3 条对角线）定位到参考序列上，并以无空位汉明距离打分；错配比例超过 `max_mismatch_fraction` 的读段不计入。结果给出总体与逐位置错配率、替换计数，以及（有质量值时）每个报告 Q 值对应的经验错误率。不对插入缺失建模，因此富含插入缺失的读段会被排除或抬高替换率。没有参考序列时，模型为 `calculate_error_profile` 的质量值推导曲线（须有质量值）。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.estimate_read_error_model(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列）。不接受读段对映射（抛出 `INVALID_INPUT`）；双端数据请对每条读段分别调用一次。选项：`reference`（序列或 FASTA 文本，最长 5 Mb；不给时模型为质量值推导曲线）、`seed_length`（默认 11）、`max_mismatch_fraction`（默认 0.2）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `model`：给出参考序列时含 `source`（`ungapped_reference_alignment`）、`reads`、`reads_aligned`、`bases_compared`、`mismatches`、`overall_mismatch_rate`、`per_position_mismatch_rate`、`substitution_counts`（如 `A>T`）、`empirical_error_by_reported_q` 和 `definition`；没有参考序列时为 `calculate_error_profile` 的字段加 `note`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import estimate_read_error_model

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

host = "ACGTTGCAAGGCTTAACCGGATCGATTCGAGGCTATTGCCAGTCAAGGTACCTGATGCCGTAAGCTTGGAACCTAGGTCCA"
other = "TTGACCGGTAACCTTGGAACCGTTAGGCCATTGACCAAGGTTCCAAGGTTCCGGAATTCCGGAACCTTGGAAGGCC"
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

mutate = {"A": "C", "C": "G", "G": "T", "T": "A"}
read1 = host[:10] + mutate[host[10]] + host[11:30]
read2 = host[20:50]
original = reverse_complement(host[40:70])
read3 = mutate[original[0]] + original[1:]
reads = fastq(("r1", read1, "I" * 30), ("r2", read2, "5" * 30), ("r3", read3, "I" * 30))
result = estimate_read_error_model(reads, options={"reference": host})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "estimate_read_error_model",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 204,
    "sequence_bases": 90,
    "lengths": [30, 30, 30],
    "sha256": "7a26ebdc2713a3ab7bca25d5b393a92b5aa85a9975f53b35f4538570801cd5bc"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": {
      "source": "ungapped_reference_alignment",
      "reads": 3,
      "reads_aligned": 3,
      "bases_compared": 90,
      "mismatches": 2,
      "overall_mismatch_rate": 0.022222,
      "per_position_mismatch_rate": [
        {
          "position": 1,
          "bases": 3,
          "mismatch_rate": 0.333333
        },
        {
          "position": 2,
          "bases": 3,
          "mismatch_rate": 0.0
        },
        "..."
      ],
      "substitution_counts": {
        "A>T": 1,
        "G>T": 1
      },
      "empirical_error_by_reported_q": {
        "20": {
          "bases": 30,
          "empirical_error_rate": 0.0,
          "reported_error_rate": 0.01
        },
        "40": {
          "bases": 60,
          "empirical_error_rate": 0.033333,
          "reported_error_rate": 0.0001
        }
      },
      "definition": "Each read is placed on the reference by k-mer seed voting (both strands, top 3 diagonals) and scored by Hamming distance without gaps; reads above max_mismatch_..."
    }
  }
}
```

## 41) extract_umi_sequences

- **作用：** 从读段中提取 UMI 序列并追加到读段名（umi_tools 约定）。

- **处理方法：** 从每条读段取前（`5prime`）或后（`3prime`）`umi_length` 个碱基，以 `<name><umi_separator><UMI>` 的形式追加到读段名。给出 `pattern`（仅限 5' 端）时按碱基逐位指定：`N` 为 UMI 碱基，`C` 为细胞条形码碱基，`X` 为丢弃碱基，如 `CCCCCCNNNNXX`；读段其余部分保留。返回剪除后的序列和质量值。长度不超过提取区域的读段会被拒绝。纯 Python 实现。

- **API：** `dnakit.sequencing_qc.extract_umi_sequences(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列）。不接受读段对映射（抛出 `INVALID_INPUT`）；双端数据请只传入带 UMI 的那条读段（通常是 read 1）。选项：`umi_length`（默认 8）与 `umi_location`（`5prime` 默认或 `3prime`），或 `pattern`（仅 5'）、`umi_separator`（默认 `_`）、`append_to_name`（默认 `True`）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`umi_location`、`umi_region_length`、`distinct_umis`、`records`（每条读段的 `id`、`umi`、`cell_barcode`、`sequence`、`quality`）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import extract_umi_sequences

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(("r1", "ACGTAAAATTTT", "ABCDEFGHIJKL"))
result = extract_umi_sequences(reads, options={"umi_length": 4})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "extract_umi_sequences",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 32,
    "sequence_bases": 12,
    "lengths": [12],
    "sha256": "b3a517a420b416cdb6f02aafbee98222c72950c3c2b06697d60edf9b11a7e7fd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 1,
    "umi_location": "5prime",
    "umi_region_length": 4,
    "distinct_umis": 1,
    "records": [
      {
        "id": "r1_ACGT",
        "umi": "ACGT",
        "cell_barcode": null,
        "sequence": "AAAATTTT",
        "quality": "EFGHIJKL"
      }
    ],
    "definition": "The first (5prime) or last (3prime) umi_length bases, or the N/C/X pattern positions at the 5' end, are removed from each read and appended to its name as '<nam..."
  }
}
```

## 42) generate_qc_report

- **作用：** 生成 FastQC 风格的质控报告：产量、长度、质量、GC、接头、重复和平台各部分，并带有 pass/warn/fail 标志。

- **处理方法：** `reads` 部分以精简形式汇集 `calculate_yield`、`calculate_read_length`、`detect_platform`、`calculate_gc_bias`、`calculate_adapter_content`、`calculate_duplication`、`calculate_qscore_summary` 和 `calculate_quality` 的结果。六个标志遵循结果中 `flag_definitions` 列出的 FastQC 阈值：逐碱基质量、逐序列质量、逐碱基 N 含量、接头含量、序列重复和 GC 含量分布形状（无法计算时为 `not_available`）。同时生成纯文本报告（`report_text`），仅在给出 `output_path` 时写入文件。纯 Python 实现，不调用 FastQC。

- **API：** `dnakit.sequencing_qc.generate_qc_report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或一对读段 `{"read1": ..., "read2": ...}`（按每条读段分别分析）。选项：`adapters`（同 `calculate_adapter_content`）、`output_path`（可选的 `.txt` 或 `.md` 文件；已存在的文件除非 `overwrite=True` 否则不会被覆盖）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `paired`、`sections`（`reads`；读段对时为 `read1` 与 `read2`；每部分含 `yield`、`read_length`、`platform`、`gc_content`、`adapter_content`、`duplication`、`max_n_percent_at_a_position`、`qscore`、`per_base_quality`、`per_read_mean_quality` 和 `flags`）、`pair_duplication`（仅读段对）、`flag_definitions`、`report_text` 以及 `output_path`（未写文件时为 `null`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import generate_qc_report

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = generate_qc_report(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "generate_qc_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "paired": false,
    "sections": {
      "reads": {
        "yield": {
          "reads": 4,
          "bases": 40,
          "mean_read_length": 10.0,
          "n_bases": 2,
          "gc_fraction": 0.5,
          "q20_bases": 38,
          "q30_bases": 28,
          "q20_fraction": 0.95,
          "q30_fraction": 0.7
        },
        "read_length": {
          "reads": 4,
          "bases": 40,
          "min": 10,
          "max": 10,
          "mean": 10.0,
          "median": 10.0,
          "stdev": 0.0,
          "n50": 10,
          "n90": 10,
          "l50": 2,
          "bin_width": 1
        },
        "platform": {
          "platform": "unknown",
          "confidence_fraction": 1.0,
          "read_name_votes": {
            "...": "1 keys"
          },
          "read_length": {
            "...": "4 keys"
          },
          "evidence": {
            "...": "0 keys"
          }
        },
        "gc_content": {
          "mean_gc_percent": 50.0,
          "stdev_gc_percent": 0.0,
          "mode_gc_percent": 50,
          "deviation_from_normal_percent": null
        },
        "adapter_content": {
          "fraction_reads_with_any_adapter": 0.0,
          "max_cumulative_percent": {
            "...": "5 keys"
          }
        },
        "duplication": {
          "distinct_sequences": 3,
          "percent_remaining_after_deduplication": 75.0
        },
        "max_n_percent_at_a_position": 25.0,
        "qscore": {
          "bases": 40,
          "mean_q": 33.0,
          "q_of_mean_error_probability": 12.79262,
          "median_q": 40.0,
          "min_q": 0,
          "max_q": 40,
          "fraction_bases_at_least": {
            "...": "4 keys"
          },
          "read_q_method": "probability",
          "reads": 4,
          "mean_read_q": 26.746991,
          "median_read_q": 30.0,
          "reads_with_mean_q_at_least": {
            "...": "5 keys"
          }
        },
        "per_base_quality": {
          "positions": 10,
          "min_median": 30.0,
          "min_lower_quartile": 15.0
        },
        "per_read_mean_quality": {
          "mean": 33.0,
          "median": 36.0,
          "mode": 40
        },
        "flags": {
          "per_base_quality": "pass",
          "per_sequence_quality": "pass",
          "per_base_n_content": "fail",
          "adapter_content": "pass",
          "sequence_duplication": "warn",
          "gc_content_shape": "not_available"
        }
      }
    },
    "flag_definitions": "FastQC thresholds: per-base quality warn if lower quartile < 10 or median < 25 (fail < 5 / < 20) at any position; per-sequence quality warn if modal read mean Q...",
    "report_text": "DNAKit sequencing QC report\n===========================\n\n[reads]\nreads: 4  bases: 40  mean length: 10.0\nQ20 fraction: 0.95  Q30 fraction: 0.7\nmean GC%: 50.0  pl...",
    "output_path": null
  }
}
```

## 43) generate_qc_report_tenx

- **作用：** 为 10x Chromium 读段生成质控报告：各区域的产量与碱基识别、条形码/UMI 有效性、重复率替代指标以及 RNA 读段的 FastQC 风格标志。

- **处理方法：** 条形码、UMI 和 RNA 读段区域给出 `calculate_yield` 与 `calculate_basecalling_qc` 的产量和碱基识别指标；RNA 读段还给出 `generate_qc_report` 的精简 FastQC 风格部分和标志（有 Phred 质量值时才有标志）。报告另含 `calculate_barcode_quality_tenx` 的条形码/UMI 报告、`calculate_duplication_tenx` 的重复率替代指标以及 `calculate_yield_tenx` 的产量汇总。返回文本报告，仅在给出 `output_path` 时写入文件。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.generate_qc_report_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有选项：`whitelist`、`adapters`、`output_path`（可选的 `.txt` 或 `.md`；除非 `overwrite=True` 否则不覆盖）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`sections`（`barcode`、`umi`、`rna_read`，含 `yield` 与 `basecalling`；`rna_read` 另有 `compact`）、`barcode_and_umi`、`duplication`、`yield`、`report_text` 以及 `output_path`（未写文件时为 `null`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import generate_qc_report_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = generate_qc_report_tenx(tenx)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "generate_qc_report_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "sections": {
      "barcode": {
        "yield": {
          "...": "9 keys"
        },
        "basecalling": {
          "...": "12 keys"
        }
      },
      "umi": {
        "yield": {
          "...": "9 keys"
        },
        "basecalling": {
          "...": "12 keys"
        }
      },
      "rna_read": {
        "yield": {
          "...": "9 keys"
        },
        "basecalling": {
          "...": "12 keys"
        },
        "compact": {
          "...": "11 keys"
        }
      }
    },
    "barcode_and_umi": {
      "chemistry": "v3",
      "barcode_length": 16,
      "umi_length": 12,
      "read_pairs": 3,
      "cell_barcode": {
        "reads": 3,
        "barcode_source": "read_start",
        "distinct_barcodes": 2,
        "reads_with_N_in_barcode": 0,
        "top_barcodes": [
          {
            "...": "2 keys"
          },
          {
            "...": "2 keys"
          }
        ],
        "barcode_entropy_bits": 0.918296,
        "mean_quality": 40.0,
        "q30_fraction": 1.0,
        "mean_quality_by_position": [
          40.0,
          40.0,
          "..."
        ]
      },
      "umi": {
        "distinct_umis": 2,
        "reads_with_N_in_umi": 0,
        "mean_quality": 33.333333,
        "q30_fraction": 0.666667,
        "homopolymer_umi_reads": 3
      }
    },
    "duplication": {
      "distinct_sequences": 2,
      "definition": "Reads are duplicates when cell barcode, UMI and the first 25 bases of read 2 are identical. This is a proxy: Cell Ranger's sequencing saturation uses (barcode, ...",
      "duplicate_rate": 0.333333
    },
    "yield": {
      "chemistry": "v3",
      "barcode_length": 16,
      "umi_length": 12,
      "read_pairs": 3,
      "total_bases_read1_read2": 174,
      "distinct_barcodes": 2,
      "median_reads_per_barcode": 1.5
    },
    "report_text": "DNAKit 10x Chromium QC report\n=============================\nchemistry v3: 16 bp barcode + 12 bp UMI; read pairs: 3\nQ30 fraction (barcode / UMI / RNA read): 1.0 ...",
    "output_path": null
  }
}
```

## 44) plot_per_base_quality

- **作用：** 绘制逐碱基质量箱线图（FastQC 风格），输出 SVG 字符串。

- **处理方法：** 对每个读段位置（读段长于 `max_columns` 时按组合并位置），箱体表示 Phred Q 的下四分位数到上四分位数，红线为中位数，须线为第 10 到第 90 百分位数，蓝线为均值。SVG 以文本返回，仅在给出 `output_path` 时写入文件。纯 Python 实现，不依赖绘图库。

- **API：** `dnakit.sequencing_qc.plot_per_base_quality(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列），或读段对（此时由 `mate` 选择读段）。必须带有 Phred 质量值。选项：`mate`（`read1` 默认或 `read2`）、`max_columns`（默认 150）、`output_path`（可选 `.svg`；除非 `overwrite=True` 否则不覆盖）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `svg`（SVG 文本）、`output_path`（未写文件时为 `null`）、`positions`、`positions_per_box` 以及 `legend`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import plot_per_base_quality

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = plot_per_base_quality(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "plot_per_base_quality",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"370\" height=\"378\" viewBox=\"0 0 370 378\" font-family=\"sans-serif\" font-size=\"11\">\n<title>Quality scores per base (...",
    "output_path": null,
    "positions": 10,
    "positions_per_box": 1,
    "legend": "box = lower to upper quartile, red line = median, whiskers = 10th-90th percentile, blue line = mean; background: green Q>=28, yellow 20-28, red <20."
  }
}
```

## 45) plot_qc_summary

- **作用：** 绘制关键质控百分比（Q20、Q30、GC、N、重复、接头）的柱状图，输出 SVG 字符串。

- **处理方法：** 各柱分别为 Q20 及以上碱基的百分比、Q30 及以上碱基的百分比、GC 含量、N 碱基、重复读段（完全序列重复）以及命中接头的读段（`calculate_adapter_content` 的 k-mer 规则）。SVG 以文本返回并附各指标数值，仅在给出 `output_path` 时写入文件。纯 Python 实现，不依赖绘图库。

- **API：** `dnakit.sequencing_qc.plot_qc_summary(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为FASTQ 文本或文件路径（可为 `.gz`）、读段列表（(id, sequence, quality) 元组、`{sequence, quality}` 映射、FASTQ 记录对象或纯序列）。选项：`adapters`（同 `calculate_adapter_content`）、`output_path`（可选 `.svg`；除非 `overwrite=True` 否则不覆盖）、`phred_offset`（33 或 64，默认 33）、`max_reads`（默认 1000000）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `svg`（SVG 文本）、`output_path`（未写文件时为 `null`）以及 `metrics`（`bases >= Q20`、`bases >= Q30`、`GC content`、`N bases`、`duplicated reads`、`reads with adapter` 的百分比）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import plot_qc_summary

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

reads = fastq(
    ("r1", "ACGTACGTAC", "I" * 10),
    ("r2", "ACGTACGTAC", "I" * 10),
    ("r3", "TTTTTGGGGG", "5" * 10),
    ("r4", "NNACGTACGT", "!!" + "I" * 8),
)
result = plot_qc_summary(reads)
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "plot_qc_summary",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 112,
    "sequence_bases": 40,
    "lengths": [10, 10, 10, 10],
    "sha256": "74e1d5397e820db4ae7e5d14647765715a6c1e7a262dff425226e4bce9a7169a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"660\" height=\"204\" viewBox=\"0 0 660 204\" font-family=\"sans-serif\" font-size=\"11\">\n<title>Sequencing QC summary</ti...",
    "output_path": null,
    "metrics": {
      "bases >= Q20": 95.0,
      "bases >= Q30": 70.0,
      "GC content": 50.0,
      "N bases": 5.0,
      "duplicated reads": 25.0,
      "reads with adapter": 0.0
    }
  }
}
```

## 46) plot_qc_summary_tenx

- **作用：** 绘制 10x 质控百分比（条形码/UMI/RNA 读段的 Q30、接头、重复率替代指标、GC）的柱状图，输出 SVG 字符串。

- **处理方法：** 各柱分别为条形码、UMI 和 RNA 读段的 Q30 比例，含接头、模板转换寡核苷酸或 poly-A 命中的 RNA 读段比例，`calculate_duplication_tenx` 的重复率替代指标，RNA 读段的 GC 含量，以及（给出 `whitelist` 时）有效条形码比例。SVG 以文本返回并附各指标数值，仅在给出 `output_path` 时写入文件。 不自动判断 10x 化学版本：由 `chemistry`（或显式给出的长度）决定如何把 read 1 拆分为细胞条形码和 UMI。

- **API：** `dnakit.sequencing_qc.plot_qc_summary_tenx(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为一对读段 `{"read1": <条形码+UMI 读段>, "read2": <cDNA 读段>}`（FASTQ 文本、路径或读段列表）；只给单个读段集会报错，且 read 1 长度不得小于 `barcode_length` + `umi_length`。选项：`chemistry`（`v1` 为 14+10、`v2` 为 16+10、`v3` 为 16+12 bp 条形码+UMI，默认 `v3`）、`barcode_length`、`umi_length`、`phred_offset`（33）、`max_reads`（1000000）。另有选项：`whitelist`（增加有效条形码柱）、`adapters`、`output_path`（可选 `.svg`；除非 `overwrite=True` 否则不覆盖）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `chemistry`、`barcode_length`、`umi_length`、`read_pairs`、`svg`（SVG 文本）、`output_path`（未写文件时为 `null`）以及 `metrics`（各项百分比，如 `Q30 bases in barcode`、`Q30 bases in UMI`、`Q30 bases in RNA read`、`RNA reads with adapter/TSO/polyA`、`duplicate rate (proxy)`、`GC content of RNA read`，给出白名单时还有 `valid barcodes`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.sequencing_qc import plot_qc_summary_tenx

def fastq(*records):
    return "".join(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)

barcode = "ACGTACGTACGTACGT"
tenx = {
    "read1": fastq(
        ("p1", barcode + "A" * 12, "I" * 28),
        ("p2", barcode + "A" * 12, "I" * 28),
        ("p3", "T" * 16 + "C" * 12, "I" * 16 + "5" * 12),
    ),
    "read2": fastq(
        ("p1", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p2", "GATTACAGATCGGAAGAGCGGTTCAGCAGG", "I" * 30),
        ("p3", "C" * 30, "5" * 30),
    ),
}
result = plot_qc_summary_tenx(tenx, options={"whitelist": [barcode]})
print(result)
```

- **示例输出：**

```json
{
  "module": "sequencing_qc",
  "operation": "plot_qc_summary_tenx",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 446,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c2bd8879256656753e784d4d4fe4f7d2de6add20603c60c539095f494a8d364c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "chemistry": "v3",
    "barcode_length": 16,
    "umi_length": 12,
    "read_pairs": 3,
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"660\" height=\"228\" viewBox=\"0 0 660 228\" font-family=\"sans-serif\" font-size=\"11\">\n<title>10x Chromium QC summary (...",
    "output_path": null,
    "metrics": {
      "Q30 bases in barcode": 100.0,
      "Q30 bases in UMI": 66.6667,
      "Q30 bases in RNA read": 66.6667,
      "RNA reads with adapter/TSO/polyA": 66.6667,
      "duplicate rate (proxy)": 33.3333,
      "GC content of RNA read": 68.8889,
      "valid barcodes": 66.6667
    }
  }
}
```

## 47) run_cutadapt

- **作用：** 为“测序质量控制”类分析提供可复用的 `run_cutadapt` 功能。候选编号 CAP-0683，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验本地路径和带类型的选项，不经 shell 调用用户自行安装的 Cutadapt 可执行程序，限制运行时间和输出大小，检查声明的修剪后 FASTQ 文件和可选的 JSON 报告，并返回结构化的来源信息。固定的选项白名单是 DNAKit 设定的边界，并不表示开放了该后端的完整命令行。封装记录了 Martin (2011), DOI 10.14806/ej.17.1.200 的引用信息和 MIT 许可证证据。

- **API：** `dnakit.sequencing_qc.run_cutadapt(input, *, options=None, verbose=True)`（正式 DNAKit 命名空间；Agent/MCP 中也注册了同样受限的外部进程调用）。

- **输入：** 必填 `input`：一个 FASTQ 路径或 `{read1, read2}` 映射。`options` 中必须给出用户自行安装程序的 `executable_path` 和输出位置；`verbose` 可选。只接受文档列出的白名单选项，且输出路径必须事先不存在。

- **输出：** 封装返回可 JSON 序列化的完成封装，包含 `status`、`backend`、实际的 argv 参数、声明的 `outputs`、截断后的后端输出和 `provenance`；后端生成声明的修剪后 FASTQ 文件和可选的 JSON 报告。可执行程序缺失、非零退出、产物缺失、JSON 报告无效以及超出输出限制，均以结构化错误返回。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_cutadapt

result = run_cutadapt(
    input={"read1": Path("reads_R1.fastq"), "read2": Path("reads_R2.fastq")},
    options={
        "executable_path": "/usr/local/bin/cutadapt",
        "output1": Path("results/trimmed_R1.fastq"),
        "output2": Path("results/trimmed_R2.fastq"),
        "adapters": ["AGATCGGAAG"],
        "json_report": Path("results/cutadapt.json"),
    },
    verbose=False,
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "backend": "cutadapt",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Martin (2011), DOI 10.14806/ej.17.1.200",
    "license": "MIT"
  }
}
```

## 48) run_fastp

- **作用：** 为“测序质量控制”类分析提供可复用的 `run_fastp` 功能。候选编号 CAP-0684，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验本地路径和带类型的选项，不经 shell 调用用户自行安装的 fastp 可执行程序，限制运行时间和输出大小，检查声明的过滤后 FASTQ 文件以及 JSON 和 HTML 报告，并返回结构化的来源信息。固定的选项白名单是 DNAKit 设定的边界，并不表示开放了该后端的完整命令行。封装记录了 Chen (2025), DOI 10.1002/imt2.70078 的引用信息和 MIT 许可证证据。

- **API：** `dnakit.sequencing_qc.run_fastp(input, *, options=None, verbose=True)`（正式 DNAKit 命名空间；Agent/MCP 中也注册了同样受限的外部进程调用）。

- **输入：** 必填 `input`：一个 FASTQ 路径或 `{read1, read2}` 映射。`options` 中必须给出用户自行安装程序的 `executable_path` 和输出位置；`verbose` 可选。只接受文档列出的白名单选项，且输出路径必须事先不存在。

- **输出：** 封装返回可 JSON 序列化的完成封装，包含 `status`、`backend`、实际的 argv 参数、声明的 `outputs`、截断后的后端输出和 `provenance`；后端生成声明的过滤后 FASTQ 文件以及 JSON 和 HTML 报告。可执行程序缺失、非零退出、产物缺失、JSON 报告无效以及超出输出限制，均以结构化错误返回。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_fastp

result = run_fastp(
    input=Path("reads.fastq"),
    options={
        "executable_path": "/usr/local/bin/fastp",
        "output": Path("results/filtered.fastq"),
        "json_report": Path("results/fastp.json"),
        "html_report": Path("results/fastp.html"),
    },
    verbose=False,
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "backend": "fastp",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Chen (2025), DOI 10.1002/imt2.70078",
    "license": "MIT"
  }
}
```

## 49) run_fastqc

- **作用：** 为“测序质量控制”类分析提供可复用的 `run_fastqc` 功能。候选编号 CAP-0685，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验本地路径和带类型的选项，不经 shell 调用用户自行安装的 FastQC 0.12.x 可执行程序，限制运行时间和输出大小，检查已有输出目录中的 FastQC HTML 和 ZIP 报告，并返回结构化的来源信息。固定的选项白名单是 DNAKit 设定的边界，并不表示开放了该后端的完整命令行。封装记录了 FastQC 项目文档作为引用（调研来源未提供 DOI），以及 GPL-3.0 许可证证据。

- **API：** `dnakit.sequencing_qc.run_fastqc(input, *, options=None, verbose=True)`（正式 DNAKit 命名空间；Agent/MCP 中也注册了同样受限的外部进程调用）。

- **输入：** 必填 `input`：1 到 8 个本地 FASTQ/BAM/SAM 类读段路径。`options` 中必须给出用户自行安装程序的 `executable_path` 和输出位置；`verbose` 可选。只接受文档列出的白名单选项，且输出路径必须事先不存在。

- **输出：** 封装返回可 JSON 序列化的完成封装，包含 `status`、`backend`、实际的 argv 参数、声明的 `outputs`、截断后的后端输出和 `provenance`；后端生成已有输出目录中的 FastQC HTML 和 ZIP 报告。可执行程序缺失、非零退出、产物缺失、JSON 报告无效以及超出输出限制，均以结构化错误返回。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_fastqc

result = run_fastqc(
    input=Path("reads.fastq.gz"),
    options={
        "executable_path": "/opt/FastQC/fastqc",
        "output_directory": Path("results/fastqc"),
        "kmer": 7,
    },
    verbose=False,
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "backend": "fastqc",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "FastQC project documentation; no DOI supplied by the survey",
    "license": "GPL-3.0"
  }
}
```

## 50) run_multiqc

- **作用：** 为“测序质量控制”类分析提供可复用的 `run_multiqc` 功能。候选编号 CAP-0686，来自 1 条原始记录，当前标记为 `new_candidate`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验本地路径和带类型的选项，不经 shell 调用用户自行安装的 MultiQC 可执行程序，限制运行时间和输出大小，检查空输出目录中声明的 MultiQC HTML 报告，并返回结构化的来源信息。固定的选项白名单是 DNAKit 设定的边界，并不表示开放了该后端的完整命令行。封装记录了 Ewels et al. (2016), DOI 10.1093/bioinformatics/btw354 的引用信息和 GPL-3.0-or-later 许可证证据。

- **API：** `dnakit.sequencing_qc.run_multiqc(input, *, options=None, verbose=True)`（正式 DNAKit 命名空间；Agent/MCP 中也注册了同样受限的外部进程调用）。

- **输入：** 必填 `input`：一个或多个包含质控输入的本地文件或目录。`options` 中必须给出用户自行安装程序的 `executable_path` 和输出位置；`verbose` 可选。只接受文档列出的白名单选项，且输出路径必须事先不存在。

- **输出：** 封装返回可 JSON 序列化的完成封装，包含 `status`、`backend`、实际的 argv 参数、声明的 `outputs`、截断后的后端输出和 `provenance`；后端生成空输出目录中声明的 MultiQC HTML 报告。可执行程序缺失、非零退出、产物缺失、JSON 报告无效以及超出输出限制，均以结构化错误返回。

- **示例代码：**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_multiqc

result = run_multiqc(
    input=[Path("results/fastqc")],
    options={
        "executable_path": "/opt/multiqc/bin/multiqc",
        "output_directory": Path("results/multiqc"),
        "filename": "multiqc_report.html",
    },
    verbose=False,
)
print(result)
```

- **示例输出：**

```json
{
  "status": "completed",
  "backend": "multiqc",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Ewels et al. (2016), DOI 10.1093/bioinformatics/btw354",
    "license": "GPL-3.0-or-later"
  }
}
```
