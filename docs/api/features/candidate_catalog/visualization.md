# 报告与图片导出

本类收录 `visualization` 模块的 **7** 个函数。其中 5 个会计算出真实结果：`plot` 把序列图绘制为独立的 SVG（共七种），`generate_analysis_summary` 计算序列统计量，`generate_qc_html_report` 与 `generate_html_genome_report` 生成自包含的 HTML 报告，`export_figure` 保存 SVG（转换为 PNG/JPG/TIFF/PDF 还需可选库 `cairosvg` 和 `Pillow`，缺少时返回 `not_implemented`）。`generate_clinical_pdf_report` 始终返回 `status: "not_implemented"`，因为 DNAKit 无法产生经验证的临床解读。`save_figure` 是已有的 Matplotlib 图保存函数。示例输出均为真实运行结果；`input_summary` 简写为 `{...}`。

## 1) export_figure

- **作用：** 把 SVG 图写入文件，或转换为 PNG/JPG/TIFF/PDF。

- **处理方法：** SVG 必须是格式良好的 XML，根元素为 `<svg>` 且有数值型 width/height（或 viewBox），并且不得含脚本、事件处理属性、`foreignObject` 或 `javascript:`（`ConfigurationError`，`UNSAFE_SVG`）；DOCTYPE/ENTITY 声明会被 SVG 写出器拒绝。`.svg` 目标按原样写出字节。其他后缀（`.png`、`.jpg`/`.jpeg`、`.tif`/`.tiff`、`.pdf`）借助可选库 `cairosvg` 和 `Pillow` 转换；缺少时返回 `status: "not_implemented"`，`required_dependency` 为 `"pip install cairosvg pillow"`。其他后缀会抛 `ConfigurationError`。

- **API：** `dnakit.visualization.export_figure(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：SVG 文本、`plot` 的结果映射（含 `svg`）、完整的 `plot` 封装，或 `SVGArtifact`。`options`：`output`（必填，目标路径）、`overwrite`（默认 false；文件已存在会抛 `ConfigurationError`）、`dpi`（栅格导出，默认 96）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `path`、`format`、`bytes` 和已写文件的 `sha256`。

- **示例代码：**

```python
import json
from dnakit.visualization import export_figure, plot

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
figure = plot(fasta, options={"kind": "gc_distribution", "bins": 4})
result = export_figure(figure, options={"output": "gc.svg"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "visualization",
  "operation": "export_figure",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "path": "gc.svg",
    "format": "svg",
    "bytes": 3778,
    "sha256": "d9141462188a87dd8ebbaf0526446a1e705268ff8ed3d744b18c190e019186cb"
  }
}
```

## 2) generate_analysis_summary

- **作用：** 把一组 DNA 序列的计数、长度统计、GC 含量和碱基计数汇总为数据和 Markdown 文本。

- **处理方法：** `gc_percent` 为全部序列上的 `100 * (G+C) / (A+C+G+T)`（没有无歧义碱基时为 `null`）。`length` 含 min、max、mean、median 和 N50（满足长度 >= L 的序列包含至少一半全部碱基的最小 L）。`non_acgt_bases` 统计其余所有符号，`n_fraction` 为 `N` 所占比例。`per_sequence` 列出至多前 1000 条序列的 id、长度和 GC 百分比。

- **API：** `dnakit.visualization.generate_analysis_summary(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA 文本或 FASTA 文件路径、单条序列字符串、序列列表、`(id, sequence)` 对或带 `id`/`sequence` 的映射，或映射 `{id: sequence}`（带 `records` 的映射会被解包）。序列会转大写、去除空白、把 `U` 读作 `T`，只接受字母和 `-*.`（否则抛 `InputFormatError`）；最多 100,000 条序列、500,000,000 个碱基。`options`：`title`（默认 `"Sequence analysis summary"`）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `sequence_count`、`total_bases`、`length`、`gc_percent`、`base_counts`、`non_acgt_bases`、`n_fraction`、`per_sequence` 和 `text`（Markdown 摘要）。

- **示例代码：**

```python
import json
from dnakit.visualization import generate_analysis_summary

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
result = generate_analysis_summary(fasta, options={"title": "Demo"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "visualization",
  "operation": "generate_analysis_summary",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "sequence_count": 3,
    "total_bases": 30,
    "length": {
      "min": 8,
      "max": 12,
      "mean": 10.0,
      "median": 10.0,
      "n50": 10
    },
    "gc_percent": 46.42857142857143,
    "base_counts": {
      "A": 9,
      "C": 7,
      "G": 6,
      "N": 2,
      "T": 6
    },
    "non_acgt_bases": 2,
    "n_fraction": 0.06666666666666667,
    "per_sequence": [
      {
        "id": "a",
        "length": 10,
        "gc_percent": 50.0
      },
      {
        "id": "b",
        "length": 12,
        "gc_percent": 80.0
      },
      {
        "id": "c",
        "length": 8,
        "gc_percent": 0.0
      }
    ],
    "text": "# Demo\n\n- sequences: 3\n- total bases: 30\n- length min / median / mean / max: 8 / 10 / 10 / 12\n- N50: 10\n- GC content: 46.43 %\n- non-ACGT bases: 2\n"
  }
}
```

## 3) generate_clinical_pdf_report

- **作用：** 未实现：临床报告需要经过验证的临床解读数据和经认证的 PDF 报告模板。

- **处理方法：** 不做任何计算。DNAKit 无法在此产生经验证的临床解读，因此不生成临床报告。非临床摘要请用 `generate_analysis_summary`、`generate_qc_html_report` 或 `generate_html_genome_report`。

- **API：** `dnakit.visualization.generate_clinical_pdf_report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅在 `input_summary` 中汇总）；`options` 不起作用。

- **输出：** 始终返回 `status: "not_implemented"` 的封装，含 `required_dependency`（“经验证的临床变异解读（ClinVar/ACMG 规则）和经认证的报告模板”）和 `message`，没有 `result`。

- **示例代码：**

```python
import json
from dnakit.visualization import generate_clinical_pdf_report

result = generate_clinical_pdf_report({"variants": [{"gene": "BRCA1", "hgvs": "c.68_69del"}]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "visualization",
  "operation": "generate_clinical_pdf_report",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "validated clinical variant interpretation (ClinVar/ACMG rules) and a certified report template",
  "message": "Clinical reports must contain validated interpretations, which DNAKit cannot produce here."
}
```

## 4) generate_html_genome_report

- **作用：** 为基因组或组装序列生成自包含的 HTML 报告，含组装统计和可选的注释摘要。

- **处理方法：** 报告由本包的 HTML 报告构建器根据序列和 `generate_analysis_summary` 的统计量（数量、总碱基数、N50、GC）生成。若给出 `options["annotation"]`，则对其调用 `dnakit.annotation.summarize` 并把结果放入 `annotation`。HTML 放在 `html` 中返回，仅在给出 `output` 时写文件。

- **API：** `dnakit.visualization.generate_html_genome_report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA 文本或 FASTA 文件路径、单条序列字符串、序列列表、`(id, sequence)` 对或带 `id`/`sequence` 的映射，或映射 `{id: sequence}`（带 `records` 的映射会被解包）。序列会转大写、去除空白、把 `U` 读作 `T`，只接受字母和 `-*.`（否则抛 `InputFormatError`）；最多 100,000 条序列、500,000,000 个碱基（contig 或染色体）。`options`：`annotation`（GFF3/GTF 文本、路径或 features）、`title`、`max_records`（默认 1000）、`max_total_sequence_symbols`（默认 10,000,000）、`output`（`.html` 路径；不给则不写文件）、`overwrite`（默认 false；文件已存在会抛 `ConfigurationError`）、`include_html`（默认 true；false 则省略 `html` 字段）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `record_count`、`result_names`、`sha256`（HTML 的摘要）、`path`、`summary`、`annotation`（仅在给出注释时）和 `html`（`include_html` 为 false 时省略；示例中已省略）。

- **示例代码：**

```python
import json
from dnakit.visualization import generate_html_genome_report

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
gff3 = "##gff-version 3\nchr1\tsrc\tgene\t1000\t5000\t.\t+\t.\tID=g1;Name=GENE1\n"
result = generate_html_genome_report(fasta, options={"annotation": gff3, "include_html": False})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "visualization",
  "operation": "generate_html_genome_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "record_count": 3,
    "result_names": [
      "annotation",
      "summary"
    ],
    "sha256": "c4d6a5fbf9302ce2e5d1ec158cfa3072ca43160d5192cd45a483bcd83d12f1a6",
    "path": null,
    "summary": {
      "sequence_count": 3,
      "total_bases": 30,
      "length": {
        "min": 8,
        "max": 12,
        "mean": 10.0,
        "median": 10.0,
        "n50": 10
      },
      "gc_percent": 46.42857142857143,
      "base_counts": {
        "A": 9,
        "C": 7,
        "G": 6,
        "N": 2,
        "T": 6
      },
      "non_acgt_bases": 2,
      "n_fraction": 0.06666666666666667
    },
    "annotation": {
      "format": "gff3",
      "feature_count": 1,
      "types": {
        "gene": {
          "count": 1,
          "total_bp": 4001,
          "mean_length": 4001.0
        }
      },
      "seqids": {
        "chr1": 1
      },
      "strands": {
        "+": 1
      },
      "gene_count": 1,
      "transcript_count": 0,
      "exon_count": 0,
      "intron_count": 0,
      "coding_transcripts": 0,
      "multi_exon_transcripts": 0,
      "transcripts_per_gene": {
        "mean": 0.0,
        "max": 0
      },
      "exons_per_transcript": {
        "mean": 0.0,
        "max": 0
      }
    }
  }
}
```

## 5) generate_qc_html_report

- **作用：** 为 DNA 序列生成自包含的 HTML 质控报告，逐条序列给出通过/未通过标记。

- **处理方法：** 序列在以下情况判为未通过：长度小于 `min_length`；N 比例（`N` 个数 / 长度）高于 `max_n_fraction`；GC 百分比（`100 * (G+C) / (A+C+G+T)`）不在 [`min_gc_percent`, `max_gc_percent`] 内；不含 A/C/G/T 的序列没有 GC 值，会在任何非默认范围下判为未通过。每个失败项都列出原因。报告包含汇总统计和 `qc` 结论。阈值只是用户设定的简单截断值，并非经过验证的标准。

- **API：** `dnakit.visualization.generate_qc_html_report(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA 文本或 FASTA 文件路径、单条序列字符串、序列列表、`(id, sequence)` 对或带 `id`/`sequence` 的映射，或映射 `{id: sequence}`（带 `records` 的映射会被解包）。序列会转大写、去除空白、把 `U` 读作 `T`，只接受字母和 `-*.`（否则抛 `InputFormatError`）；最多 100,000 条序列、500,000,000 个碱基。`options`：`min_length`（默认 1）、`max_n_fraction`（默认 0.1）、`min_gc_percent`（默认 0）、`max_gc_percent`（默认 100）、`title`、`max_records`（默认 1000）、`max_total_sequence_symbols`（默认 10,000,000）、`output`（`.html` 路径；不给则不写文件）、`overwrite`（默认 false；文件已存在会抛 `ConfigurationError`）、`include_html`（默认 true；false 则省略 `html` 字段）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 含 `record_count`、`result_names`、`sha256`、`path`、`summary`、`qc`（`passed`、`failed_count`、含 `id` 和 `reasons` 的 `failed`、`thresholds`）和 `html`（`include_html` 为 false 时省略；示例中已省略）。

- **示例代码：**

```python
import json
from dnakit.visualization import generate_qc_html_report

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
result = generate_qc_html_report(fasta, options={"min_length": 9, "include_html": False})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "visualization",
  "operation": "generate_qc_html_report",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "record_count": 3,
    "result_names": [
      "qc",
      "summary"
    ],
    "sha256": "38211562920fb60a5494c1b34c47b0000a3d83a29af1d6bce34dea2d8c9681f6",
    "path": null,
    "summary": {
      "sequence_count": 3,
      "total_bases": 30,
      "length": {
        "min": 8,
        "max": 12,
        "mean": 10.0,
        "median": 10.0,
        "n50": 10
      },
      "gc_percent": 46.42857142857143,
      "base_counts": {
        "A": 9,
        "C": 7,
        "G": 6,
        "N": 2,
        "T": 6
      },
      "non_acgt_bases": 2,
      "n_fraction": 0.06666666666666667
    },
    "qc": {
      "passed": 1,
      "failed_count": 2,
      "failed": [
        {
          "id": "b",
          "reasons": [
            "N fraction 0.167 > 0.1"
          ]
        },
        {
          "id": "c",
          "reasons": [
            "length 8 < 9"
          ]
        }
      ],
      "thresholds": {
        "min_length": 9,
        "max_n_fraction": 0.1,
        "min_gc_percent": 0.0,
        "max_gc_percent": 100.0
      }
    }
  }
}
```

## 6) save_figure

- **作用：** 为“报告与图片导出”类分析提供可复用的 `save_figure` 功能。候选编号 CAP-0828，来自 180 条原始记录，当前标记为 `review_needed`；输入类型、后端、格式、平台、提供方、检测方法或生物学对象上的差异都放在同一个带类型的接口里，而不是为每个变体各建一个函数。

- **处理方法：** 校验带类型的输入封装，规范化记录与元数据，执行 `save_figure` 操作，并返回附带来源与诊断信息的结构化结果。可选参数 `plot_kind` 用于选择一种支持的变体（取值见源目录记录）。正式包复用现有的 Matplotlib 图片保存接口；图表类型和输出格式的差异仍作为该函数的选项。审核依据：`generic_operation_name`。已有同名函数证据：`dnakit.visualization.save_figure|save_figure`。

- **API：** `dnakit.visualization.save_figure(figure[required], target[required], overwrite[optional], tight_layout[optional], close[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 必填参数：`figure`、`target`；可选参数：`overwrite`、`tight_layout`、`close`。IO 参数按函数签名接受文本、字节、映射、路径或已读入的记录。

- **输出：** 现有的正式函数 `save_figure` 返回写出文件的 `Path`。它接受 Matplotlib 图对象，并根据目标文件后缀选择输出格式。

- **示例代码：**

```python
from pathlib import Path

import matplotlib.pyplot as plt

from dnakit.visualization import save_figure

figure, axis = plt.subplots()
axis.plot([0, 1], [0, 1])
result = save_figure(figure, Path("results/save_figure.svg"), overwrite=True)
print(result)
```

- **示例输出：**

```json
{
  "status": "formal_existing",
  "capability_id": "CAP-0828",
  "operation": "save_figure",
  "result_type": "PosixPath",
  "result_preview": "results/save_figure.svg",
  "provenance": {
    "schema": "formal_dnakit",
    "implementation": "existing_visualization_api"
  }
}
```

## 7) plot

- **作用：** 把基于序列的图绘制为独立的 SVG（合法 XML、无脚本），并连同所绘数值一并返回。

- **处理方法：** `kind` 选择图的类型，由本包无依赖的 SVG 渲染器绘制：`length_histogram`（默认；序列长度，单位 nt）、`gc_distribution`（每条序列的 `100 * (G+C) / (A+C+G+T)`，在 0-100 上分成 `bins` 个等宽区间；没有无歧义碱基的序列被跳过并计入 `skipped`）、`base_composition`（所有序列中各符号的比例，A、C、G、T 在前）、`sequence`（第一条序列或 `id` 指定序列的符号）、`linear_map` 与 `circular_map`（第一条序列或 `id` 指定的序列，`features` 为 0 起点半开坐标的 `{type, start, end}`；环形图总是把序列视为环状：无论 `topology` 为何值都会强制使用环状拓扑，因此无需另行指定环状序列）以及 `label_comparison`（`true` 与 `predicted` 类别标签的成对比例条形图）。`data` 保存所绘制的精确计数、边界或数值。除非给出 `output`，否则不写文件。

- **API：** `dnakit.visualization.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：FASTA 文本或 FASTA 文件路径、单条序列字符串、序列列表、`(id, sequence)` 对或带 `id`/`sequence` 的映射，或映射 `{id: sequence}`（带 `records` 的映射会被解包）。序列会转大写、去除空白、把 `U` 读作 `T`，只接受字母和 `-*.`（否则抛 `InputFormatError`）；最多 100,000 条序列、500,000,000 个碱基；`label_comparison` 则为映射 `{"true": [...], "predicted": [...]}`。`options`：`kind`（默认 `"length_histogram"`；未知类型会抛 `ConfigurationError`）、`bins`（长度默认 30，GC 默认 20，范围 1 到 1000）、`id`、`features`、`topology`（默认 `"linear"`，用于 `sequence`/`linear_map`）、`title`、`output`（要写入的 `.svg` 路径）、`overwrite`（默认 false）。

- **输出：** 封装包含 `status: "ok"`、`module`、`operation`、`provenance`、`input_summary`、`warnings` 和 `result`；`result` 是含 `kind`、`svg`（SVG 文本）、`width`、`height`、`sha256`（SVG 的摘要）、`data`（所绘数值，键随 `kind` 而异）和 `path`（已写文件或 `null`）的字典。示例中冗长的 `svg` 文本简写为 `<svg ... </svg>`。

- **示例代码：**

```python
import json
from dnakit.visualization import plot

fasta = ">a\nACGTACGTAC\n>b\nGGGGCCCCAANN\n>c\nATATATAT\n"
result = plot(fasta, options={"kind": "gc_distribution", "bins": 4, "title": "GC content"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "visualization",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {...},
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "gc_distribution",
    "svg": "<svg ... </svg>",
    "width": 720,
    "height": 720,
    "sha256": "78bc7424b731541585975c3e833d831387d12f7cabb3de97d995009360bc6861",
    "data": {
      "count": 3,
      "skipped": 0,
      "gc_percent": [
        50.0,
        80.0,
        0.0
      ],
      "bin_edges": [
        0.0,
        25.0,
        50.0,
        75.0,
        100.0
      ],
      "counts": [
        1,
        0,
        1,
        1
      ]
    },
    "path": null
  }
}
```
