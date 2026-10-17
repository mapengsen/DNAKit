# 分子实验

本类收录 **10** 项功能。10 项均在本地实际计算并返回 `status: "ok"`（纯 Python；Tm 与消光系数来自 `dnakit.thermodynamics`），没有返回 `not_implemented` 的函数，也都不需要外部程序或可选库。它们是实验台计算器、方案规划工具和小型序列工具（Sanger 峰图碱基识别与读段组装、引物分析、优化与检验），并不能替代专门的组装软件或引物设计软件。这些函数通过正式的 DNAKit 命名空间导出，JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。下面的示例输出均为真实运行结果，较长的数组和字符串以 `...` 省略。

## 1) analyze

- **作用：** 分析寡核苷酸序列：碱基组成、GC 含量、反向互补序列、分子量、A260 消光系数和熔解温度。

- **处理方法：** 逐条序列统计碱基数，计算 GC 含量、反向互补序列、最长同聚物以及 3' 端是否有 G/C 夹。分子量为带 5'-OH 的无水单链分子量，260 nm 消光系数为 `dnakit.thermodynamics` 给出的最近邻理论值（含减色效应），并据此给出每 OD260 的 nmol 与 ug。Tm 默认采用 SantaLucia (1998) 最近邻模型（`sodium_mm`、`primer_nm`）；`wallace`（每个 A/T 记 2，每个 G/C 记 4）适用于 2-13 nt，`basic`（64.9 + 41 (G+C-16.4)/N）适用于不短于 14 nt 的序列，否则 Tm 为 `null` 并附 `tm_note`。纯 Python 实现，不依赖外部程序。

- **API：** `dnakit.molecular_lab.analyze(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 DNA 字符串（仅限 A、C、G、T）、FASTA 文本、序列或 `{id, sequence}` 映射的列表，或名称 -> 序列的映射（最多 1000 条、每条 100 kb）。选项：`tm_method`（`nearest_neighbor` 默认、`wallace`、`basic`）、`sodium_mm`（默认 50）、`primer_nm`（默认 250）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `sequences`（每条序列含 `length`、`composition`、`gc_fraction`、`reverse_complement`、`is_self_complementary`、`longest_homopolymer`、`ends_with_gc_clamp`、`molecular_weight_da`、`extinction_coefficient_260nm`、`nmol_per_od260`、`ug_per_od260`、`tm_celsius`、`tm_method`）、`count` 以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import analyze

result = analyze({"primer": "CGTTGACCGTAGCTAGGCTA"}, options={"sodium_mm": 50, "primer_nm": 250})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "analyze",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "874a8555eeae8766308c0b460b1ba35a68a63691bbf0f8bb6ce8790a59570a39"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "sequences": {
      "primer": {
        "length": 20,
        "composition": {
          "A": 4,
          "C": 5,
          "G": 6,
          "T": 5
        },
        "gc_fraction": 0.55,
        "reverse_complement": "TAGCCTAGCTACGGTCAACG",
        "is_self_complementary": false,
        "longest_homopolymer": 2,
        "ends_with_gc_clamp": false,
        "molecular_weight_da": 6133.04,
        "extinction_coefficient_260nm": 191900.0,
        "nmol_per_od260": 5.21105,
        "ug_per_od260": 31.9596,
        "tm_celsius": 55.432,
        "tm_method": "nearest_neighbor"
      }
    },
    "count": 1,
    "definition": "Composition, GC fraction, reverse complement, anhydrous single-strand molecular weight, theoretical 260 nm extinction coefficient (nearest-neighbour hypochromic..."
  }
}
```

## 2) assemble_sanger_reads

- **作用：** 通过贪心的后缀/前缀重叠合并，将相互重叠的 Sanger 读段组装成 contig，并同时考虑两条链。

- **处理方法：** 贪心的重叠-布局组装：反复合并后缀/前缀重叠最长、重叠不少于 `min_overlap` 个碱基且无空位一致性不低于 `min_identity` 的一对 contig（或读段）。启用 `both_strands` 时，读段可取反向互补后再拼接。重叠区的一致序列保留较早 contig 的碱基，不做质量加权，也不处理空位，因此读段应先做质量修剪。被包含的读段会并入；无法连接的读段作为 singlet 返回。纯 Python 实现，不依赖外部组装软件。

- **API：** `dnakit.molecular_lab.assemble_sanger_reads(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 FASTA 文本、序列或 `{id, sequence}` 映射的列表，或名称 -> 序列的映射（最多 100 条、每条 20 kb；允许 IUPAC 字母）。选项：`min_overlap`（默认 20 个碱基）、`min_identity`（默认 0.95，无空位）、`both_strands`（默认 `True`）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `reads`、`contigs`（每个 contig 含 `id`、`length`、`sequence` 以及 `reads`——其 `id`、`strand`、`offset`、`length`）、`singlets`、所用选项和 `definition` 说明。contig 的方向是任意的，可能是原模板的反向互补序列。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import assemble_sanger_reads

template = "ATGGCTAGCTAGGATCCGATTACAGGCTTAACGGATCGCTAGGCATCGATTCGGCAAGTCTTGACCATGCAGG"
reverse_complement = template[40:].translate(str.maketrans("ACGT", "TGCA"))[::-1]
reads = {"f1": template[:30], "f2": template[20:50], "r3": reverse_complement}
result = assemble_sanger_reads(reads, options={"min_overlap": 10})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "assemble_sanger_reads",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 123,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5b7ac51ed6b5dfc15af47ad4da1d341477c59246037d2cf25852f57ad54a4d0b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "reads": 3,
    "contigs": [
      {
        "id": "contig1",
        "length": 73,
        "sequence": "CCTGCATGGTCAAGACTTGCCGAATCGATGCCTAGCGATCCGTTAAGCCTGTAATCGGATCCTAGCTAGCCAT",
        "reads": [
          {
            "id": "r3",
            "strand": "+",
            "offset": 0,
            "length": 33
          },
          {
            "id": "f2",
            "strand": "-",
            "offset": 23,
            "length": 30
          },
          {
            "id": "f1",
            "strand": "-",
            "offset": 43,
            "length": 30
          }
        ]
      }
    ],
    "singlets": [],
    "min_overlap": 10,
    "min_identity": 0.95,
    "both_strands": true,
    "definition": "Greedy overlap-layout: repeatedly merge the pair of contigs with the longest suffix/prefix overlap (>= min_overlap bases, identity >= min_identity, reads or the..."
  }
}
```

## 3) calculate

- **作用：** 按 `options["calculation"]` 选择并执行一项实验台计算：稀释、摩尔浓度、拷贝数、连接、A260、Tm、PCR 产量或凝胶片段大小。

- **处理方法：** 参数取自输入映射和 `options`（以 options 为准）。`dilution`：c1、v1、c2、v2 中给出三个（C1 V1 = C2 V2），求第四个及稀释剂体积。`molarity`：按每 bp 660 g/mol（`ds`）或每 nt 330 g/mol（`ss`）在质量浓度与摩尔浓度间换算，给出 `sequence` 时用精确分子量。`copy_number`：拷贝数 = 质量 / 分子量 x 6.02214076e23。`ligation`：插入片段质量 = 比例 x 载体质量 x 插入 bp / 载体 bp，给出浓度时同时算出体积。`a260`：按每 A260 单位 50（dsDNA）、33（ssDNA）或 40（RNA）ng/uL 换算浓度，给出 `sequence` 时使用寡核苷酸消光系数。`tm`：最近邻（默认）、`wallace` 或 `basic`。`pcr_yield`：N0 x (1 + E)^循环数。`gel_size`：对 log10(片段大小) 与迁移距离做最小二乘拟合。纯 Python 实现；Tm 与消光系数来自 `dnakit.thermodynamics`。

- **API：** `dnakit.molecular_lab.calculate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为参数映射；必须给出 `options["calculation"]` 以选择计算类型。`dilution`：`c1`、`v1`、`c2`、`v2` 中的三个。`molarity`：`length`（bp）或 `sequence`、`strand`（`ds` 默认或 `ss`），以及 `mass_concentration_ng_per_ul` 或 `molar_concentration_nm` 之一。`copy_number`：`mass_ng`、`length` 或 `sequence`、`strand`。`ligation`：`vector_bp`、`insert_bp`、`vector_ng`、`insert_to_vector_ratio`（默认 3），可选 `vector_conc_ng_per_ul` 与 `insert_conc_ng_per_ul`。`a260`：`a260`，`nucleic_acid`（`dsDNA`、`ssDNA`、`RNA`）或 `sequence`，`dilution_factor`、`path_length_cm`。`tm`：`sequence`、`method`（`nearest_neighbor` 默认、`wallace`、`basic`）、`sodium_mm`（50）、`primer_nm`（250）。`pcr_yield`：`template_copies`、`cycles`、`efficiency`（默认 1.0）。`gel_size`：`ladder`（`{bp, distance}` 列表）和 `distances`。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `calculation`、`formula` 说明以及随计算类型而异的字段：`dilution` 为四个量及 `solved`、`diluent_volume`、`dilution_factor`；`molarity` 为 `molecular_weight`、`mw_basis`、`mass_concentration_ng_per_ul`、`molar_concentration_nm`、`pmol_per_ul`；`copy_number` 为 `molecular_weight`、`mass_ng`、`copies`；`ligation` 为 `insert_ng`、`vector_ng`、`insert_to_vector_ratio`、`vector_pmol`、`insert_pmol`（给出浓度时还有体积）；`a260` 为 `method`、`mass_concentration_ng_per_ul`；`tm` 为 `tm_celsius`、`gc_fraction`、`conditions`、`model`；`pcr_yield` 为 `product_copies`、`fold_amplification`；`gel_size` 为 `slope_log10bp_per_distance`、`intercept_log10bp`、`r_squared`、`estimated_bp`、`extrapolated_points`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import calculate

parameters = {"vector_bp": 3000, "insert_bp": 1000, "vector_ng": 50, "vector_conc_ng_per_ul": 25, "insert_conc_ng_per_ul": 10}
result = calculate(parameters, options={"calculation": "ligation"})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 113,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "50280e11e9bf70b18a72becc28c313b3cb8c918d937adb6ffd3f004a5dd1f027"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "calculation": "ligation",
    "insert_ng": 50.0,
    "vector_ng": 50.0,
    "insert_to_vector_ratio": 3.0,
    "vector_pmol": 0.0252525,
    "insert_pmol": 0.0757576,
    "formula": "insert_ng = ratio * vector_ng * insert_bp / vector_bp (equal ratio of ends, 660 g/mol per bp)",
    "vector_volume_ul": 2.0,
    "insert_volume_ul": 5.0
  }
}
```

## 4) call_bases_from_sanger

- **作用：** 从 Sanger 色谱峰图中逐个峰位选取信号最强的通道来识别碱基，混合峰用 IUPAC 简并码表示。

- **处理方法：** 在每个峰位上取信号最强的通道对应的碱基；若次强通道不低于最强通道的 `mixed_peak_ratio`，则报告双碱基 IUPAC 简并码，并列入 `mixed_positions`。对 ABIF `.ab1` 文件，四个通道为按 FWO_1 排序的 DATA9-12，峰位取自 PLOC2，结果同时给出与仪器自身碱基识别的一致率。未提供峰位时，自动检测四个通道逐点最大值曲线上高于 `min_peak_height`、间隔不小于 `min_peak_spacing` 个采样点的局部极大值。不估计 Phred 式质量值。纯 Python 实现，不依赖外部碱基识别程序。

- **API：** `dnakit.molecular_lab.call_bases_from_sanger(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为 ABIF `.ab1` 文件路径，或映射 `{"traces": {"A": [...], "C": [...], "G": [...], "T": [...]}, "peak_locations": [...]}`，其中 `peak_locations` 可省略。选项：`mixed_peak_ratio`（默认 0.5）、`min_peak_height`（默认为全局最大值的 0.1，检测峰位时使用）、`min_peak_spacing`（默认 3 个采样点）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `source`、`samples`、`peaks`、`peaks_detected_automatically`、`sequence`、`secondary_to_primary_peak_ratio`（每个峰）、`mixed_positions`（`index`、`sample`、`bases`、`code`、`secondary_to_primary`）、`peak_locations` 以及 `definition` 说明；`.ab1` 文件还含通道顺序和 `agreement_with_instrument_basecalls`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import call_bases_from_sanger

import math

def peak(n, centres, sigma=2.0):
    return [sum(h * math.exp(-((i - c) ** 2) / (2 * sigma**2)) for c, h in centres.items()) for i in range(n)]

traces = {
    "A": peak(100, {10: 100, 90: 100, 50: 60}),
    "C": peak(100, {30: 100}),
    "G": peak(100, {50: 100}),
    "T": peak(100, {70: 100}),
}
result = call_bases_from_sanger({"traces": traces, "peak_locations": [10, 30, 50, 70, 90]})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "call_bases_from_sanger",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 9056,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8570c3e0ee1c1dfb6b9f2ea92f12e1b21c7946c975fc500d682a7d0a69c963e9"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "source": "mapping",
    "samples": 100,
    "peaks": 5,
    "peaks_detected_automatically": false,
    "sequence": "ACRTA",
    "secondary_to_primary_peak_ratio": [
      0.0,
      0.0,
      0.6,
      0.0,
      "..."
    ],
    "mixed_positions": [
      {
        "index": 2,
        "sample": 50,
        "bases": "GA",
        "code": "R",
        "secondary_to_primary": 0.6
      }
    ],
    "peak_locations": [
      10,
      30,
      50,
      70,
      "..."
    ],
    "definition": "At each peak location the base is the channel with the highest signal; if the second-highest channel is at least mixed_peak_ratio (default 0.5) of it, the IUPAC..."
  }
}
```

## 5) export

- **作用：** 将实验记录（样品、引物、序列）序列化为 CSV、TSV、JSON 或 FASTA 文本，可选择写入文件。

- **处理方法：** 每条记录转为一行（CSV/TSV，标准引号规则）、一个 JSON 对象，或一条 FASTA 记录（标题取记录的 `id` 或 `name`，序列按 `line_width` 折行）。列默认按字段首次出现的顺序排列。仅在给出 `output_path` 时写文件，已存在的文件除非 `overwrite=True` 否则不会被覆盖。纯 Python 实现。

- **API：** `dnakit.molecular_lab.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为记录映射（或序列字符串）的列表、单条记录，或 id -> 序列的映射。选项：`format`（`csv` 默认、`tsv`、`json`、`fasta`）、`columns`（字段顺序，默认按字段出现顺序）、`line_width`（FASTA，默认 80）、`output_path`（可选文件）、`overwrite`（默认 `False`）。FASTA 需要 `sequence`（或 `seq`）字段。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `format`、`records`、`columns`、`text`（序列化后的文本）和 `output_path`（未写文件时为 `null`）。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import export

records = [{"id": "p1", "sequence": "ACGT", "conc": 5}, {"id": "p2", "sequence": "GG", "conc": None}]
result = export(records, options={"format": "fasta"})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 91,
    "sequence_bases": 6,
    "lengths": [
      4,
      2
    ],
    "sha256": "540186960d98208a3e6d57ad53aa6c394164e5559c0e9ac130c9bdb7d63aa38e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "fasta",
    "records": 2,
    "columns": [
      "id",
      "sequence",
      "conc"
    ],
    "text": ">p1\nACGT\n>p2\nGG\n",
    "output_path": null
  }
}
```

## 6) optimize

- **作用：** 通过调整引物长度使其熔解温度接近目标值。

- **处理方法：** 引物取模板前 n 个碱基（`forward`）或模板末尾 n 个碱基的反向互补序列（`reverse`），n 从 `min_length` 到 `max_length`。对每个候选计算 Tm（采用 `dnakit.thermodynamics` 中 SantaLucia 1998 最近邻模型，50 mM Na+、250 nM 引物，或 `basic` 公式），选取 Tm 最接近 `target_tm` 的长度，并列时取较短者。启用 `require_gc_clamp` 时，只考虑 3' 端碱基为 G 或 C 的候选（若没有则报错）。不评估二级结构和特异性。

- **API：** `dnakit.molecular_lab.optimize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为模板 DNA 字符串。选项：`target_tm`（默认 60）、`min_length`（18）、`max_length`（30）、`direction`（`forward` 默认或 `reverse`）、`require_gc_clamp`（默认 `False`）、`tm_method`（`nearest_neighbor` 默认或 `basic`；前者使用 `sodium_mm` 50 与 `primer_nm` 250）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `direction`、`target_tm`、`tm_method`、`best`（所选候选：`length`、`primer`、`tm_celsius`、`tm_error`、`gc_fraction`、`gc_clamp`）、`candidates`（每个尝试长度的相同字段）以及 `definition` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import optimize

template = "CGTTGACCGTAGCTAGGCTAATCGGATTCGCAGGT"
result = optimize(template, options={"target_tm": 60})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "optimize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 35,
    "sequence_bases": 35,
    "lengths": [
      35
    ],
    "sha256": "3de911088164c6583899211652ff4c24cfac73b79ceab8df82c00b738ee6f019"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "direction": "forward",
    "target_tm": 60.0,
    "tm_method": "nearest_neighbor",
    "best": {
      "length": 25,
      "primer": "CGTTGACCGTAGCTAGGCTAATCGG",
      "tm_celsius": 60.835,
      "tm_error": 0.83537,
      "gc_fraction": 0.56,
      "gc_clamp": true
    },
    "candidates": [
      {
        "length": 18,
        "primer": "CGTTGACCGTAGCTAGGC",
        "tm_celsius": 54.425,
        "tm_error": 5.5755,
        "gc_fraction": 0.611111,
        "gc_clamp": true
      },
      {
        "length": 19,
        "primer": "CGTTGACCGTAGCTAGGCT",
        "tm_celsius": 55.786,
        "tm_error": 4.2139,
        "gc_fraction": 0.578947,
        "gc_clamp": false
      },
      {
        "length": 20,
        "primer": "CGTTGACCGTAGCTAGGCTA",
        "tm_celsius": 55.432,
        "tm_error": 4.5675,
        "gc_fraction": 0.55,
        "gc_clamp": false
      },
      "..."
    ],
    "definition": "Primers are the first (forward) or reverse complement of the last (reverse) min_length..max_length bases of the template; the length whose Tm is closest to targ..."
  }
}
```

## 7) parse

- **作用：** 将实验室量值字符串（`50 ng/uL`、`10 mM`、`200 uL`）解析为 SI 数值，或展开孔板孔位范围（`A1:B3`）。

- **处理方法：** 量值格式为 `<数值> <单位>`，单位可为 g、L、mol、M，前缀可为 k、m、u（或微符号）、n、p、f，也可为质量/体积或物质的量/体积的比值，`X` 表示倍数浓度；数值被换算为 g、L、mol、g/L 或 mol/L，并给出量纲。`kind` 为 `wells` 时，把 96 孔（8x12）或 384 孔（16x24）板上的单孔和矩形范围按行优先或列优先展开。纯 Python 实现。

- **API：** `dnakit.molecular_lab.parse(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为量值字符串或其列表（`kind` 为 `quantity`），或形如 `A1:B3, H12` 的孔位文本（`kind` 为 `wells`）。选项：`kind`（`quantity` 默认或 `wells`）；对孔位另有 `plate`（96 默认或 384）和 `order`（`row` 默认或 `column`）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 对量值：`kind`、`quantities`（每项含 `text`、`value`、`unit`、`dimension`、`si_value`、`si_unit`）和 `definition` 说明；对孔位：`kind`、`plate`、`order`、`wells`、`count` 和 `definition`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import parse

result = parse(["50 ng/uL", "1.5 mM", "200 \u00b5L", "5X"])
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "parse",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 38,
    "sequence_bases": 8,
    "lengths": [4, 2, 1, 1],
    "sha256": "825a57d9c055f721e499e4afc93c58da3805cc480a7ed180994258ebc4f722f0"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "quantity",
    "quantities": [
      {
        "text": "50 ng/uL",
        "value": 50.0,
        "unit": "ng/uL",
        "dimension": "mass_concentration",
        "si_value": 0.05,
        "si_unit": "g/L"
      },
      {
        "text": "1.5 mM",
        "value": 1.5,
        "unit": "mM",
        "dimension": "molar_concentration",
        "si_value": 0.0015,
        "si_unit": "mol/L"
      },
      {
        "text": "200 µL",
        "value": 200.0,
        "unit": "µL",
        "dimension": "volume",
        "si_value": 0.0002,
        "si_unit": "L"
      },
      "..."
    ],
    "definition": "'<number> <unit>' with units g, L, mol, M and prefixes k, m, u (or the micro sign), n, p, f; ratios mass/volume and amount/volume; 'X' for fold concentration. S..."
  }
}
```

## 8) plan

- **作用：** 规划带精确体积的实验方案：PCR 预混液、梯度稀释或连接反应体系。

- **处理方法：** 参数取自输入映射和 `options`（以 options 为准），由 `protocol` 选择方案。`pcr`：各组分体积 = 终浓度 / 母液浓度 x 反应体积，其余用水补足，预混液（除模板外的全部组分）按 反应数 x (1 + 余量) 放大。`serial_dilution`：每管由上一管按 `dilution_factor` 稀释到 `final_volume_ul`，给出每管的转移体积与稀释剂体积。`ligation`：按摩尔比（同 `calculate`）求插入片段质量，再给出载体、插入片段、10x 连接酶缓冲液、连接酶和水的体积。组分总体积超出反应体积的方案会被拒绝。纯 Python 实现。

- **API：** `dnakit.molecular_lab.plan(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为参数映射（可为空）；`options["protocol"]` 取 `pcr`（默认）、`serial_dilution` 或 `ligation`。`pcr`：`reaction_volume_ul`（25）、`n_reactions`（1）、`overage_fraction`（0.1）、`template_volume_ul`（1）、`components`（`{name, stock, final, unit}` 列表；默认 5x 母液配 1x 缓冲液、0.2 mM dNTP、0.5 uM 引物、0.025 U/uL 聚合酶）。`serial_dilution`：`start_concentration`、`dilution_factor`（10）、`steps`、`final_volume_ul`（100）。`ligation`：`vector_bp`、`insert_bp`、`vector_ng`、`insert_to_vector_ratio`（3），可选浓度，`reaction_volume_ul`（20）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `protocol`、`definition` 说明以及随方案而异的字段：`pcr` 为 `reaction_volume_ul`、`n_reactions`、`overage_fraction`、`master_mix_for_reactions`、`components`（每个组分的 `stock`、`final`、`unit`、`volume_per_reaction_ul`、`master_mix_volume_ul`，含水）、`master_mix_total_ul`、`master_mix_aliquot_per_tube_ul`、`template_per_tube_ul`；`serial_dilution` 为 `dilution_factor`、`steps`、`tubes`、`total_diluent_ul`；`ligation` 为 `calculate` 的连接用量加 `reaction_volume_ul`、`ligase_buffer_10x_ul`、`ligase_ul`，给出浓度时还有载体、插入片段和水的体积。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import plan

result = plan({}, options={"protocol": "pcr", "n_reactions": 10})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "plan",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 2,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "protocol": "pcr",
    "reaction_volume_ul": 25.0,
    "n_reactions": 10,
    "overage_fraction": 0.1,
    "master_mix_for_reactions": 11.0,
    "components": [
      {
        "component": "reaction buffer",
        "stock": 5.0,
        "final": 1.0,
        "unit": "X",
        "volume_per_reaction_ul": 5.0,
        "master_mix_volume_ul": 55.0
      },
      {
        "component": "dNTP mix",
        "stock": 10.0,
        "final": 0.2,
        "unit": "mM each",
        "volume_per_reaction_ul": 0.5,
        "master_mix_volume_ul": 5.5
      },
      {
        "component": "forward primer",
        "stock": 10.0,
        "final": 0.5,
        "unit": "uM",
        "volume_per_reaction_ul": 1.25,
        "master_mix_volume_ul": 13.75
      },
      "..."
    ],
    "master_mix_total_ul": 264.0,
    "master_mix_aliquot_per_tube_ul": 24.0,
    "template_per_tube_ul": 1.0,
    "definition": "volume = final / stock * reaction volume per component; water makes up the rest; the master mix (everything but the template) is scaled by n_reactions * (1 + ov..."
  }
}
```

## 9) plot

- **作用：** 绘制示意性的虚拟琼脂糖凝胶图，输出 SVG 字符串。

- **处理方法：** 每条泳道是一组以 bp 计的片段大小。迁移距离在 `min_bp` 与 `max_bp` 之间与 log10(片段大小) 呈线性（琼脂糖常见的半对数关系），条带亮度随片段质量增加。这只是示意图，并非对具体凝胶、缓冲液或电压的模拟。SVG 以文本返回，仅在给出 `output_path` 时写入文件。纯 Python 实现，不依赖绘图库。

- **API：** `dnakit.molecular_lab.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为泳道名 -> 片段大小（bp）列表的映射（或此类列表的列表），例如 `{"ladder": [1000, 500, 100], "PCR": [650]}`。选项：`min_bp`（默认 100）、`max_bp`（默认 10000）、`output_path`（可选 `.svg`；已存在的文件除非 `overwrite=True` 否则不会被覆盖）。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `svg`（SVG 文本）、`output_path`（未写文件时为 `null`）、`lanes`、`min_bp`、`max_bp` 以及 `model` 说明。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import plot

result = plot({"ladder": [1000, 500, 100], "PCR": [650]})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1363816d06183a2cf60426b59ffcdf3e55f048b4076c468c821cc930afe14db7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"188\" height=\"400\" viewBox=\"0 0 188 400\" font-family=\"sans-serif\" font-size=\"10\">\n<title>Virtual agarose gel</titl...",
    "output_path": null,
    "lanes": 2,
    "min_bp": 100.0,
    "max_bp": 10000.0,
    "model": "Schematic gel: migration distance is linear in log10(size) between min_bp and max_bp (the usual semilog behaviour of agarose); band brightness grows with fragme..."
  }
}
```

## 10) validate

- **作用：** 依据经验性的设计准则检验 PCR 引物。

- **处理方法：** 对每条引物检查长度、GC 含量、3' G/C 夹（3' 端碱基为 G 或 C，且最后 5 个碱基中 G/C 不超过 3 个）、同聚物连续长度、同一二核苷酸单元重复 4 次及以上，以及 Tm（默认 SantaLucia 1998 最近邻）；对正向/反向引物对还检查 Tm 差。每项检查都给出 `passed` 与 `detail`。通过检查并不能预测扩增是否成功：二级结构、引物二聚体和脱靶结合均未评估。纯 Python 实现。

- **API：** `dnakit.molecular_lab.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input` 为单条引物序列、引物列表或映射，或引物对 `{"forward": ..., "reverse": ...}`。选项：`min_length`（18）、`max_length`（30）、`min_gc`（0.40）、`max_gc`（0.60）、`min_tm`（52）、`max_tm`（65）、`max_homopolymer`（4）、`max_tm_difference`（引物对，5）、`tm_method`（默认 `nearest_neighbor`）、`sodium_mm`、`primer_nm`。

- **输出：** 返回可 JSON 序列化的封装，包含 `status`（`"ok"`）、`module`、`operation`、`provenance`（`implementation: "native"`、`backend: "local"`）、`input_summary`、`warnings` 与 `result`。`result` 包含 `kind`（正向/反向引物对为 `primer_pair`，否则为 `primers`）、总体 `passed`、各引物的报告（引物对为 `forward` 与 `reverse`，否则为名称 -> 报告的 `primers` 映射；每份报告含 `sequence`、`tm_celsius`、`passed` 与 `checks`，每项检查含 `check`、`passed`、`detail`），引物对另有 `pair_checks`（Tm 差）以及说明所用阈值的 `rules`。输入或选项无效时抛出 `InputFormatError` 或 `ConfigurationError`。

- **示例代码：**

```python
from dnakit.molecular_lab import validate

result = validate({"forward": "CGTTGACCGTAGCTAGGCTA", "reverse": "GCGCGCGCGCGCGCGCGCGC"})
print(result)
```

- **示例输出：**

```json
{
  "module": "molecular_lab",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 70,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6a15be3b2d065bb6549f70ec94b2bba749a5f7594bfe14e633230b4ed427261e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "primer_pair",
    "passed": false,
    "forward": {
      "sequence": "CGTTGACCGTAGCTAGGCTA",
      "tm_celsius": 55.432,
      "passed": false,
      "checks": [
        {
          "check": "length",
          "passed": true,
          "detail": "20 nt (allowed 18-30)"
        },
        {
          "check": "gc_content",
          "passed": true,
          "detail": "55.0% (allowed 40-60%)"
        },
        {
          "check": "three_prime_gc_clamp",
          "passed": false,
          "detail": "last 5 nt GGCTA has 3 G/C; 3' base A (want 3' G/C and <= 3 G/C in the last 5)"
        },
        "..."
      ]
    },
    "reverse": {
      "sequence": "GCGCGCGCGCGCGCGCGCGC",
      "tm_celsius": 80.078,
      "passed": false,
      "checks": [
        {
          "check": "length",
          "passed": true,
          "detail": "20 nt (allowed 18-30)"
        },
        {
          "check": "gc_content",
          "passed": false,
          "detail": "100.0% (allowed 40-60%)"
        },
        {
          "check": "three_prime_gc_clamp",
          "passed": false,
          "detail": "last 5 nt CGCGC has 5 G/C; 3' base C (want 3' G/C and <= 3 G/C in the last 5)"
        },
        "..."
      ]
    },
    "pair_checks": [
      {
        "check": "tm_difference",
        "passed": false,
        "detail": "|Tm1 - Tm2| = 24.6 C (allowed <= 5)"
      }
    ],
    "rules": "Rule-of-thumb PCR primer criteria (adjustable): length 18-30 nt, GC 40-60%, Tm 52-65 C (nearest-neighbour, 50 mM Na+, 250 nM), 3' terminal G/C with at most 3 G/..."
  }
}
```
