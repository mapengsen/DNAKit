# 单细胞与空间组学

本类收录 **54** 项针对单细胞与空间斑点计数矩阵的有文档的函数。其中 52 项对传入的数据计算出真实结果：43 项为纯 Python，7 项在调用时导入 numpy（`cluster`、`cluster_spatial`、`infer_trajectory`、`infer_trajectory_spatial`、`map_cells`、`map_cells_spatial`、`map_single_cell_to_reference`），2 项导入 numpy 与 scipy（`deconvolute_spots`、`deconvolute_spots_spatial`）；所需软件包缺失时它们改为返回 `status: "not_implemented"`。另有 2 项，`remove_doublets` 与 `remove_doublets_spatial`，因未内置双细胞检测模型而始终返回 `status: "not_implemented"`。计数矩阵为基因×细胞（列为细胞或斑点）；`_spatial` 变体把斑点按细胞处理，并在方法需要时接收 `coordinates`。输入格式错误或缺少必填选项时抛出 `InputFormatError` / `ConfigurationError`。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) analyze_crispr_perturbations

- **作用：** 汇总 CRISPR 扰动筛选：每个 guide 的细胞数，以及各 guide 的靶基因相对对照细胞的敲低程度。

- **处理方法：** 表达量先按文库大小归一化（`count / 文库大小 * target_sum`，默认 10000）。对每个非对照标签，比较其 guide 细胞与对照细胞中靶基因（由 `targets` 指定，默认即 guide 标签本身）的平均归一化表达：`log2fc = log2(mean_guide + 1) - log2(mean_control + 1)`，`percent_knockdown = 100 * (1 - mean_guide / mean_control)`（对照均值为 0 时为 `null`）。不做显著性检验。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.analyze_crispr_perturbations(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`perturbations`（必填；`{cell: guide 标签}` 或按列顺序的列表）、`control`（`"control"`；对照细胞的标签）、`targets`（`{guide: 靶基因}`；默认 guide 标签即基因名）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`control`、`n_control_cells`、`cells_per_label`、`target_sum`、`knockdown_definition` 与 `guides`；`guides` 为 guide -> `n_cells`、`target_gene`、`target_in_matrix`，以及（基因在矩阵中时）`mean_normalized_guide`、`mean_normalized_control`、`log2fc`、`percent_knockdown`。靶基因不在矩阵中时 `target_in_matrix` 为 false 并给出警告。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import analyze_crispr_perturbations

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
perturbations = {"c1": "control", "c2": "control", "c3": "control", "c4": "CD3D", "c5": "CD3D", "c6": "CD3D"}

result = analyze_crispr_perturbations(counts, options={"perturbations": perturbations})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "analyze_crispr_perturbations",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "control": "control",
    "n_control_cells": 3,
    "cells_per_label": {
      "CD3D": 3,
      "control": 3
    },
    "target_sum": 10000.0,
    "knockdown_definition": "100 * (1 - mean CP10k(guide cells) / mean CP10k(control cells)) of the target gene",
    "guides": {
      "CD3D": {
        "n_cells": 3,
        "target_gene": "CD3D",
        "target_in_matrix": true,
        "mean_normalized_guide": 61.7283950617284,
        "mean_normalized_control": 1363.1378477316855,
        "log2fc": -4.442726976505825,
        "percent_knockdown": 95.47159554226693
      }
    }
  }
}
```

## 2) annotate_cells

- **作用：** 用用户提供的标志基因打分，为细胞标注细胞类型（不使用训练模型）。

- **处理方法：** 取 ln(1 + CP10k) 值。对每个标志基因，在所有细胞间做 z 分数（总体标准差；常数基因得 0）。某细胞类型在某细胞中的得分为该类型标志基因 z 分数的均值，细胞被赋予得分最高的类型；最高分不大于 `min_score` 时为 `"Unassigned"`。矩阵中缺失的标志基因会给出警告（某类型的标志基因全部缺失则报错）。不使用训练模型。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.annotate_cells(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`markers`（必填；`{细胞类型: [标志基因]}`）、`min_score`（0.0）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`assignments`（{cell: 类型}）、`best_score`（{cell: 得分}）、`counts`（每种类型的细胞数）与 `scores`（{类型: {cell: 得分}}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import annotate_cells

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
markers = {"T cell": ["CD3D", "CD3E"], "B cell": ["CD79A", "MS4A1"]}

result = annotate_cells(counts, options={"markers": markers})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "annotate_cells",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score(type, cell) = mean over the type's marker genes of the per-gene z-score of ln(1 + CP10k); the highest score wins if > min_score",
    "assignments": {
      "c1": "T cell",
      "c2": "T cell",
      "c3": "T cell",
      "c4": "B cell",
      "c5": "B cell",
      "c6": "B cell"
    },
    "best_score": {
      "c1": 0.8356997436691734,
      "c2": 0.8567110675970135,
      "c3": 0.8259936626183231,
      "c4": 0.8350643975343444,
      "c5": 0.8410178115070588,
      "c6": 0.8413837216201503
    },
    "counts": {
      "B cell": 3,
      "T cell": 3
    },
    "scores": {
      "T cell": {
        "c1": 0.8356997436691734,
        "c2": 0.8567110675970135,
        "c3": 0.8259936626183231,
        "c4": -0.5599901482596032,
        "c5": -0.5763799279190517,
        "c6": -1.3820343977058553
      },
      "B cell": {
        "c1": -0.5632664728908073,
        "c2": -0.5716923966489433,
        "c3": -1.3825070611218027,
        "c4": 0.8350643975343444,
        "c5": 0.8410178115070588,
        "c6": 0.8413837216201503
      }
    }
  }
}
```

## 3) annotate_cells_spatial

- **作用：** 用用户提供的标志基因打分，为空间斑点标注细胞类型（不使用训练模型）。

- **处理方法：** 取 ln(1 + CP10k) 值。对每个标志基因，在所有细胞间做 z 分数（总体标准差；常数基因得 0）。某细胞类型在某细胞中的得分为该类型标志基因 z 分数的均值，细胞被赋予得分最高的类型；最高分不大于 `min_score` 时为 `"Unassigned"`。矩阵中缺失的标志基因会给出警告（某类型的标志基因全部缺失则报错）。不使用训练模型。纯 Python 实现。 定义与 `annotate_cells` 相同，斑点相当于细胞。

- **API：** `dnakit.singlecell_spatial.annotate_cells_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`markers`（必填；`{细胞类型: [标志基因]}`）、`min_score`（0.0）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`assignments`（{cell: 类型}）、`best_score`（{cell: 得分}）、`counts`（每种类型的细胞数）与 `scores`（{类型: {cell: 得分}}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import annotate_cells_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
markers = {"T cell": ["CD3D", "CD3E"], "B cell": ["CD79A", "MS4A1"]}

result = annotate_cells_spatial(counts, options={"markers": markers})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "annotate_cells_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score(type, cell) = mean over the type's marker genes of the per-gene z-score of ln(1 + CP10k); the highest score wins if > min_score",
    "assignments": {
      "c1": "T cell",
      "c2": "T cell",
      "c3": "T cell",
      "c4": "B cell",
      "c5": "B cell",
      "c6": "B cell"
    },
    "best_score": {
      "c1": 0.8356997436691734,
      "c2": 0.8567110675970135,
      "c3": 0.8259936626183231,
      "c4": 0.8350643975343444,
      "c5": 0.8410178115070588,
      "c6": 0.8413837216201503
    },
    "counts": {
      "B cell": 3,
      "T cell": 3
    },
    "scores": {
      "T cell": {
        "c1": 0.8356997436691734,
        "c2": 0.8567110675970135,
        "c3": 0.8259936626183231,
        "c4": -0.5599901482596032,
        "c5": -0.5763799279190517,
        "c6": -1.3820343977058553
      },
      "B cell": {
        "c1": -0.5632664728908073,
        "c2": -0.5716923966489433,
        "c3": -1.3825070611218027,
        "c4": 0.8350643975343444,
        "c5": 0.8410178115070588,
        "c6": 0.8413837216201503
      }
    }
  }
}
```

## 4) build_cell_gene_matrix

- **作用：** 由逐读段或逐 UMI 的（细胞，基因）归属记录构建基因×细胞计数矩阵。

- **处理方法：** 记录带 UMI 时，（细胞，基因）的计数为不同 UMI 的个数（可先用 `correct_umis` 合并测序错误）；否则为记录条数，或可选 `count` 列之和。基因与细胞按字母序排列。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.build_cell_gene_matrix(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** TSV/CSV 文本或文件路径形式的记录（每行一条：cell, gene[, umi]；可选表头，列名为 `cell`/`barcode`、`gene`/`feature`、`umi`、`count`）、`[cell, gene[, umi]]` 列表的列表，或带这些键的映射列表。选项：`max_matrix_elements`（5000000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_records`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import build_cell_gene_matrix

records = """cell\tgene\tumi
c1\tCD3D\tAAAA
c1\tCD3D\tCCCC
c1\tCD3D\tCCCC
c1\tACTB\tGGGG
c2\tCD3D\tAAAA
c2\tCD79A\tTTTT
"""

result = build_cell_gene_matrix(records)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "build_cell_gene_matrix",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 7,
    "text_length": 93,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1a65abffca3e87366af0ad395ad70b1ccd7cd20630b9717e39096fb6da41429d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "counts = number of distinct UMIs per (cell, gene)",
    "n_records": 6,
    "n_genes": 3,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "ACTB",
        "CD3D",
        "CD79A"
      ],
      "cells": [
        "c1",
        "c2"
      ],
      "values": [
        [
          1.0,
          0.0
        ],
        [
          2.0,
          1.0
        ],
        [
          0.0,
          1.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 5) build_cell_peak_matrix

- **作用：** 由 scATAC 片段和峰列表构建峰×细胞计数矩阵。

- **处理方法：** 片段为带细胞条形码的类 BED 区间（0 起始、半开）。`count: "fragments"` 时，一个片段与某个峰至少重叠 1 bp 就对该峰计一次（按可选的第 5 列加权）；`"cut_sites"` 时，片段的首、末碱基落入峰内才计数。峰为矩阵的行，命名为 `chrom:start-end`。位于没有峰的染色体上的片段只会把其条形码加入细胞列表。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.build_cell_peak_matrix(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** TSV 文本或文件路径形式的片段（`chrom, start, end, barcode[, count]`；忽略 `#` 开头的行），或此类列表的列表。选项：`peaks`（必填；`"chr1:100-600"` 字符串、`[chrom, start, end]` 列表或 BED 文本）、`count`（默认 `"fragments"`，或 `"cut_sites"`）、`max_matrix_elements`（5000000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_fragments`、`mode`、`n_peaks`（峰的个数，即矩阵行数；此字段原名 `n_genes`）、`n_cells`、`matrix`（`orientation`、`genes`（存放峰名）、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import build_cell_peak_matrix

fragments = """chr1\t100\t300\tAAAC\t1
chr1\t150\t250\tAAAC\t1
chr1\t900\t1000\tAAAC\t1
chr1\t120\t220\tTTTG\t1
chr2\t50\t150\tTTTG\t1
"""
peaks = ["chr1:100-600", "chr1:800-1200", "chr2:0-100"]

result = build_cell_peak_matrix(fragments, options={"peaks": peaks})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "build_cell_peak_matrix",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4f988916a04d2fe97bb9a93fe5194995673409cb28d6eef8bb12999615188d02"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "counts = fragments (weighted by the optional 5th column) overlapping each peak by >= 1 bp; each fragment counts once per peak",
    "n_fragments": 5,
    "mode": "fragments",
    "n_peaks": 3,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "chr1:100-600",
        "chr1:800-1200",
        "chr2:0-100"
      ],
      "cells": [
        "AAAC",
        "TTTG"
      ],
      "values": [
        [
          2.0,
          1.0
        ],
        [
          1.0,
          0.0
        ],
        [
          0.0,
          1.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 6) calculate_cell_cycle_score

- **作用：** 用用户提供的基因集为每个细胞计算 S 期与 G2/M 期得分，并判定 G1、S 或 G2M。

- **处理方法：** 取 ln(1 + CP10k) 值。基因集得分 = 集合内基因的平均表达 - 表达量匹配的对照基因的平均表达：对照基因取自平均表达秩相同的分箱（`n_bins`，优先取秩最近的，每个集合基因最多 `ctrl_size` 个）。对照基因的选取是确定性的，因此数值与 scanpy 的随机对照略有差异。两个得分都小于 0 时为 `G1`，否则取得分较大者（`S` 或 `G2M`）。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.calculate_cell_cycle_score(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`s_genes`（必填）、`g2m_genes`（必填）、`n_bins`（25）、`ctrl_size`（50）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`s_score` 与 `g2m_score`（{cell: 得分}）、`phase`（{cell: `G1`/`S`/`G2M`}）和 `phase_counts`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import calculate_cell_cycle_score

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
options = {"s_genes": ["CD3D", "CD3E"], "g2m_genes": ["CD79A", "MS4A1"], "n_bins": 3, "ctrl_size": 2}

result = calculate_cell_cycle_score(counts, options=options)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_cell_cycle_score",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score = mean ln(1+CP10k) of the gene set minus mean of expression-matched control genes (genes in th...",
    "s_score": {
      "c1": 4.309330241732457,
      "c2": 4.2525698094105415,
      "c3": 4.305041081459723,
      "c4": -2.6471780733568737,
      "c5": -3.891133860838857,
      "c6": -5.096677528135468
    },
    "g2m_score": {
      "c1": -2.4136427067829653,
      "c2": -2.537043500514193,
      "c3": -6.331498665182803,
      "c4": 4.257662662904112,
      "c5": 3.040942474661665,
      "c6": 5.755455155631239
    },
    "phase": {
      "c1": "S",
      "c2": "S",
      "c3": "S",
      "c4": "G2M",
      "c5": "G2M",
      "c6": "G2M"
    },
    "phase_counts": {
      "G1": 0,
      "S": 3,
      "G2M": 3
    }
  }
}
```

## 7) calculate_cell_interactions

- **作用：** 根据平均表达量为细胞类型之间的配体-受体相互作用打分（不做显著性检验）。

- **处理方法：** 取 ln(1 + CP10k) 值。对每一个有序细胞类型对（发送者，接收者）和每个给定的配体-受体对，`score = 配体在发送细胞中的均值 * 受体在接收细胞中的均值`；`min_pct` 用于剔除配体在发送细胞中、或受体在接收细胞中表达比例低于该值的相互作用。含缺失基因的配对会被跳过并给出警告。不内置配体-受体数据库，必须自行提供配对，也不做置换等显著性检验。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.calculate_cell_interactions(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`cell_types`（必填；`{cell: 类型}` 或按列顺序的列表）、`pairs`（必填；`[[配体, 受体], ...]`）、`min_pct`（0.0）、`top_n`（100）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`no_significance_test`（true）、`n_interactions`（`top_n` 截断之前的个数）与 `interactions`；`interactions` 按 `score` 降序排列（最多 `top_n` 项），每项含 `ligand`、`receptor`、`sender`、`receiver`、`ligand_mean`、`receptor_mean`、`score`、`pct_ligand`、`pct_receptor`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import calculate_cell_interactions

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
cell_types = ["T", "T", "T", "B", "B", "B"]
pairs = [["CD3D", "MS4A1"], ["CD79A", "CD3E"]]

result = calculate_cell_interactions(counts, options={"cell_types": cell_types, "pairs": pairs})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_cell_interactions",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score = mean ln(1+CP10k) of the ligand in sender cells x mean of the receptor in receiver cells",
    "no_significance_test": true,
    "n_interactions": 8,
    "interactions": [
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "T",
        "receiver": "B",
        "ligand_mean": 7.212575686641653,
        "receptor_mean": 7.161175211635025,
        "score": 51.65051821921968,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0
      },
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "B",
        "receiver": "T",
        "ligand_mean": 7.296230322614124,
        "receptor_mean": 7.056039166585788,
        "score": 51.48248692479612,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0
      },
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "B",
        "receiver": "B",
        "ligand_mean": 7.296230322614124,
        "receptor_mean": 1.7422472658960821,
        "score": 12.711837330922547,
        "pct_ligand": 1.0,
        "pct_receptor": 0.3333333333333333
      },
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "T",
        "receiver": "T",
        "ligand_mean": 7.212575686641653,
        "receptor_mean": 1.7612005927668049,
        "score": 12.702792574688724,
        "pct_ligand": 1.0,
        "pct_receptor": 0.3333333333333333
      },
      "..."
    ]
  }
}
```

## 8) calculate_cell_interactions_spatial

- **作用：** 为距离在给定半径内的相邻斑点（或细胞）之间的配体-受体相互作用打分。

- **处理方法：** 取 ln(1 + CP10k) 值。对彼此距离不超过 `radius` 的每一个有序斑点对 (a, b)（a、b 不同），计算 `配体在 a 中的值 * 受体在 b 中的值`；`score` 为按（发送者类型，接收者类型）分组后这些乘积的均值，`n_cell_pairs` 为对应的斑点对数。含缺失基因的配对会被跳过并给出警告；没有相邻斑点对的类型对不产生结果行。不内置配体-受体数据库，也不做显著性检验。纯 Python 实现；最多 5000 个斑点。

- **API：** `dnakit.singlecell_spatial.calculate_cell_interactions_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`cell_types`（必填）、`pairs`（必填；`[[配体, 受体], ...]`）、`radius`（必填；发送者与接收者的最大距离，单位同坐标）、`coordinates`（必填；`{spot: [x, y]}`、按列顺序的 `[x, y]` 列表，或 CSV/TSV/Visium `tissue_positions` 文件）、`min_pct`（0.0）、`top_n`（100）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`no_significance_test`（true）、`n_interactions`（`top_n` 截断之前的个数）与 `interactions`；`interactions` 按 `score` 降序排列（最多 `top_n` 项），每项含 `ligand`、`receptor`、`sender`、`receiver`、`ligand_mean`、`receptor_mean`、`score`、`pct_ligand`、`pct_receptor`与 `n_cell_pairs`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import calculate_cell_interactions_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}
cell_types = ["T", "T", "T", "B", "B", "B"]
pairs = [["CD3D", "MS4A1"], ["CD79A", "CD3E"]]

result = calculate_cell_interactions_spatial(
    counts,
    options={"cell_types": cell_types, "pairs": pairs, "radius": 1.5, "coordinates": coordinates},
)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_cell_interactions_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score = mean over sender-receiver cell pairs within `radius` of ln(1+CP10k)(ligand in sender) x ln(1+CP10k)(receptor in receiver)",
    "no_significance_test": true,
    "n_interactions": 8,
    "interactions": [
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "B",
        "receiver": "T",
        "ligand_mean": 7.296979225009493,
        "receptor_mean": 7.085909282507124,
        "score": 51.70338903539261,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0,
        "n_cell_pairs": 7
      },
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "T",
        "receiver": "B",
        "ligand_mean": 7.198078570902774,
        "receptor_mean": 7.1621555536358175,
        "score": 51.54990974557287,
        "pct_ligand": 1.0,
        "pct_receptor": 1.0,
        "n_cell_pairs": 7
      },
      {
        "ligand": "CD3D",
        "receptor": "MS4A1",
        "sender": "B",
        "receiver": "B",
        "ligand_mean": 2.6133708988441233,
        "receptor_mean": 7.16289081013641,
        "score": 18.705839909215634,
        "pct_ligand": 0.5,
        "pct_receptor": 1.0,
        "n_cell_pairs": 4
      },
      {
        "ligand": "CD79A",
        "receptor": "CD3E",
        "sender": "T",
        "receiver": "T",
        "ligand_mean": 2.661704053969166,
        "receptor_mean": 7.108311869448126,
        "score": 18.502819124425383,
        "pct_ligand": 0.5,
        "pct_receptor": 1.0,
        "n_cell_pairs": 4
      },
      "..."
    ]
  }
}
```

## 9) calculate_markers

- **作用：** 对每个分组（聚类或细胞类型），通过与其余所有细胞比较来寻找标志基因。

- **处理方法：** 每个分组作为第一组，其余所有细胞作为第二组。对每个基因：`norm = count / 文库大小 * target_sum`（`normalization: "none"` 时直接使用输入值），检验在 `log2(norm + pseudocount)` 上进行，`log2fc = log2(mean(norm, 第一组) + pseudocount) - log2(mean(norm, 第二组) + pseudocount)`，`pct_*` 为原始计数 > 0 的细胞比例，`p_value` 来自 Welch 不等方差 t 检验（默认）或 Mann-Whitney U / Wilcoxon 秩和检验（正态近似，含并列校正与连续性校正），`q_value` 为对所有检验基因做的 Benjamini-Hochberg 校正。基因按 p 值升序排列（并列时 |log2fc| 大者在前）。`only_positive`（默认 true）时只保留 `log2fc > 0` 的基因。使用 `welch` 时，细胞数少于 2 的分组（或其余细胞少于 2 个）会被跳过并给出警告（警告中指明分组）；只有当没有任何分组可检验时才抛出 `GROUPS_TOO_SMALL`，只有一个标签时抛出 `GROUPS_TOO_FEW`。`find_differential_features` 与 `transcriptomics.calculate_differential_expression` 的行为不同：遇到过小的分组会抛出 `InputFormatError`（`GROUP_TOO_SMALL`）。每个基因行中，`mean_group1`/`pct_group1` 描述被检验的分组，`mean_group2`/`pct_group2` 描述所有其他细胞。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.calculate_markers(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`groups`（必填；`{cell: 标签}` 或按列顺序的列表）、`only_positive`（true）、`method`（默认 `"welch"`，或 `"wilcoxon"`）、`normalization`（默认 `"cpm"` = 将每个细胞缩放到 `target_sum`，或 `"none"`）、`target_sum`（10000）、`pseudocount`（1.0）、`top_n`（100）、`min_pct`（0.0；保留任一组中 pct 不低于该值的基因）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`method`、`only_positive`、`comparison` 与 `markers`；`markers` 为 分组 -> `n_cells` 与 `genes`（按排名取 top_n 个基因；每个基因行含 `gene`、`mean_group1`、`mean_group2`（归一化值的均值）、`log2fc`、`pct_group1`、`pct_group2`、`statistic`、`p_value` 与 `q_value`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import calculate_markers

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
groups = ["T", "T", "T", "B", "B", "B"]

result = calculate_markers(counts, options={"groups": groups, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_markers",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "welch",
    "only_positive": true,
    "comparison": "each group versus all other cells",
    "markers": {
      "B": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD79A",
            "mean_group1": 1481.4814814814815,
            "mean_group2": 68.02721088435375,
            "log2fc": 4.424705081255279,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.1091841711508663,
            "p_value": 0.0893622027560642,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "MS4A1",
            "mean_group1": 1296.2962962962963,
            "mean_group2": 65.35947712418302,
            "log2fc": 4.289061692972825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.062672412366318,
            "p_value": 0.09161271184455913,
            "q_value": 0.18438518211850805
          }
        ]
      },
      "T": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD3D",
            "mean_group1": 1363.1378477316855,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.442726976505825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.136919192646547,
            "p_value": 0.08794875669662708,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "CD3E",
            "mean_group1": 1173.3582321817617,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.226610101334473,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.0438371186837414,
            "p_value": 0.09219259105925402,
            "q_value": 0.18438518211850805
          }
        ]
      }
    }
  }
}
```

## 10) calculate_markers_spatial

- **作用：** 对每个空间区域（或斑点聚类），通过与其余所有斑点比较来寻找标志基因。

- **处理方法：** 方法与 `calculate_markers` 相同，斑点相当于细胞。对每个基因：`norm = count / 文库大小 * target_sum`（`normalization: "none"` 时直接使用输入值），检验在 `log2(norm + pseudocount)` 上进行，`log2fc = log2(mean(norm, 第一组) + pseudocount) - log2(mean(norm, 第二组) + pseudocount)`，`pct_*` 为原始计数 > 0 的细胞比例，`p_value` 来自 Welch 不等方差 t 检验（默认）或 Mann-Whitney U / Wilcoxon 秩和检验（正态近似，含并列校正与连续性校正），`q_value` 为对所有检验基因做的 Benjamini-Hochberg 校正。基因按 p 值升序排列（并列时 |log2fc| 大者在前）。`only_positive`（默认 true）时只保留 `log2fc > 0` 的基因。与 `calculate_markers` 相同，使用 `welch` 时，斑点数少于 2 的区域（或其余斑点少于 2 个）会被跳过并给出警告，只有当没有任何区域可检验时才抛出 `GROUPS_TOO_SMALL`；而 `find_differential_features_spatial` 会抛出 `GROUP_TOO_SMALL`。`mean_group1`/`pct_group1` 描述被检验的区域，`mean_group2`/`pct_group2` 描述所有其他斑点。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.calculate_markers_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`groups`（必填；`{spot: 区域标签}` 或按列顺序的列表）、`only_positive`（true）、`method`（默认 `"welch"`，或 `"wilcoxon"`）、`normalization`（默认 `"cpm"` = 将每个细胞缩放到 `target_sum`，或 `"none"`）、`target_sum`（10000）、`pseudocount`（1.0）、`top_n`（100）、`min_pct`（0.0；保留任一组中 pct 不低于该值的基因）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`method`、`only_positive`、`comparison` 与 `markers`；`markers` 为 区域 -> `n_cells` 与 `genes`（按排名取 top_n 个基因；每个基因行含 `gene`、`mean_group1`、`mean_group2`（归一化值的均值）、`log2fc`、`pct_group1`、`pct_group2`、`statistic`、`p_value` 与 `q_value`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import calculate_markers_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
domains = ["tissue1", "tissue1", "tissue1", "tissue2", "tissue2", "tissue2"]

result = calculate_markers_spatial(counts, options={"groups": domains, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_markers_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "welch",
    "only_positive": true,
    "comparison": "each group versus all other cells",
    "markers": {
      "tissue1": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD3D",
            "mean_group1": 1363.1378477316855,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.442726976505825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.136919192646547,
            "p_value": 0.08794875669662708,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "CD3E",
            "mean_group1": 1173.3582321817617,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.226610101334473,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.0438371186837414,
            "p_value": 0.09219259105925402,
            "q_value": 0.18438518211850805
          }
        ]
      },
      "tissue2": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD79A",
            "mean_group1": 1481.4814814814815,
            "mean_group2": 68.02721088435375,
            "log2fc": 4.424705081255279,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.1091841711508663,
            "p_value": 0.0893622027560642,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "MS4A1",
            "mean_group1": 1296.2962962962963,
            "mean_group2": 65.35947712418302,
            "log2fc": 4.289061692972825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.062672412366318,
            "p_value": 0.09161271184455913,
            "q_value": 0.18438518211850805
          }
        ]
      }
    }
  }
}
```

## 11) calculate_perturbation_effect

- **作用：** 对每种扰动，检验其相对对照细胞的全转录组效应。

- **处理方法：** 每个非对照标签逐基因地与对照细胞比较，所用双组检验与 `find_differential_features` 相同（扰动为第一组，对照为第二组）。对每个基因：`norm = count / 文库大小 * target_sum`（`normalization: "none"` 时直接使用输入值），检验在 `log2(norm + pseudocount)` 上进行，`log2fc = log2(mean(norm, 第一组) + pseudocount) - log2(mean(norm, 第二组) + pseudocount)`，`pct_*` 为原始计数 > 0 的细胞比例，`p_value` 来自 Welch 不等方差 t 检验（默认）或 Mann-Whitney U / Wilcoxon 秩和检验（正态近似，含并列校正与连续性校正），`q_value` 为对所有检验基因做的 Benjamini-Hochberg 校正。基因按 p 值升序排列（并列时 |log2fc| 大者在前）。使用 `welch` 时，细胞数少于 2 的扰动（或对照细胞少于 2 个）会被跳过并给出警告；若没有任何扰动可检验则报错。每个基因行中，`_group1` 字段（`mean_group1`、`pct_group1`）描述被扰动的细胞，`_group2` 字段描述对照细胞；`log2fc` 为正表示在被扰动细胞中更高（扰动对对照，而不是某一组对所有其他细胞）。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.calculate_perturbation_effect(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`perturbations`（必填；`{cell: 扰动标签}` 或列表）、`control`（`"control"`）、`method`（默认 `"welch"`，或 `"wilcoxon"`）、`normalization`（默认 `"cpm"` = 将每个细胞缩放到 `target_sum`，或 `"none"`）、`target_sum`（10000）、`pseudocount`（1.0）、`top_n`（100）、`min_pct`（0.0；保留任一组中 pct 不低于该值的基因）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`control`、`n_control_cells`、`comparison`、`method` 与 `effects`；`effects` 为 扰动 -> `n_cells` 与 `genes`（按排名取 top_n 个基因；每个基因行含 `gene`、`mean_group1`、`mean_group2`（归一化值的均值）、`log2fc`、`pct_group1`、`pct_group2`、`statistic`、`p_value` 与 `q_value`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import calculate_perturbation_effect

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
perturbations = ["KO", "KO", "KO", "control", "control", "control"]

result = calculate_perturbation_effect(counts, options={"perturbations": perturbations, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "calculate_perturbation_effect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "control": "control",
    "n_control_cells": 3,
    "comparison": "each perturbation versus control cells",
    "method": "welch",
    "effects": {
      "KO": {
        "n_cells": 3,
        "genes": [
          {
            "gene": "CD3D",
            "mean_group1": 1363.1378477316855,
            "mean_group2": 61.7283950617284,
            "log2fc": 4.442726976505825,
            "pct_group1": 1.0,
            "pct_group2": 0.3333333333333333,
            "statistic": 3.136919192646547,
            "p_value": 0.08794875669662708,
            "q_value": 0.18438518211850805
          },
          {
            "gene": "CD79A",
            "mean_group1": 68.02721088435375,
            "mean_group2": 1481.4814814814815,
            "log2fc": -4.424705081255279,
            "pct_group1": 0.3333333333333333,
            "pct_group2": 1.0,
            "statistic": -3.1091841711508663,
            "p_value": 0.0893622027560642,
            "q_value": 0.18438518211850805
          }
        ]
      }
    }
  }
}
```

## 12) cluster

- **作用：** 对对数归一化表达量的前若干主成分做 k-means，对细胞聚类。

- **处理方法：** 表达量取 ln(1 + CP10k)（每个细胞缩放到 `target_sum` = 10000 个计数，再取 ln(1 + x)）（取 ln(1 + CP10k) 值并对每个基因中心化（可选缩放到单位方差，并限于方差最高的 `n_top_genes` 个基因）；主成分来自该矩阵的 SVD，并固定每个主成分的符号使其绝对值最大的载荷为正）。k-means 采用确定性的最远点初始化和 Lloyd 迭代；标签按首次出现顺序编号。`n_clusters: "auto"` 时在 2 到 `max_clusters` 之间选平均轮廓系数最高的 k。需要 numpy（缺失时返回 `not_implemented`）；最多 20000 个细胞（`"auto"` 时最多 5000 个）。

- **API：** `dnakit.singlecell_spatial.cluster(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`n_clusters`（默认 `"auto"`，或整数）、`max_clusters`（10）、`n_pcs`（10；0 表示不做 PCA）、`n_top_genes`（默认全部，否则取方差最高的基因）、`scale`（false）、`target_sum`（10000）、`return_embedding`（false）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_clusters`、`labels`（{cell: 字符串形式的聚类编号，按首次出现顺序编号}）、`cluster_sizes`、`inertia`、`explained_variance_ratio`（`n_pcs` 为 0 时为 `null`）、`silhouette_by_k`（`n_clusters` 为 `"auto"` 时给出各候选 k 的平均轮廓系数，否则为 `null`），以及在 `return_embedding` 时的 `embedding`（{cell: 坐标}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import cluster

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""

result = cluster(counts, options={"n_clusters": 2, "n_pcs": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "cluster",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-means (deterministic farthest-point start, Lloyd iterations) on the top PCs of ln(1+CP10k) values; k from options or the best mean silhouette",
    "n_clusters": 2,
    "labels": {
      "c1": "0",
      "c2": "0",
      "c3": "0",
      "c4": "1",
      "c5": "1",
      "c6": "1"
    },
    "cluster_sizes": {
      "0": 3,
      "1": 3
    },
    "inertia": 61.986183555554646,
    "explained_variance_ratio": [
      0.650159103044713,
      0.17804888640013142
    ],
    "silhouette_by_k": null
  }
}
```

## 13) cluster_spatial

- **作用：** 对经空间平滑的主成分做 k-means，把斑点聚成空间区域。

- **处理方法：** 启发式方法，并非已发表的方法：先按 `cluster` 的方式计算主成分，再把每个斑点的主成分向量替换为 `(1 - alpha) * 自身 + alpha * 其 n_neighbors 个空间最近邻的均值`（按坐标的欧氏距离取最近邻），然后在平滑后的向量上做 k-means（`n_clusters`/轮廓系数规则同 `cluster`）。需要 numpy（缺失时返回 `not_implemented`）；最多 5000 个斑点。

- **API：** `dnakit.singlecell_spatial.cluster_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`coordinates`（必填；`{spot: [x, y]}`、按列顺序的 `[x, y]` 列表，或 CSV/TSV/Visium `tissue_positions` 文件）、`n_neighbors`（6；须小于斑点数）、`alpha`（0.5）、`n_clusters`（`"auto"` 或整数）、`max_clusters`（10）、`n_pcs`（10；0 表示不做 PCA）、`n_top_genes`（默认全部，否则取方差最高的基因）、`scale`（false）、`target_sum`（10000）、`return_embedding`（false）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_clusters`、`labels`（{cell: 字符串形式的聚类编号，按首次出现顺序编号}）、`cluster_sizes`、`inertia`、`explained_variance_ratio`（`n_pcs` 为 0 时为 `null`）、`silhouette_by_k`（`n_clusters` 为 `"auto"` 时给出各候选 k 的平均轮廓系数，否则为 `null`），以及在 `return_embedding` 时的 `embedding`（{cell: 坐标}）（`inertia` 与 `embedding` 针对平滑后的向量）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import cluster_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}

result = cluster_spatial(counts, options={"coordinates": coordinates, "n_clusters": 2, "n_pcs": 2, "n_neighbors": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "cluster_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "spatial-domain heuristic: each spot's PC vector is replaced by (1-alpha)*own + alpha*mean of its 2 n...",
    "n_clusters": 2,
    "labels": {
      "c1": "0",
      "c2": "0",
      "c3": "0",
      "c4": "1",
      "c5": "1",
      "c6": "1"
    },
    "cluster_sizes": {
      "0": 3,
      "1": 3
    },
    "inertia": 11.41227429529395,
    "explained_variance_ratio": [
      0.650159103044713,
      0.17804888640013142
    ],
    "silhouette_by_k": null
  }
}
```

## 14) correct_umis

- **作用：** 围绕高丰度 UMI 合并 UMI 组，以校正 UMI 中的测序错误（UMI-tools 的 directional、adjacency、cluster 算法）。

- **处理方法：** 在汉明距离（至多 `threshold` 个错配）上实现 UMI-tools 的网络方法。`directional`：从最丰富的 UMI 开始，吸收满足 `count_a >= 2 * count_b - 1` 的邻居 `b`，并可传递；`adjacency`：最丰富的未分配 UMI 一步吸收其所有未分配邻居；`cluster`：连通分量，以最丰富的 UMI 为代表；`unique`：不合并。所有 UMI 长度必须相同，且只能含 A、C、G、T、N。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.correct_umis(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** UMI 字符串列表（每次出现计 1 条读段）、`{UMI: 读段数}`，或每行一个 `UMI` 或 `UMI<TAB>计数` 的文本/文件。选项：`method`（默认 `"directional"`，另有 `"adjacency"`、`"cluster"`、`"unique"`）、`threshold`（1；最大汉明距离，1 到 4）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`method`、`threshold`、`n_reads`、`n_unique_before`、`n_unique_after`、`corrected`（{被合并的 UMI: 代表 UMI}）与 `groups`（每组的代表、成员和读段数）及 `groups_truncated`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import correct_umis

umis = ["AAAA", "AAAA", "AAAA", "AAAA", "AAAT", "CCCC", "CCCC", "GGGG"]

result = correct_umis(umis)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "correct_umis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 8,
    "text_length": 64,
    "sequence_bases": 32,
    "lengths": [4, 4, 4, 4, 4, 4, 4, 4],
    "sha256": "7051d7f210a596683b53aab27d96e3402faa226fa42320cbb7d67085c757c7ec"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "UMI-tools directional: starting from the most abundant UMI, absorb neighbours within the Hamming thr...",
    "method": "directional",
    "threshold": 1,
    "n_reads": 8,
    "n_unique_before": 4,
    "n_unique_after": 3,
    "corrected": {
      "AAAT": "AAAA"
    },
    "groups": [
      {
        "representative": "AAAA",
        "members": [
          "AAAA",
          "AAAT"
        ],
        "reads": 5
      },
      {
        "representative": "CCCC",
        "members": [
          "CCCC"
        ],
        "reads": 2
      },
      {
        "representative": "GGGG",
        "members": [
          "GGGG"
        ],
        "reads": 1
      }
    ],
    "groups_truncated": false
  }
}
```

## 15) correct_umis_spatial

- **作用：** 围绕高丰度 UMI 合并 UMI 组，以校正空间条形码文库中 UMI 的测序错误（UMI-tools 算法）。

- **处理方法：** 在汉明距离（至多 `threshold` 个错配）上实现 UMI-tools 的网络方法。`directional`：从最丰富的 UMI 开始，吸收满足 `count_a >= 2 * count_b - 1` 的邻居 `b`，并可传递；`adjacency`：最丰富的未分配 UMI 一步吸收其所有未分配邻居；`cluster`：连通分量，以最丰富的 UMI 为代表；`unique`：不合并。所有 UMI 长度必须相同，且只能含 A、C、G、T、N。纯 Python 实现。 算法与 `correct_umis` 相同；应对同一斑点、同一基因的 UMI 调用。

- **API：** `dnakit.singlecell_spatial.correct_umis_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** UMI 字符串列表、`{UMI: 读段数}`，或每行一个 `UMI` 或 `UMI<TAB>计数` 的文本/文件（同一斑点、同一基因的 UMI）。选项：`method`（默认 `"directional"`，另有 `"adjacency"`、`"cluster"`、`"unique"`）、`threshold`（1；最大汉明距离，1 到 4）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`method`、`threshold`、`n_reads`、`n_unique_before`、`n_unique_after`、`corrected`（{被合并的 UMI: 代表 UMI}）与 `groups`（每组的代表、成员和读段数）及 `groups_truncated`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import correct_umis_spatial

umi_counts = {"ACGTAC": 10, "ACGTAA": 1, "TTTTTT": 3}

result = correct_umis_spatial(umi_counts)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "correct_umis_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 40,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5b6db5fa8fe58914529df0bf23f8e0e2578f99c5ee8538e1b252c9973fa8c775"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "UMI-tools directional: starting from the most abundant UMI, absorb neighbours within the Hamming thr...",
    "method": "directional",
    "threshold": 1,
    "n_reads": 14,
    "n_unique_before": 3,
    "n_unique_after": 2,
    "corrected": {
      "ACGTAA": "ACGTAC"
    },
    "groups": [
      {
        "representative": "ACGTAC",
        "members": [
          "ACGTAC",
          "ACGTAA"
        ],
        "reads": 11
      },
      {
        "representative": "TTTTTT",
        "members": [
          "TTTTTT"
        ],
        "reads": 3
      }
    ],
    "groups_truncated": false
  }
}
```

## 16) deconvolute_spots

- **作用：** 在签名矩阵上用非负最小二乘估计每个斑点的细胞类型比例。

- **处理方法：** 在与签名矩阵共有的基因上，用非负最小二乘（`scipy.optimize.nnls`）把每个斑点的表达谱拟合为签名列的非负组合；系数除以其和即得细胞类型比例，`rmse` 为拟合残差。`normalize`（默认 true）时先把每一列缩放到总和为 1。若给出 `reference` 与 `reference_labels` 而非 `signature`，则签名为每种细胞类型的平均文库大小归一化表达谱。需要 numpy 和 scipy（缺失时返回 `not_implemented`）。它不是 RCTD 或 cell2location 之类基于模型的方法。

- **API：** `dnakit.singlecell_spatial.deconvolute_spots(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 待解卷积的斑点：基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`signature`（基因×细胞类型矩阵，可用任一矩阵形式）或 `reference`（带注释细胞的计数矩阵）加 `reference_labels`（`{cell: 类型}` 或列表）、`normalize`（true）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`cell_types`、`n_shared_genes`、`fractions`（{斑点: {细胞类型: 比例}}）与 `rmse`（{斑点: 残差}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import deconvolute_spots

spots = """gene,s1,s2
A,100,10
B,10,100
C,50,50
"""
signature = """gene,typeX,typeY
A,10,1
B,1,10
C,5,5
"""

result = deconvolute_spots(spots, options={"signature": signature})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "deconvolute_spots",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 37,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d0743d02cb300ebb3b2966930e2cef0eac47e89410bcaa86785506703f7c27f3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "non-negative least squares of each mixture column on the signature columns (scipy.optimize.nnls); fr...",
    "cell_types": [
      "typeX",
      "typeY"
    ],
    "n_shared_genes": 3,
    "fractions": {
      "s1": {
        "typeX": 1.0,
        "typeY": 4.636616585542447e-17
      },
      "s2": {
        "typeX": 0.0,
        "typeY": 1.0
      }
    },
    "rmse": {
      "s1": 3.0913372807800604e-17,
      "s2": 2.1817011766850644e-17
    }
  }
}
```

## 17) deconvolute_spots_spatial

- **作用：** 在签名矩阵上用非负最小二乘估计每个空间斑点的细胞类型比例。

- **处理方法：** 在与签名矩阵共有的基因上，用非负最小二乘（`scipy.optimize.nnls`）把每个斑点的表达谱拟合为签名列的非负组合；系数除以其和即得细胞类型比例，`rmse` 为拟合残差。`normalize`（默认 true）时先把每一列缩放到总和为 1。若给出 `reference` 与 `reference_labels` 而非 `signature`，则签名为每种细胞类型的平均文库大小归一化表达谱。需要 numpy 和 scipy（缺失时返回 `not_implemented`）。它不是 RCTD 或 cell2location 之类基于模型的方法。 方法与 `deconvolute_spots` 相同。

- **API：** `dnakit.singlecell_spatial.deconvolute_spots_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 待解卷积的斑点：基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`signature`（基因×细胞类型矩阵，可用任一矩阵形式）或 `reference`（带注释细胞的计数矩阵）加 `reference_labels`（`{cell: 类型}` 或列表）、`normalize`（true）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`cell_types`、`n_shared_genes`、`fractions`（{斑点: {细胞类型: 比例}}）与 `rmse`（{斑点: 残差}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import deconvolute_spots_spatial

spots = """gene,s1,s2
A,100,10
B,10,100
C,50,50
"""
signature = """gene,typeX,typeY
A,10,1
B,1,10
C,5,5
"""

result = deconvolute_spots_spatial(spots, options={"signature": signature})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "deconvolute_spots_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 37,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d0743d02cb300ebb3b2966930e2cef0eac47e89410bcaa86785506703f7c27f3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "non-negative least squares of each mixture column on the signature columns (scipy.optimize.nnls); fr...",
    "cell_types": [
      "typeX",
      "typeY"
    ],
    "n_shared_genes": 3,
    "fractions": {
      "s1": {
        "typeX": 1.0,
        "typeY": 4.636616585542447e-17
      },
      "s2": {
        "typeX": 0.0,
        "typeY": 1.0
      }
    },
    "rmse": {
      "s1": 3.0913372807800604e-17,
      "s2": 2.1817011766850644e-17
    }
  }
}
```

## 18) demultiplex

- **作用：** 用简单的计数阈值规则，根据标签寡核苷酸（HTO）计数把细胞分配到样本。

- **处理方法：** 基于标签寡核苷酸（HTO）计数的阈值启发式规则，并非 HTODemux 的混合模型拟合。设最高的标签计数为 `top`、次高为 `second`：`top < min_counts` 时为 `Negative`；`second >= min_counts` 且 `top < ratio * second` 时为 `Doublet`；否则为 `Singlet`，赋予最高标签的名称。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.demultiplex(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 标签×细胞计数矩阵（行为标签，列为细胞），可用任一矩阵形式：映射 `{hashtag: [各细胞计数]}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为标签名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`min_counts`（5.0）、`ratio`（3.0）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`min_counts`、`ratio`、`assignments`（{cell: 标签名、`Doublet` 或 `Negative`}）与 `counts`（每种归属的细胞数）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import demultiplex

hashtags = """hashtag,c1,c2,c3,c4
HTO1,50,2,30,1
HTO2,3,40,28,0
"""

result = demultiplex(hashtags)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "demultiplex",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 50,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "209261238cd9bd6c030720432ba3df609b7a02316b34df54fa4d4bc42bd778f7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "Negative if top hashtag count < 5.0; Doublet if second count >= 5.0 and top < 3.0 x second; otherwis...",
    "min_counts": 5.0,
    "ratio": 3.0,
    "assignments": {
      "c1": "HTO1",
      "c2": "HTO2",
      "c3": "Doublet",
      "c4": "Negative"
    },
    "counts": {
      "Doublet": 1,
      "HTO1": 1,
      "HTO2": 1,
      "Negative": 1
    }
  }
}
```

## 19) demultiplex_spatial

- **作用：** 用简单的计数阈值规则，根据标签寡核苷酸计数把斑点分配到样本。

- **处理方法：** 基于标签寡核苷酸（HTO）计数的阈值启发式规则，并非 HTODemux 的混合模型拟合。设最高的标签计数为 `top`、次高为 `second`：`top < min_counts` 时为 `Negative`；`second >= min_counts` 且 `top < ratio * second` 时为 `Doublet`；否则为 `Singlet`，赋予最高标签的名称。纯 Python 实现。 规则与 `demultiplex` 相同，斑点相当于细胞。

- **API：** `dnakit.singlecell_spatial.demultiplex_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 标签×斑点计数矩阵（行为标签，列为斑点），可用任一矩阵形式：映射 `{hashtag: [各斑点计数]}`、行列表、CSV/TSV 文本或文件路径（表头为斑点名，首列为标签名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`min_counts`（5.0）、`ratio`（3.0）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`min_counts`、`ratio`、`assignments`（{cell: 标签名、`Doublet` 或 `Negative`}）与 `counts`（每种归属的细胞数）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import demultiplex_spatial

hashtags = """hashtag,s1,s2,s3,s4
HTO1,50,2,30,1
HTO2,3,40,28,0
"""

result = demultiplex_spatial(hashtags)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "demultiplex_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 50,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3613a5996b191d009f5788ae9d7e173227292c4dae93e4de467f7299851a7c9a"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "Negative if top hashtag count < 5.0; Doublet if second count >= 5.0 and top < 3.0 x second; otherwis...",
    "min_counts": 5.0,
    "ratio": 3.0,
    "assignments": {
      "s1": "HTO1",
      "s2": "HTO2",
      "s3": "Doublet",
      "s4": "Negative"
    },
    "counts": {
      "Doublet": 1,
      "HTO1": 1,
      "HTO2": 1,
      "Negative": 1
    }
  }
}
```

## 20) extract_barcodes

- **作用：** 从读段中提取细胞条形码（固定的读段窗口）并统计各条形码，可选与白名单比对。

- **处理方法：** 条形码 = `read[barcode_start : barcode_start + barcode_length]`（从 0 起算）。短于窗口的读段计入 `n_too_short` 并被跳过（`n_reads = n_extracted + n_unmatched + n_too_short`）；含 N 的条形码计入 `n_barcode_with_N`，且无法被校正。给定白名单时，条形码与白名单完全一致，或与唯一一个白名单条形码的汉明距离为 1（此时校正为该条形码）才被接受（白名单条目会转为大写并去掉 `-1` 之类的后缀）。每条列出的读段带有 `corrected` 字段（无法匹配时为 `null`），未匹配的读段不计入 `barcode_counts`，`whitelist` 给出 `size`、`exact`、`corrected` 与 `unmatched`。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.extract_barcodes(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** FASTQ/FASTA 文本或文件路径、序列列表，或 `{读段名: 序列}`（液滴方案的 read 1）。选项：`barcode_start`（0）、`barcode_length`（16）、`whitelist`（有效条形码的列表、文本或路径）、`top_n`（100）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_reads`（全部读段）、`n_extracted`（条形码被计入的读段：无白名单时为所有长度足够的读段，有白名单时仅为精确或校正匹配）、`n_unmatched`（条形码在白名单中找不到不超过 1 个错配的条目的读段；无白名单时恒为 0）、`n_too_short`、`n_barcode_with_N`、`n_unique_barcodes`、`barcode_counts`（{条形码: 读段数}，按频率降序，最多 `top_n` 项）与 `reads`（每条读段的条形码，最多列出 1000 条）及 `reads_truncated`；给定白名单时还有 `whitelist`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import extract_barcodes

reads = """@r1
ACGTTTGGAACC
+
IIIIIIIIIIII
@r2
ACGATTGGAACC
+
IIIIIIIIIIII
@r3
TTTTTTGGAACC
+
IIIIIIIIIIII
@r4
AC
+
II
@r5
GGGGTTGGAACC
+
IIIIIIIIIIII
"""

result = extract_barcodes(reads, options={"barcode_start": 0, "barcode_length": 4, "whitelist": ["ACGT", "TTTT"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_barcodes",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 140,
    "sequence_bases": 50,
    "lengths": [12, 12, 12, 2, 12],
    "sha256": "c0075b27af6dfe28d96714bc947d410ffc32d803fd54970c92ac95498e32e0a3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "barcode = read[0:4] (0-based, 4 bases) matched to the whitelist exactly or by a unique 1-mismatch correction",
    "n_reads": 5,
    "n_extracted": 3,
    "n_unmatched": 1,
    "n_too_short": 1,
    "n_barcode_with_N": 0,
    "n_unique_barcodes": 2,
    "barcode_counts": {
      "ACGT": 2,
      "TTTT": 1
    },
    "reads": [
      {
        "read": "r1",
        "barcode": "ACGT",
        "corrected": "ACGT"
      },
      {
        "read": "r2",
        "barcode": "ACGA",
        "corrected": "ACGT"
      },
      {
        "read": "r3",
        "barcode": "TTTT",
        "corrected": "TTTT"
      },
      {
        "read": "r5",
        "barcode": "GGGG",
        "corrected": null
      }
    ],
    "reads_truncated": false,
    "whitelist": {
      "size": 2,
      "exact": 2,
      "corrected": 1,
      "unmatched": 1
    }
  }
}
```

## 21) extract_barcodes_spatial

- **作用：** 从读段中提取空间（斑点）条形码并统计各条形码，可选与白名单比对。

- **处理方法：** 条形码 = `read[barcode_start : barcode_start + barcode_length]`（从 0 起算），与 `extract_barcodes` 相同。短于窗口的读段计入 `n_too_short` 并被跳过（`n_reads = n_extracted + n_unmatched + n_too_short`）。给定白名单时，条形码与白名单完全一致，或与唯一一个白名单条形码的汉明距离为 1（此时校正为该条形码）才被接受（白名单条目会转为大写并去掉 `-1` 之类的后缀）。每条列出的读段带有 `corrected` 字段（无法匹配时为 `null`），未匹配的读段不计入 `barcode_counts`，`whitelist` 给出 `size`、`exact`、`corrected` 与 `unmatched`。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.extract_barcodes_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** FASTQ/FASTA 文本或文件路径、序列列表，或 `{读段名: 序列}`（例如 Visium 的 read 1）。选项：`barcode_start`（0）、`barcode_length`（16）、`whitelist`（有效条形码的列表、文本或路径）、`top_n`（100）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_reads`（全部读段）、`n_extracted`（条形码被计入的读段：无白名单时为所有长度足够的读段，有白名单时仅为精确或校正匹配）、`n_unmatched`（条形码在白名单中找不到不超过 1 个错配的条目的读段；无白名单时恒为 0）、`n_too_short`、`n_barcode_with_N`、`n_unique_barcodes`、`barcode_counts`（{条形码: 读段数}，按频率降序，最多 `top_n` 项）与 `reads`（每条读段的条形码，最多列出 1000 条）及 `reads_truncated`；给定白名单时还有 `whitelist`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import extract_barcodes_spatial

reads = {"r1": "ACGTTTGGAACC", "r2": "ACGTTTGGAAGG", "r3": "TTTTCCGGAACC"}

result = extract_barcodes_spatial(reads, options={"barcode_length": 4})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_barcodes_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 66,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d8f2418172f6179af1f8f4228df8b4c23252a51e8cb59cf96a82640b296f3947"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "barcode = read[0:4] (0-based, 4 bases)",
    "n_reads": 3,
    "n_extracted": 3,
    "n_unmatched": 0,
    "n_too_short": 0,
    "n_barcode_with_N": 0,
    "n_unique_barcodes": 2,
    "barcode_counts": {
      "ACGT": 2,
      "TTTT": 1
    },
    "reads": [
      {
        "read": "r1",
        "barcode": "ACGT"
      },
      {
        "read": "r2",
        "barcode": "ACGT"
      },
      {
        "read": "r3",
        "barcode": "TTTT"
      }
    ],
    "reads_truncated": false
  }
}
```

## 22) extract_umis

- **作用：** 从读段中提取 UMI（固定的读段窗口）并统计各 UMI 的数量。

- **处理方法：** UMI = `read[umi_start : umi_start + umi_length]`（从 0 起算）；默认值对应 10x Chromium v3 的结构（16 bp 条形码后接 12 bp UMI）。短于窗口的读段计入 `n_too_short` 并被跳过；含 N 的 UMI 计入 `n_umi_with_N`。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.extract_umis(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** FASTQ/FASTA 文本或文件路径、序列列表，或 `{读段名: 序列}`（液滴方案的 read 1）。选项：`umi_start`（16）、`umi_length`（12）、`top_n`（100）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_reads`、`n_extracted`、`n_too_short`、`n_umi_with_N`、`n_unique_umis`、`umi_counts`（{UMI: 读段数}，按频率降序，最多 `top_n` 项）与 `reads`（每条读段的 UMI，最多列出 1000 条）及 `reads_truncated`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import extract_umis

reads = """@r1
ACGTTTGGAACC
+
IIIIIIIIIIII
@r2
ACGTTTGGAACC
+
IIIIIIIIIIII
@r3
ACGTCCCCAAAA
+
IIIIIIIIIIII
"""

result = extract_umis(reads, options={"umi_start": 4, "umi_length": 4})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_umis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 96,
    "sequence_bases": 36,
    "lengths": [12, 12, 12],
    "sha256": "872d4f5e9ac65abcc42c51edf93e26d961e0279f6e8b0fa5ebab99dd354d59ed"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "umi = read[4:8] (0-based, 4 bases)",
    "n_reads": 3,
    "n_extracted": 3,
    "n_too_short": 0,
    "n_umi_with_N": 0,
    "n_unique_umis": 2,
    "umi_counts": {
      "TTGG": 2,
      "CCCC": 1
    },
    "reads": [
      {
        "read": "r1",
        "umi": "TTGG"
      },
      {
        "read": "r2",
        "umi": "TTGG"
      },
      {
        "read": "r3",
        "umi": "CCCC"
      }
    ],
    "reads_truncated": false
  }
}
```

## 23) extract_umis_spatial

- **作用：** 从空间转录组读段中提取 UMI（固定的读段窗口）并统计各 UMI 的数量。

- **处理方法：** UMI = `read[umi_start : umi_start + umi_length]`（从 0 起算），与 `extract_umis` 相同。短于窗口的读段计入 `n_too_short` 并被跳过；含 N 的 UMI 计入 `n_umi_with_N`。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.extract_umis_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** FASTQ/FASTA 文本或文件路径、序列列表，或 `{读段名: 序列}`（例如 Visium 的 read 1）。选项：`umi_start`（16）、`umi_length`（12）、`top_n`（100）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_reads`、`n_extracted`、`n_too_short`、`n_umi_with_N`、`n_unique_umis`、`umi_counts`（{UMI: 读段数}，按频率降序，最多 `top_n` 项）与 `reads`（每条读段的 UMI，最多列出 1000 条）及 `reads_truncated`。以普通列表给出的读段依次命名为 `read_1`、`read_2`、……。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import extract_umis_spatial

reads = ["ACGTTTGGAACC", "ACGTTTGGAACC", "ACGTCCCCAAAA"]

result = extract_umis_spatial(reads, options={"umi_start": 4, "umi_length": 4})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "extract_umis_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 48,
    "sequence_bases": 36,
    "lengths": [
      12,
      12,
      12
    ],
    "sha256": "ea8969ac1fdce4d906a276683a791ef5ae593aed7e80670ff58a88fe43c9eaa4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "umi = read[4:8] (0-based, 4 bases)",
    "n_reads": 3,
    "n_extracted": 3,
    "n_too_short": 0,
    "n_umi_with_N": 0,
    "n_unique_umis": 2,
    "umi_counts": {
      "TTGG": 2,
      "CCCC": 1
    },
    "reads": [
      {
        "read": "read_1",
        "umi": "TTGG"
      },
      {
        "read": "read_2",
        "umi": "TTGG"
      },
      {
        "read": "read_3",
        "umi": "CCCC"
      }
    ],
    "reads_truncated": false
  }
}
```

## 24) find_differential_features

- **作用：** 检验两组细胞之间的差异表达（Welch t 检验或 Wilcoxon 秩和检验）。

- **处理方法：** 对每个基因：`norm = count / 文库大小 * target_sum`（`normalization: "none"` 时直接使用输入值），检验在 `log2(norm + pseudocount)` 上进行，`log2fc = log2(mean(norm, 第一组) + pseudocount) - log2(mean(norm, 第二组) + pseudocount)`，`pct_*` 为原始计数 > 0 的细胞比例，`p_value` 来自 Welch 不等方差 t 检验（默认）或 Mann-Whitney U / Wilcoxon 秩和检验（正态近似，含并列校正与连续性校正），`q_value` 为对所有检验基因做的 Benjamini-Hochberg 校正。基因按 p 值升序排列（并列时 |log2fc| 大者在前）。`group1` 为第一组，`group2` 为第二组。Welch 检验要求每组至少 2 个细胞（否则抛出 `InputFormatError`，错误码 `GROUP_TOO_SMALL`；可改用 `wilcoxon`）。这与 `calculate_markers` 不同，后者对过小的分组只给出警告并跳过。比较的是 `group1` 对 `group2`，而不是某一组对所有其他细胞；`mean_group1`/`pct_group1` 描述 `group1`，`log2fc` 为正表示在 `group1` 中更高。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.find_differential_features(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`groups`（必填；`{cell: 标签}` 或列表）、`group1` 与 `group2`（要比较的标签；恰有两个标签时可省略，此时按排序后的标签取）、`method`（默认 `"welch"`，或 `"wilcoxon"`）、`normalization`（默认 `"cpm"` = 将每个细胞缩放到 `target_sum`，或 `"none"`）、`target_sum`（10000）、`pseudocount`（1.0）、`top_n`（100）、`min_pct`（0.0；保留任一组中 pct 不低于该值的基因）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`group1`、`group2`、`n_group1`、`n_group2`、`method`、`normalization`、`target_sum`、`pseudocount`、`n_genes_tested` 与 `genes`（按排名取 top_n 个基因；每个基因行含 `gene`、`mean_group1`、`mean_group2`（归一化值的均值）、`log2fc`、`pct_group1`、`pct_group2`、`statistic`、`p_value` 与 `q_value`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import find_differential_features

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
groups = ["T", "T", "T", "B", "B", "B"]

result = find_differential_features(counts, options={"groups": groups, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "find_differential_features",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "B",
    "group2": "T",
    "n_group1": 3,
    "n_group2": 3,
    "method": "welch",
    "normalization": "cpm",
    "target_sum": 10000.0,
    "pseudocount": 1.0,
    "n_genes_tested": 8,
    "genes": [
      {
        "gene": "CD3D",
        "mean_group1": 61.7283950617284,
        "mean_group2": 1363.1378477316855,
        "log2fc": -4.442726976505825,
        "pct_group1": 0.3333333333333333,
        "pct_group2": 1.0,
        "statistic": -3.136919192646547,
        "p_value": 0.08794875669662708,
        "q_value": 0.18438518211850805
      },
      {
        "gene": "CD79A",
        "mean_group1": 1481.4814814814815,
        "mean_group2": 68.02721088435375,
        "log2fc": 4.424705081255279,
        "pct_group1": 1.0,
        "pct_group2": 0.3333333333333333,
        "statistic": 3.1091841711508663,
        "p_value": 0.0893622027560642,
        "q_value": 0.18438518211850805
      }
    ]
  }
}
```

## 25) find_differential_features_spatial

- **作用：** 检验两组斑点（或两个空间区域）之间的差异表达。

- **处理方法：** 对每个基因：`norm = count / 文库大小 * target_sum`（`normalization: "none"` 时直接使用输入值），检验在 `log2(norm + pseudocount)` 上进行，`log2fc = log2(mean(norm, 第一组) + pseudocount) - log2(mean(norm, 第二组) + pseudocount)`，`pct_*` 为原始计数 > 0 的细胞比例，`p_value` 来自 Welch 不等方差 t 检验（默认）或 Mann-Whitney U / Wilcoxon 秩和检验（正态近似，含并列校正与连续性校正），`q_value` 为对所有检验基因做的 Benjamini-Hochberg 校正。基因按 p 值升序排列（并列时 |log2fc| 大者在前）。`group1` 为第一组，`group2` 为第二组，斑点相当于细胞。Welch 检验要求每组至少 2 个斑点（否则抛出 `InputFormatError`，错误码 `GROUP_TOO_SMALL`；可改用 `wilcoxon`）。这与 `calculate_markers_spatial` 不同，后者对过小的区域只给出警告并跳过。比较的是 `group1` 对 `group2`，而不是某一组对所有其他斑点；`log2fc` 为正表示在 `group1` 中更高。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.find_differential_features_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`groups`（必填）、`group1` 与 `group2`（恰有两个标签时可省略）、`method`（默认 `"welch"`，或 `"wilcoxon"`）、`normalization`（默认 `"cpm"` = 将每个细胞缩放到 `target_sum`，或 `"none"`）、`target_sum`（10000）、`pseudocount`（1.0）、`top_n`（100）、`min_pct`（0.0；保留任一组中 pct 不低于该值的基因）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`group1`、`group2`、`n_group1`、`n_group2`、`method`、`normalization`、`target_sum`、`pseudocount`、`n_genes_tested` 与 `genes`（按排名取 top_n 个基因；每个基因行含 `gene`、`mean_group1`、`mean_group2`（归一化值的均值）、`log2fc`、`pct_group1`、`pct_group2`、`statistic`、`p_value` 与 `q_value`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import find_differential_features_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
domains = ["tissue1", "tissue1", "tissue1", "tissue2", "tissue2", "tissue2"]

result = find_differential_features_spatial(counts, options={"groups": domains, "top_n": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "find_differential_features_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "tissue1",
    "group2": "tissue2",
    "n_group1": 3,
    "n_group2": 3,
    "method": "welch",
    "normalization": "cpm",
    "target_sum": 10000.0,
    "pseudocount": 1.0,
    "n_genes_tested": 8,
    "genes": [
      {
        "gene": "CD3D",
        "mean_group1": 1363.1378477316855,
        "mean_group2": 61.7283950617284,
        "log2fc": 4.442726976505825,
        "pct_group1": 1.0,
        "pct_group2": 0.3333333333333333,
        "statistic": 3.136919192646547,
        "p_value": 0.08794875669662708,
        "q_value": 0.18438518211850805
      },
      {
        "gene": "CD79A",
        "mean_group1": 68.02721088435375,
        "mean_group2": 1481.4814814814815,
        "log2fc": -4.424705081255279,
        "pct_group1": 0.3333333333333333,
        "pct_group2": 1.0,
        "statistic": -3.1091841711508663,
        "p_value": 0.0893622027560642,
        "q_value": 0.18438518211850805
      }
    ]
  }
}
```

## 26) infer_cell_type_from_markers

- **作用：** 按标志基因得分排序，为每个细胞推断细胞类型（标志基因由用户提供）。

- **处理方法：** 取 ln(1 + CP10k) 值。对每个标志基因，在所有细胞间做 z 分数（总体标准差；常数基因得 0）。某细胞类型在某细胞中的得分为该类型标志基因 z 分数的均值，细胞被赋予得分最高的类型；最高分不大于 `min_score` 时为 `"Unassigned"`。矩阵中缺失的标志基因会给出警告（某类型的标志基因全部缺失则报错）。不使用训练模型。纯 Python 实现。 计算方式与 `annotate_cells` 相同。

- **API：** `dnakit.singlecell_spatial.infer_cell_type_from_markers(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`markers`（必填；`{细胞类型: [标志基因]}`）、`min_score`（0.0）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`assignments`（{cell: 类型}）、`best_score`（{cell: 得分}）、`counts`（每种类型的细胞数）与 `scores`（{类型: {cell: 得分}}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import infer_cell_type_from_markers

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
markers = {"T cell": ["CD3D", "CD3E"], "B cell": ["CD79A", "MS4A1"]}

result = infer_cell_type_from_markers(counts, options={"markers": markers})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_cell_type_from_markers",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "score(type, cell) = mean over the type's marker genes of the per-gene z-score of ln(1 + CP10k); the highest score wins if > min_score",
    "assignments": {
      "c1": "T cell",
      "c2": "T cell",
      "c3": "T cell",
      "c4": "B cell",
      "c5": "B cell",
      "c6": "B cell"
    },
    "best_score": {
      "c1": 0.8356997436691734,
      "c2": 0.8567110675970135,
      "c3": 0.8259936626183231,
      "c4": 0.8350643975343444,
      "c5": 0.8410178115070588,
      "c6": 0.8413837216201503
    },
    "counts": {
      "B cell": 3,
      "T cell": 3
    },
    "scores": {
      "T cell": {
        "c1": 0.8356997436691734,
        "c2": 0.8567110675970135,
        "c3": 0.8259936626183231,
        "c4": -0.5599901482596032,
        "c5": -0.5763799279190517,
        "c6": -1.3820343977058553
      },
      "B cell": {
        "c1": -0.5632664728908073,
        "c2": -0.5716923966489433,
        "c3": -1.3825070611218027,
        "c4": 0.8350643975343444,
        "c5": 0.8410178115070588,
        "c6": 0.8413837216201503
      }
    }
  }
}
```

## 27) infer_cnv

- **作用：** 相对参考细胞，沿染色体推断类似拷贝数的表达谱（inferCNV 风格，不含 HMM）。

- **处理方法：** inferCNV 风格的表达谱，不含 HMM：把计数转为 `log2(CPM + 1)`（`target_sum` 为 1000000），逐基因减去参考细胞的均值，在每条染色体内按位置排序后用 `window` 个基因的滑动平均平滑，最后对每个细胞做中位数居中。在 `gene_positions` 中没有记录、或平均原始计数低于 `min_mean` 的基因被剔除；染色体按自然顺序、基因按起点排序，滑动平均窗口在染色体两端会收缩。细胞的 `cnv_score` 为其修正后表达值的平均绝对值。不输出 CNV 区段，也不判定增益/缺失。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.infer_cnv(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`gene_positions`（必填；`{gene: [染色体, 起点]}`）、`groups`（必填；`{cell: 标签}` 或列表）、`reference_group`（必填；正常细胞的标签）、`window`（101 个基因）、`min_mean`（0.1）、`target_sum`（1000000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`window`、`reference_group`、`n_reference_cells`、`chromosome_of_gene`、`cnv_score`（{cell: 得分}）、`cnv_score_definition`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（平滑并居中后的表达谱）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import infer_cnv

counts = """gene,ref1,ref2,t1,t2
g1,10,12,10,11
g2,10,11,10,12
g3,10,10,20,22
g4,10,12,21,20
g5,10,11,10,12
g6,10,10,11,10
"""
options = {
    "gene_positions": {"g1": ["chr1", 100], "g2": ["chr1", 200], "g3": ["chr1", 300],
                       "g4": ["chr1", 400], "g5": ["chr2", 100], "g6": ["chr2", 200]},
    "groups": ["ref", "ref", "tumor", "tumor"],
    "reference_group": "ref",
    "window": 2,
}

result = infer_cnv(counts, options=options)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_cnv",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 7,
    "text_length": 111,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fde96bb232f2a9eec82e6b8896860c6131be1a2adbd848214e2b2f93e05d43ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "inferCNV-style expression profile: log2(CPM+1) minus the mean of the reference cells per gene, movin...",
    "window": 2,
    "reference_group": "ref",
    "n_reference_cells": 2,
    "chromosome_of_gene": {
      "g1": "chr1",
      "g2": "chr1",
      "g3": "chr1",
      "g4": "chr1",
      "g5": "chr2",
      "g6": "chr2"
    },
    "cnv_score": {
      "ref1": 0.016522693522969274,
      "ref2": 0.01652269352296878,
      "t1": 0.31512640468826447,
      "t2": 0.31465201165528606
    },
    "cnv_score_definition": "mean absolute value of the modified expression of the cell",
    "n_genes": 6,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2",
        "g3",
        "g4",
        "g5",
        "g6"
      ],
      "cells": [
        "ref1",
        "ref2",
        "t1",
        "t2"
      ],
      "values": [
        [
          -0.03387681047634672,
          0.03387681047634524,
          -0.2506093734737208,
          -0.17289386641218188
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.1160999446066177,
          0.1728938664121819
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.47289432167668544,
          0.4603909802784264
        ],
        [
          0.0004988540582449019,
          -0.0004988540582446062,
          0.8189548991593275,
          0.7299596507252003
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ]
      ]
    },
    "output_path": null
  }
}
```

## 28) infer_cnv_spatial

- **作用：** 相对参考斑点，沿染色体推断斑点的类似拷贝数的表达谱（inferCNV 风格）。

- **处理方法：** inferCNV 风格的表达谱，不含 HMM：把计数转为 `log2(CPM + 1)`（`target_sum` 为 1000000），逐基因减去参考细胞的均值，在每条染色体内按位置排序后用 `window` 个基因的滑动平均平滑，最后对每个细胞做中位数居中。在 `gene_positions` 中没有记录、或平均原始计数低于 `min_mean` 的基因被剔除；染色体按自然顺序、基因按起点排序，滑动平均窗口在染色体两端会收缩。细胞的 `cnv_score` 为其修正后表达值的平均绝对值。不输出 CNV 区段，也不判定增益/缺失。纯 Python 实现。 定义与 `infer_cnv` 相同，斑点相当于细胞。

- **API：** `dnakit.singlecell_spatial.infer_cnv_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`gene_positions`（必填；`{gene: [染色体, 起点]}`）、`groups`（必填；`{spot: 标签}` 或列表）、`reference_group`（必填；正常细胞的标签）、`window`（101 个基因）、`min_mean`（0.1）、`target_sum`（1000000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`window`、`reference_group`、`n_reference_cells`、`chromosome_of_gene`、`cnv_score`（{cell: 得分}）、`cnv_score_definition`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（平滑并居中后的表达谱）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import infer_cnv_spatial

counts = """gene,ref1,ref2,t1,t2
g1,10,12,10,11
g2,10,11,10,12
g3,10,10,20,22
g4,10,12,21,20
g5,10,11,10,12
g6,10,10,11,10
"""
options = {
    "gene_positions": {"g1": ["chr1", 100], "g2": ["chr1", 200], "g3": ["chr1", 300],
                       "g4": ["chr1", 400], "g5": ["chr2", 100], "g6": ["chr2", 200]},
    "groups": ["ref", "ref", "tumor", "tumor"],
    "reference_group": "ref",
    "window": 2,
}

result = infer_cnv_spatial(counts, options=options)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_cnv_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 7,
    "text_length": 111,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fde96bb232f2a9eec82e6b8896860c6131be1a2adbd848214e2b2f93e05d43ad"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "inferCNV-style expression profile: log2(CPM+1) minus the mean of the reference cells per gene, movin...",
    "window": 2,
    "reference_group": "ref",
    "n_reference_cells": 2,
    "chromosome_of_gene": {
      "g1": "chr1",
      "g2": "chr1",
      "g3": "chr1",
      "g4": "chr1",
      "g5": "chr2",
      "g6": "chr2"
    },
    "cnv_score": {
      "ref1": 0.016522693522969274,
      "ref2": 0.01652269352296878,
      "t1": 0.31512640468826447,
      "t2": 0.31465201165528606
    },
    "cnv_score_definition": "mean absolute value of the modified expression of the cell",
    "n_genes": 6,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2",
        "g3",
        "g4",
        "g5",
        "g6"
      ],
      "cells": [
        "ref1",
        "ref2",
        "t1",
        "t2"
      ],
      "values": [
        [
          -0.03387681047634672,
          0.03387681047634524,
          -0.2506093734737208,
          -0.17289386641218188
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.1160999446066177,
          0.1728938664121819
        ],
        [
          -0.0004988540582449019,
          0.0004988540582446057,
          0.47289432167668544,
          0.4603909802784264
        ],
        [
          0.0004988540582449019,
          -0.0004988540582446062,
          0.8189548991593275,
          0.7299596507252003
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ],
        [
          0.03188139424336711,
          -0.03188139424336681,
          -0.11609994460661768,
          -0.17588685305186283
        ]
      ]
    },
    "output_path": null
  }
}
```

## 29) infer_trajectory

- **作用：** 沿表达空间的 k 近邻图，计算从根细胞出发的测地线伪时间。

- **处理方法：** 细胞用中心化后的 ln(1 + CP10k) 向量表示（`n_pcs` > 0 时用主成分）。构建以欧氏距离为边权的对称 k 近邻图，`pseudotime` 为从根细胞出发的最短路径距离（Dijkstra）除以最大的有限距离。与根不连通的细胞为 `null` 并给出警告。这只是基于图距离的伪时间，不是带分支或谱系的轨迹推断方法。需要 numpy（缺失时返回 `not_implemented`）；最多 5000 个细胞。

- **API：** `dnakit.singlecell_spatial.infer_trajectory(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`root`（必填；细胞名）、`n_neighbors`（15）、`n_pcs`（0 = 使用中心化向量；> 0 = 先做 PCA）、`n_top_genes`、`scale`、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`root`、`n_neighbors` 与 `pseudotime`（{cell: 0 到 1 的值，与根不连通时为 `null`}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import infer_trajectory

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""

result = infer_trajectory(counts, options={"root": "c1", "n_neighbors": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_trajectory",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "shortest-path (geodesic) distance from the root cell in a symmetric kNN graph of the cells' expressi...",
    "root": "c1",
    "n_neighbors": 3,
    "pseudotime": {
      "c1": 0.0,
      "c2": 0.600610963142177,
      "c3": 0.5933128032311472,
      "c4": 0.8570303182661878,
      "c5": 0.931692072807665,
      "c6": 1.0
    }
  }
}
```

## 30) infer_trajectory_spatial

- **作用：** 计算斑点沿组织的伪时间：在空间邻接图上从根斑点出发的最短路径。

- **处理方法：** 在斑点坐标上构建对称 k 近邻图（`n_neighbors`，默认 6），边权为两个斑点中心化 ln(1 + CP10k) 表达向量（`n_pcs` > 0 时为主成分）之间的欧氏距离。`pseudotime` 为从根斑点出发的最短路径长度除以最大值。与根不连通的斑点为 `null` 并给出警告。需要 numpy（缺失时返回 `not_implemented`）；最多 5000 个斑点。

- **API：** `dnakit.singlecell_spatial.infer_trajectory_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`root`（必填；斑点名）、`coordinates`（必填；`{spot: [x, y]}`、按列顺序的 `[x, y]` 列表，或 CSV/TSV/Visium `tissue_positions` 文件）、`n_neighbors`（6；空间近邻数）、`n_pcs`（0）、`n_top_genes`、`scale`、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`root`、`n_neighbors` 与 `pseudotime`（{spot: 0 到 1 的值，与根不连通时为 `null`}）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import infer_trajectory_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}

result = infer_trajectory_spatial(counts, options={"root": "c1", "coordinates": coordinates, "n_neighbors": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_trajectory_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "shortest-path (geodesic) distance from the root cell in a symmetric kNN graph built on spatial coord...",
    "root": "c1",
    "n_neighbors": 2,
    "pseudotime": {
      "c1": 0.0,
      "c2": 0.28163935878865004,
      "c3": 0.5616005800468636,
      "c4": 0.40187989249504213,
      "c5": 0.7227424793087153,
      "c6": 1.0
    }
  }
}
```

## 31) infer_velocity

- **作用：** 用稳态模型，由剪接与未剪接计数估计每个基因、每个细胞的 RNA 速率。

- **处理方法：** La Manno 等（2018）的稳态模型：对每个基因，在该基因剪接（`s`）或未剪接（`u`）表达位于前 `quantile` 的细胞上，用过原点的最小二乘拟合降解比 `gamma = sum(u * s) / sum(s^2)`，`velocity = u - gamma * s`。应传入已归一化的数据层（或设 `normalize: true`，把每层缩放到每细胞 10000）。不把速率投影到嵌入图上。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.infer_velocity(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 剪接计数：基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`unspliced`（必填；基因和细胞与剪接矩阵相同且顺序一致的矩阵）、`quantile`（0.05）、`normalize`（false）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`quantile`、`gamma`（{gene: gamma}）、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（速率矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import infer_velocity

spliced = """gene,c1,c2,c3,c4
g1,1,2,3,4
g2,4,3,2,1
"""
unspliced = """gene,c1,c2,c3,c4
g1,2,3,5,6
g2,1,1,1,1
"""

result = infer_velocity(spliced, options={"unspliced": unspliced, "quantile": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_velocity",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 39,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fad0b78cf99c22d8ad13bf8a0c0cc01340461e0a4940257b6b5569cf483f5dba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "steady-state model (La Manno et al. 2018): gamma_g = sum(u*s)/sum(s^2) over the cells whose s or u i...",
    "quantile": 0.5,
    "gamma": {
      "g1": 1.56,
      "g2": 0.3333333333333333
    },
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2",
        "c3",
        "c4"
      ],
      "values": [
        [
          0.43999999999999995,
          -0.1200000000000001,
          0.3200000000000003,
          -0.2400000000000002
        ],
        [
          -0.33333333333333326,
          0.0,
          0.33333333333333337,
          0.6666666666666667
        ]
      ]
    },
    "output_path": null
  }
}
```

## 32) infer_velocity_spatial

- **作用：** 用稳态模型，由剪接与未剪接计数估计每个基因、每个斑点的 RNA 速率。

- **处理方法：** La Manno 等（2018）的稳态模型：对每个基因，在该基因剪接（`s`）或未剪接（`u`）表达位于前 `quantile` 的细胞上，用过原点的最小二乘拟合降解比 `gamma = sum(u * s) / sum(s^2)`，`velocity = u - gamma * s`。应传入已归一化的数据层（或设 `normalize: true`，把每层缩放到每细胞 10000）。不把速率投影到嵌入图上。纯 Python 实现。 模型与 `infer_velocity` 相同。

- **API：** `dnakit.singlecell_spatial.infer_velocity_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 剪接计数：基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`unspliced`（必填；基因和细胞与剪接矩阵相同且顺序一致的矩阵）、`quantile`（0.05）、`normalize`（false）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`quantile`、`gamma`（{gene: gamma}）、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（速率矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import infer_velocity_spatial

spliced = """gene,s1,s2,s3,s4
g1,1,2,3,4
g2,4,3,2,1
"""
unspliced = """gene,s1,s2,s3,s4
g1,2,3,5,6
g2,1,1,1,1
"""

result = infer_velocity_spatial(spliced, options={"unspliced": unspliced, "quantile": 0.5})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "infer_velocity_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 39,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "41d0822c1a43c1e406482fe3ce6fa42ad899944a4ee67f6be5a69ff110be1048"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "steady-state model (La Manno et al. 2018): gamma_g = sum(u*s)/sum(s^2) over the cells whose s or u i...",
    "quantile": 0.5,
    "gamma": {
      "g1": 1.56,
      "g2": 0.3333333333333333
    },
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "s1",
        "s2",
        "s3",
        "s4"
      ],
      "values": [
        [
          0.43999999999999995,
          -0.1200000000000001,
          0.3200000000000003,
          -0.2400000000000002
        ],
        [
          -0.33333333333333326,
          0.0,
          0.33333333333333337,
          0.6666666666666667
        ]
      ]
    },
    "output_path": null
  }
}
```

## 33) integrate

- **作用：** 对对数归一化表达量做逐基因的批次均值中心化，整合多个批次的细胞。

- **处理方法：** 逐基因的批次均值中心化：`x' = x - 该基因在所属批次的均值 + 该基因在所有细胞的均值`，作用于 ln(1 + CP10k) 值（`log_normalize` 为 true，默认）或输入值。它只对齐批次均值，不校正批次特异的方差或细胞类型构成，也不是 Harmony、scVI 或 MNN。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.integrate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`batches`（必填；`{cell: 批次标签}` 或按列顺序的列表）、`log_normalize`（true）、`target_sum`（10000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`batches`（{批次: 细胞数}）、`scale`、`note`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（校正后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import integrate

counts = """gene,a1,a2,b1,b2
g1,10,12,20,22
g2,5,6,5,7
"""

result = integrate(counts, options={"batches": ["A", "A", "B", "B"], "log_normalize": False})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "integrate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5f82895f5b4d8d8706f997e655451c6413349069dbf3487fe01fb0ba40faf310"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "per-gene batch mean-centring on the supplied values",
    "batches": {
      "A": 2,
      "B": 2
    },
    "scale": "as supplied",
    "note": "Aligns batch means only; it does not correct batch-specific variance or cell-type composition differences (not Harmony/scVI/MNN).",
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "a1",
        "a2",
        "b1",
        "b2"
      ],
      "values": [
        [
          15.0,
          17.0,
          15.0,
          17.0
        ],
        [
          5.25,
          6.25,
          4.75,
          6.75
        ]
      ]
    },
    "output_path": null
  }
}
```

## 34) integrate_single_cell_batches

- **作用：** 对对数归一化表达量做逐基因的批次均值中心化，整合多个单细胞批次。

- **处理方法：** 逐基因的批次均值中心化：`x' = x - 该基因在所属批次的均值 + 该基因在所有细胞的均值`，作用于 ln(1 + CP10k) 值（`log_normalize` 为 true，默认）或输入值。它只对齐批次均值，不校正批次特异的方差或细胞类型构成，也不是 Harmony、scVI 或 MNN。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.integrate_single_cell_batches(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`batches`（必填；`{cell: 批次标签}` 或按列顺序的列表）、`log_normalize`（true）、`target_sum`（10000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`batches`（{批次: 细胞数}）、`scale`、`note`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（校正后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import integrate_single_cell_batches

counts = """gene,a1,a2,b1,b2
g1,10,12,20,22
g2,5,6,5,7
"""

result = integrate_single_cell_batches(counts, options={"batches": ["A", "A", "B", "B"], "log_normalize": False})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "integrate_single_cell_batches",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5f82895f5b4d8d8706f997e655451c6413349069dbf3487fe01fb0ba40faf310"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "per-gene batch mean-centring on the supplied values",
    "batches": {
      "A": 2,
      "B": 2
    },
    "scale": "as supplied",
    "note": "Aligns batch means only; it does not correct batch-specific variance or cell-type composition differences (not Harmony/scVI/MNN).",
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "a1",
        "a2",
        "b1",
        "b2"
      ],
      "values": [
        [
          15.0,
          17.0,
          15.0,
          17.0
        ],
        [
          5.25,
          6.25,
          4.75,
          6.75
        ]
      ]
    },
    "output_path": null
  }
}
```

## 35) integrate_spatial

- **作用：** 对对数归一化表达量做逐基因的批次均值中心化，整合多张空间切片。

- **处理方法：** 逐基因的批次均值中心化：`x' = x - 该基因在所属批次的均值 + 该基因在所有细胞的均值`，作用于 ln(1 + CP10k) 值（`log_normalize` 为 true，默认）或输入值。它只对齐批次均值，不校正批次特异的方差或细胞类型构成，也不是 Harmony、scVI 或 MNN。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.integrate_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`batches`（必填；`{spot: 批次标签}` 或按列顺序的列表）、`log_normalize`（true）、`target_sum`（10000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`batches`（{批次: 细胞数}）、`scale`、`note`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（校正后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import integrate_spatial

counts = """gene,a1,a2,b1,b2
g1,10,12,20,22
g2,5,6,5,7
"""

result = integrate_spatial(counts, options={"batches": ["slice1", "slice1", "slice2", "slice2"], "log_normalize": False})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "integrate_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 43,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "5f82895f5b4d8d8706f997e655451c6413349069dbf3487fe01fb0ba40faf310"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "per-gene batch mean-centring on the supplied values",
    "batches": {
      "slice1": 2,
      "slice2": 2
    },
    "scale": "as supplied",
    "note": "Aligns batch means only; it does not correct batch-specific variance or cell-type composition differences (not Harmony/scVI/MNN).",
    "n_genes": 2,
    "n_cells": 4,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "a1",
        "a2",
        "b1",
        "b2"
      ],
      "values": [
        [
          15.0,
          17.0,
          15.0,
          17.0
        ],
        [
          5.25,
          6.25,
          4.75,
          6.75
        ]
      ]
    },
    "output_path": null
  }
}
```

## 36) map_cells

- **作用：** 通过 k 近邻标签转移，把查询细胞映射到带注释的参考上。

- **处理方法：** 查询数据与参考数据只保留共有基因（至少 2 个），并都转换为 ln(1 + CP10k) 尺度。每个查询细胞取欧氏距离最近的 `n_neighbors` 个参考细胞中的多数标签（票数相同时取最近的并列邻居的标签）；`confidence = 得票数 / k`。这是朴素的 k 近邻标签转移，不是 Seurat anchors 或 scANVI。需要 numpy（缺失时返回 `not_implemented`）；查询与参考各最多 5000 个细胞。

- **API：** `dnakit.singlecell_spatial.map_cells(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 查询：基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`reference`（必填；带注释的参考计数矩阵，可用任一矩阵形式）、`reference_labels`（必填；`{cell: 标签}` 或按参考矩阵列顺序的列表）、`n_neighbors`（5）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_shared_genes` 与 `assignments`；`assignments` 为 查询细胞 -> `label`、`confidence`、`nearest_reference_cell`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import map_cells

reference = """gene,r1,r2,r3,r4
CD3D,9,8,0,1
CD79A,0,1,9,8
ACTB,20,19,21,20
"""
query = """gene,q1,q2
CD3D,7,0
CD79A,1,10
ACTB,20,22
"""

result = map_cells(query, options={"reference": reference, "reference_labels": ["T", "T", "B", "B"], "n_neighbors": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "map_cells",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1c42bc93471c1856987fb453c0089d4e8835c4e8d79ab188f49b6b39dd729a51"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-nearest-neighbour label transfer (k=3): Euclidean distance on ln(1+CP10k) values of the 3 shared g...",
    "n_shared_genes": 3,
    "assignments": {
      "q1": {
        "label": "T",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r2"
      },
      "q2": {
        "label": "B",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r3"
      }
    }
  }
}
```

## 37) map_cells_spatial

- **作用：** 通过 k 近邻标签转移，把空间斑点映射到带注释的单细胞参考上。

- **处理方法：** 查询数据与参考数据只保留共有基因（至少 2 个），并都转换为 ln(1 + CP10k) 尺度。每个查询细胞取欧氏距离最近的 `n_neighbors` 个参考细胞中的多数标签（票数相同时取最近的并列邻居的标签）；`confidence = 得票数 / k`。这是朴素的 k 近邻标签转移，不是 Seurat anchors 或 scANVI。需要 numpy（缺失时返回 `not_implemented`）；查询与参考各最多 5000 个细胞。 查询中的斑点相当于查询细胞。

- **API：** `dnakit.singlecell_spatial.map_cells_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 查询斑点：基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`reference`（必填；带注释的参考计数矩阵，可用任一矩阵形式）、`reference_labels`（必填；`{cell: 标签}` 或按参考矩阵列顺序的列表）、`n_neighbors`（5）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_shared_genes` 与 `assignments`；`assignments` 为 查询细胞 -> `label`、`confidence`、`nearest_reference_cell`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import map_cells_spatial

reference = """gene,r1,r2,r3,r4
CD3D,9,8,0,1
CD79A,0,1,9,8
ACTB,20,19,21,20
"""
spots = """gene,s1,s2
CD3D,7,0
CD79A,1,10
ACTB,20,22
"""

result = map_cells_spatial(spots, options={"reference": reference, "reference_labels": ["T", "T", "B", "B"], "n_neighbors": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "map_cells_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "27a116a639b79546f99b5ef15311b99dbe8acc233d94e46701791788b26a3da5"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-nearest-neighbour label transfer (k=3): Euclidean distance on ln(1+CP10k) values of the 3 shared g...",
    "n_shared_genes": 3,
    "assignments": {
      "s1": {
        "label": "T",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r2"
      },
      "s2": {
        "label": "B",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r3"
      }
    }
  }
}
```

## 38) map_single_cell_to_reference

- **作用：** 通过 k 近邻标签转移，把单细胞映射到带注释的参考图谱上。

- **处理方法：** 查询数据与参考数据只保留共有基因（至少 2 个），并都转换为 ln(1 + CP10k) 尺度。每个查询细胞取欧氏距离最近的 `n_neighbors` 个参考细胞中的多数标签（票数相同时取最近的并列邻居的标签）；`confidence = 得票数 / k`。这是朴素的 k 近邻标签转移，不是 Seurat anchors 或 scANVI。需要 numpy（缺失时返回 `not_implemented`）；查询与参考各最多 5000 个细胞。

- **API：** `dnakit.singlecell_spatial.map_single_cell_to_reference(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 查询：基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`reference`（必填；带注释的参考计数矩阵，可用任一矩阵形式）、`reference_labels`（必填；`{cell: 标签}` 或按参考矩阵列顺序的列表）、`n_neighbors`（5）、`target_sum`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_shared_genes` 与 `assignments`；`assignments` 为 查询细胞 -> `label`、`confidence`、`nearest_reference_cell`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import map_single_cell_to_reference

reference = """gene,r1,r2,r3,r4
CD3D,9,8,0,1
CD79A,0,1,9,8
ACTB,20,19,21,20
"""
query = """gene,q1,q2
CD3D,7,0
CD79A,1,10
ACTB,20,22
"""

result = map_single_cell_to_reference(
    query, options={"reference": reference, "reference_labels": ["T", "T", "B", "B"], "n_neighbors": 3}
)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "map_single_cell_to_reference",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1c42bc93471c1856987fb453c0089d4e8835c4e8d79ab188f49b6b39dd729a51"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "k-nearest-neighbour label transfer (k=3): Euclidean distance on ln(1+CP10k) values of the 3 shared g...",
    "n_shared_genes": 3,
    "assignments": {
      "q1": {
        "label": "T",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r2"
      },
      "q2": {
        "label": "B",
        "confidence": 0.6666666666666666,
        "nearest_reference_cell": "r3"
      }
    }
  }
}
```

## 39) normalize

- **作用：** 按文库大小归一化计数矩阵，可选再取 ln(1 + x)。

- **处理方法：** 每个细胞（列）按文库大小缩放到 `target_sum`；`method: "log1p"`（默认）随后再取 ln(1 + x)，`"total"` 只缩放，`"cpm"` 则缩放到 1000000。总计数为零的细胞保持全零。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.normalize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`method`（默认 `"log1p"`，另有 `"total"`、`"cpm"`）、`target_sum`（10000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`method`、`target_sum`、`library_sizes`（{cell: 总计数}）、`log_base`（此处为 `null`）、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（归一化后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import normalize

counts = """gene,c1,c2
g1,1,10
g2,3,30
"""

result = normalize(counts, options={"method": "total", "target_sum": 100})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 27,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "94e87f829256b4e9de053af8b28b46acfc7b746db0aa5b9f2626d508d2123e26"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "total",
    "target_sum": 100.0,
    "library_sizes": {
      "c1": 4.0,
      "c2": 40.0
    },
    "log_base": null,
    "n_genes": 2,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2"
      ],
      "values": [
        [
          25.0,
          25.0
        ],
        [
          75.0,
          75.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 40) normalize_spatial

- **作用：** 按文库大小归一化斑点计数矩阵，可选再取 ln(1 + x)。

- **处理方法：** 每个细胞（列）按文库大小缩放到 `target_sum`；`method: "log1p"`（默认）随后再取 ln(1 + x)，`"total"` 只缩放，`"cpm"` 则缩放到 1000000。总计数为零的细胞保持全零。纯 Python 实现。 定义与 `normalize` 相同，斑点相当于细胞。

- **API：** `dnakit.singlecell_spatial.normalize_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`method`（默认 `"log1p"`，另有 `"total"`、`"cpm"`）、`target_sum`（10000）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`method`、`target_sum`、`library_sizes`（{cell: 总计数}）、`log_base`（此处为 `null`）、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（归一化后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import normalize_spatial

counts = """gene,s1,s2
g1,1,10
g2,3,30
"""

result = normalize_spatial(counts, options={"method": "total", "target_sum": 100})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "normalize_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 27,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ee2c9bba2f5b8bbaa511f940886f514dbe9fa98bbca3baef34e40ab8a12dce54"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "total",
    "target_sum": 100.0,
    "library_sizes": {
      "s1": 4.0,
      "s2": 40.0
    },
    "log_base": null,
    "n_genes": 2,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "s1",
        "s2"
      ],
      "values": [
        [
          25.0,
          25.0
        ],
        [
          75.0,
          75.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 41) plot_embedding

- **作用：** 把二维嵌入（UMAP、t-SNE 或 PCA 坐标）画成 SVG 散点图。

- **处理方法：** 每个细胞按给定坐标画一个点。给出 `labels` 时按类别着色并带图例（图例最多列出 20 个标签；`labels` 优先于 `values`）；给出 `values` 时使用连续色带；两者都没有时所有点同色。SVG 以文本方式直接生成，不依赖绘图库；给出 `output_path` 时也会写出（文件名须以 `.svg` 结尾；文件已存在则抛出 `ConfigurationError`）。

- **API：** `dnakit.singlecell_spatial.plot_embedding(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 每个细胞一行的坐标：`{cell: [x, y]}`、`[x, y]` 列表、CSV/TSV 文本或文件（表头为维度名，首列为细胞名），或 numpy 数组。选项：`labels`（`{cell: 标签}` 或列表；按类别着色并带图例）、`values`（`{cell: 数值}` 或列表；连续色带）、`output_path`（可选，要写出的 `.svg` 文件；不会覆盖已有文件）、`title`。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`format`（`"svg"`）、`svg`（SVG 文本）、`output_path`（写出的 `.svg` 文件，否则为 `null`）、`n_points`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import plot_embedding

embedding = {"c1": [0.0, 0.1], "c2": [0.2, 0.0], "c3": [1.0, 1.1], "c4": [1.2, 0.9]}

result = plot_embedding(embedding, options={"labels": ["T", "T", "B", "B"], "title": "UMAP"})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_embedding",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 72,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "29c10aeb08c17f75aa107412cf896fb1568539d195e831ede3be015149e5ab12"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "n_points": 4,
    "output_path": null
  }
}
```

## 42) plot_embedding_spatial

- **作用：** 把空间数据集的斑点坐标画成按标签或数值着色的 SVG 散点图。

- **处理方法：** 每个斑点按其组织坐标画一个点（y 轴向下，与图像一致）。给出 `labels` 时按类别着色并带图例（`labels` 优先于 `values`）；给出 `values` 时使用连续色带。SVG 以文本方式直接生成，不依赖绘图库；给出 `output_path` 时也会写出（文件名须以 `.svg` 结尾；文件已存在则抛出 `ConfigurationError`）。

- **API：** `dnakit.singlecell_spatial.plot_embedding_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 每个斑点一行的坐标：`{spot: [x, y]}`、`[x, y]` 列表、CSV/TSV 文本或文件，或 numpy 数组。选项：`labels`、`values`、`output_path`（可选，要写出的 `.svg` 文件；不会覆盖已有文件）、`title`。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`format`（`"svg"`）、`svg`（SVG 文本）、`output_path`（写出的 `.svg` 文件，否则为 `null`）、`n_points`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import plot_embedding_spatial

spots = {"s1": [0, 0], "s2": [1, 0], "s3": [0, 1], "s4": [1, 1]}

result = plot_embedding_spatial(spots, options={"values": [1.0, 2.0, 3.0, 4.0], "title": "Spots"})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_embedding_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "314985b9143128b56c67c5110e2e2e99c929a010e0138a149feb01de78c4bcd6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "n_points": 4,
    "output_path": null
  }
}
```

## 43) plot_heatmap

- **作用：** 把基因×细胞表达热图画成 SVG。

- **处理方法：** 绘制在所示细胞中方差最高的 `max_genes` 个基因和前 `max_cells` 个细胞；默认取 `log1p(x)`（`log1p`），可选对每个基因做 z 分数（`row_zscore`）。SVG 以文本方式直接生成，不依赖绘图库；给出 `output_path` 时也会写出（文件名须以 `.svg` 结尾；文件已存在则抛出 `ConfigurationError`）。

- **API：** `dnakit.singlecell_spatial.plot_heatmap(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`max_genes`（50）、`max_cells`（200）、`log1p`（true）、`row_zscore`（false）、`output_path`（可选，要写出的 `.svg` 文件；不会覆盖已有文件）、`title`、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`format`（`"svg"`）、`svg`（SVG 文本）、`output_path`（写出的 `.svg` 文件，否则为 `null`）、`genes_shown`、`n_cells_shown` 与 `selection`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import plot_heatmap

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""

result = plot_heatmap(counts, options={"max_genes": 4, "title": "Expression"})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_heatmap",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"230\" height=\"140\" viewBox=\"0 0 230 140\"><rect width=\"...",
    "genes_shown": [
      "CD79A",
      "CD3D",
      "MS4A1",
      "CD3E"
    ],
    "n_cells_shown": 6,
    "selection": "genes with the highest variance across the shown cells",
    "output_path": null
  }
}
```

## 44) plot_heatmap_spatial

- **作用：** 把单个基因的空间表达图（斑点按 ln(1 + CP10k) 着色）画成 SVG。

- **处理方法：** 每个斑点按其坐标绘制，并按该基因的 ln(1 + CP10k) 值（`target_sum` 为 10000）用连续色带着色。SVG 以文本方式直接生成，不依赖绘图库；给出 `output_path` 时也会写出（文件名须以 `.svg` 结尾；文件已存在则抛出 `ConfigurationError`）。

- **API：** `dnakit.singlecell_spatial.plot_heatmap_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`gene`（必填）、`coordinates`（必填；`{spot: [x, y]}`、按列顺序的 `[x, y]` 列表，或 CSV/TSV/Visium `tissue_positions` 文件）、`target_sum`（10000）、`output_path`（可选，要写出的 `.svg` 文件；不会覆盖已有文件）、`title`、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`format`（`"svg"`）、`svg`（SVG 文本）、`output_path`（写出的 `.svg` 文件，否则为 `null`）、`gene`、`n_spots`、`min_value` 与 `max_value`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import plot_heatmap_spatial

counts = """gene,c1,c2,c3,c4,c5,c6
CD3D,8,6,7,0,1,0
CD3E,5,7,6,1,0,0
CD79A,0,1,0,9,8,7
MS4A1,1,0,0,6,7,8
ACTB,20,18,22,19,21,20
GAPDH,15,14,16,15,13,17
MT-CO1,2,3,2,4,3,2
LYZ,0,0,1,0,1,0
"""
coordinates = {"c1": [0, 0], "c2": [1, 0], "c3": [2, 0], "c4": [0, 1], "c5": [1, 1], "c6": [2, 1]}

result = plot_heatmap_spatial(counts, options={"gene": "CD3D", "coordinates": coordinates})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "plot_heatmap_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 9,
    "text_length": 175,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "547af3b969e0a9463bbcf4f67f3895fd1eb999f74e7bcaf007bf628988f8507b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "gene": "CD3D",
    "n_spots": 6,
    "min_value": 0.0,
    "max_value": 7.358593577814888,
    "output_path": null
  }
}
```

## 45) qc

- **作用：** 计算计数矩阵的逐细胞、逐基因质控指标，并标出未通过阈值的细胞。

- **处理方法：** 每个细胞：总计数、检出基因数（计数 > 0）和线粒体计数百分比（名称以 `mito_prefix` 开头的基因的计数除以总计数）。未通过任一给定阈值的细胞连同未通过的检查项一并列出，检出细胞数少于 `min_cells` 的基因也会列出。未给定的阈值只报告、不过滤。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.qc(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`min_genes`、`min_counts`、`max_counts`、`max_pct_mito`、`min_cells`（均默认 `null`，即只报告）、`mito_prefix`（`"MT-"` 与 `"mt-"`）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`n_cells`、`n_genes`、`n_mito_genes`、`median_total_counts`、`median_genes_detected`、`median_pct_mito`、`genes_detected_in_any_cell`、`thresholds`、`n_cells_failing`、`failing_cells`、`n_genes_below_min_cells`、`genes_below_min_cells` 与 `per_cell`（`cell`、`total_counts`、`n_genes_detected`、`pct_mito`、`fails`）及 `per_cell_truncated`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import qc

matrix = {"MT-A": [2, 0], "B": [4, 5], "C": [4, 5]}

result = qc(matrix, options={"cells": ["x", "y"], "max_pct_mito": 10})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "edd1c156fa59a3ff619c27e7b9283d63d5a49d9dfb445ea4617ec6867f80b124"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_cells": 2,
    "n_genes": 3,
    "n_mito_genes": 1,
    "median_total_counts": 10.0,
    "median_genes_detected": 2.5,
    "median_pct_mito": 10.0,
    "genes_detected_in_any_cell": 3,
    "thresholds": {
      "min_genes": null,
      "min_counts": null,
      "max_counts": null,
      "max_pct_mito": 10.0,
      "min_cells": null
    },
    "n_cells_failing": 1,
    "failing_cells": [
      "x"
    ],
    "n_genes_below_min_cells": 0,
    "genes_below_min_cells": [],
    "per_cell": [
      {
        "cell": "x",
        "total_counts": 10.0,
        "n_genes_detected": 3,
        "pct_mito": 20.0,
        "fails": [
          "max_pct_mito"
        ]
      },
      {
        "cell": "y",
        "total_counts": 10.0,
        "n_genes_detected": 2,
        "pct_mito": 0.0,
        "fails": []
      }
    ],
    "per_cell_truncated": false
  }
}
```

## 46) qc_spatial

- **作用：** 计算空间计数矩阵的逐斑点、逐基因质控指标，并标出未通过阈值的斑点。

- **处理方法：** 每个细胞：总计数、检出基因数（计数 > 0）和线粒体计数百分比（名称以 `mito_prefix` 开头的基因的计数除以总计数）。未通过任一给定阈值的细胞连同未通过的检查项一并列出，检出细胞数少于 `min_cells` 的基因也会列出。未给定的阈值只报告、不过滤。纯 Python 实现。 可选的 `coordinates` 会额外给出有坐标的斑点数和坐标范围。

- **API：** `dnakit.singlecell_spatial.qc_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`min_genes`、`min_counts`、`max_counts`、`max_pct_mito`、`min_cells`（均默认 `null`，即只报告）、`mito_prefix`（`"MT-"` 与 `"mt-"`）、`coordinates`（可选；`{spot: [x, y]}`、按列顺序的 `[x, y]` 列表，或 CSV/TSV/Visium `tissue_positions` 文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`n_cells`、`n_genes`、`n_mito_genes`、`median_total_counts`、`median_genes_detected`、`median_pct_mito`、`genes_detected_in_any_cell`、`thresholds`、`n_cells_failing`、`failing_cells`、`n_genes_below_min_cells`、`genes_below_min_cells` 与 `per_cell`（`cell`、`total_counts`、`n_genes_detected`、`pct_mito`、`fails`）及 `per_cell_truncated`；给出 `coordinates` 时还有 `n_spots_with_coordinates`、`x_range`、`y_range`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import qc_spatial

matrix = {"MT-A": [2, 0], "B": [4, 5], "C": [4, 5]}
spots = {"x": [0, 0], "y": [3, 4]}

result = qc_spatial(matrix, options={"cells": ["x", "y"], "max_pct_mito": 10, "coordinates": spots})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "qc_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "edd1c156fa59a3ff619c27e7b9283d63d5a49d9dfb445ea4617ec6867f80b124"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_cells": 2,
    "n_genes": 3,
    "n_mito_genes": 1,
    "median_total_counts": 10.0,
    "median_genes_detected": 2.5,
    "median_pct_mito": 10.0,
    "genes_detected_in_any_cell": 3,
    "thresholds": {
      "min_genes": null,
      "min_counts": null,
      "max_counts": null,
      "max_pct_mito": 10.0,
      "min_cells": null
    },
    "n_cells_failing": 1,
    "failing_cells": [
      "x"
    ],
    "n_genes_below_min_cells": 0,
    "genes_below_min_cells": [],
    "per_cell": [
      {
        "cell": "x",
        "total_counts": 10.0,
        "n_genes_detected": 3,
        "pct_mito": 20.0,
        "fails": [
          "max_pct_mito"
        ]
      },
      {
        "cell": "y",
        "total_counts": 10.0,
        "n_genes_detected": 2,
        "pct_mito": 0.0,
        "fails": []
      }
    ],
    "per_cell_truncated": false,
    "n_spots_with_coordinates": 2,
    "x_range": [
      0.0,
      3.0
    ],
    "y_range": [
      0.0,
      4.0
    ]
  }
}
```

## 47) read

- **作用：** 读取计数矩阵（映射、行列表、CSV/TSV、数组或 10x 目录），并汇总其形状与稀疏度。

- **处理方法：** 共用的矩阵解析器接受所列的各种输入形式，并对每种形式得到同样的基因×细胞矩阵。`sparsity = 1 - 非零元素数 / (基因数 * 细胞数)`。矩阵元素个数不超过 `max_return_elements` 时返回矩阵数值，否则省略并给出警告。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`max_return_elements`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`n_genes`、`n_cells`、`total_counts`、`nonzero_entries`、`sparsity`、`genes_preview`、`cells_preview`（前 10 个名称）、`orientation`，以及（矩阵足够小时）`matrix`（`orientation`、`genes`、`cells`、`values`）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import read

matrix = {"g1": [1, 0, 3], "g2": [0, 0, 2]}

result = read(matrix, options={"cells": ["c1", "c2", "c3"]})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ac109c8082501eb269579a666d0f9049f1507060d522b056e4cd9c762cdf5f4e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 3,
    "total_counts": 6.0,
    "nonzero_entries": 3,
    "sparsity": 0.5,
    "genes_preview": [
      "g1",
      "g2"
    ],
    "cells_preview": [
      "c1",
      "c2",
      "c3"
    ],
    "orientation": "genes_by_cells",
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2",
        "c3"
      ],
      "values": [
        [
          1.0,
          0.0,
          3.0
        ],
        [
          0.0,
          0.0,
          2.0
        ]
      ]
    }
  }
}
```

## 48) read_spatial

- **作用：** 读取空间计数矩阵及斑点坐标，并汇总。

- **处理方法：** 计数矩阵的解析方式与 `read` 相同。坐标按名称与斑点对应（列表则按列顺序）；可识别无表头的 Visium v1 `tissue_positions_list` 文件（x = 像素列，y = 像素行）。没有坐标时只读取计数矩阵并给出警告。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.read_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`coordinates`（可选；`{spot: [x, y]}`、`[x, y]` 列表，或 CSV/TSV/Visium `tissue_positions` 文件）、`x_column` 与 `y_column`（表头名称）、`max_return_elements`（10000）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`read` 的各字段，以及（有坐标时）`coordinates`（{spot: [x, y]}）、`n_spots_with_coordinates`、`x_range`、`y_range`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import read_spatial

matrix = {"g1": [1, 0, 3], "g2": [0, 0, 2]}
spots = {"c1": [0, 0], "c2": [10, 0], "c3": [10, 5]}

result = read_spatial(matrix, options={"cells": ["c1", "c2", "c3"], "coordinates": spots})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "read_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ac109c8082501eb269579a666d0f9049f1507060d522b056e4cd9c762cdf5f4e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 3,
    "total_counts": 6.0,
    "nonzero_entries": 3,
    "sparsity": 0.5,
    "genes_preview": [
      "g1",
      "g2"
    ],
    "cells_preview": [
      "c1",
      "c2",
      "c3"
    ],
    "orientation": "genes_by_cells",
    "coordinates": {
      "c1": [
        0.0,
        0.0
      ],
      "c2": [
        10.0,
        0.0
      ],
      "c3": [
        10.0,
        5.0
      ]
    },
    "n_spots_with_coordinates": 3,
    "x_range": [
      0.0,
      10.0
    ],
    "y_range": [
      0.0,
      5.0
    ],
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2",
        "c3"
      ],
      "values": [
        [
          1.0,
          0.0,
          3.0
        ],
        [
          0.0,
          0.0,
          2.0
        ]
      ]
    }
  }
}
```

## 49) remove_ambient

- **作用：** 利用环境 RNA 谱和污染比例，从液滴计数矩阵中扣除环境 RNA 计数。

- **处理方法：** 简化的减法式 SoupX 风格校正，污染比例由用户给定：`校正值 = max(0, count - contamination * 细胞总计数 * 该基因的环境占比)`。环境谱可直接给出（`ambient_profile`，会归一化到总和为 1），或取 `empty_cells` 各列计数之和，或取总计数不超过 `empty_max_counts` 的各列之和；按此优先级使用，且使用后两个空液滴选项时这些空液滴会从输出中移除。不估计污染比例，也不把移除的计数重新分配为整数（`round` 仅对校正值取整）。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.remove_ambient(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`contamination`（必填；每个细胞计数中属于环境 RNA 的比例，取值于 [0, 1)）、`ambient_profile`（`{gene: 权重}`）或 `empty_cells`（空液滴列的列表）或 `empty_max_counts`（总计数不超过该值的列视为空液滴，并从输出中移除）、`round`（false；对校正值取整）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`contamination`、`ambient_source`、`ambient_profile`（{gene: 占比}）、`counts_removed`、`excluded_empty_droplets`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（校正后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import remove_ambient

counts = """gene,c1,c2
g1,10,0
g2,90,100
"""

result = remove_ambient(counts, options={"contamination": 0.1, "ambient_profile": {"g1": 1, "g2": 9}})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_ambient",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fb38f8f48b60feca759bbab925123eb87142fcf8fa7b1b5f22563bed973635e4"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "corrected = max(0, count - contamination * cell_total * ambient_fraction(gene)); simplified subtract...",
    "contamination": 0.1,
    "ambient_source": "ambient_profile",
    "ambient_profile": {
      "g1": 0.1,
      "g2": 0.9
    },
    "counts_removed": 19.0,
    "excluded_empty_droplets": [],
    "n_genes": 2,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "c1",
        "c2"
      ],
      "values": [
        [
          9.0,
          0.0
        ],
        [
          81.0,
          91.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 50) remove_ambient_spatial

- **作用：** 利用环境 RNA 谱和污染比例，从空间计数矩阵中扣除环境 RNA 计数。

- **处理方法：** 简化的减法式 SoupX 风格校正，污染比例由用户给定：`校正值 = max(0, count - contamination * 细胞总计数 * 该基因的环境占比)`。环境谱可直接给出（`ambient_profile`，会归一化到总和为 1），或取 `empty_cells` 各列计数之和，或取总计数不超过 `empty_max_counts` 的各列之和；按此优先级使用，且使用后两个空液滴选项时这些空液滴会从输出中移除。不估计污染比例，也不把移除的计数重新分配为整数（`round` 仅对校正值取整）。纯 Python 实现。 校正方式与 `remove_ambient` 相同，斑点相当于细胞。

- **API：** `dnakit.singlecell_spatial.remove_ambient_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`contamination`（必填；每个细胞计数中属于环境 RNA 的比例，取值于 [0, 1)）、`ambient_profile`（`{gene: 权重}`）或 `empty_cells`（空液滴列的列表）或 `empty_max_counts`（总计数不超过该值的列视为空液滴，并从输出中移除）、`round`（false；对校正值取整）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`contamination`、`ambient_source`、`ambient_profile`（{gene: 占比}）、`counts_removed`、`excluded_empty_droplets`、`n_genes`、`n_cells`、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（校正后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import remove_ambient_spatial

counts = """gene,s1,s2
g1,10,0
g2,90,100
"""

result = remove_ambient_spatial(counts, options={"contamination": 0.1, "ambient_profile": {"g1": 1, "g2": 9}})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_ambient_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8bd3aebc702036c3c3c2c8bce05e801f0123aa42609892abb4c79e24c6682b4d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "corrected = max(0, count - contamination * cell_total * ambient_fraction(gene)); simplified subtract...",
    "contamination": 0.1,
    "ambient_source": "ambient_profile",
    "ambient_profile": {
      "g1": 0.1,
      "g2": 0.9
    },
    "counts_removed": 19.0,
    "excluded_empty_droplets": [],
    "n_genes": 2,
    "n_cells": 2,
    "matrix": {
      "orientation": "genes_by_cells",
      "genes": [
        "g1",
        "g2"
      ],
      "cells": [
        "s1",
        "s2"
      ],
      "values": [
        [
          9.0,
          0.0
        ],
        [
          81.0,
          91.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 51) remove_doublets

- **作用：** 原本用于从计数矩阵中去除双细胞；目前未实现，始终返回 `not_implemented`。

- **处理方法：** 不进行任何计算。双细胞检测需要“模拟加分类器”的实现（例如 Scrublet 或 DoubletFinder）：模拟人工双细胞，再由训练好的分类器或近邻模型给每个细胞打分。DNAKit 不会调用任何此类软件包，所以安装它们不会改变结果，本库也没有内置确定性的参考定义，因此该函数仅用于保持名称和调用签名稳定。

- **API：** `dnakit.singlecell_spatial.remove_doublets(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 接受任意输入（只在 `input_summary` 中做摘要）；预期输入为计数矩阵。`options` 可选，目前未使用。

- **输出：** 始终返回封装：`status: "not_implemented"`，含 `module`、`operation`、`provenance`、`input_summary`、`warnings`、`required_dependency`（缺少内容的描述而非安装命令：一个基于模拟和分类器的双细胞检测实现，DNAKit 并不调用它）和 `message`（未实现的原因）。没有 `result`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import remove_doublets

counts = """gene,c1,c2
g1,10,0
g2,90,100
"""

result = remove_doublets(counts)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_doublets",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "fb38f8f48b60feca759bbab925123eb87142fcf8fa7b1b5f22563bed973635e4"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a doublet-detection implementation based on simulation and a classifier (for example Scrublet or DoubletFinder); DNAKit does not call any such package",
  "message": "Doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call; no deterministic reference definition is bundled."
}
```

## 52) remove_doublets_spatial

- **作用：** 原本用于去除空间斑点数据中的多重细胞斑点；目前未实现，始终返回 `not_implemented`。

- **处理方法：** 不进行任何计算。斑点级双细胞检测需要“模拟加分类器”的实现（例如 Scrublet 或 DoubletFinder），DNAKit 不会调用它们，所以安装也不会改变结果；本库也未内置此类实现，因此该函数仅用于保持名称和调用签名稳定。

- **API：** `dnakit.singlecell_spatial.remove_doublets_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 接受任意输入（只在 `input_summary` 中做摘要）；预期输入为斑点计数矩阵。`options` 可选，目前未使用。

- **输出：** 始终返回封装：`status: "not_implemented"`，含 `module`、`operation`、`provenance`、`input_summary`、`warnings`、`required_dependency`（缺少内容的描述而非安装命令：一个基于模拟和分类器的双细胞检测实现，DNAKit 并不调用它）和 `message`（未实现的原因）。没有 `result`。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import remove_doublets_spatial

counts = """gene,s1,s2
g1,10,0
g2,90,100
"""

result = remove_doublets_spatial(counts)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "remove_doublets_spatial",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 29,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "8bd3aebc702036c3c3c2c8bce05e801f0123aa42609892abb4c79e24c6682b4d"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a doublet-detection implementation based on simulation and a classifier (for example Scrublet or DoubletFinder); DNAKit does not call any such package",
  "message": "Spot-level doublet detection needs a simulation-and-classifier implementation (for example Scrublet or DoubletFinder) which DNAKit does not call; none is bundled."
}
```

## 53) validate

- **作用：** 校验计数矩阵：形状、数值、非负性、整数计数、名称唯一性以及空细胞/空基因。

- **处理方法：** 错误（使 `valid` 为 false）：矩阵无法解析、基因或细胞名称重复、存在负值，以及（`require_integer` 为 true，默认）存在非整数值。警告：总计数为零的细胞和在所有细胞中均为零的基因。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.validate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）。选项：`require_integer`（true；非整数值视为错误）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`n_genes`、`n_cells`、`empty_cells`、`empty_genes`、`valid`、`errors` 与 `warnings`（消息列表）；无法解析的矩阵不会抛出异常，而是返回 `valid: false` 并附上解析错误。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import validate

result = validate("gene,c1,c2\na,1,2\na,3,4\n")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 23,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f6eb24e753294ea677e9b11c564784e44b668c8f1c28050fd650d961bb38203e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 2,
    "empty_cells": 0,
    "empty_genes": 0,
    "valid": false,
    "errors": [
      "duplicate gene names: ['a']"
    ],
    "warnings": []
  }
}
```

## 54) validate_spatial

- **作用：** 校验空间计数矩阵及其斑点坐标（每个斑点都必须有有限的坐标）。

- **处理方法：** 包含 `validate` 的各项检查，并要求给出 `coordinates` 且每个斑点都有有限坐标，否则添加一条错误。纯 Python 实现。

- **API：** `dnakit.singlecell_spatial.validate_spatial(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×细胞计数矩阵（列为细胞）：映射 `{gene: [各细胞计数]}` 或 `{gene: {cell: count}}`、行列表、CSV/TSV 文本或文件路径（表头为细胞名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录（`matrix.mtx`、`features.tsv`、`barcodes.tsv`）；斑点（spot）按细胞处理，即列为斑点。选项：`coordinates`（必填；`{spot: [x, y]}`、按列顺序的 `[x, y]` 列表，或 CSV/TSV/Visium `tissue_positions` 文件）、`require_integer`（true）、`orientation`（默认 `"genes_by_cells"`；行为细胞时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`validate` 的各字段，以及（坐标有效时）`n_spots_with_coordinates`、`x_range`、`y_range`；坐标缺失会作为错误报告。

- **示例代码：**

```python
import json

from dnakit.singlecell_spatial import validate_spatial

matrix = {"g1": [1, 0, 3], "g2": [0, 0, 2]}
spots = {"c1": [0, 0], "c2": [1, 1]}

result = validate_spatial(matrix, options={"cells": ["c1", "c2", "c3"], "coordinates": spots})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "singlecell_spatial",
  "operation": "validate_spatial",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 34,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ac109c8082501eb269579a666d0f9049f1507060d522b056e4cd9c762cdf5f4e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 3,
    "empty_cells": 1,
    "empty_genes": 0,
    "valid": false,
    "errors": [
      "coordinates: [MISSING_COORDINATES] 1 cells have no coordinates (first: 'c3'). Hint: Provide coordinates for every cell/spot in the matrix."
    ],
    "warnings": [
      "1 cells have zero total counts (first: 'c2')"
    ]
  }
}
```
