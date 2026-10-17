# 转录组学

本类收录 **18** 项针对 bulk RNA-seq 的有文档的函数。其中 11 项对传入的数据计算出真实结果（10 项为纯 Python；`deconvolute` 在调用时导入 numpy 和 scipy，任一缺失则返回 `status: "not_implemented"`）。5 项封装外部程序（`align` 运行 STAR 或 HISAT2，`quantify` 调度后面三个函数，`run_featurecounts`、`run_kallisto_quantification`、`run_salmon_quantification` 分别运行 featureCounts、kallisto、salmon）：它们会校验输入，并在所需程序不在 PATH 中时返回 `status: "not_implemented"`。另有 2 项，`predict_protein_function` 与 `predict_protein_localization`，因未内置训练模型而始终返回 `status: "not_implemented"`。计数矩阵为基因×样本。输入格式错误或缺少必填选项时抛出 `InputFormatError` / `ConfigurationError`。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) align

- **作用：** 通过运行 STAR 或 HISAT2（外部程序）把 RNA-seq 读段比对到基因组。

- **处理方法：** 该函数本身不做比对。它先检查输入，再运行 PATH 中找到的 `STAR`（`--genomeDir <index> --readFilesIn ... --outSAMtype BAM SortedByCoordinate`，`.gz` 输入时用 `zcat`）或 `hisat2`（`-x <index>`，配合 `-U` 或 `-1/-2 ... -S`）（`aligner: "auto"` 优先 STAR）。PATH 中两者都没有时返回 `status: "not_implemented"`，而不是报错。本文档未实际运行 STAR/HISAT2（此处均未安装）。

- **API：** `dnakit.transcriptomics.align(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 一个 FASTQ 路径（单端）或两个 FASTQ 路径的列表（双端）；STAR 可用 `.gz`。选项：`index`（必填；STAR 基因组目录或 HISAT2 索引前缀）、`output_prefix`（必填；该文件不能已存在，其所在目录必须存在）、`aligner`（默认 `"auto"`，或 `"star"`、`"hisat2"`）、`threads`（1）、`timeout_seconds`（3600）。

- **输出：** PATH 中有 STAR 或 HISAT2 时：返回 `status: "ok"` 的封装，`result` 含 `program`、`aligner`、`command`、`alignment_file`（STAR：`<prefix>Aligned.sortedByCoord.out.bam`；HISAT2：`<prefix>.sam`）、`alignment_bytes` 与 `stderr_tail`。外部程序出错退出时抛出 `BackendExecutionError`，超时则抛出 `BackendTimeoutError`。路径必须存在且不得以 `-` 开头；输出路径不能已存在。两者都没有时：返回 `status: "not_implemented"`，含 `required_dependency`（`"STAR or hisat2 on PATH"`）和 `message`（PATH 中未找到该程序），如示例所示。

- **示例代码：**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import align

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
Path("genome_index").mkdir()
Path("out").mkdir()
result = align("reads.fastq", options={"index": "genome_index", "output_prefix": "out/sample_"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "align",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "STAR or hisat2 on PATH",
  "message": "transcriptomics.align runs the external program STAR / hisat2, which was not found on PATH."
}
```

## 2) annotate

- **作用：** 用 GTF 文件为基因标识符（Ensembl ID 或基因符号）注释坐标、链方向和生物类型。

- **处理方法：** 只扫描 GTF（纯文本或 `.gz`，最大 2 GiB）中的 `gene` 特征。标识符与 `gene_id`（含或不含 `.version` 后缀）或 `gene_name` 匹配；取第一个匹配的基因特征。坐标为 1 起始闭区间，与 GTF 一致；`biotype` 取自 `gene_biotype` 或 `gene_type`。无匹配的标识符列在 `not_found` 中，注释为 `null`。纯 Python 实现。

- **API：** `dnakit.transcriptomics.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因 ID 或符号的列表、每行一个的文本、文件路径，或以基因为键的映射/计数表。选项：`gtf`（必填；已存在的 `.gtf` 或 `.gtf.gz` 路径，须含 `gene` 特征）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`gtf`、`n_requested`、`n_found`、`not_found` 与 `annotations`（{标识符: {`gene_id`、`gene_name`、`chromosome`、`start`、`end`、`strand`、`biotype`、`source`} 或 `null`}）。

- **示例代码：**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import annotate

gtf = """chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG1.5"; gene_name "TP53"; gene_biotype "protein_coding";
chr2\tsrc\tgene\t1000\t2000\t.\t-\t.\tgene_id "ENSG2.1"; gene_name "MYC"; gene_biotype "lncRNA";
"""
Path("genes.gtf").write_text(gtf)

result = annotate(["TP53", "ENSG2", "UNKNOWN"], options={"gtf": "genes.gtf"})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "transcriptomics",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 28,
    "sequence_bases": 13,
    "lengths": [
      2,
      4,
      7
    ],
    "sha256": "d35e660823da88a58853e330de3187d8719d657015cb617d8028c633735ff4fa"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "lookup of each identifier (gene_id with or without version suffix, or gene_name) among the 'gene' fe...",
    "gtf": "genes.gtf",
    "n_requested": 3,
    "n_found": 2,
    "not_found": [
      "UNKNOWN"
    ],
    "annotations": {
      "TP53": {
        "gene_id": "ENSG1.5",
        "gene_name": "TP53",
        "chromosome": "chr1",
        "start": 100,
        "end": 500,
        "strand": "+",
        "biotype": "protein_coding",
        "source": "src"
      },
      "ENSG2": {
        "gene_id": "ENSG2.1",
        "gene_name": "MYC",
        "chromosome": "chr2",
        "start": 1000,
        "end": 2000,
        "strand": "-",
        "biotype": "lncRNA",
        "source": "src"
      },
      "UNKNOWN": null
    }
  }
}
```

## 3) calculate_differential_expression

- **作用：** 检验两组样本之间的基因差异表达（对 log2 CPM 做 Welch t 检验或 Wilcoxon 秩和检验）。

- **处理方法：** 对每个基因：`norm = count / 文库大小 * target_sum`（`normalization: "none"` 时直接使用输入值），检验在 `log2(norm + pseudocount)` 上进行，`log2fc = log2(mean(norm, 第一组) + pseudocount) - log2(mean(norm, 第二组) + pseudocount)`，`pct_*` 为原始计数 > 0 的细胞比例，`p_value` 来自 Welch 不等方差 t 检验（默认）或 Mann-Whitney U / Wilcoxon 秩和检验（正态近似，含并列校正与连续性校正），`q_value` 为对所有检验基因做的 Benjamini-Hochberg 校正。基因按 p 值升序排列（并列时 |log2fc| 大者在前）。`group1` 为第一组，`group2` 为第二组；`target_sum` 默认为 1000000，因此数值为 CPM。该检验只是逐基因的 t 检验或秩和检验，不是 DESeq2、edgeR 那样的负二项模型；Welch 检验要求每组至少 2 个样本（否则抛出 `InputFormatError`，错误码 `GROUP_TOO_SMALL`；可改用 `wilcoxon`）。与 `singlecell_spatial.calculate_markers` 对过小分组只给出警告并跳过不同，本函数会抛出异常。比较的是 `group1` 对 `group2`。纯 Python 实现。

- **API：** `dnakit.transcriptomics.calculate_differential_expression(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`groups`（必填；`{sample: 条件}` 或按列顺序的列表）、`group1` 与 `group2`（标签；恰有两个标签时可省略，此时按排序后的标签取）、`method`（默认 `"welch"`，或 `"wilcoxon"`）、`normalization`（默认 `"cpm"` = 将每个样本缩放到 `target_sum`，或 `"none"` = 直接检验输入值）、`target_sum`（1000000）、`pseudocount`（1.0）、`top_n`（100）、`min_pct`（0.0；任一组中计数 > 0 的样本比例的下限）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`group1`、`group2`、`n_group1`、`n_group2`、`method`、`normalization`、`target_sum`、`pseudocount`、`n_genes_tested`、`definition` 与 `genes`（按 p 值升序取 top_n 个基因；每个基因行含 `gene`、`mean_group1`、`mean_group2`（归一化值的均值）、`log2fc`、`pct_group1`、`pct_group2`、`statistic`、`p_value` 与 `q_value`，其中 `pct_*` 为计数 > 0 的样本比例）。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import calculate_differential_expression

counts = """gene,ctrl1,ctrl2,ctrl3,treat1,treat2,treat3
G1,100,120,110,400,380,420
G2,50,55,45,52,48,50
G3,200,210,190,60,70,65
G4,10,12,8,11,9,10
G5,500,520,480,510,490,505
"""
groups = ["ctrl", "ctrl", "ctrl", "treat", "treat", "treat"]

result = calculate_differential_expression(counts, options={"groups": groups, "top_n": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "transcriptomics",
  "operation": "calculate_differential_expression",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 162,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "6bd60df05ff99d244ba3316a53a33fb300eb356cb150f138b797452f18d5f8ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "ctrl",
    "group2": "treat",
    "n_group1": 3,
    "n_group2": 3,
    "method": "welch",
    "normalization": "cpm",
    "target_sum": 1000000.0,
    "pseudocount": 1.0,
    "n_genes_tested": 5,
    "definition": "value = count / library size * target_sum; test on log2(value + pseudocount); log2fc = log2(mean(gro...",
    "genes": [
      {
        "gene": "G5",
        "mean_group1": 574897.4540996943,
        "mean_group2": 488711.4839501525,
        "log2fc": 0.23432119327432943,
        "pct_group1": 1.0,
        "pct_group2": 1.0,
        "statistic": 14.920024280439547,
        "p_value": 0.000125519935865829,
        "q_value": 0.000627599679329145
      },
      {
        "gene": "G1",
        "mean_group1": 126397.79860106656,
        "mean_group2": 389455.03823506663,
        "log2fc": -1.6234777351714555,
        "pct_group1": 1.0,
        "pct_group2": 1.0,
        "statistic": -25.979707994055687,
        "p_value": 0.000443874559514615,
        "q_value": 0.0011096863987865376
      },
      {
        "gene": "G3",
        "mean_group1": 229885.66987242262,
        "mean_group2": 63399.54882086602,
        "log2fc": 1.8583555787894248,
        "pct_group1": 1.0,
        "pct_group2": 1.0,
        "statistic": 23.047641606359313,
        "p_value": 0.0016789159197535598,
        "q_value": 0.002798193199589266
      }
    ]
  }
}
```

## 4) calculate_differential_splicing

- **作用：** 根据剪接位点读段计数，检验两组之间外显子包含率（PSI，percent spliced in）的差异。

- **处理方法：** 对每个事件，逐样本计算 `PSI = (inc / inclusion_length) / (inc / inclusion_length + exc / exclusion_length)`（无读段的样本被忽略）；`psi_group` 为该组样本的均值，`delta_psi = psi_group1 - psi_group2`。`p_value` 为对两组合并的包含/排除读段数做的双侧 Fisher 精确检验；它忽略重复间的变异，因而偏不保守（anti-conservative）。`q_value` 为对被检验事件做的 Benjamini-Hochberg 校正。任一组合并读段数少于 `min_reads` 的事件会列出但不检验（`delta_psi`、`p_value`、`q_value` 为 `null`）。它不是 rMATS，也不是 beta-二项模型。纯 Python 实现。

- **API：** `dnakit.transcriptomics.calculate_differential_splicing(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 包含型剪接位点计数，事件×样本：每行一个剪接事件的包含型剪接位点计数矩阵，可用任一矩阵形式（映射 `{event: [各样本计数]}`、行列表、CSV/TSV 文本或文件路径、numpy 数组或 DataFrame，或 10x 目录）；计数必须为整数。选项：`exclusion`（必填；排除型剪接位点计数，事件和样本及其顺序须与输入一致）、`groups`（必填；`{sample: 条件}` 或列表）、`group1` 与 `group2`（恰有两个标签时可省略）、`inclusion_length`（1.0）与 `exclusion_length`（1.0；用于归一化读段数的有效长度）、`min_reads`（10；每组需要的合并读段数）、`top_n`（100）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`group1`、`group2`、`n_events`、`n_tested`、`definition`、`min_reads` 与 `events`（按 p 值升序取 top_n 个事件，未检验者排在最后），每项含 `event`、`psi_group1`、`psi_group2`、`inclusion_group1`、`exclusion_group1`、`inclusion_group2`、`exclusion_group2`（合并的读段数）、`delta_psi`、`p_value` 与 `q_value`。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import calculate_differential_splicing

inclusion = """event,a1,a2,b1,b2
exon1,80,90,20,25
exon2,50,55,52,48
"""
exclusion = """event,a1,a2,b1,b2
exon1,20,10,80,75
exon2,50,45,48,52
"""

result = calculate_differential_splicing(
    inclusion, options={"exclusion": exclusion, "groups": ["A", "A", "B", "B"]}
)
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "transcriptomics",
  "operation": "calculate_differential_splicing",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 54,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "611a697ca02148bba4442fbcd806ec930e8f3c16da9b6d0a0aa3311d22f8b3c3"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group1": "A",
    "group2": "B",
    "n_events": 2,
    "n_tested": 2,
    "definition": "PSI per sample = (inc/inclusion_length) / (inc/inclusion_length + exc/exclusion_length); psi_group =...",
    "min_reads": 10,
    "events": [
      {
        "event": "exon1",
        "psi_group1": 0.8500000000000001,
        "psi_group2": 0.225,
        "inclusion_group1": 170,
        "exclusion_group1": 30,
        "inclusion_group2": 45,
        "exclusion_group2": 155,
        "delta_psi": 0.6250000000000001,
        "p_value": 3.5477193042748607e-38,
        "q_value": 7.095438608549721e-38
      },
      {
        "event": "exon2",
        "psi_group1": 0.525,
        "psi_group2": 0.5,
        "inclusion_group1": 105,
        "exclusion_group1": 95,
        "inclusion_group2": 100,
        "exclusion_group2": 100,
        "delta_psi": 0.025000000000000022,
        "p_value": 0.6891172486827962,
        "q_value": 0.6891172486827962
      }
    ]
  }
}
```

## 5) calculate_translation_efficiency

- **作用：** 逐基因计算翻译效率（核糖体印记丰度除以 RNA-seq 丰度）。

- **处理方法：** 对每个共有基因和每个样本，`TE = 核糖体印记丰度 / RNA 丰度`，两种丰度均为各自文库的每百万计数（给出 `gene_lengths` 时为 RPKM；同一长度同时除以核糖体印记丰度和 RNA 丰度，在比值中相互抵消：在默认 `pseudocount` 为 0 时，`gene_lengths` 对 TE 值没有影响，只改变结果 `definition` 中的单位名称）。设置 `pseudocount` 时两种丰度在除以长度之后都加上它（因此伪计数大于 0 时，长度会轻微影响 TE），设置 `log2` 时报告 `log2(TE)`。RNA 丰度为零时 `TE` 为 `null`。`mean_te` 为 TE 有定义的样本上的均值。核糖体印记与 RNA-seq 样本按相同名称配对，否则在样本数相同时按位置配对。纯 Python 实现。

- **API：** `dnakit.transcriptomics.calculate_translation_efficiency(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 核糖体印记计数，基因×样本：核糖体印记计数矩阵，可用任一矩阵形式（映射 `{gene: [各样本计数]}`、行列表、CSV/TSV 文本或文件路径、numpy 数组或 DataFrame，或 10x 目录）。选项：`rna`（必填；同一批基因的 RNA-seq 计数矩阵；样本按相同名称配对，否则按位置配对并沿用核糖体印记的样本名）、`gene_lengths`（`{gene: bp}`；仅把丰度的单位标成 RPKM 而不是 CPM，不改变 TE，见处理方法）、`log2`（false）、`pseudocount`（0.0）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`n_genes`、`samples` 与 `genes`（{gene: {`te`: {样本: TE 或 `null`}, `mean_te`}}）。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import calculate_translation_efficiency

ribo = """gene,s1,s2
g1,100,200
g2,300,400
"""
rna = """gene,s1,s2
g1,50,100
g2,600,800
"""

result = calculate_translation_efficiency(ribo, options={"rna": rna})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "calculate_translation_efficiency",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 33,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "18571d7b56efc4b7e5fbfa750d75cb89f9eec543617624660f0c7d8783b9dffd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "TE = ribo abundance / RNA abundance, each as CPM (counts per million) of its own library; null when the RNA abundance is zero",
    "n_genes": 2,
    "samples": [
      "s1",
      "s2"
    ],
    "genes": {
      "g1": {
        "te": {
          "s1": 3.25,
          "s2": 3.0
        },
        "mean_te": 3.125
      },
      "g2": {
        "te": {
          "s1": 0.8125,
          "s2": 0.75
        },
        "mean_te": 0.78125
      }
    }
  }
}
```

## 6) compare

- **作用：** 在共有基因上用 Pearson 或 Spearman 相关系数比较表达谱。

- **处理方法：** 在两个矩阵共有的基因（至少 3 个）上，计算输入的每个样本与 `other` 的每个样本之间的相关系数（未给出 `other` 时计算输入内所有样本对）。`log2`（默认 true）时先取 `log2(value + 1)`；`normalization: "cpm"` 时先把每个样本缩放到 1000000。某个表达谱为常数时相关系数为 `null`。纯 Python 实现。

- **API：** `dnakit.transcriptomics.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`other`（可选的第二个矩阵；缺省时比较输入内所有样本对）、`method`（默认 `"pearson"`，或 `"spearman"`）、`log2`（true；对 `log2(value + 1)` 求相关）、`normalization`（默认 `"none"`，或 `"cpm"`）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`method`、`n_shared_genes`、`mode`（`"all pairs within the input"` 或 `"input versus other"`）与 `correlation`（{输入样本: {other 样本: r}}）。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import compare

counts = """gene,s1,s2,s3
g1,1,2,10
g2,2,4,5
g3,3,6,1
g4,4,8,0
"""

result = compare(counts, options={"log2": False})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 51,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0b35d186e0b0240699683fa5bfc3679b87c7da2fc9ae6d575703917756be5849"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "definition": "pearson correlation across 4 shared genes of the values",
    "method": "pearson",
    "n_shared_genes": 4,
    "mode": "all pairs within the input",
    "correlation": {
      "s1": {
        "s1": 1.0,
        "s2": 1.0,
        "s3": -0.9655351182200101
      },
      "s2": {
        "s1": 1.0,
        "s2": 1.0,
        "s3": -0.9655351182200101
      },
      "s3": {
        "s1": -0.9655351182200101,
        "s2": -0.9655351182200101,
        "s3": 1.0
      }
    }
  }
}
```

## 7) deconvolute

- **作用：** 在签名矩阵上用非负最小二乘估计 bulk RNA-seq 样本的细胞类型比例。

- **处理方法：** 在与签名矩阵共有的基因上，用非负最小二乘（`scipy.optimize.nnls`）把每个样本的表达谱拟合为签名列的非负组合；系数除以其和即得比例，`rmse` 为拟合残差。`normalize`（默认 true）时先把每一列缩放到总和为 1。若给出 `reference` 与 `reference_labels` 而非 `signature`，则签名为每种细胞类型的平均文库大小归一化表达谱。需要 numpy 和 scipy（缺失时返回 `not_implemented`）。它不是 CIBERSORT（不含 nu-SVR），也不是基于模型的方法。

- **API：** `dnakit.transcriptomics.deconvolute(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`signature`（基因×细胞类型矩阵，可用任一矩阵形式）或 `reference`（带注释细胞的计数矩阵）加 `reference_labels`（`{cell: 类型}` 或列表）、`normalize`（true）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`definition`、`cell_types`、`n_shared_genes`、`fractions`（{样本: {细胞类型: 比例}}）与 `rmse`（{样本: 残差}）。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import deconvolute

bulk = """gene,mix1,mix2
A,100,10
B,10,100
C,50,50
"""
signature = """gene,typeX,typeY
A,10,1
B,1,10
C,5,5
"""

result = deconvolute(bulk, options={"signature": signature})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "transcriptomics",
  "operation": "deconvolute",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 41,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "de454e1e3f3d2c243a5fd4e7192b152fcd71e2561fa478617bee4a934cc5c305"
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
      "mix1": {
        "typeX": 1.0,
        "typeY": 4.636616585542447e-17
      },
      "mix2": {
        "typeX": 0.0,
        "typeY": 1.0
      }
    },
    "rmse": {
      "mix1": 3.0913372807800604e-17,
      "mix2": 2.1817011766850644e-17
    }
  }
}
```

## 8) detect

- **作用：** 用读段计数（或 CPM）阈值检测 RNA-seq 样本中表达的基因。

- **处理方法：** 当 `count >= min_count`（若给出 `min_cpm`，则改为 `count / 文库大小 * 1e6 >= min_cpm`，并取代计数规则）时，基因在该样本中被检出。在至少 `min_samples` 个样本中被检出的基因视为整体表达。这只是简单的阈值规则，不是统计检验。纯 Python 实现。

- **API：** `dnakit.transcriptomics.detect(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`min_count`（1）、`min_cpm`（`null`；设置后改用 CPM >= `min_cpm`）、`min_samples`（1；需检出该基因的样本数）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`rule`（所用规则的文字描述）、`n_genes`、`n_detected`、`detected_genes`（最多列出 5000 个）及 `detected_truncated`、`per_sample_detected`（{样本: 检出基因数}）。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import detect

counts = """gene,s1,s2,s3
g1,0,0,0
g2,5,0,0
g3,5,7,9
"""

result = detect(counts, options={"min_count": 5, "min_samples": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "detect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 41,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e0dd490a7d446a9a0f50bf0dd5e73f4ae92c1b6c496f7e711eb6246ee72eb249"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "rule": "a gene is detected in a sample when count >= 5.0; it is expressed overall when detected in >= 2 samples",
    "n_genes": 3,
    "n_detected": 1,
    "detected_genes": [
      "g3"
    ],
    "detected_truncated": false,
    "per_sample_detected": {
      "s1": 2,
      "s2": 1,
      "s3": 1
    }
  }
}
```

## 9) normalize

- **作用：** 把 RNA-seq 计数矩阵归一化为 CPM、TPM、FPKM、median-of-ratios 或上四分位尺度。

- **处理方法：** `cpm`：`count / 文库大小 * 1e6`。`fpkm`/`rpkm`：`count / (基因长度(kb) * 文库大小(百万))`。`tpm`：`(count / 长度(kb)) / 所有基因 (count / 长度(kb)) 之和 * 1e6`。`median_of_ratios`：样本的 DESeq2 size factor 为，在所有样本计数均为正的基因上，`count / 样本间几何平均` 的中位数，数值为 `count / size factor`。`upper_quartile`：size factor 为该样本非零计数的 75 分位数除以各样本此分位数的均值，数值为 `count / size factor`。TPM 与 FPKM 需要为每个基因提供 `gene_lengths`。设置 `log2` 时结果为 `log2(value + pseudocount)`。cpm/fpkm/tpm 遇到总计数为零的样本会报错。纯 Python 实现。

- **API：** `dnakit.transcriptomics.normalize(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`method`（默认 `"cpm"`；另有 `"tpm"`、`"fpkm"`/`"rpkm"`、`"median_of_ratios"`、`"upper_quartile"`）、`gene_lengths`（`{gene: bp}`；tpm/fpkm 必填）、`log2`（false）、`pseudocount`（1.0）、`return_matrix`（true）、`output_path`（可选，要写出的 TSV 文件；不会覆盖已有文件）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`method`、`definition`、`size_factors`（{样本: 因子}；仅 `median_of_ratios` 与 `upper_quartile` 有）、`pseudocount`（仅 `log2` 时有）、`n_genes`、`n_cells`（即样本数）、`matrix`（`orientation`、`genes`、`cells`、`values`；`return_matrix` 为 false 时省略）与 `output_path`（写出的 TSV 文件，否则为 `null`）（归一化后的矩阵）。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import normalize

counts = """gene,s1,s2
g1,10,20
g2,90,180
"""

result = normalize(counts, options={"method": "cpm"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "normalize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 30,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7346eedbe68369b8e087ddc644b6a50133cb1b06e9fd39344d46c93c7998ca09"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "cpm",
    "definition": "CPM = count / library size * 1e6",
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
          100000.0,
          100000.0
        ],
        [
          900000.0,
          900000.0
        ]
      ]
    },
    "output_path": null
  }
}
```

## 10) plot

- **作用：** 把 RNA-seq 表达热图或样本对样本的散点图画成 SVG。

- **处理方法：** `kind: "heatmap"` 绘制样本间方差最高的 `max_genes` 个基因；`kind: "scatter"` 对样本 `x` 与 `y` 每个基因画一个点，并报告二者的 Pearson 相关系数 `pearson_r`。`log2`（默认 true）时绘制的数值为 `log2(value + 1)`。SVG 以文本方式直接生成，不依赖绘图库；给出 `output_path` 时也会写出（文件名须以 `.svg` 结尾；文件已存在则抛出 `ConfigurationError`）。

- **API：** `dnakit.transcriptomics.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`kind`（默认 `"heatmap"`，或 `"scatter"`）、`x` 与 `y`（样本名，`"scatter"` 必填；点为基因）、`log2`（true）、`max_genes`（50；热图展示方差最大的基因）、`output_path`（可选，要写出的 `.svg` 文件；不会覆盖已有文件）、`title`、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`format`（`"svg"`）、`svg`（SVG 文本）、`output_path`（写出的 `.svg` 文件，否则为 `null`）、`kind`；热图还有 `genes_shown`；散点图还有 `x`、`y` 与 `pearson_r`。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import plot

counts = """gene,s1,s2,s3
g1,1,2,10
g2,2,4,5
g3,3,6,1
"""

result = plot(counts, options={"kind": "scatter", "x": "s1", "y": "s2"})
print(json.dumps(result, indent=2))
```

- **示例输出：**（过长的列表和字符串用 `...` 缩略）

```json
{
  "module": "transcriptomics",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "eab73a065bf739aaf7971512f189b5aca923aad5062f3d120371b0fc77bb69ae"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "svg",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"480\" viewBox=\"0 0 720 480\"><rect width=\"...",
    "kind": "scatter",
    "output_path": null,
    "x": "s1",
    "y": "s2",
    "pearson_r": 0.9997907600601286
  }
}
```

## 11) predict_protein_function

- **作用：** 原本用于预测蛋白质功能；目前未实现，始终返回 `not_implemented`。

- **处理方法：** 不进行任何计算。蛋白质功能预测需要训练好的模型或功能注释数据库（例如 InterPro 或 DeepGOPlus 权重），而本库没有内置，因此该函数仅用于保持名称和调用签名稳定。

- **API：** `dnakit.transcriptomics.predict_protein_function(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 接受任意输入（只在 `input_summary` 中做摘要）；预期输入为蛋白质序列。`options` 可选，目前未使用。

- **输出：** 始终返回封装：`status: "not_implemented"`，含 `module`、`operation`、`provenance`、`input_summary`、`warnings`、`required_dependency` 与 `message`，没有 `result`。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import predict_protein_function

result = predict_protein_function("MKTAYIAKQR")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "predict_protein_function",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0b2b6b8411d220f8648b67b5f839a4336470f4fc93d9604009a1cb3a093ae29e"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained protein-function model or annotation database",
  "message": "Protein function prediction needs a trained model or functional annotation database (for example InterPro/DeepGOPlus weights) that is not bundled."
}
```

## 12) predict_protein_localization

- **作用：** 原本用于预测蛋白质亚细胞定位；目前未实现，始终返回 `not_implemented`。

- **处理方法：** 不进行任何计算。定位预测需要训练好的预测器（例如 DeepLoc 或 SignalP 权重），而本库没有内置，因此该函数仅用于保持名称和调用签名稳定。

- **API：** `dnakit.transcriptomics.predict_protein_localization(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 接受任意输入（只在 `input_summary` 中做摘要）；预期输入为蛋白质序列。`options` 可选，目前未使用。

- **输出：** 始终返回封装：`status: "not_implemented"`，含 `module`、`operation`、`provenance`、`input_summary`、`warnings`、`required_dependency` 与 `message`，没有 `result`。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import predict_protein_localization

result = predict_protein_localization("MKTAYIAKQR")
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "predict_protein_localization",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "0b2b6b8411d220f8648b67b5f839a4336470f4fc93d9604009a1cb3a093ae29e"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "trained subcellular-localisation model",
  "message": "Protein localisation prediction needs a trained predictor (for example DeepLoc or SignalP weights) that is not bundled."
}
```

## 13) qc

- **作用：** 计算每个样本的 RNA-seq 质控指标：文库大小、检出基因数、头部基因占比和线粒体百分比。

- **处理方法：** 每个样本：`library_size`（总计数）、`genes_detected`（计数 >= `min_count` 的基因数）、`fraction_zero_genes`、`top_gene_fraction`（最丰富基因的占比）、`top10_genes_fraction`（最丰富的 10 个基因的占比）和 `pct_mito`（名称以 `mito_prefix` 开头的基因的计数除以总计数再乘 100）。未通过任一给定阈值的样本连同未通过的检查项一并列出；未给定的阈值只报告、不过滤。纯 Python 实现。

- **API：** `dnakit.transcriptomics.qc(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`min_count`（1；判定“检出”的计数阈值）、`min_library_size`、`min_genes`、`max_pct_mito`（均默认 `null`，即只报告）、`mito_prefix`（`"MT-"` 与 `"mt-"`）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`n_samples`、`n_genes`、`min_count`、`n_mito_genes`、`median_library_size`、`samples`（{样本: 上述指标及 `fails`}）、`failing_samples` 与 `thresholds`。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import qc

counts = """gene,s1,s2
MT-ND1,5,0
A,45,50
B,50,50
"""

result = qc(counts, options={"max_pct_mito": 3})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 4,
    "text_length": 38,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "c3047d0b68648cb9963418c78c82ee1c83131f4a520b723a7153e31e7acf07bc"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_samples": 2,
    "n_genes": 3,
    "min_count": 1.0,
    "n_mito_genes": 1,
    "median_library_size": 100.0,
    "samples": {
      "s1": {
        "library_size": 100.0,
        "genes_detected": 3,
        "fraction_zero_genes": 0.0,
        "top_gene_fraction": 0.5,
        "top10_genes_fraction": 1.0,
        "pct_mito": 5.0,
        "fails": [
          "max_pct_mito"
        ]
      },
      "s2": {
        "library_size": 100.0,
        "genes_detected": 2,
        "fraction_zero_genes": 0.33333333333333337,
        "top_gene_fraction": 0.5,
        "top10_genes_fraction": 1.0,
        "pct_mito": 0.0,
        "fails": []
      }
    },
    "failing_samples": [
      "s1"
    ],
    "thresholds": {
      "min_library_size": null,
      "min_genes": null,
      "max_pct_mito": 3.0
    }
  }
}
```

## 14) quantify

- **作用：** 通过运行 salmon 或 kallisto（FASTQ 输入）或 featureCounts（BAM/SAM 输入）来定量转录本或基因表达。

- **处理方法：** 这是 `run_salmon_quantification`、`run_kallisto_quantification` 与 `run_featurecounts` 的调度函数，本身不做定量。以 `.bam` 或 `.sam` 结尾的输入交给 featureCounts，其他输入按 FASTQ 处理，交给 PATH 中找到的 salmon、kallisto 中的第一个（`method: "auto"`）。PATH 中没有合适的程序时返回 `status: "not_implemented"`。本文档未实际运行外部程序（此处均未安装）。

- **API：** `dnakit.transcriptomics.quantify(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** salmon/kallisto 用的 FASTQ 路径（单端一个、双端两个），或 featureCounts 用的 BAM/SAM 路径。选项：`method`（默认 `"auto"`；也可为 `"salmon"`、`"kallisto"`、`"featurecounts"`，须与输入类型匹配），以及所选执行函数（`run_salmon_quantification`、`run_kallisto_quantification` 或 `run_featurecounts`）的选项。

- **输出：** PATH 中有合适的程序时：`status: "ok"`，`result` 为所选执行函数的结果再加 `method`。没有时：`status: "not_implemented"`，含 `required_dependency`（例如 `"salmon or kallisto on PATH"`）和 `message`，如示例所示。`method` 与输入类型不匹配时抛出 `ConfigurationError`。

- **示例代码：**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import quantify

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
result = quantify("reads.fastq", options={"index": "transcripts_index", "output_dir": "quant_out"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "quantify",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "salmon or kallisto on PATH",
  "message": "transcriptomics.quantify runs the external program salmon / kallisto, which was not found on PATH."
}
```

## 15) read

- **作用：** 读取表达计数矩阵（映射、行列表、CSV/TSV、数组或 10x 目录），并汇总基因数、样本数与稀疏度。

- **处理方法：** 共用的矩阵解析器接受所列的各种输入形式，并对每种形式得到同样的基因×样本矩阵。`sparsity = 1 - 非零元素数 / (基因数 * 样本数)`。矩阵元素个数不超过 `max_return_elements` 时返回矩阵数值，否则省略并给出警告。纯 Python 实现。

- **API：** `dnakit.transcriptomics.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 基因×样本计数矩阵（列为样本）：映射 `{gene: [各样本计数]}` 或 `{gene: {sample: count}}`、行列表、CSV/TSV 文本或文件路径（表头为样本名，首列为基因名）、numpy 数组或 pandas DataFrame，或 10x 目录。选项：`max_return_elements`（10000）、`orientation`（默认 `"genes_by_cells"`，行为基因；行为样本时用 `"cells_by_genes"`）、`genes` 与 `cells`（无名称输入形式下的名称列表；`cells` 即样本名）、`max_matrix_elements`（5000000）。

- **输出：** 返回的封装包含 `status`（`"ok"`）、`module`、`operation`、`provenance`、`input_summary`、`warnings` 与 `result`。输入格式错误或缺少必填选项时不返回封装，而是抛出 `InputFormatError` / `ConfigurationError`。`result`：`n_genes`、`n_cells`（即样本数）、`total_counts`、`nonzero_entries`、`sparsity`、`genes_preview`、`cells_preview`（前 10 个名称）、`orientation`，以及（矩阵足够小时）`matrix`（`orientation`、`genes`、`cells`、`values`）。

- **示例代码：**

```python
import json

from dnakit.transcriptomics import read

counts = """gene,s1,s2
g1,5,0
g2,3,8
"""

result = read(counts)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 25,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e9a9a4e08438f684e4fe883bbb5c626311d31f99954fd46f604ece73b69a4b4d"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_genes": 2,
    "n_cells": 2,
    "total_counts": 16.0,
    "nonzero_entries": 3,
    "sparsity": 0.25,
    "genes_preview": [
      "g1",
      "g2"
    ],
    "cells_preview": [
      "s1",
      "s2"
    ],
    "orientation": "genes_by_cells",
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
          5.0,
          0.0
        ],
        [
          3.0,
          8.0
        ]
      ]
    }
  }
}
```

## 16) run_featurecounts

- **作用：** 通过运行 Subread 的 featureCounts（外部程序），对比对文件按基因计数。

- **处理方法：** 该函数本身不做计数。它先检查输入，再对比对文件运行 `featureCounts -a <annotation> -o <output> -T <threads> -s <strand> -t <feature_type> -g <attribute>`（`paired` 时再加 `-p --countReadPairs`），并汇总输出表。PATH 中没有 featureCounts 时返回 `status: "not_implemented"`。本文档未实际运行 featureCounts（此处未安装）。

- **API：** `dnakit.transcriptomics.run_featurecounts(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 一个 BAM/SAM 路径或路径列表。选项：`annotation`（必填；已存在的 GTF/GFF 路径）、`output`（必填；要创建的计数文件，不能已存在）、`strand`（0 无链特异、1 正链特异、2 反链特异；默认 0）、`feature_type`（`"exon"`）、`attribute`（`"gene_id"`）、`paired`（false；按读段对计数）、`threads`（1）、`timeout_seconds`（3600）。

- **输出：** PATH 中有 featureCounts 时：`status: "ok"`，`result` 含 `program`、`command`、`counts_file`、`summary_file`、`n_genes`、`samples`、`assigned_counts`（{样本: 分配到基因的计数总和}）与 `stderr_tail`。外部程序出错退出时抛出 `BackendExecutionError`，超时则抛出 `BackendTimeoutError`。路径必须存在且不得以 `-` 开头；输出路径不能已存在。没有 featureCounts 时：`status: "not_implemented"`，含 `required_dependency`（`"featureCounts (subread) on PATH"`）和 `message`，如示例所示。

- **示例代码：**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import run_featurecounts

Path("sample.bam").write_bytes(b"BAM\x01")
Path("genes.gtf").write_text("")
result = run_featurecounts("sample.bam", options={"annotation": "genes.gtf", "output": "counts.txt"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "run_featurecounts",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 10,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "46082487f59dd44b1aabcca96ce194fd906f7cd4b90ecb2c50b79aa9e1cd0637"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "featureCounts (subread) on PATH",
  "message": "transcriptomics.run_featurecounts runs the external program featureCounts (subread), which was not found on PATH."
}
```

## 17) run_kallisto_quantification

- **作用：** 通过对 FASTQ 读段运行 `kallisto quant`（外部程序）来定量转录本丰度。

- **处理方法：** 该函数本身不做定量。它先检查输入，再运行 `kallisto quant -i <index> -o <output_dir> -t <threads>`（单端读段再加 `--single -l <fragment_length> -s <fragment_sd>`），并汇总 `abundance.tsv`。PATH 中没有 kallisto 时返回 `status: "not_implemented"`。本文档未实际运行 kallisto（此处未安装）。

- **API：** `dnakit.transcriptomics.run_kallisto_quantification(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 一个 FASTQ 路径（单端）或两个 FASTQ 路径的列表（双端）。选项：`index`（必填；已存在的 kallisto 索引文件）、`output_dir`（必填；不能已存在或必须为空，其父目录必须存在）、`fragment_length` 与 `fragment_sd`（单端读段必填）、`threads`（1）、`timeout_seconds`（3600）。

- **输出：** PATH 中有 kallisto 时：`status: "ok"`，`result` 含 `program`、`command`、`abundance_file`、`n_transcripts`、`total_est_counts`、`total_tpm`、`top_transcripts`（TPM 最高的 20 个转录本）与 `stderr_tail`。外部程序出错退出时抛出 `BackendExecutionError`，超时则抛出 `BackendTimeoutError`。路径必须存在且不得以 `-` 开头；输出路径不能已存在。没有 kallisto 时：`status: "not_implemented"`，含 `required_dependency`（`"kallisto on PATH"`）和 `message`，如示例所示。

- **示例代码：**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import run_kallisto_quantification

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
Path("transcripts.idx").write_bytes(b"")
result = run_kallisto_quantification("reads.fastq", options={"index": "transcripts.idx", "output_dir": "kallisto_out", "fragment_length": 200, "fragment_sd": 20})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "run_kallisto_quantification",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "kallisto on PATH",
  "message": "transcriptomics.run_kallisto_quantification runs the external program kallisto, which was not found on PATH."
}
```

## 18) run_salmon_quantification

- **作用：** 通过对 FASTQ 读段运行 `salmon quant`（外部程序）来定量转录本丰度。

- **处理方法：** 该函数本身不做定量。它先检查输入，再运行 `salmon quant -i <index> -l <library_type> -o <output_dir> -p <threads> --validateMappings`（单端用 `-r`，双端用 `-1/-2`），并汇总 `quant.sf`。PATH 中没有 salmon 时返回 `status: "not_implemented"`。本文档未实际运行 salmon（此处未安装）。

- **API：** `dnakit.transcriptomics.run_salmon_quantification(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** 一个 FASTQ 路径（单端）或两个 FASTQ 路径的列表（双端）。选项：`index`（必填；已存在的 salmon 索引目录）、`output_dir`（必填；不能已存在或必须为空，其父目录必须存在）、`library_type`（`"A"` = 自动）、`threads`（1）、`timeout_seconds`（3600）。

- **输出：** PATH 中有 salmon 时：`status: "ok"`，`result` 含 `program`、`command`、`quant_file`、`n_transcripts`、`total_reads_assigned`、`total_tpm`、`top_transcripts`（TPM 最高的 20 个转录本）与 `stderr_tail`。外部程序出错退出时抛出 `BackendExecutionError`，超时则抛出 `BackendTimeoutError`。路径必须存在且不得以 `-` 开头；输出路径不能已存在。没有 salmon 时：`status: "not_implemented"`，含 `required_dependency`（`"salmon on PATH"`）和 `message`，如示例所示。

- **示例代码：**

```python
import json
from pathlib import Path

from dnakit.transcriptomics import run_salmon_quantification

Path("reads.fastq").write_text("@r1\nACGTACGT\n+\nIIIIIIII\n")
Path("salmon_index").mkdir()
result = run_salmon_quantification("reads.fastq", options={"index": "salmon_index", "output_dir": "salmon_out"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "transcriptomics",
  "operation": "run_salmon_quantification",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 11,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d9d1590a8b2978eb1aa73ef39296ead72a0a304c3b4f3d46e696cca37a3baae9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "salmon on PATH",
  "message": "transcriptomics.run_salmon_quantification runs the external program salmon, which was not found on PATH."
}
```
