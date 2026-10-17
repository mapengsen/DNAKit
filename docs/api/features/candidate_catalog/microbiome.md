# 微生物组

本类收录 **21** 项有文档的候选功能，其中 12 项在本地原生计算出真实结果（丰度表处理、多样性、k-mer ANI 估计、共现网络、差异丰度、去冗余、绘图与导出），9 项始终返回 `status: "not_implemented"`（`assemble`、`bin`、`build_microbial_gene_catalog`、`build_pathogen_transmission_tree`、`calculate_digital_ddh`、`classify`、`infer_viral_quasispecies`、`predict`、`track`），因为它们需要外部程序、参考数据库或已训练的模型。`calculate_differential_abundance` 在未安装 scipy 时也会返回 `not_implemented`。有两处结果明确是启发式：`calculate_ani`（以及 `ani` 模式的 `dereplicate`）是基于精确 k-mer 集合的 Mash 公式估计，不是基于比对的 ANI；共现网络只使用相对丰度的普通相关系数。这些函数都通过正式的 DNAKit 命名空间导出；JSON 调用方案兼容的函数也会注册到 Agent/MCP 工具库。

## 1) annotate

- **作用：** 把分类学谱系字符串解析成分级字段（域到种），并给出最低的已分配等级。

- **处理方法：** 纯文本解析：谱系按 `;` 或 `|` 切分。若每一级都带等级前缀（`k__`/`d__`、`p__`、`c__`、`o__`、`f__`、`g__`、`s__`；`k__` 与 `d__` 都对应 `domain`），由前缀决定等级；否则按位置依次读作域、门、纲、目、科、属、种（无前缀的级别超过 7 个会报错）。空名称或占位名称（`__`、`Unassigned`、`unclassified`、`unknown`、`uncultured`、`ambiguous`、`incertae sedis` 等）变为 `null`。它不会对序列做分类，只是解释已经分配好的谱系文本。

- **API：** `dnakit.microbiome.annotate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：单个谱系字符串、谱系字符串列表（编号为 `lineage1` 等）或映射 `{feature_id: 谱系}`，至多 1,000,000 条。`options`：`entity`（只支持 `"taxonomy"`，即默认值）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 是列表，每条谱系一项：`id`、`lineage`、`ranks`（`domain`、`phylum`、`class`、`order`、`family`、`genus`、`species` 到名称或 `null` 的映射）、`lowest_rank` 和 `lowest_taxon`。

- **示例代码：**

```python
import json
from dnakit.microbiome import annotate

lineages = {
    "otu1": "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__",
    "otu2": "Archaea;Euryarchaeota;Methanobacteria",
}
result = annotate(lineages)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "annotate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 157,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "1c91805a12665e845e3dc84ecad93538c3da24850a4389154d1366268c535f9e"
  },
  "warnings": [],
  "status": "ok",
  "result": [
    {
      "id": "otu1",
      "lineage": "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__",
      "ranks": {
        "domain": "Bacteria",
        "phylum": "Firmicutes",
        "class": "Bacilli",
        "order": "Lactobacillales",
        "family": "Lactobacillaceae",
        "genus": "Lactobacillus",
        "species": null
      },
      "lowest_rank": "genus",
      "lowest_taxon": "Lactobacillus"
    },
    {
      "id": "otu2",
      "lineage": "Archaea;Euryarchaeota;Methanobacteria",
      "ranks": {
        "domain": "Archaea",
        "phylum": "Euryarchaeota",
        "class": "Methanobacteria",
        "order": null,
        "family": null,
        "genus": null,
        "species": null
      },
      "lowest_rank": "class",
      "lowest_taxon": "Methanobacteria"
    }
  ]
}
```

## 2) assemble

- **作用：** 把宏基因组读段组装成 contig。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算，也没有接入任何组装程序；需要另行运行 MEGAHIT、metaSPAdes 之类的宏基因组组装程序。不会编造结果。

- **API：** `dnakit.microbiome.assemble(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a metagenome assembler (MEGAHIT, metaSPAdes) on PATH`，`message` 为 `Metagenome assembly is not implemented natively; no assembler is wired into this function.`

- **示例代码：**

```python
import json
from dnakit.microbiome import assemble

result = assemble(["ACGTACGTACGGT", "CGTACGGTTTACA"])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "assemble",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 34,
    "sequence_bases": 26,
    "lengths": [
      13,
      13
    ],
    "sha256": "b1c3a29937b1ddad6d5b46dca97b726278cd13bba542ce3c40feb88148f58dd9"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a metagenome assembler (MEGAHIT, metaSPAdes) on PATH",
  "message": "Metagenome assembly is not implemented natively; no assembler is wired into this function."
}
```

## 3) bin

- **作用：** 把宏基因组 contig 分箱为基因组 bin。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。分箱需要读段覆盖度数据以及 MetaBAT2、MaxBin2、CONCOCT 之类的分箱程序，本函数均未接入。不会编造结果。

- **API：** `dnakit.microbiome.bin(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a contig binner (MetaBAT2, MaxBin2, CONCOCT) plus read-coverage data`，`message` 为 `Contig binning is not implemented natively.`

- **示例代码：**

```python
import json
from dnakit.microbiome import bin as bin_contigs

result = bin_contigs({"contig1": "ACGTACGTACGGT", "contig2": "CGTACGGTTTACA"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "bin",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 56,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "d7a452d267f1cfc5abdacd676e9400fb869a1b9218a885ce2183a8ad6016f4a7"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a contig binner (MetaBAT2, MaxBin2, CONCOCT) plus read-coverage data",
  "message": "Contig binning is not implemented natively."
}
```

## 4) build_cooccurrence_network

- **作用：** 根据各样本中相对丰度的两两相关性构建物种共现网络。

- **处理方法：** 先把计数转换为每个样本内的相对丰度，在少于 `min_samples_present` 个样本中非零的物种被丢弃。对剩余物种两两计算相关系数：`method="pearson"` 时为相对丰度的 Pearson r，`"spearman"` 时为其平均秩的 Pearson r，即 Spearman rho；|r| >= `min_abs_correlation` 的物种对成为边。不做显著性检验或多重检验校正，相对丰度的相关性存在组成性偏差（调用会附带警告）。至少需要 3 个样本，保留的物种至多 2000 个；含常数物种（r 无定义）的物种对被跳过。

- **API：** `dnakit.microbiome.build_cooccurrence_network(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格。`options`：`method`（`"spearman"` 默认 | `"pearson"`）、`min_abs_correlation`（默认 0.6）、`min_samples_present`（默认 3）、文本表的 `orientation`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `method`、`min_abs_correlation`、`nodes`（保留的物种）、`edges`（`source`、`target`、`correlation`，按 |r| 降序）和 `n_samples`。

- **示例代码：**

```python
import json
from dnakit.microbiome import build_cooccurrence_network

table = {
    "s1": {"a": 10, "b": 20, "c": 30, "d": 1},
    "s2": {"a": 20, "b": 41, "c": 5, "d": 2},
    "s3": {"a": 30, "b": 58, "c": 25, "d": 3},
    "s4": {"a": 40, "b": 80, "c": 9, "d": 4},
}
result = build_cooccurrence_network(table, options={"method": "spearman", "min_abs_correlation": 0.9})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "build_cooccurrence_network",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 170,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7db91f0cb61b29d7f403f76b490c90844fc8908f32ac9d7676109e0c58527b77"
  },
  "warnings": [
    "Correlations of relative abundances are compositional; interpret with care."
  ],
  "status": "ok",
  "result": {
    "method": "spearman",
    "min_abs_correlation": 0.9,
    "nodes": [
      "a",
      "b",
      "c",
      "d"
    ],
    "edges": [
      {
        "source": "a",
        "target": "c",
        "correlation": -1.0
      },
      {
        "source": "a",
        "target": "d",
        "correlation": 1.0
      },
      {
        "source": "c",
        "target": "d",
        "correlation": -1.0
      }
    ],
    "n_samples": 4
  }
}
```

## 5) build_microbial_gene_catalog

- **作用：** 构建非冗余的微生物基因目录。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。基因目录需要基因预测（如 Prodigal）和按相似度聚类（CD-HIT 或 MMseqs2），本函数均未原生实现。不会编造结果。

- **API：** `dnakit.microbiome.build_microbial_gene_catalog(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a gene predictor (Prodigal) and a sequence clusterer (CD-HIT, MMseqs2) on PATH`，`message` 为 `Gene-catalog construction needs gene calling and identity clustering, which are not implemented natively.`

- **示例代码：**

```python
import json
from dnakit.microbiome import build_microbial_gene_catalog

result = build_microbial_gene_catalog({"contig1": "ATGAAACCCGGGTTTTAA"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "build_microbial_gene_catalog",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 33,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "546f6096d459247e3dae06ce404ad533735c07131858cf5b58f6222ad31b5654"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a gene predictor (Prodigal) and a sequence clusterer (CD-HIT, MMseqs2) on PATH",
  "message": "Gene-catalog construction needs gene calling and identity clustering, which are not implemented natively."
}
```

## 6) build_pathogen_transmission_tree

- **作用：** 推断病原体传播树。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。传播树推断需要 outbreaker2、TransPhylo 之类的统计（流行病学/系统动力学）模型和采样日期，本函数均未实现。不会编造结果。

- **API：** `dnakit.microbiome.build_pathogen_transmission_tree(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a transmission inference model (outbreaker2, TransPhylo) with sampling dates`，`message` 为 `Transmission-tree inference requires a statistical model and dated samples; none is implemented.`

- **示例代码：**

```python
import json
from dnakit.microbiome import build_pathogen_transmission_tree

result = build_pathogen_transmission_tree({"patientA": "ACGTACGT", "patientB": "ACGTACGA"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "build_pathogen_transmission_tree",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e245b9769d512ad6d4af24d7600f9e72ec40ec8f60b79441a98398ef3dd81411"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a transmission inference model (outbreaker2, TransPhylo) with sampling dates",
  "message": "Transmission-tree inference requires a statistical model and dated samples; none is implemented."
}
```

## 7) calculate_ani

- **作用：** 用 Mash 公式根据共有 k-mer 估计基因组间的平均核苷酸一致性；这是估计值，不是基于比对的 ANI。

- **处理方法：** 对每一对序列构建精确的 k-mer 集合 A 和 B（含非 ACGT 符号的 k-mer 被跳过；`canonical` 为真时 k-mer 与其反向互补视为同一个词）。由 Jaccard 指数 j = |A&B| / |A|B|，Mash 距离 D = -ln(2j / (1 + j)) / k，`ani_estimate` = 1 - D（没有共有 k-mer 时为 `null`）；同时报告包含度 |A&B|/|A| 和 |A&B|/|B|。与 Mash 本身不同，这里不使用 MinHash 草图，而是比较完整的 k-mer 集合。它不是 FastANI 或 ANIb（没有片段比对），估计值只在基因组大小相近且 ANI 高于约 0.8 时可靠；此处没有拿它与基于比对的 ANI 做验证。调用会附带相应警告。

- **API：** `dnakit.microbiome.calculate_ani(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：DNA 序列：字符串列表、映射 `{id: 序列}` 或 FASTA 文本/路径（ACGT 加 IUPAC 简并码；U 视为 T）；2 到 50 条序列，总长至多 20 Mb。`options`：`k`（默认 21，范围 2-32）、`canonical`（默认 `true`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `method`（`"mash_kmer_jaccard"`）和 `pairs`（每个无序序列对一项：`a`、`b`、`k`、`shared_kmers`、`jaccard`、`containment_a_in_b`、`containment_b_in_a`、`mash_distance`、`ani_estimate`）。

- **示例代码：**

```python
import json
import random
from dnakit.microbiome import calculate_ani

rng = random.Random(1)
genome_a = "".join(rng.choice("ACGT") for _ in range(3000))
genome_b = "".join(rng.choice("ACGT".replace(c, "")) if i % 50 == 0 else c   # change every 50th base
                   for i, c in enumerate(genome_a))
result = calculate_ani({"A": genome_a, "B": genome_b}, options={"k": 21})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "calculate_ani",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 6018,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "b5509626509603ca23fb0e0c995cf060a54a720b39e6e6a3a06732cf2ab35eaf"
  },
  "warnings": [
    "k-mer set (Mash-formula) ANI estimate, not an alignment-based ANI."
  ],
  "status": "ok",
  "result": {
    "method": "mash_kmer_jaccard",
    "pairs": [
      {
        "a": "A",
        "b": "B",
        "k": 21,
        "shared_kmers": 1740.0,
        "jaccard": 0.41232227488151657,
        "containment_a_in_b": 0.5838926174496645,
        "containment_b_in_a": 0.5838926174496645,
        "mash_distance": 0.025620866061470255,
        "ani_estimate": 0.9743791339385297
      }
    ]
  }
}
```

## 8) calculate_differential_abundance

- **作用：** 检验每个物种在两组样本间的相对丰度差异（Mann-Whitney U 检验，Benjamini-Hochberg 校正）。

- **处理方法：** 先把计数转换为每个样本内的相对丰度，出现样本比例低于 `min_prevalence` 的物种被跳过。对其余每个物种计算两组均值，`log2_fold_change` = log2((mean_b + pseudocount) / (mean_a + pseudocount))，并用 `scipy.stats.mannwhitneyu` 做 A 组对 B 组的双侧 Mann-Whitney U 检验（报告 A 组的 U；所有值相同时 U = n_a*n_b/2、p = 1）。所有被检验物种的 p 值经 Benjamini-Hochberg 校正得到 `q_value`，物种按 p 值排序。这是对组成性相对丰度的秩检验：除相对丰度外没有其他归一化，不含协变量，也没有零膨胀模型。需要 scipy：缺少时函数返回 `status: "not_implemented"`，`required_dependency` 为 `pip install scipy`。

- **API：** `dnakit.microbiome.calculate_differential_abundance(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格。`options`：`groups`（必填：按样本顺序排列的分组标签列表，或映射 样本 -> 标签）；`group_a`/`group_b`（默认取排序后的两个标签，此时标签必须恰好两个）；`pseudocount`（默认 1e-6）；`min_prevalence`（默认 0.1）。每组至少需要 2 个样本。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `group_a`、`group_b`、`n_a`、`n_b`、`test`（`"mann_whitney_u"`）和 `taxa`（每个物种有 `taxon`、`mean_a`、`mean_b`、`log2_fold_change`、`u_statistic`、`p_value`、`q_value`）。

- **示例代码：**

```python
import json
from dnakit.microbiome import calculate_differential_abundance

table = {
    "c1": {"a": 90, "b": 10, "c": 50}, "c2": {"a": 80, "b": 20, "c": 55},
    "c3": {"a": 85, "b": 15, "c": 45},
    "t1": {"a": 20, "b": 80, "c": 50}, "t2": {"a": 30, "b": 70, "c": 52},
    "t3": {"a": 25, "b": 75, "c": 48},
}
groups = {"c1": "control", "c2": "control", "c3": "control", "t1": "treated", "t2": "treated", "t3": "treated"}
result = calculate_differential_abundance(table, options={"groups": groups})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "calculate_differential_abundance",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 210,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3d66b69cc2ad0a8cda972aa64d7d905e73c3c44f81aff7edb99364aaa3b4ad65"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "group_a": "control",
    "group_b": "treated",
    "n_a": 3,
    "n_b": 3,
    "taxa": [
      {
        "taxon": "a",
        "mean_a": 0.5674453096032629,
        "mean_b": 0.16654022443496128,
        "log2_fold_change": -1.7686045585338197,
        "u_statistic": 9.0,
        "p_value": 0.1,
        "q_value": 0.15000000000000002
      },
      {
        "taxon": "b",
        "mean_a": 0.09971573353108391,
        "mean_b": 0.5002054686265213,
        "log2_fold_change": 2.326616184784275,
        "u_statistic": 0.0,
        "p_value": 0.1,
        "q_value": 0.15000000000000002
      },
      {
        "taxon": "c",
        "mean_a": 0.3328389568656532,
        "mean_b": 0.3332543069385175,
        "log2_fold_change": 0.0017992128170873996,
        "u_statistic": 4.5,
        "p_value": 1.0,
        "q_value": 1.0
      }
    ],
    "test": "mann_whitney_u"
  }
}
```

## 9) calculate_digital_ddh

- **作用：** 计算数字 DNA-DNA 杂交值（dDDH）。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。dDDH 需要基因组比对（BLAST、BLAT 或 MUMmer 的 HSP）和 GGDC 回归公式，均未实现。如需基于 k-mer 的一致性估计，请用 `calculate_ani`。

- **API：** `dnakit.microbiome.calculate_digital_ddh(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `genome alignments (BLAST/BLAT/MUMmer HSPs) and the GGDC regression formulas`，`message` 为 `dDDH is not implemented; use calculate_ani for a k-mer identity estimate instead.`

- **示例代码：**

```python
import json
from dnakit.microbiome import calculate_digital_ddh

result = calculate_digital_ddh({"A": "ACGTACGTACGT", "B": "ACGTACGTACGA"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "calculate_digital_ddh",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 42,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "18186543adfa540845c61ad49fbb8968a47e742c04b3030c93180cc48804768a"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "genome alignments (BLAST/BLAT/MUMmer HSPs) and the GGDC regression formulas",
  "message": "dDDH is not implemented; use calculate_ani for a k-mer identity estimate instead."
}
```

## 10) classify

- **作用：** 为读段或 contig 分配分类。没有原生实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 没有参考数据库，因此不做任何计算。返回信息指向同一模块中的 `classify_*` 封装（如 `classify_kraken2`、`classify_kaiju`、`classify_centrifuge`、`classify_gtdbtk`），它们调用已安装的分类器和数据库。不会编造结果。

- **API：** `dnakit.microbiome.classify(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a taxonomic classifier database (Kraken2, Kaiju, Centrifuge, GTDB-Tk; see classify_kraken2 etc.)`，`message` 为 `Native taxonomic classification is not implemented; use the classify_* wrappers with an installed database.`

- **示例代码：**

```python
import json
from dnakit.microbiome import classify

result = classify(["ACGTACGTACGTACGT"])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "classify",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 20,
    "sequence_bases": 16,
    "lengths": [
      16
    ],
    "sha256": "87cd4b17356c3dee4b13d54fba3beefea5678a7b651616a23b27c782da631af2"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a taxonomic classifier database (Kraken2, Kaiju, Centrifuge, GTDB-Tk; see classify_kraken2 etc.)",
  "message": "Native taxonomic classification is not implemented; use the classify_* wrappers with an installed database."
}
```

## 11) compare

- **作用：** 根据丰度表计算样本间的 beta 多样性距离矩阵。

- **处理方法：** 对每一对样本 a、b：`braycurtis` = sum|a-b| / sum(a+b)（用原始或相对计数）；`jaccard` = 1 - 共有/并集（按计数 > 0 判断存在）；`aitchison` = 计数 + `pseudocount` 后的中心对数比（CLR）谱之间的欧氏距离；`euclidean` = 数值的欧氏距离。`relative` 为真时先把每个样本转换为相对丰度。矩阵对称、对角线为 0；不做稀释，也不使用系统发育距离（UniFrac）。

- **API：** `dnakit.microbiome.compare(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格；至少 2 个、至多 2000 个样本。`options`：`entity`（只支持 `"beta_diversity"`）；`metric`（`"braycurtis"` 默认 | `"jaccard"` | `"aitchison"` | `"euclidean"`）；`pseudocount`（默认 1.0，供 `aitchison` 使用）；`relative`（默认 `false`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `metric`、`samples` 和对称的 `distance_matrix`（行列顺序与 `samples` 一致）。

- **示例代码：**

```python
import json
from dnakit.microbiome import compare

table = {"s1": {"a": 5, "b": 5, "c": 5, "d": 5},
         "s2": {"a": 1, "b": 1, "c": 2, "d": 6},
         "s3": {"a": 10, "b": 0, "c": 3, "d": 0}}
result = compare(table, options={"metric": "braycurtis"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "compare",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 121,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "4193a3bec6078b3f816617bfa2ee0d04d6170d8243e9c6557184867a5162df81"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "metric": "braycurtis",
    "samples": [
      "s1",
      "s2",
      "s3"
    ],
    "distance_matrix": [
      [
        0.0,
        0.4,
        0.5151515151515151
      ],
      [
        0.4,
        0.0,
        0.7391304347826086
      ],
      [
        0.5151515151515151,
        0.7391304347826086,
        0.0
      ]
    ]
  }
}
```

## 12) dereplicate

- **作用：** 把冗余的基因组/序列记录合并成簇，每簇一个代表。

- **处理方法：** `method="exact"` 把完全相同的序列分为一组（不区分大小写；`canonical` 为真时，互为反向互补的序列也视为相同）。`method="ani"` 为贪心：序列按长度从长到短（再按 id）排序，每条序列加入第一个其代表序列与之 k-mer ANI 估计值 >= `min_ani` 的已有簇，否则自成新簇。ANI 估计值是 `calculate_ani` 在精确 k-mer 集合上的 Mash 公式估计（此处 `k` 默认 16），所以 `"ani"` 模式是 k-mer 估计，不是 dRep 那样基于比对的聚类；并且只与各簇的代表序列比较。

- **API：** `dnakit.microbiome.dereplicate(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：DNA 序列：字符串列表、映射 `{id: 序列}` 或 FASTA 文本/路径（ACGT 加 IUPAC 简并码；U 视为 T）；1 到 2000 条序列，总长至多 50 Mb。`options`：`method`（`"exact"` 默认 | `"ani"`）、`canonical`（默认 `true`）、`k`（默认 16，范围 2-32，用于 `ani` 模式）、`min_ani`（默认 0.95，用于 `ani` 模式）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `method`、`input_count`、`cluster_count` 和 `clusters`（每个簇含 `representative` 和 `members`）。

- **示例代码：**

```python
import json
import random
from dnakit.microbiome import dereplicate

rng = random.Random(2)
g1 = "".join(rng.choice("ACGT") for _ in range(2000))
g2 = "".join(rng.choice("ACGT".replace(c, "")) if i % 200 == 0 else c for i, c in enumerate(g1))   # ~0.5 % changed
g3 = "".join(rng.choice("ACGT") for _ in range(1800))                                              # unrelated
result = dereplicate({"g1": g1, "g2": g2, "g3": g3}, options={"method": "ani", "k": 16, "min_ani": 0.95})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "dereplicate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 5830,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "340d4746b48c4e13c53b16f2046187b727af0cd4ee1b052ad4cea41fe5176b8e"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "ani",
    "input_count": 3,
    "cluster_count": 2,
    "clusters": [
      {
        "representative": "g1",
        "members": [
          "g1",
          "g2"
        ]
      },
      {
        "representative": "g3",
        "members": [
          "g3"
        ]
      }
    ]
  }
}
```

## 13) detect

- **作用：** 检测物种的流行率和核心微生物组（在足够多样本中高于检测水平的物种）。

- **处理方法：** 物种在样本中的相对丰度为 计数 / 样本总数。相对丰度严格大于 `detection` 时认为该物种存在于该样本；流行率为其存在的样本比例。核心物种的流行率 >= `prevalence`。物种按流行率、再按平均相对丰度排序。不做统计检验，也不使用按丰度加权的核心定义。

- **API：** `dnakit.microbiome.detect(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格；每个样本的总数必须为正。`options`：`entity`（`"core_microbiome"` 默认，会给出核心物种 | `"prevalence"`，只给出表）；`detection`（默认 0.0）；`prevalence`（默认 0.5）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `detection`、`n_samples` 和 `prevalence_table`（每个物种有 `taxon`、`samples_present`、`prevalence`、`mean_relative_abundance`）；`core_microbiome` 时还有 `prevalence_threshold` 和 `core_taxa`。

- **示例代码：**

```python
import json
from dnakit.microbiome import detect

table = {"s1": {"a": 50, "b": 10, "c": 0}, "s2": {"a": 40, "b": 0, "c": 0}, "s3": {"a": 30, "b": 5, "c": 1}}
result = detect(table, options={"entity": "core_microbiome", "detection": 0.0, "prevalence": 0.6})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "detect",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 100,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "417c10ac425501d532651729f25cccaad64f71d0508e4209810a058e285d982c"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "detection": 0.0,
    "n_samples": 3,
    "prevalence_table": [
      {
        "taxon": "a",
        "samples_present": 3,
        "prevalence": 1.0,
        "mean_relative_abundance": 0.888888888888889
      },
      {
        "taxon": "b",
        "samples_present": 2,
        "prevalence": 0.6666666666666666,
        "mean_relative_abundance": 0.10185185185185186
      },
      {
        "taxon": "c",
        "samples_present": 1,
        "prevalence": 0.3333333333333333,
        "mean_relative_abundance": 0.009259259259259259
      }
    ],
    "prevalence_threshold": 0.6,
    "core_taxa": [
      "a",
      "b"
    ]
  }
}
```

## 14) export

- **作用：** 把丰度表序列化为 QIIME 风格的 TSV（物种为行）或 JSON，可选写入文件。

- **处理方法：** 先解析并校验表，再输出为 TSV（表头为 `#OTU ID` 加各样本 id，每个物种一行，整数不带小数点）或包含 `samples`、`taxa`、`counts` 的 JSON（键已排序）。只有 `overwrite` 为真时才覆盖已有输出文件。文本同时返回在 `result` 中。

- **API：** `dnakit.microbiome.export(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格。`options`：`format`（`"tsv"` 默认 | `"json"`）；`output`（可选文件路径）；`overwrite`（默认 `false`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `format`、`text` 和 `output_path`（已写入的文件，未写入为 `null`）。

- **示例代码：**

```python
import json
from dnakit.microbiome import export

table = {"s1": {"a": 5, "b": 0}, "s2": {"a": 1, "b": 7}}
result = export(table, options={"format": "tsv"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "826172e6734fde444cb3c4c03ccd74ba350458e29235dc64e5827464038ab612"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "text": "#OTU ID\ts1\ts2\na\t5\t1\nb\t0\t7\n",
    "output_path": null
  }
}
```

## 15) infer_viral_quasispecies

- **作用：** 重建病毒准种单倍型。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。单倍型重建需要已比对的读段以及 QuRe、ShoRAH、PredictHaplo 之类的程序，本函数均未接入。不会编造结果。

- **API：** `dnakit.microbiome.infer_viral_quasispecies(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a viral haplotype reconstruction tool (QuRe, ShoRAH, PredictHaplo) and aligned reads`，`message` 为 `Quasispecies inference is not implemented natively.`

- **示例代码：**

```python
import json
from dnakit.microbiome import infer_viral_quasispecies

result = infer_viral_quasispecies(["ACGTACGTAC", "ACGTACGTAA"])
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "infer_viral_quasispecies",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 28,
    "sequence_bases": 20,
    "lengths": [
      10,
      10
    ],
    "sha256": "87931cb5c372cda553027980d7b4eb79ea0a83bd1432f85ad02bba23469f8907"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a viral haplotype reconstruction tool (QuRe, ShoRAH, PredictHaplo) and aligned reads",
  "message": "Quasispecies inference is not implemented natively."
}
```

## 16) plot

- **作用：** 绘制堆叠的相对丰度条形图（每个样本一根柱），输出 SVG 字符串。

- **处理方法：** 相对丰度以堆叠方块的形式画在手写的 SVG 中（不依赖绘图库）。按所有样本总计数最大的 `top_n` 个物种各取固定 10 色调色板中的一种颜色，其余物种合并为灰色的 `Other`；图例列出主要物种。每个方块是带 `data-sample`、`data-taxon`、`data-fraction` 属性的 `<rect class="bar">`。

- **API：** `dnakit.microbiome.plot(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格；每个样本的总数必须为正。`options`：`top_n`（默认 8，范围 1-10）；`output`（可选 `.svg` 路径）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `svg`（SVG 文本）和 `output_path`（已写入的文件，未写入为 `null`）。

- **示例代码：**

```python
import json
from dnakit.microbiome import plot

table = {"s1": {"a": 8, "b": 2}, "s2": {"a": 3, "b": 7}}
result = plot(table, options={"top_n": 2})
print(json.dumps(result, indent=2))
```

- **示例输出：** 输出中较长的数组和字符串用 `...` 缩短。

```json
{
  "module": "microbiome",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ff58a93da8baaf87e71302a304526230cd06541032b0298e374387dc98fd9d4b"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"312\" height=\"280\"><rect class=\"bar\" data-sample=\"s1\" data-taxon=\"a\" data-fraction=\"...",
    "output_path": null
  }
}
```

## 17) predict

- **作用：** 预测微生物组的表型或功能。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。由丰度谱预测功能或表型需要 PICRUSt2、BugBase 之类的已训练参考数据，DNAKit 并未内置。不会编造结果。

- **API：** `dnakit.microbiome.predict(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `a trained model or reference database (e.g. PICRUSt2, BugBase)`，`message` 为 `Prediction from a microbiome profile needs trained reference data that is not bundled.`

- **示例代码：**

```python
import json
from dnakit.microbiome import predict

result = predict({"s1": {"a": 8, "b": 2}})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "predict",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 24,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "beafe4e6da1b369b3fbf3192253bb6e2e35f1331de917590e4047cfcb11ee9ce"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "a trained model or reference database (e.g. PICRUSt2, BugBase)",
  "message": "Prediction from a microbiome profile needs trained reference data that is not bundled."
}
```

## 18) qc

- **作用：** 检查丰度表的质量：每个样本的深度、丰富度、零值比例以及低深度标记。

- **处理方法：** 对每个样本：`total`（计数之和）、`observed`（计数 > 0 的物种数）、`zero_fraction` = 1 - observed / 物种数，总数低于 `min_depth` 时 `low_depth` 为真。在所有样本中都为 0 的物种列为 `empty_taxa`。深度阈值是用户设定的经验值，并非由数据推出。

- **API：** `dnakit.microbiome.qc(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格。`options`：`entity`（只支持 `"abundance_table"`）；`min_depth`（默认 1000）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `n_samples`、`n_taxa`、`min_depth`、`samples`（每个样本有 `sample`、`total`、`observed`、`zero_fraction`、`low_depth`）、`empty_taxa` 和 `flagged_samples`。

- **示例代码：**

```python
import json
from dnakit.microbiome import qc

table = {"s1": {"a": 900, "b": 100, "c": 0}, "s2": {"a": 30, "b": 0, "c": 0}}
result = qc(table, options={"min_depth": 1000})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "qc",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 69,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3f254bebd28a310387a75982cd4bf7766e8f0eb10be64aee81359716e2eb18ba"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n_samples": 2,
    "n_taxa": 3,
    "min_depth": 1000.0,
    "samples": [
      {
        "sample": "s1",
        "total": 1000.0,
        "observed": 2,
        "zero_fraction": 0.33333333333333337,
        "low_depth": false
      },
      {
        "sample": "s2",
        "total": 30.0,
        "observed": 1,
        "zero_fraction": 0.6666666666666667,
        "low_depth": true
      }
    ],
    "empty_taxa": [
      "c"
    ],
    "flagged_samples": [
      "s2"
    ]
  }
}
```

## 19) quantify

- **作用：** 对丰度表做定量：相对丰度、alpha 多样性、稀释或中心对数比（CLR）。这里不做读段比对，计数必须已经存在。

- **处理方法：** `relative_abundance`（别名 `abundance`）：计数 / 样本总数。`alpha_diversity`，对每个样本：观测物种数、Shannon H = -sum p ln p（默认自然对数，可用 `log_base` 改底）、Simpson（Gini-Simpson）= 1 - sum p^2、逆 Simpson = 1 / sum p^2、Pielou 均匀度 = H / ln(观测数)（只有一个物种时为 `null`）、偏差校正的 Chao1 = S_obs + F1 (F1 - 1) / (2 (F2 + 1)) 和 Good 覆盖度 = 1 - F1 / N，其中 F1/F2 为单例/双例物种数，N 为总计数（Chao1 与覆盖度只对整数计数计算，否则为 `null`）。`rarefaction`：用 `random.Random(seed)` 对每个样本无放回抽样到 `depth` 条读段；读段不足的样本被丢弃并列出（`drop_short` 为假时改为报错）。`clr`：log(计数 + `pseudocount`) 减去其在各物种上的均值。

- **API：** `dnakit.microbiome.quantify(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格；除 `rarefaction` 外每个样本的总数必须为正。`options`：`entity`（`"relative_abundance"` 默认 | `"abundance"` | `"alpha_diversity"` | `"rarefaction"` | `"clr"`）；`log_base`（默认 e，用于 `alpha_diversity`）；`depth`（`rarefaction` 必填，仅限整数计数）、`seed`（默认 0）、`drop_short`（默认 `true`）；`pseudocount`（默认 1.0，用于 `clr`）。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 取决于 `entity`：`relative_abundance` 给出 `samples`、`taxa`、`relative_abundance`；`alpha_diversity` 给出 `samples` 和 `metrics`（每个样本有 `observed`、`shannon`、`simpson`、`inverse_simpson`、`pielou_evenness`、`chao1`、`goods_coverage`）；`rarefaction` 给出 `depth`、`dropped_samples`、`samples`、`taxa`、`counts`；`clr` 给出 `samples`、`taxa`、`clr`。

- **示例代码：**

```python
import json
from dnakit.microbiome import quantify

table = {"s1": {"a": 5, "b": 5, "c": 5, "d": 5}, "s2": {"a": 1, "b": 1, "c": 2, "d": 6}}
result = quantify(table, options={"entity": "alpha_diversity"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "quantify",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 80,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ea761b00be122b589f5902ee4cb3fb88280b57aad31592d0352b9e8c8102d34f"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [
      "s1",
      "s2"
    ],
    "metrics": {
      "s1": {
        "observed": 4.0,
        "shannon": 1.3862943611198906,
        "simpson": 0.75,
        "inverse_simpson": 4.0,
        "pielou_evenness": 1.0,
        "chao1": 4.0,
        "goods_coverage": 1.0
      },
      "s2": {
        "observed": 4.0,
        "shannon": 1.0888999753452238,
        "simpson": 0.5800000000000001,
        "inverse_simpson": 2.380952380952381,
        "pielou_evenness": 0.7854752972273344,
        "chao1": 4.5,
        "goods_coverage": 0.8
      }
    }
  }
}
```

## 20) read

- **作用：** 把丰度表（QIIME/BIOM-TSV 或 CSV 文本、文件路径或内存结构）读入为 样本 x 物种 的计数矩阵。

- **处理方法：** 文本按制表符（或逗号）切分；第一条不是 `#` 注释的行（以 `#OTU`、`#Feature` 或 `#ID` 开头的表头视为表头）给出样本 id，第一列给出物种 id。`taxa_rows` 时表会被转置，使结果的行为样本。数值必须是有限的非负数，id 必须唯一，各行宽度必须一致，否则抛出输入错误。不做任何归一化。

- **API：** `dnakit.microbiome.read(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：丰度表：映射 `{样本: {物种: 计数}}`、`{"samples": [...], "taxa": [...], "counts": [[...]]}`（行为样本）、带 `options["samples"]`/`options["taxa"]` 的计数行列表，或 QIIME/BIOM-TSV、CSV 文本/路径（`options["orientation"]`：`"taxa_rows"` 默认 | `"samples_rows"`）；数值须为有限的非负数，至多 5000 个样本、200000 个物种、2000 万个单元格。`options`：文本表的 `orientation`。

- **输出：** 返回封装包含 `status`（`"ok"`）、`module`（`"microbiome"`）、`operation`、`provenance`（`implementation`、`backend`）、`input_summary`（通用的文本级摘要：记录数、文本长度、碱基数、SHA-256；映射类型的输入按 1 条记录计）、`warnings` 和 `result`。`result` 含 `samples`、`taxa` 和 `counts`（行为样本，数值为浮点数）。

- **示例代码：**

```python
import json
from dnakit.microbiome import read

text = "#OTU ID\ts1\ts2\nEscherichia\t10\t0\nBacteroides\t5\t20\n"
result = read(text)
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a3f580da762132a4f472fe7656ecdcf87576a39c444c4a62471110133decc416"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "samples": [
      "s1",
      "s2"
    ],
    "taxa": [
      "Escherichia",
      "Bacteroides"
    ],
    "counts": [
      [
        10.0,
        5.0
      ],
      [
        0.0,
        20.0
      ]
    ]
  }
}
```

## 21) track

- **作用：** 跨样本追踪菌株或病原体。尚未实现：该函数始终返回 `status: "not_implemented"`。

- **处理方法：** 不做任何计算。菌株追踪需要菌株水平的标记基因或 SNP 参考（如 StrainPhlAn、inStrain）和流行病学元数据，本函数无法原生提供。不会编造结果。

- **API：** `dnakit.microbiome.track(input[required], options[optional])`（正式 DNAKit 命名空间；能否通过 Agent/MCP 调用取决于 JSON 调用方案）。

- **输入：** `input`：任意值（仅做摘要）。`options`：被忽略。

- **输出：** 调用返回 `status: "not_implemented"`（`provenance.implementation` 也是 `"not_implemented"`），含常规封装字段，没有 `result`。`required_dependency` 为 `strain-level marker/SNP references (e.g. StrainPhlAn, inStrain) and epidemiological metadata`，`message` 为 `Strain/pathogen tracking is not implemented natively.`

- **示例代码：**

```python
import json
from dnakit.microbiome import track

result = track({"isolate1": "ACGTACGT", "isolate2": "ACGTACGA"})
print(json.dumps(result, indent=2))
```

- **示例输出：**

```json
{
  "module": "microbiome",
  "operation": "track",
  "provenance": {
    "implementation": "not_implemented",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 48,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "130efd49bcf3ffab1415417f4a23967523f5c89a74650caf0c6f17f20e766ba1"
  },
  "warnings": [],
  "status": "not_implemented",
  "required_dependency": "strain-level marker/SNP references (e.g. StrainPhlAn, inStrain) and epidemiological metadata",
  "message": "Strain/pathogen tracking is not implemented natively."
}
```
