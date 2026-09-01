# 可选生物信息功能

本页收录 DNAKit 原来没有、现在通过可选科学后端直接调用的功能。公开 API 按功能放在
`annotation`、`comparative` 和 `molbio` 中，
不需要使用 `integrations.tooluniverse` 一类入口。

安装可选依赖：

```bash
python -m pip install "dnakit[external-tools]"
```

这些方法不会训练或微调模型。Ensembl、ClinVar、dbSNP 和 gnomAD 需要联网；dN/dS 与
Golden Gate 在本地执行。每次调用只加载对应的一个
后端工具。统一返回 `ProviderResult`：主要结果位于 `result.data`，调用参数、后端版本和
来源分别位于 `result.parameters`、`result.metadata` 和 `result.provenance`。

## 1) HGVS 变异 VEP 注释

- **作用：** 通过 Ensembl VEP 返回 HGVS 变异的转录本后果和影响信息。
- **API：** `dnakit.annotation.annotate_variant_vep(hgvs_notation[必须], species[可选])`
- **输入：** 一个 HGVS 表达式；默认物种为 `human`。
- **示例代码：**

```python
from dnakit.annotation import annotate_variant_vep

result = annotate_variant_vep("NM_000546.6:c.215C>G")
print(result.data)
```

- **示例结果：** Ensembl VEP 返回的结构化 consequence 记录。

## 2) rsID 的 VEP 注释

- **作用：** 使用 Ensembl VEP 注释一个 dbSNP rsID。
- **API：** `dnakit.annotation.annotate_rsid_vep(rsid[必须], species[可选])`
- **输入：** 一个 rsID，例如 `rs28934578`。
- **示例代码：**

```python
from dnakit.annotation import annotate_rsid_vep

result = annotate_rsid_vep("rs28934578")
print(result.data)
```

- **示例结果：** 与该 rsID 对应的 VEP consequence 记录。

## 3) 变异标识转换

- **作用：** 将 rsID、HGVS 等变异标识转换为 Ensembl 可提供的等价表示。
- **API：** `dnakit.annotation.recode_variant(variant_id[必须], species[可选])`
- **输入：** 一个变异标识和可选物种。
- **示例代码：**

```python
from dnakit.annotation import recode_variant

result = recode_variant("rs28934578")
print(result.data)
```

- **示例结果：** 对应的 HGVS、基因组位置和其他标识。

## 4) ClinVar 变异搜索

- **作用：** 按基因、疾病、ClinVar ID 或 HGVS/蛋白变化搜索 ClinVar。
- **API：** `dnakit.annotation.search_clinvar_variants(gene[可选], condition[可选], variant_id[可选], variant_name[可选], clinical_significance[可选], limit[可选])`
- **输入：** 至少提供一个查询字段；`limit` 为 1–100。
- **示例代码：**

```python
from dnakit.annotation import search_clinvar_variants

result = search_clinvar_variants(gene="TP53", variant_name="c.215C>G", limit=10)
print(result.data)
```

- **示例结果：** ClinVar 变异标识和摘要记录。

## 5) ClinVar 变异详情

- **作用：** 按 ClinVar variation ID 获取 accession、基因、坐标和审核状态。
- **API：** `dnakit.annotation.get_clinvar_variant(variant_id[必须])`
- **输入：** ClinVar variation ID，例如 `12345`。
- **示例代码：**

```python
from dnakit.annotation import get_clinvar_variant

result = get_clinvar_variant("12345")
print(result.data)
```

- **示例结果：** ClinVar 变异详情。

## 6) ClinVar 临床意义

- **作用：** 获取指定 ClinVar variation ID 的临床意义和解释记录。
- **API：** `dnakit.annotation.get_clinvar_significance(variant_id[必须])`
- **输入：** ClinVar variation ID。
- **示例代码：**

```python
from dnakit.annotation import get_clinvar_significance

result = get_clinvar_significance("12345")
print(result.data)
```

- **示例结果：** pathogenicity 分类和临床解释。

## 7) dbSNP 变异信息

- **作用：** 按 rsID 查询 dbSNP 等位基因及 GRCh37/GRCh38 坐标。
- **API：** `dnakit.annotation.get_dbsnp_variant(rsid[必须])`
- **输入：** 带或不带 `rs` 前缀的 rsID。
- **示例代码：**

```python
from dnakit.annotation import get_dbsnp_variant

result = get_dbsnp_variant("rs28934578")
print(result.data)
```

- **示例结果：** dbSNP 坐标、assembly 和等位基因记录。

## 8) dbSNP 等位基因频率

- **作用：** 获取指定 rsID 的群体等位基因频率记录。
- **API：** `dnakit.annotation.get_dbsnp_frequencies(rsid[必须])`
- **输入：** 一个 rsID。
- **示例代码：**

```python
from dnakit.annotation import get_dbsnp_frequencies

result = get_dbsnp_frequencies("rs28934578")
print(result.data)
```

- **示例结果：** 按研究或群体整理的 allele count/频率数据。

## 9) gnomAD 变异搜索

- **作用：** 使用 rsID 或变异文本查找 gnomAD 的标准 `chrom-pos-ref-alt` 标识。
- **API：** `dnakit.annotation.search_gnomad_variants(query[必须], dataset[可选])`
- **输入：** 查询文本；默认统一使用 `gnomad_r4`。
- **示例代码：**

```python
from dnakit.annotation import search_gnomad_variants

result = search_gnomad_variants("rs28934578")
print(result.data)
```

- **示例结果：** 匹配的 gnomAD variant ID 和基础信息。

## 10) gnomAD 变异汇总

- **作用：** 获取一个 gnomAD variant ID 的总体等位基因频率和元数据。
- **API：** `dnakit.annotation.get_gnomad_variant(variant_id[必须], dataset[可选])`
- **输入：** `chrom-pos-ref-alt` 标识；默认 `gnomad_r4`。
- **示例代码：**

```python
from dnakit.annotation import get_gnomad_variant

result = get_gnomad_variant("17-7675088-C-G")
print(result.data)
```

- **示例结果：** genome/exome 总体 allele count、allele number 和频率。

## 11) gnomAD 分群频率

- **作用：** 获取一个变异在不同 ancestry 和性别分组中的等位基因频率。
- **API：** `dnakit.annotation.get_gnomad_population_frequencies(variant_id[必须], dataset[可选])`
- **输入：** gnomAD variant ID；默认 `gnomad_r4`。
- **示例代码：**

```python
from dnakit.annotation import get_gnomad_population_frequencies

result = get_gnomad_population_frequencies("17-7675088-C-G")
print(result.data)
```

- **示例结果：** 按 ancestry/sex 分层的 AC、AN 和 AF。

## 12) dN/dS 选择分析

- **作用：** 对两条密码子对齐 CDS 计算 Nei–Gojobori dN、dS 和 dN/dS，并进行 Jukes–Cantor 修正。
- **API：** `dnakit.comparative.calculate_dn_ds(sequence_a[必须], sequence_b[必须])`
- **输入：** 等长、线性、canonical A/C/G/T、长度为 3 倍数的 `DNASequence`/`DNARecord`/`DNA`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.comparative import calculate_dn_ds

result = calculate_dn_ds(DNASequence("ATGGCTGAA"), DNASequence("ATGGCCGAG"))
print(result.data["dN"], result.data["dS"], result.data["dN_dS"])
```

- **示例结果：** `dN`、`dS`、`dN_dS`、位点/差异计数和解释；`dS=0` 时比值为 `None`。

## 13) Golden Gate 部件设计

- **作用：** 为多个 DNA 部件分配唯一非回文连接 overhang，并添加 BsaI/BbsI 位点。
- **API：** `dnakit.molbio.design_golden_gate(parts[必须], enzyme[可选])`
- **输入：** 至少两条 canonical DNA；设计阶段支持 `BsaI` 或 `BbsI`。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.molbio import design_golden_gate

result = design_golden_gate(
    (DNASequence("ATGGCTGAA"), DNASequence("GCCAAATAA")),
    enzyme="BsaI",
)
print(result.data["parts_with_overhangs"])
```

- **示例结果：** 每个部件的左右 overhang、完整合成序列和实验说明。

## 14) Golden Gate 反应组装

- **作用：** 模拟 Type IIS 酶切、按 overhang 匹配连接，并返回组装产物和片段顺序。
- **API：** `dnakit.molbio.assemble_golden_gate(fragments[必须], enzyme[可选], circular[可选], labels[可选])`
- **输入：** 至少两条带相应 Type IIS 位点的 canonical DNA；支持 BsaI、BbsI、Esp3I/BsmBI 和 SapI。
- **示例代码：**

```python
from dnakit import DNASequence
from dnakit.molbio import assemble_golden_gate

result = assemble_golden_gate(
    (
        DNASequence("GGTCTCAAAACATGGAGGAGCCGCAGTCAGATCCTAGCGTTGAATTCGGATCCCTTTTGAGACC"),
        DNASequence("GGTCTCAAAAGGGATCCAAGCTTACGTACGTATGGAGGAGCCGCAGTCAGATATTTTGAGACC"),
        DNASequence("GGTCTCAAAATCCTAGCGTTGAATTCGGATCCAAGCTTACGTACGTATGGAGCGTTTGAGACC"),
        DNASequence("GGTCTCAAACGATCGATCGATCGATCGATCGGTTTTGAGACC"),
    ),
    labels=("vector", "insert-1", "insert-2", "insert-3"),
)
print(result.data["product_sequence"])
```

- **示例结果：** 产物序列、长度、assembly order、junction overhang 和未使用片段。

## 使用边界

- 未安装可选依赖时会抛出 `BackendUnavailableError`，并提示安装 `dnakit[external-tools]`。
- 远程提供方可能限流、维护或改变数据库版本；请保存 `result.provenance`。
- ClinVar、VEP、dbSNP 和 gnomAD 输出仅用于研究注释，不能直接作为临床诊断。
- DNAKit 不公开任意 ToolUniverse 工具执行入口，只允许调用本页列出的功能。
