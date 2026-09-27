# 变异效应评分与实验数据

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: variant-scores -->

<span id="variant-scores"></span>

## 1) 查询整条蛋白质的 AlphaMissense 致病性分类 {#dnakit.databases.alpha_missense_get_protein_scores}

- **作用：** 按 UniProt 编号汇总蛋白质各位点的替换分类、位点平均分数与阈值；可用 max_residues 限制查询位点数。
- **处理方法：** 通过 `AlphaMissense` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.alpha_missense_get_protein_scores(*, uniprot_id, max_residues=None, verbose=True)`
- **输入：** 必须提供 `uniprot_id`。
- **输出：** 返回字典，包含 `uniprot_id`, `protein_length`, `scores`, `n_positions_returned` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.alpha_missense_get_protein_scores(
    uniprot_id='P01116',
    max_residues=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "uniprot_id": "<string>",
  "protein_length": "<integer>",
  "scores": [
    "<object>"
  ],
  "n_positions_returned": "<integer>",
  "n_positions_attempted": "<integer>",
  "max_residues_cap": "<integer>",
  "thresholds": "<object>",
  "pdb_download": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.alpha_missense_get_protein_scores-parameters)

## 2) 查询蛋白质指定位点的替换分类与平均评分 {#dnakit.databases.alpha_missense_get_residue_scores}

- **作用：** 查看指定位点的良性、不确定和致病替换分类，以及位点平均分数；返回值不是每种替换各自的数值分数。
- **处理方法：** 通过 `AlphaMissense` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.alpha_missense_get_residue_scores(*, uniprot_id, position, verbose=True)`
- **输入：** 必须提供 `uniprot_id`, `position`。
- **输出：** 返回字典，包含 `status`, `data`, `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.alpha_missense_get_residue_scores(
    uniprot_id='P15056',
    position=600,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "data": {
    "uniprot_id": "<string>",
    "position": "<integer>",
    "scores": "<object>",
    "thresholds": "<object>"
  },
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.alpha_missense_get_residue_scores-parameters)

## 3) 查询单个错义变异的 AlphaMissense 分类 {#dnakit.databases.alpha_missense_get_variant_score}

- **作用：** 查询单个错义变异的 AlphaMissense 分类。上游不提供该替换的独立数值分数，pathogenicity_score 为 null，score_available 为 false。
- **处理方法：** 通过 `AlphaMissense` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.alpha_missense_get_variant_score(*, uniprot_id, variant, verbose=True)`
- **输入：** 必须提供 `uniprot_id`, `variant`。
- **输出：** 返回字典，包含 `status`, `data`, `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.alpha_missense_get_variant_score(
    uniprot_id='P15056',
    variant='V600E',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "data": {
    "uniprot_id": "<string>",
    "variant": "<string>",
    "position": "<integer>",
    "reference_aa": "<string>",
    "variant_aa": "<string>",
    "pathogenicity_score": "<number>",
    "score_available": "<boolean>",
    "classification": "<string>"
  },
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.alpha_missense_get_variant_score-parameters)

## 4) 查询基因组指定位点的全部替换 CADD 分数 {#dnakit.databases.cadd_get_position_scores}

- **作用：** 查询基因组指定位点的全部替换 CADD 分数。
- **处理方法：** 通过 `CADD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cadd_get_position_scores(*, chrom, pos, version=None, verbose=True)`
- **输入：** 必须提供 `chrom`, `pos`。
- **输出：** 返回字典，包含 `chrom`, `pos`, `variants`, `version` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cadd_get_position_scores(
    chrom='7',
    pos=140753336,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "chrom": "<string>",
  "pos": "<integer>",
  "variants": [
    "<value>"
  ],
  "version": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_position_scores-parameters)

## 5) 查询基因组区间内的 CADD 分数 {#dnakit.databases.cadd_get_range_scores}

- **作用：** 查询基因组区间内的 CADD 分数。
- **处理方法：** 通过 `CADD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cadd_get_range_scores(*, chrom, start, end, version=None, verbose=True)`
- **输入：** 必须提供 `chrom`, `start`, `end`。
- **输出：** 返回字典，包含 `chrom`, `start`, `end`, `variants` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cadd_get_range_scores(
    chrom='7',
    start=140753330,
    end=140753340,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "chrom": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "variants": [
    "<value>"
  ],
  "count": "<integer>",
  "version": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_range_scores-parameters)

## 6) 通过 MyVariant 查询已收录的 CADD 注释 {#dnakit.databases.cadd_get_scores}

- **作用：** 通过 MyVariant 查询已收录的 CADD 注释。
- **处理方法：** 通过 `API` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cadd_get_scores(*, variant, assembly='hg38', verbose=True)`
- **输入：** 必须提供 `variant`。
- **输出：** 返回 MyVariant 检索结果；字典结果另含 dnakit_score_source，表明分数来自已索引的 CADD 注释。缺失分数保持缺失，不执行新的模型预测。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cadd_get_scores(
    variant='chr7:g.140753336A>T',
    assembly='hg38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 MyVariant 检索结果；字典结果另含 dnakit_score_source，表明分数来自已索引的 CADD 注释。缺失分数保持缺失，不执行新的模型预测。

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_scores-parameters)

## 7) 查询单个变异的 CADD 有害性分数 {#dnakit.databases.cadd_get_variant_score}

- **作用：** 查询单个变异的 CADD 有害性分数。
- **处理方法：** 通过 `CADD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cadd_get_variant_score(*, chrom, pos, ref, alt, version=None, include_annotations=None, verbose=True)`
- **输入：** 必须提供 `chrom`, `pos`, `ref`, `alt`。
- **输出：** 返回字典，包含 `chrom`, `pos`, `ref`, `alt` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cadd_get_variant_score(
    chrom='7',
    pos=140753336,
    ref='A',
    alt='T',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "chrom": "<string>",
  "pos": "<integer>",
  "ref": "<string>",
  "alt": "<string>",
  "phred_score": "<number>",
  "raw_score": "<number>",
  "interpretation": "<string>",
  "version": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.cadd_get_variant_score-parameters)

## 8) 查询基因的 EVE 评分覆盖情况 {#dnakit.databases.eve_get_gene_info}

- **作用：** 查询基因的 EVE 评分覆盖情况。
- **处理方法：** 通过 `EVE` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.eve_get_gene_info(*, gene_symbol, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** 返回字典，包含 `status`, `data`, `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.eve_get_gene_info(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "data": {
    "gene_symbol": "<string>",
    "ensembl_id": "<string>",
    "eve_website": "<string>"
  },
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.eve_get_gene_info-parameters)

## 9) 查询变异的 EVE 致病性分数 {#dnakit.databases.eve_get_variant_score}

- **作用：** 查询变异的 EVE 致病性分数。
- **处理方法：** 通过 `EVE` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.eve_get_variant_score(*, variant=None, chrom=None, pos=None, ref=None, alt=None, species=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data`, `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.eve_get_variant_score(
    chrom='17',
    pos=7674220,
    ref='G',
    alt='A',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "data": {
    "variant": "<string>",
    "eve_scores": [
      "<value>"
    ],
    "threshold": "<object>"
  },
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.eve_get_variant_score-parameters)

## 10) 将深度突变扫描数据整理为氨基酸效应矩阵 {#dnakit.databases.mave_db_get_effect_matrix}

- **作用：** 将深度突变扫描数据整理为氨基酸效应矩阵。
- **处理方法：** 通过 `MaveDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mave_db_get_effect_matrix(*, urn, uniprot_id=None, score_field=None, verbose=True)`
- **输入：** 必须提供 `urn`。
- **输出：** 返回字典，包含 `urn`, `matrix`, `positions`, `amino_acid_order` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mave_db_get_effect_matrix(
    urn='urn:mavedb:00000115-a-7',
    uniprot_id='P01116',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "urn": "<string>",
  "matrix": [
    [
      "<number>"
    ]
  ],
  "positions": [
    "<integer>"
  ],
  "amino_acid_order": "<string>",
  "shape": [
    "<integer>"
  ],
  "n_parsed_single_missense": "<integer>",
  "n_dropped": "<integer>",
  "n_total_rows": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.mave_db_get_effect_matrix-parameters)

## 11) 汇总变异的多来源致病性预测分数 {#dnakit.databases.my_variant_get_pathogenicity_scores}

- **作用：** 汇总变异的多来源致病性预测分数。
- **处理方法：** 通过 `MyVariant` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.my_variant_get_pathogenicity_scores(*, variant_id, fields=None, assembly=None, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回 dbNSFP 已有变异评分，包括数据源提供的 REVEL、CADD、SIFT 等；缺失评分不会补算。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.my_variant_get_pathogenicity_scores(
    variant_id='rs45478192',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 dbNSFP 已有变异评分，包括数据源提供的 REVEL、CADD、SIFT 等；缺失评分不会补算。

[完整参数与使用说明](../../faq/api/36_variant_effect_scores.md#dnakit.databases.my_variant_get_pathogenicity_scores-parameters)

<!-- END GENERATED API: variant-scores -->
