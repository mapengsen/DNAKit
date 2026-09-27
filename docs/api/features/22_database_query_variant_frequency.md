# 变异记录与群体频率查询

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

## 1) dbSNP 变异信息 {#7-dbsnp}

- **作用：** 按 rsID 查询 dbSNP 等位基因及 GRCh37/GRCh38 坐标。
- **处理方法：** 按 rsID 查询 dbSNP，并保留等位基因与组装坐标。
- **API：** `dnakit.annotation.get_dbsnp_variant(rsid[必须])`
- **输入：** 带或不带 `rs` 前缀的 rsID。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import get_dbsnp_variant

result = get_dbsnp_variant("rs28934578")
print(result.data)
```

- **示例结果：** dbSNP 坐标、assembly 和等位基因记录。

## 2) dbSNP 等位基因频率 {#8-dbsnp}

- **作用：** 获取指定 rsID 的群体等位基因频率记录。
- **处理方法：** 读取 dbSNP 各研究或群体的等位基因计数与频率。
- **API：** `dnakit.annotation.get_dbsnp_frequencies(rsid[必须])`
- **输入：** 一个 rsID。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import get_dbsnp_frequencies

result = get_dbsnp_frequencies("rs28934578")
print(result.data)
```

- **示例结果：** 按研究或群体整理的 allele count/频率数据。

## 3) gnomAD 变异搜索 {#9-gnomad}

- **作用：** 使用 rsID 或变异文本查找 gnomAD 的标准 `chrom-pos-ref-alt` 标识。
- **处理方法：** 在指定 gnomAD 数据集中检索匹配的变异 ID。
- **API：** `dnakit.annotation.search_gnomad_variants(query[必须], dataset[可选])`
- **输入：** 查询文本；默认统一使用 `gnomad_r4`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import search_gnomad_variants

result = search_gnomad_variants("rs28934578")
print(result.data)
```

- **示例结果：** 匹配的 gnomAD variant ID 和基础信息。

## 4) gnomAD 变异汇总 {#10-gnomad}

- **作用：** 获取一个 gnomAD variant ID 的总体等位基因频率和元数据。
- **处理方法：** 按 variant ID 汇总 gnomAD 的 AC、AN 与 AF。
- **API：** `dnakit.annotation.get_gnomad_variant(variant_id[必须], dataset[可选])`
- **输入：** `chrom-pos-ref-alt` 标识；默认 `gnomad_r4`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import get_gnomad_variant

result = get_gnomad_variant("17-7675088-C-G")
print(result.data)
```

- **示例结果：** genome/exome 总体 allele count、allele number 和频率。

## 5) gnomAD 分群频率 {#11-gnomad}

- **作用：** 获取一个变异在不同 ancestry 和性别分组中的等位基因频率。
- **处理方法：** 从指定 gnomAD 记录提取 ancestry/sex 分层频率。
- **API：** `dnakit.annotation.get_gnomad_population_frequencies(variant_id[必须], dataset[可选])`
- **输入：** gnomAD variant ID；默认 `gnomad_r4`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import get_gnomad_population_frequencies

result = get_gnomad_population_frequencies("17-7675088-C-G")
print(result.data)
```

- **示例结果：** 按 ancestry/sex 分层的 AC、AN 和 AF。

<!-- BEGIN GENERATED API: variant-records -->

<span id="variant-records"></span>

## 6) 查询 dbSNP 等位基因频率 {#dnakit.databases.dbsnp_get_frequencies}

- **作用：** 查询 dbSNP 等位基因频率。
- **处理方法：** 通过 `dbSNPGetFrequencies` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dbsnp_get_frequencies(*, rsid, verbose=True)`
- **输入：** 必须提供 `rsid`。
- **输出：** 返回字典，包含 `status`, `data`, `rsid`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dbsnp_get_frequencies(
    rsid='rs429358',
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
    "frequencies": [
      "<object>"
    ],
    "total_studies": "<integer>"
  },
  "rsid": "<string>",
  "url": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.dbsnp_get_frequencies-parameters)

## 7) 按 rsID 查询 dbSNP 变异坐标与等位基因 {#dnakit.databases.dbsnp_get_variant_by_rsid}

- **作用：** 按 rsID 查询 dbSNP 变异坐标与等位基因。
- **处理方法：** 通过 `dbSNPGetVariantByRsID` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dbsnp_get_variant_by_rsid(*, rsid, verbose=True)`
- **输入：** 必须提供 `rsid`。
- **输出：** 返回字典，包含 `status`, `data`, `rsid`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dbsnp_get_variant_by_rsid(
    rsid='rs12345',
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
    "refsnp_id": "<string>",
    "snp_id": "<integer>",
    "chromosome": "<string>",
    "position": "<string>",
    "assembly": "<string>",
    "position_grch38": "<string>",
    "position_grch37": "<string>",
    "allele": "<string>"
  },
  "rsid": "<string>",
  "url": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.dbsnp_get_variant_by_rsid-parameters)

## 8) 查询 gnomAD 基因变异约束指标 {#dnakit.databases.gnomad_get_gene_constraints}

- **作用：** 查询 gnomAD 基因变异约束指标。
- **处理方法：** 通过 `gnomADGetGeneConstraints` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gnomad_get_gene_constraints(*, gene_symbol, reference_genome=None, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** 返回字典，包含 `gene`, `note` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gnomad_get_gene_constraints(
    gene_symbol='BRCA1',
    reference_genome='GRCh38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene": {
    "symbol": "<string>",
    "gene_id": "<string>",
    "exac_constraint": "<object>",
    "gnomad_constraint": {
      "exp_lof": "<number>",
      "obs_lof": "<integer>",
      "oe_lof": "<number>",
      "pLI": "<number>",
      "exp_mis": "<number>",
      "obs_mis": "<integer>",
      "oe_mis": "<number>",
      "exp_syn": "<number>"
    },
    "mitochondrial_constraint": "<object>"
  },
  "note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_get_gene_constraints-parameters)

## 9) 查询 gnomAD 变异基本信息 {#dnakit.databases.gnomad_get_variant}

- **作用：** 查询 gnomAD 变异基本信息。
- **处理方法：** 通过 `gnomADGraphQLQuery` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gnomad_get_variant(*, variant_id, dataset=None, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回字典，包含 `variant`, `dataset`, `reference_genome`, `dataset_note` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gnomad_get_variant(
    variant_id='19-44908822-C-T',
    dataset='gnomad_r3',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": {
    "variant_id": "<string>",
    "chrom": "<string>",
    "pos": "<integer>",
    "ref": "<string>",
    "alt": "<string>",
    "rsid": "<string>",
    "genome": {
      "af": "<number>",
      "an": "<integer>",
      "ac": "<integer>"
    },
    "exome": {
      "af": "<number>",
      "an": "<integer>",
      "ac": "<integer>"
    }
  },
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_get_variant-parameters)

## 10) 查询 gnomAD 变异的分群等位基因频率 {#dnakit.databases.gnomad_get_variant_populations}

- **作用：** 查询 gnomAD 变异的分群等位基因频率。
- **处理方法：** 通过 `gnomADGetVariantPopulations` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gnomad_get_variant_populations(*, variant_id, dataset=None, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回字典，包含 `variant_id`, `chrom`, `pos`, `ref` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gnomad_get_variant_populations(
    variant_id='1-55051215-G-GA',
    dataset='gnomad_r4',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant_id": "<string>",
  "chrom": "<string>",
  "pos": "<integer>",
  "ref": "<string>",
  "alt": "<string>",
  "rsid": "<string>",
  "dataset": "<string>",
  "genome": {
    "ac": "<integer>",
    "an": "<integer>",
    "af": "<number>",
    "populations": [
      "<object>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_get_variant_populations-parameters)

## 11) 检索 gnomAD 变异编号 {#dnakit.databases.gnomad_search_variants}

- **作用：** 检索 gnomAD 变异编号。
- **处理方法：** 通过 `gnomADGraphQLQuery` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gnomad_search_variants(*, query=None, dataset=None, variant=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `variant_search`, `dataset`, `reference_genome`, `dataset_note` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gnomad_search_variants(
    query='rs7412',
    dataset='gnomad_r3',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant_search": [
    {
      "variant_id": "<string>"
    }
  ],
  "dataset": "<string>",
  "reference_genome": "<string>",
  "dataset_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_variant_frequency.md#dnakit.databases.gnomad_search_variants-parameters)

<!-- END GENERATED API: variant-records -->
