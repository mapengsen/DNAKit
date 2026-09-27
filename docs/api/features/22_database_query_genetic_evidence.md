# 遗传病关联、剂量与可干预性

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: genetic-evidence -->

<span id="genetic-evidence"></span>

## 1) 按基因查询 ClinGen 剂量敏感性审定结果 {#dnakit.databases.clin_gen_dosage_by_gene}

- **作用：** 按基因查询 ClinGen 剂量敏感性审定结果。
- **处理方法：** 通过 `ClinGenDosage` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_dosage_by_gene(*, gene=None, gene_symbol=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `symbol`, `hgnc_id`, `location`, `grch37` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_dosage_by_gene(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "symbol": "<string>",
    "hgnc_id": "<string>",
    "location": "<string>",
    "grch37": "<string>",
    "grch38": "<string>",
    "haplo_assertion": "<integer>",
    "triplo_assertion": "<integer>",
    "pli": "<number>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_dosage_by_gene-parameters)

## 2) 按区间检索 ClinGen 剂量敏感性记录 {#dnakit.databases.clin_gen_dosage_region_search}

- **作用：** 按区间检索 ClinGen 剂量敏感性记录。
- **处理方法：** 通过 `ClinGenDosage` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_dosage_region_search(*, chromosome, start, end, assembly=None, verbose=True)`
- **输入：** 必须提供 `chromosome`, `start`, `end`。
- **输出：** 返回记录列表；每条记录包含 `symbol`, `hgnc_id`, `type`, `location` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_dosage_region_search(
    chromosome='17',
    start=43044295,
    end=43125370,
    assembly='GRCh38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "symbol": "<string>",
    "hgnc_id": "<string>",
    "type": "<integer>",
    "location": "<string>",
    "grch37": "<string>",
    "grch38": "<string>",
    "haplo_assertion": "<integer>",
    "triplo_assertion": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_dosage_region_search-parameters)

## 3) 查询成人遗传发现的临床可干预性证据 {#dnakit.databases.clin_gen_get_actionability_adult}

- **作用：** 查询成人遗传发现的临床可干预性证据。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_get_actionability_adult(*, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `docId`, `topicIri`, `curationType`, `latestSearchDate` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_get_actionability_adult(
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "docId": "<string>",
    "topicIri": "<string>",
    "curationType": "<string>",
    "latestSearchDate": "<string>",
    "lastUpdated": "<string>",
    "lastAuthor": "<string>",
    "context": "<string>",
    "contextIri": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_actionability_adult-parameters)

## 4) 查询儿童遗传发现的临床可干预性证据 {#dnakit.databases.clin_gen_get_actionability_pediatric}

- **作用：** 查询儿童遗传发现的临床可干预性证据。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_get_actionability_pediatric(*, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `docId`, `topicIri`, `curationType`, `latestSearchDate` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_get_actionability_pediatric(
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "docId": "<string>",
    "topicIri": "<string>",
    "curationType": "<string>",
    "latestSearchDate": "<string>",
    "lastUpdated": "<string>",
    "lastAuthor": "<string>",
    "context": "<string>",
    "contextIri": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_actionability_pediatric-parameters)

## 5) 获取 ClinGen 剂量敏感性审定记录 {#dnakit.databases.clin_gen_get_dosage_sensitivity}

- **作用：** 获取 ClinGen 剂量敏感性审定记录。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_get_dosage_sensitivity(*, gene=None, include_regions=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回“获取 ClinGen 剂量敏感性审定记录”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_get_dosage_sensitivity(
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<object>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_dosage_sensitivity-parameters)

## 6) 获取 ClinGen 基因与疾病关联有效性审定记录 {#dnakit.databases.clin_gen_get_gene_validity}

- **作用：** 获取 ClinGen 基因与疾病关联有效性审定记录。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_get_gene_validity(*, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回“获取 ClinGen 基因与疾病关联有效性审定记录”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_get_gene_validity(
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<object>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_gene_validity-parameters)

## 7) 查询 ClinGen 变异致病性审定结果 {#dnakit.databases.clin_gen_get_variant_classifications}

- **作用：** 查询 ClinGen 变异致病性审定结果。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_get_variant_classifications(*, gene=None, variant=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `Variation`, `ClinVar Variation Id`, `HGNC Gene Symbol`, `Disease` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_get_variant_classifications(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "Variation": "<string>",
    "ClinVar Variation Id": "<string>",
    "HGNC Gene Symbol": "<string>",
    "Disease": "<string>",
    "Mondo Id": "<string>",
    "Assertion": "<string>",
    "Expert Panel": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_get_variant_classifications-parameters)

## 8) 按基因检索成人与儿童临床可干预性 {#dnakit.databases.clin_gen_search_actionability}

- **作用：** 按基因检索成人与儿童临床可干预性。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_search_actionability(*, gene, verbose=True)`
- **输入：** 必须提供 `gene`。
- **输出：** 返回字典，包含 `data`, `gene_searched`, `adult_count`, `pediatric_count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_search_actionability(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "Adult": [
      "<value>"
    ],
    "Pediatric": [
      "<value>"
    ]
  },
  "gene_searched": "<string>",
  "adult_count": "<integer>",
  "pediatric_count": "<integer>",
  "failed_contexts": "<object>",
  "note": "<string>",
  "source": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_search_actionability-parameters)

## 9) 按基因名称检索剂量敏感性评分 {#dnakit.databases.clin_gen_search_dosage_sensitivity}

- **作用：** 按基因名称检索剂量敏感性评分。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_search_dosage_sensitivity(*, gene, verbose=True)`
- **输入：** 必须提供 `gene`。
- **输出：** 返回记录列表；每条记录包含 `GENE SYMBOL`, `HGNC ID`, `HAPLOINSUFFICIENCY`, `TRIPLOSENSITIVITY` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_search_dosage_sensitivity(
    gene='MECP2',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "GENE SYMBOL": "<string>",
    "HGNC ID": "<string>",
    "HAPLOINSUFFICIENCY": "<string>",
    "TRIPLOSENSITIVITY": "<string>",
    "ONLINE REPORT": "<string>",
    "DATE": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_search_dosage_sensitivity-parameters)

## 10) 按基因名称检索疾病关联有效性 {#dnakit.databases.clin_gen_search_gene_validity}

- **作用：** 按基因名称检索疾病关联有效性。
- **处理方法：** 通过 `ClinGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_gen_search_gene_validity(*, gene, verbose=True)`
- **输入：** 必须提供 `gene`。
- **输出：** 返回记录列表；每条记录包含 `GENE SYMBOL`, `GENE ID (HGNC)`, `DISEASE LABEL`, `DISEASE ID (MONDO)` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_gen_search_gene_validity(
    gene='BRCA1',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "GENE SYMBOL": "<string>",
    "GENE ID (HGNC)": "<string>",
    "DISEASE LABEL": "<string>",
    "DISEASE ID (MONDO)": "<string>",
    "CLASSIFICATION": "<string>",
    "MOI": "<string>",
    "SOP": "<string>",
    "GCEP": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_genetic_evidence.md#dnakit.databases.clin_gen_search_gene_validity-parameters)

<!-- END GENERATED API: genetic-evidence -->
