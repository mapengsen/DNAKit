# 疾病、表型与本体查询

查询疾病、表型及相关本体条目。

<!-- BEGIN GENERATED API: phenotypes -->

<span id="phenotypes"></span>

## 1) 查询疾病相关基因及关联评分 {#dnakit.databases.dis_ge_net_get_disease_genes}

- **作用：** 查询疾病相关基因及关联评分。
- **处理方法：** 通过 `DisGeNET` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dis_ge_net_get_disease_genes(*, operation=None, disease, min_score=None, limit=None, verbose=True)`
- **输入：** 必须提供 `disease`。 需配置：`DISGENET_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dis_ge_net_get_disease_genes(
    disease='C0006142',
    min_score=0.3,
    limit=30,
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
    "disease": "<string>",
    "genes": [
      "<object>"
    ],
    "gene_count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_get_disease_genes-parameters)

## 2) 按来源与评分筛选基因疾病关联 {#dnakit.databases.dis_ge_net_get_gda}

- **作用：** 按来源与评分筛选基因疾病关联。
- **处理方法：** 通过 `DisGeNET` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dis_ge_net_get_gda(*, operation=None, gene=None, disease=None, source=None, min_score=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需配置：`DISGENET_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dis_ge_net_get_gda(
    gene='TP53',
    min_score=0.3,
    limit=20,
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
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_get_gda-parameters)

## 3) 查询变异与疾病的关联证据 {#dnakit.databases.dis_ge_net_get_vda}

- **作用：** 查询变异与疾病的关联证据。
- **处理方法：** 通过 `DisGeNET` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dis_ge_net_get_vda(*, operation=None, variant=None, gene=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需配置：`DISGENET_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dis_ge_net_get_vda(
    gene='BRCA1',
    limit=20,
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
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_get_vda-parameters)

## 4) 按疾病概念编号检索相关基因 {#dnakit.databases.dis_ge_net_search_disease}

- **作用：** 按疾病概念编号检索相关基因。
- **处理方法：** 通过 `DisGeNET` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dis_ge_net_search_disease(*, operation=None, disease, limit=None, verbose=True)`
- **输入：** 必须提供 `disease`。 需配置：`DISGENET_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dis_ge_net_search_disease(
    disease='C0006142',
    limit=10,
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
    "disease": "<string>",
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_search_disease-parameters)

## 5) 按基因名称检索关联疾病 {#dnakit.databases.dis_ge_net_search_gene}

- **作用：** 按基因名称检索关联疾病。
- **处理方法：** 通过 `DisGeNET` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dis_ge_net_search_gene(*, operation=None, gene, limit=None, verbose=True)`
- **输入：** 必须提供 `gene`。 需配置：`DISGENET_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dis_ge_net_search_gene(
    gene='BRCA1',
    limit=10,
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
    "gene": "<string>",
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_search_gene-parameters)

## 6) 统计 GenCC 基因疾病有效性分类 {#dnakit.databases.gen_cc_get_classifications}

- **作用：** 统计 GenCC 基因疾病有效性分类。
- **处理方法：** 通过 `GenCC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gen_cc_get_classifications(*, operation=None, gene_symbol=None, disease=None, submitter=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gen_cc_get_classifications(
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
    "classifications": [
      "<object>"
    ],
    "total_submissions": "<integer>",
    "unique_genes": "<integer>",
    "unique_diseases": "<integer>",
    "top_submitters": [
      "<object>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.gen_cc_get_classifications-parameters)

## 7) 按疾病检索 GenCC 基因关联证据 {#dnakit.databases.gen_cc_search_disease}

- **作用：** 按疾病检索 GenCC 基因关联证据。
- **处理方法：** 通过 `GenCC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gen_cc_search_disease(*, operation=None, disease=None, disease_name=None, classification=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gen_cc_search_disease(
    disease='Marfan syndrome',
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
    "disease": "<string>",
    "submissions": [
      "<object>"
    ],
    "submission_count": "<integer>",
    "unique_genes": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.gen_cc_search_disease-parameters)

## 8) 按基因查询 GenCC 疾病关联分类 {#dnakit.databases.gen_cc_search_gene}

- **作用：** 按基因查询 GenCC 疾病关联分类。
- **处理方法：** 通过 `GenCC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gen_cc_search_gene(*, operation=None, gene_symbol, classification=None, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gen_cc_search_gene(
    gene_symbol='BRCA2',
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
    "submissions": [
      "<object>"
    ],
    "submission_count": "<integer>",
    "unique_diseases": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.gen_cc_search_gene-parameters)

## 9) 按多个 HPO 表型检索关联疾病 {#dnakit.databases.get_joint_associated_diseases_by_hpo_id_list}

- **作用：** 按多个 HPO 表型检索关联疾病。
- **处理方法：** 通过 `MonarchDiseasesForMultiplePheno` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.get_joint_associated_diseases_by_hpo_id_list(*, HPO_ID_list, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `HPO_ID_list`。
- **输出：** 返回“按多个 HPO 表型检索关联疾病”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.get_joint_associated_diseases_by_hpo_id_list(
    HPO_ID_list=['HP:0001250'],
    limit=5,
    offset=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<string>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.get_joint_associated_diseases_by_hpo_id_list-parameters)

## 10) 查询 HPO 表型对应的疾病 {#dnakit.databases.hpo_get_diseases_by_phenotype}

- **作用：** 查询 HPO 表型对应的疾病。
- **处理方法：** 通过 `HPO` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpo_get_diseases_by_phenotype(*, term_id=None, limit=None, offset=None, hp_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `diseases` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpo_get_diseases_by_phenotype(
    term_id='HP:0001250',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "diseases": [
    {
      "id": "<string>",
      "name": "<string>",
      "mondo_id": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_get_diseases_by_phenotype-parameters)

## 11) 查询 HPO 表型对应的基因 {#dnakit.databases.hpo_get_genes_by_phenotype}

- **作用：** 查询 HPO 表型对应的基因。
- **处理方法：** 通过 `HPO` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpo_get_genes_by_phenotype(*, term_id, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `term_id`。
- **输出：** 返回字典，包含 `genes` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpo_get_genes_by_phenotype(
    term_id='HP:0001250',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "genes": [
    {
      "id": "<string>",
      "name": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_get_genes_by_phenotype-parameters)

## 12) 查询 HPO 术语定义与同义词 {#dnakit.databases.hpo_get_term}

- **作用：** 查询 HPO 术语定义与同义词。
- **处理方法：** 通过 `HPO` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpo_get_term(*, term_id, verbose=True)`
- **输入：** 必须提供 `term_id`。
- **输出：** 返回字典，包含 `id`, `name`, `definition`, `comment` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpo_get_term(
    term_id='HP:0001250',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "id": "<string>",
  "name": "<string>",
  "definition": "<string>",
  "comment": "<string>",
  "synonyms": [
    "<string>"
  ],
  "descendant_count": "<integer>",
  "xrefs": [
    "<string>"
  ],
  "translations": [
    "<object>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_get_term-parameters)

## 13) 按关键词检索 HPO 表型术语 {#dnakit.databases.hpo_search_terms}

- **作用：** 按关键词检索 HPO 表型术语。
- **处理方法：** 通过 `HPO` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpo_search_terms(*, query, max_results=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `id`, `name`, `definition`, `descendant_count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpo_search_terms(
    query='seizure',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "id": "<string>",
    "name": "<string>",
    "definition": "<string>",
    "descendant_count": "<integer>",
    "synonyms": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_search_terms-parameters)

## 14) 汇总人类基因标识与注释 {#dnakit.databases.marrvel_get_gene}

- **作用：** 汇总人类基因标识与注释。
- **处理方法：** 通过 `MARRVELGene` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.marrvel_get_gene(*, symbol=None, gene_symbol=None, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组至少选一组：`symbol` / `gene_symbol` / `gene`。
- **输出：** 返回基因标识交叉引用、名称、位置、别名与摘要；不包含 OMIM 编号。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.marrvel_get_gene(
    symbol='CFTR',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回基因标识交叉引用、名称、位置、别名与摘要；不包含 OMIM 编号。

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.marrvel_get_gene-parameters)

## 15) 查询基因的 OMIM 表型与遗传方式 {#dnakit.databases.marrvel_get_omim_phenotypes}

- **作用：** 查询基因的 OMIM 表型与遗传方式。
- **处理方法：** 通过 `MARRVELOmim` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.marrvel_get_omim_phenotypes(*, symbol=None, gene_symbol=None, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组至少选一组：`symbol` / `gene_symbol` / `gene`。
- **输出：** 返回“查询基因的 OMIM 表型与遗传方式”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.marrvel_get_omim_phenotypes(
    symbol='CFTR',
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

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.marrvel_get_omim_phenotypes-parameters)

## 16) 查询 Monarch 基因疾病关联 {#dnakit.databases.monarch_get_gene_diseases}

- **作用：** 查询 Monarch 基因疾病关联。
- **处理方法：** 通过 `Monarch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_get_gene_diseases(*, subject, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `subject`。
- **输出：** 返回字典，包含 `items`, `total` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_get_gene_diseases(
    subject='HGNC:11998',
    limit=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "items": [
    "<value>"
  ],
  "total": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_get_gene_diseases-parameters)

## 17) 查询 Monarch 跨物种基因表型关联 {#dnakit.databases.monarch_get_gene_phenotypes}

- **作用：** 查询 Monarch 跨物种基因表型关联。
- **处理方法：** 通过 `Monarch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_get_gene_phenotypes(*, subject, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `subject`。
- **输出：** 返回字典，包含 `items`, `total` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_get_gene_phenotypes(
    subject='HGNC:11998',
    limit=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "items": [
    "<value>"
  ],
  "total": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_get_gene_phenotypes-parameters)

## 18) 按名称检索 Monarch 生物医学实体 {#dnakit.databases.monarch_search}

- **作用：** 按名称检索 Monarch 生物医学实体。
- **处理方法：** 通过 `MonarchV3` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_search(*, query, category=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `id`, `name`, `category`, `symbol` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_search(
    query='Alzheimer',
    category='biolink:Disease',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "id": "<string>",
    "name": "<string>",
    "category": "<string>",
    "symbol": "<string>",
    "description": "<string>",
    "taxon": "<string>",
    "taxon_label": "<string>",
    "deprecated": "<boolean>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_search-parameters)

## 19) 检索 Monarch 基因与跨物种标识 {#dnakit.databases.monarch_search_gene}

- **作用：** 检索 Monarch 基因与跨物种标识。
- **处理方法：** 通过 `Monarch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_search_gene(*, query, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `items`, `total` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_search_gene(
    query='TP53',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "items": [
    "<value>"
  ],
  "total": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_search_gene-parameters)

## 20) 查询知识图谱中的跨物种实体关联 {#dnakit.databases.monarch_v3_get_associations}

- **作用：** 查询知识图谱中的跨物种实体关联。
- **处理方法：** 通过 `MonarchV3` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_v3_get_associations(*, subject=None, object=None, category, limit=None, verbose=True)`
- **输入：** 必须提供 `category`。
- **输出：** 返回“查询知识图谱中的跨物种实体关联”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_associations(
    subject='HGNC:11998',
    category='biolink:GeneToPhenotypicFeatureAssociation',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<value>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_associations-parameters)

## 21) 按 CURIE 查询生物医学实体详情 {#dnakit.databases.monarch_v3_get_entity}

- **作用：** 按 CURIE 查询生物医学实体详情。
- **处理方法：** 通过 `MonarchV3` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_v3_get_entity(*, entity_id, verbose=True)`
- **输入：** 必须提供 `entity_id`。
- **输出：** data 包含实体名称、类别、描述、同义词和交叉引用。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_entity(
    entity_id='HGNC:11998',
    verbose=False,
)
print(result)
```

- **示例结果：** data 包含实体名称、类别、描述、同义词和交叉引用。

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_entity-parameters)

## 22) 按身体系统汇总疾病表型分布 {#dnakit.databases.monarch_v3_get_histopheno}

- **作用：** 按身体系统汇总疾病表型分布。
- **处理方法：** 通过 `MonarchV3` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_v3_get_histopheno(*, entity_id, verbose=True)`
- **输入：** 必须提供 `entity_id`。
- **输出：** 返回“按身体系统汇总疾病表型分布”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_histopheno(
    entity_id='MONDO:0007947',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<value>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_histopheno-parameters)

## 23) 查询生物医学实体的跨本体标识映射 {#dnakit.databases.monarch_v3_get_mappings}

- **作用：** 查询生物医学实体的跨本体标识映射。
- **处理方法：** 通过 `MonarchV3` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_v3_get_mappings(*, entity_id, limit=None, verbose=True)`
- **输入：** 必须提供 `entity_id`。
- **输出：** 返回“查询生物医学实体的跨本体标识映射”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_mappings(
    entity_id='MONDO:0005148',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<value>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_mappings-parameters)

## 24) 按表型语义相似性检索疾病 {#dnakit.databases.monarch_v3_phenotype_similarity_search}

- **作用：** 按表型语义相似性检索疾病。
- **处理方法：** 通过 `MonarchV3` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_v3_phenotype_similarity_search(*, phenotypes, group=None, limit=None, verbose=True)`
- **输入：** 必须提供 `phenotypes`。
- **输出：** 返回“按表型语义相似性检索疾病”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_v3_phenotype_similarity_search(
    phenotypes=['HP:0001250', 'HP:0001249', 'HP:0002353'],
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<value>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_phenotype_similarity_search-parameters)

## 25) 检索 Monarch 知识图谱实体与类别 {#dnakit.databases.monarch_v3_search}

- **作用：** 检索 Monarch 知识图谱实体与类别。
- **处理方法：** 通过 `MonarchV3` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.monarch_v3_search(*, query=None, limit=None, category=None, q=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回“检索 Monarch 知识图谱实体与类别”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.monarch_v3_search(
    query='BRCA1',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<value>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_search-parameters)

## 26) 查询本体术语的上位概念 {#dnakit.databases.ols_get_term_ancestors}

- **作用：** 查询本体术语的上位概念。
- **处理方法：** 通过 `OLS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ols_get_term_ancestors(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `page`, `numElements`, `totalPages`, `totalElements` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ols_get_term_ancestors(
    term_iri='http://www.ebi.ac.uk/efo/EFO_0000384',
    ontology='efo',
    include_obsolete=False,
    size=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "page": "<integer>",
  "numElements": "<integer>",
  "totalPages": "<integer>",
  "totalElements": "<integer>",
  "elements": [
    "<value>"
  ],
  "term_iri": "<string>",
  "ontology": "<string>",
  "filters": {
    "include_obsolete": "<boolean>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_get_term_ancestors-parameters)

## 27) 查询本体术语的下位概念 {#dnakit.databases.ols_get_term_children}

- **作用：** 查询本体术语的下位概念。
- **处理方法：** 通过 `OLS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ols_get_term_children(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `page`, `numElements`, `totalPages`, `totalElements` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ols_get_term_children(
    term_iri='http://www.ebi.ac.uk/efo/EFO_0000384',
    ontology='efo',
    include_obsolete=False,
    size=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "page": "<integer>",
  "numElements": "<integer>",
  "totalPages": "<integer>",
  "totalElements": "<integer>",
  "elements": [
    "<value>"
  ],
  "term_iri": "<string>",
  "ontology": "<string>",
  "filters": {
    "include_obsolete": "<boolean>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_get_term_children-parameters)

## 28) 查询本体术语详情 {#dnakit.databases.ols_get_term_info}

- **作用：** 查询本体术语详情。
- **处理方法：** 通过 `OLS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ols_get_term_info(*, operation=None, id=None, term_id=None, term_iri=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `iri`, `ontologyName`, `shortForm`, `label` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ols_get_term_info(
    id='HP:0001903',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "iri": "<string>",
  "ontologyName": "<string>",
  "shortForm": "<string>",
  "label": "<string>",
  "oboId": "<string>",
  "isObsolete": "<boolean>",
  "description": [
    "<value>"
  ],
  "synonyms": [
    "<value>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_get_term_info-parameters)

## 29) 按关键词检索本体术语 {#dnakit.databases.ols_search_terms}

- **作用：** 按关键词检索本体术语。
- **处理方法：** 通过 `OLS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ols_search_terms(*, operation=None, query, rows=None, size=None, ontology=None, exact_match=None, include_obsolete=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `response`, `responseHeader`, `facet_counts`, `query` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ols_search_terms(
    query='gene:TP53',
    rows=10,
    ontology='efo',
    exact_match=False,
    include_obsolete=False,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "response": {
    "docs": [
      "<value>"
    ],
    "numFound": "<integer>",
    "start": "<integer>"
  },
  "responseHeader": {
    "QTime": "<integer>",
    "status": "<integer>"
  },
  "facet_counts": {
    "facet_fields": {
      "ontologyPreferredPrefix": "<array>",
      "isDefiningOntology": "<array>",
      "ontologyId": "<array>",
      "ontologyIri": "<array>",
      "isObsolete": "<array>",
      "type": "<array>"
    }
  },
  "query": "<string>",
  "filters": {
    "ontology": "<string>",
    "exact_match": "<boolean>",
    "exact_match_fields": "<string>",
    "include_obsolete": "<boolean>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_search_terms-parameters)

## 30) 获取 OMIM 疾病的临床表型摘要 {#dnakit.databases.omim_get_clinical_synopsis}

- **作用：** 获取 OMIM 疾病的临床表型摘要。
- **处理方法：** 通过 `OMIM` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.omim_get_clinical_synopsis(*, operation=None, mim_number, verbose=True)`
- **输入：** 必须提供 `mim_number`。 需配置：`OMIM_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.omim_get_clinical_synopsis(
    mim_number='219700',
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
    "mimNumber": "<integer>",
    "inheritance": "<string>",
    "neurologicCentralNervousSystem": [
      "<value>"
    ],
    "cardiovascular": [
      "<value>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_get_clinical_synopsis-parameters)

## 31) 按 MIM 编号查询 OMIM 条目 {#dnakit.databases.omim_get_entry}

- **作用：** 按 MIM 编号查询 OMIM 条目。
- **处理方法：** 通过 `OMIM` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.omim_get_entry(*, operation=None, mim_number, include=None, verbose=True)`
- **输入：** 必须提供 `mim_number`。 需配置：`OMIM_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.omim_get_entry(
    mim_number='164730',
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
    "mimNumber": "<integer>",
    "titles": "<object>",
    "textSectionList": [
      "<value>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_get_entry-parameters)

## 32) 查询 OMIM 基因疾病与染色体映射 {#dnakit.databases.omim_get_gene_map}

- **作用：** 查询 OMIM 基因疾病与染色体映射。
- **处理方法：** 通过 `OMIM` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.omim_get_gene_map(*, operation=None, mim_number=None, chromosome=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需配置：`OMIM_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.omim_get_gene_map(
    mim_number='164730',
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
    "gene_maps": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_get_gene_map-parameters)

## 33) 检索孟德尔遗传疾病与相关基因 {#dnakit.databases.omim_search}

- **作用：** 检索孟德尔遗传疾病与相关基因。
- **处理方法：** 通过 `OMIM` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.omim_search(*, operation=None, query, limit=None, start=None, verbose=True)`
- **输入：** 必须提供 `query`。 需配置：`OMIM_API_KEY`.
- **输出：** 返回字典，包含 `status`, `data`, `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.omim_search(
    query='BRCA1',
    limit=5,
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
    "total_results": "<integer>",
    "entries": [
      "<value>"
    ]
  },
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_search-parameters)

## 34) 查询罕见病分类层级 {#dnakit.databases.orphanet_get_classification}

- **作用：** 查询罕见病分类层级。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_classification(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_classification(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "parent": "<object>",
    "children": [
      "<value>"
    ],
    "member_of_classifications": "<object>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_classification-parameters)

## 35) 按 ORPHA 编号查询罕见病详情 {#dnakit.databases.orphanet_get_disease}

- **作用：** 按 ORPHA 编号查询罕见病详情。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_disease(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_disease(
    orpha_code='558',
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
    "ORPHAcode": "<string>",
    "Preferred term": "<string>",
    "Definition": "<string>",
    "Synonyms": [
      "<value>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_disease-parameters)

## 36) 按 ORPHA 编号获取已报告的基因关联 {#dnakit.databases.orphanet_get_disease_genes}

- **作用：** 按 ORPHA 编号获取已报告的基因关联。
- **处理方法：** 通过 `API` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_disease_genes(*, orpha_code, verbose=True)`
- **输入：** 必须提供 `orpha_code`。
- **输出：** 返回 Orphanet_get_genes 的原始关联记录结构，见 [orphanet_get_genes](22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes)。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_disease_genes(
    orpha_code='558',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 Orphanet_get_genes 的原始关联记录结构，见 [orphanet_get_genes](22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes)。

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_disease_genes-parameters)

## 37) 查询罕见病患病率与发病率 {#dnakit.databases.orphanet_get_epidemiology}

- **作用：** 查询罕见病患病率与发病率。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_epidemiology(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_epidemiology(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "preferred_term": "<string>",
    "prevalences": [
      "<object>"
    ],
    "prevalence_count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_epidemiology-parameters)

## 38) 按基因检索关联罕见病 {#dnakit.databases.orphanet_get_gene_diseases}

- **作用：** 按基因检索关联罕见病。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_gene_diseases(*, operation=None, gene_name=None, gene_symbol=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_gene_diseases(
    gene_name='FBN1',
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
    "gene_name": "<string>",
    "diseases": [
      "<object>"
    ],
    "disease_count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_gene_diseases-parameters)

## 39) 查询罕见病关联基因与关联类型 {#dnakit.databases.orphanet_get_genes}

- **作用：** 查询罕见病关联基因与关联类型。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_genes(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_genes(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "genes": [
      "<object>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes-parameters)

## 40) 映射罕见病与其他医学编码 {#dnakit.databases.orphanet_get_icd_mapping}

- **作用：** 映射罕见病与其他医学编码。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_icd_mapping(*, operation=None, orpha_code=None, coding_system=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_icd_mapping(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "mappings": "<object>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_icd_mapping-parameters)

## 41) 查询罕见病起病年龄与遗传方式 {#dnakit.databases.orphanet_get_natural_history}

- **作用：** 查询罕见病起病年龄与遗传方式。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_natural_history(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_natural_history(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "preferred_term": "<string>",
    "average_age_of_onset": [
      "<string>"
    ],
    "type_of_inheritance": [
      "<string>"
    ],
    "disorder_group": "<string>",
    "typology": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_natural_history-parameters)

## 42) 查询罕见病表型及出现频率 {#dnakit.databases.orphanet_get_phenotypes}

- **作用：** 查询罕见病表型及出现频率。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_get_phenotypes(*, operation=None, orpha_code=None, orpha_id=None, orphacode=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_get_phenotypes(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "preferred_term": "<string>",
    "phenotypes": [
      "<object>"
    ],
    "phenotype_count": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_phenotypes-parameters)

## 43) 按完整或部分名称检索罕见病 {#dnakit.databases.orphanet_search_by_name}

- **作用：** 按完整或部分名称检索罕见病。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_search_by_name(*, operation=None, name=None, exact=None, lang=None, limit=None, query=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_search_by_name(
    name='Marfan syndrome',
    exact=True,
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
    "results": [
      "<value>"
    ],
    "count": "<integer>",
    "exact_match": "<boolean>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_search_by_name-parameters)

## 44) 按关键词检索罕见病 {#dnakit.databases.orphanet_search_diseases}

- **作用：** 按关键词检索罕见病。
- **处理方法：** 通过 `Orphanet` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.orphanet_search_diseases(*, operation=None, query, limit=None, lang=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.orphanet_search_diseases(
    query='Marfan syndrome',
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
    "results": [
      "<value>"
    ],
    "count": "<integer>",
    "total_count": "<integer>",
    "query": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_search_diseases-parameters)

## 45) 将疾病名称解析为 EFO 术语编号 {#dnakit.databases.osl_get_efo_id_by_disease_name}

- **作用：** 将疾病名称解析为 EFO 术语编号。
- **处理方法：** 通过 `EFO` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.osl_get_efo_id_by_disease_name(*, disease, verbose=True)`
- **输入：** 必须提供 `disease`。
- **输出：** 返回字典，包含 `efo_id`, `name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.osl_get_efo_id_by_disease_name(
    disease='Diabetes Mellitus',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "efo_id": "<string>",
  "name": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.osl_get_efo_id_by_disease_name-parameters)

<!-- END GENERATED API: phenotypes -->
