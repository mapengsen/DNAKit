# ClinVar 临床意义查询

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

## 1) ClinVar 变异搜索 {#4-clinvar}

- **作用：** 按基因、疾病、ClinVar ID 或 HGVS/蛋白变化搜索 ClinVar。
- **处理方法：** 组合查询条件检索 ClinVar，并限制返回条数。
- **API：** `dnakit.annotation.search_clinvar_variants(gene[可选], condition[可选], variant_id[可选], variant_name[可选], clinical_significance[可选], limit[可选])`
- **输入：** 至少提供一个查询字段；`limit` 为 1–100。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import search_clinvar_variants

result = search_clinvar_variants(gene="TP53", variant_name="c.215C>G", limit=10)
print(result.data)
```

- **示例结果：** ClinVar 变异标识和摘要记录。

## 2) ClinVar 变异详情 {#5-clinvar}

- **作用：** 按 ClinVar variation ID 获取 accession、基因、坐标和审核状态。
- **处理方法：** 按 variation ID 查询 ClinVar 详情记录。
- **API：** `dnakit.annotation.get_clinvar_variant(variant_id[必须])`
- **输入：** ClinVar variation ID，例如 `12345`。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import get_clinvar_variant

result = get_clinvar_variant("12345")
print(result.data)
```

- **示例结果：** ClinVar 变异详情。

## 3) ClinVar 临床意义 {#6-clinvar}

- **作用：** 获取指定 ClinVar variation ID 的临床意义和解释记录。
- **处理方法：** 读取指定 ClinVar 记录中的临床意义与解释。
- **API：** `dnakit.annotation.get_clinvar_significance(variant_id[必须])`
- **输入：** ClinVar variation ID。
- **输出：** 返回 `ProviderResult`，查询内容在 `data`，并附有 `provider` 和来源信息。
- **示例代码：**

```python
from dnakit.annotation import get_clinvar_significance

result = get_clinvar_significance("12345")
print(result.data)
```

- **示例结果：** pathogenicity 分类和临床解释。

<!-- BEGIN GENERATED API: clinical-variants -->

<span id="clinical-variants"></span>

## 4) 查询 ClinVar 变异临床意义 {#dnakit.databases.clin_var_get_clinical_significance}

- **作用：** 查询 ClinVar 变异临床意义。
- **处理方法：** 通过 `ClinVarGetClinicalSignificance` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_var_get_clinical_significance(*, variant_id, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回字典，包含 `status`, `data`, `variant_id`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_var_get_clinical_significance(
    variant_id='12345',
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
    "variant_id": "<string>",
    "germline_classification": {
      "description": "<string>",
      "review_status": "<string>",
      "last_evaluated": "<string>",
      "fda_recognized": "<string>",
      "traits": "<array>"
    },
    "clinical_impact": {
      "description": "<string>",
      "review_status": "<string>",
      "last_evaluated": "<string>"
    },
    "oncogenicity": {
      "description": "<string>",
      "review_status": "<string>",
      "last_evaluated": "<string>"
    },
    "raw_data": "<object>"
  },
  "variant_id": "<string>",
  "url": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_clinical_variants.md#dnakit.databases.clin_var_get_clinical_significance-parameters)

## 5) 查询 ClinVar 变异摘要与审查状态 {#dnakit.databases.clin_var_get_variant_details}

- **作用：** 查询 ClinVar 变异摘要与审查状态。
- **处理方法：** 通过 `ClinVarGetVariantDetails` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_var_get_variant_details(*, variant_id, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回字典，包含 `status`, `data`, `variant_id`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_var_get_variant_details(
    variant_id='12345',
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
    "variant_id": "<string>",
    "accession": "<string>",
    "obj_type": "<string>",
    "chromosome": "<string>",
    "location": "<string>",
    "variation_name": "<string>",
    "title": "<string>",
    "genes": [
      "<string>"
    ]
  },
  "variant_id": "<string>",
  "url": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_clinical_variants.md#dnakit.databases.clin_var_get_variant_details-parameters)

## 6) 按基因、疾病或编号检索 ClinVar 变异 {#dnakit.databases.clin_var_search_variants}

- **作用：** 按基因、疾病或编号检索 ClinVar 变异。
- **处理方法：** 通过 `ClinVarSearchVariants` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.clin_var_search_variants(*, gene=None, condition=None, variant_id=None, variant_name=None, max_results=None, limit=None, clinical_significance=None, gene_symbol=None, significance=None, query=None, variant=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data`, `search_params`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.clin_var_search_variants(
    gene='BRCA1',
    max_results=10,
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
    "esearchresult": {
      "count": "<string>",
      "idlist": "<array>",
      "querytranslation": "<string>"
    }
  },
  "search_params": {
    "gene": "<string>",
    "condition": "<string>",
    "variant_id": "<string>"
  },
  "url": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_clinical_variants.md#dnakit.databases.clin_var_search_variants-parameters)

<!-- END GENERATED API: clinical-variants -->
