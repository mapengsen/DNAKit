# 队列表型关联与精细定位查询

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: phewas -->

<span id="phewas"></span>

## 1) 查询日本生物样本库中的变异表型关联 {#dnakit.databases.bio_bank_japan_phewas_by_variant}

- **作用：** 查询日本生物样本库中的变异表型关联。
- **处理方法：** 通过 `PheWebPheWAS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.bio_bank_japan_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `variant`, `rsids`, `nearest_genes`, `associations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.bio_bank_japan_phewas_by_variant(
    rsid='rs7903146',
    max_pval=5e-08,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": "<string>",
  "rsids": "<string>",
  "nearest_genes": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "phenostring": "<string>",
      "category": "<string>",
      "pval": "<number>",
      "beta": "<number>",
      "sebeta": "<number>",
      "af": "<number>",
      "num_cases": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phewas.md#dnakit.databases.bio_bank_japan_phewas_by_variant-parameters)

## 2) 查询 FinnGen 变异精细定位结果 {#dnakit.databases.finn_gen_get_variant_finemapping}

- **作用：** 查询 FinnGen 变异精细定位结果。
- **处理方法：** 通过 `FinnGen` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.finn_gen_get_variant_finemapping(*, variant, verbose=True)`
- **输入：** 必须提供 `variant`。
- **输出：** 返回记录列表；每条记录包含 `phenocode`, `chromosome`, `start`, `end` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.finn_gen_get_variant_finemapping(
    variant='19:44908684:T:C',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "phenocode": "<string>",
    "chromosome": "<integer>",
    "start": "<integer>",
    "end": "<integer>",
    "type": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_phewas.md#dnakit.databases.finn_gen_get_variant_finemapping-parameters)

## 3) 查询基因稀有变异负荷与表型的关联 {#dnakit.databases.genebass_gene_burden_phewas}

- **作用：** 查询基因稀有变异负荷与表型的关联。
- **处理方法：** 通过 `Genebass` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.genebass_gene_burden_phewas(*, gene, burden_set=None, limit=None, max_pval=None, verbose=True)`
- **输入：** 必须提供 `gene`。
- **输出：** 返回字典，包含 `gene_id`, `gene_symbol`, `burden_set`, `associations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.genebass_gene_burden_phewas(
    gene='TCF7L2',
    burden_set='pLoF',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_id": "<string>",
  "gene_symbol": "<string>",
  "burden_set": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "description": "<string>",
      "trait_type": "<string>",
      "pval": "<number>",
      "pval_skat": "<number>",
      "beta": "<number>",
      "total_variants": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phewas.md#dnakit.databases.genebass_gene_burden_phewas-parameters)

## 4) 查询台湾精准医疗队列的变异表型关联 {#dnakit.databases.tpmi_phewas_by_variant}

- **作用：** 查询台湾精准医疗队列的变异表型关联。
- **处理方法：** 通过 `PheWebPheWAS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.tpmi_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `variant`, `rsids`, `nearest_genes`, `associations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.tpmi_phewas_by_variant(
    rsid='rs7903146',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": "<string>",
  "rsids": "<string>",
  "nearest_genes": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "phenostring": "<string>",
      "category": "<string>",
      "pval": "<number>",
      "beta": "<number>",
      "sebeta": "<number>",
      "af": "<number>",
      "num_cases": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phewas.md#dnakit.databases.tpmi_phewas_by_variant-parameters)

## 5) 查询 UKB-TOPMed 的变异表型关联 {#dnakit.databases.ukbtopmed_phewas_by_variant}

- **作用：** 查询 UKB-TOPMed 的变异表型关联。
- **处理方法：** 通过 `PheWebPheWAS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ukbtopmed_phewas_by_variant(*, rsid=None, variant=None, limit=None, max_pval=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `variant`, `rsids`, `nearest_genes`, `associations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ukbtopmed_phewas_by_variant(
    rsid='rs7903146',
    max_pval=5e-08,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": "<string>",
  "rsids": "<string>",
  "nearest_genes": "<string>",
  "associations": [
    {
      "phenocode": "<string>",
      "phenostring": "<string>",
      "category": "<string>",
      "pval": "<number>",
      "beta": "<number>",
      "sebeta": "<number>",
      "af": "<number>",
      "num_cases": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_phewas.md#dnakit.databases.ukbtopmed_phewas_by_variant-parameters)

<!-- END GENERATED API: phewas -->
