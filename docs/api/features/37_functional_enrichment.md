# 基因集与通路富集

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: enrichment -->

<span id="enrichment"></span>

## 1) 对基因集合执行指定库的富集分析 {#dnakit.databases.enrichr_enrich}

- **作用：** 对基因集合执行指定库的富集分析。
- **处理方法：** 通过 `EnrichrExt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.enrichr_enrich(*, operation=None, gene_list, library=None, top_n=None, species=None, verbose=True)`
- **输入：** 必须提供 `gene_list`。
- **输出：** 返回富集条目、P 值、z 分数、combined score 和命中基因。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.enrichr_enrich(
    gene_list=['TP53', 'BRCA1', 'EGFR', 'KRAS', 'MYC', 'BRAF', 'PTEN'],
    library='GO_Biological_Process_2023',
    top_n=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回富集条目、P 值、z 分数、combined score 和命中基因。

[完整参数与使用说明](../../faq/api/37_functional_enrichment.md#dnakit.databases.enrichr_enrich-parameters)

## 2) 查询基因与富集通路的连接关系 {#dnakit.databases.enrichr_gene_enrichment_analysis}

- **作用：** 查询基因与富集通路的连接关系。
- **处理方法：** 通过 `Enrichr` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.enrichr_gene_enrichment_analysis(*, gene_list, libs=None, verbose=True)`
- **输入：** 必须提供 `gene_list`。
- **输出：** 返回字典，包含 `connected_paths`, `connections` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.enrichr_gene_enrichment_analysis(
    gene_list=['TP53', 'BRCA1', 'EGFR'],
    libs=['KEGG_2021_Human'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "connected_paths": "<object>",
  "connections": "<object>"
}
```

[完整参数与使用说明](../../faq/api/37_functional_enrichment.md#dnakit.databases.enrichr_gene_enrichment_analysis-parameters)

## 3) 执行 PANTHER 基因集合过度代表分析 {#dnakit.databases.panther_enrichment}

- **作用：** 执行 PANTHER 基因集合过度代表分析。
- **处理方法：** 通过 `PANTHER` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.panther_enrichment(*, gene_list, organism=None, annotation_dataset=None, verbose=True)`
- **输入：** 必须提供 `gene_list`。
- **输出：** 返回字典，包含 `gene_list`, `organism`, `annotation_dataset`, `result_count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.panther_enrichment(
    gene_list='TP53,BRCA1,EGFR,KRAS,PTEN',
    organism=9606,
    annotation_dataset='GO:0008150',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_list": "<string>",
  "organism": "<integer>",
  "annotation_dataset": "<string>",
  "result_count": "<integer>",
  "enriched_terms": [
    {
      "term_id": "<string>",
      "term_label": "<string>",
      "number_in_list": "<integer>",
      "number_in_reference": "<integer>",
      "expected": "<number>",
      "fold_enrichment": "<number>",
      "pvalue": "<number>",
      "fdr": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/37_functional_enrichment.md#dnakit.databases.panther_enrichment-parameters)

## 4) 执行 Reactome 通路富集分析 {#dnakit.databases.reactome_analysis_pathway_enrichment}

- **作用：** 执行 Reactome 通路富集分析。
- **处理方法：** 通过 `ReactomeAnalysis` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.reactome_analysis_pathway_enrichment(*, identifiers, page_size=None, include_disease=None, projection=None, verbose=True)`
- **输入：** 必须提供 `identifiers`。
- **输出：** 返回字典，包含 `token`, `analysis_type`, `projection`, `identifiers_not_found` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.reactome_analysis_pathway_enrichment(
    identifiers='TP53\nBRCA1\nCDH1\nEGFR\nKRAS',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "token": "<string>",
  "analysis_type": "<string>",
  "projection": "<boolean>",
  "identifiers_not_found": "<integer>",
  "pathways_found": "<integer>",
  "pathways": [
    {
      "pathway_id": "<string>",
      "name": "<string>",
      "species": "<string>",
      "is_disease": "<boolean>",
      "is_lowest_level": "<boolean>",
      "entities_found": "<integer>",
      "entities_total": "<integer>",
      "entities_coverage": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/37_functional_enrichment.md#dnakit.databases.reactome_analysis_pathway_enrichment-parameters)

## 5) 分析蛋白质集合的功能与通路富集 {#dnakit.databases.string_functional_enrichment}

- **作用：** 分析蛋白质集合的功能与通路富集。
- **处理方法：** 通过 `STRING` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.string_functional_enrichment(*, protein_ids, species=None, category=None, verbose=True)`
- **输入：** 必须提供 `protein_ids`。
- **输出：** 返回记录列表；每条记录包含 `category`, `term`, `number_of_genes`, `number_of_genes_in_background` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.string_functional_enrichment(
    protein_ids=['TP53', 'MDM2', 'ATM', 'CHEK2', 'BRCA1'],
    species=9606,
    category='Process',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "category": "<string>",
    "term": "<string>",
    "number_of_genes": "<integer>",
    "number_of_genes_in_background": "<integer>",
    "ncbiTaxonId": "<integer>",
    "inputGenes": "<string>",
    "preferredNames": "<string>",
    "p_value": "<number>"
  }
]
```

[完整参数与使用说明](../../faq/api/37_functional_enrichment.md#dnakit.databases.string_functional_enrichment-parameters)

## 6) 使用多行蛋白标识执行 STRING 富集分析 {#dnakit.databases.string_get_functional_enrichment}

- **作用：** 使用多行蛋白标识执行 STRING 富集分析。
- **处理方法：** 通过 `API` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.string_get_functional_enrichment(*, identifiers, species=9606, category=None, verbose=True)`
- **输入：** 必须提供 `identifiers`。
- **输出：** 返回 STRING_functional_enrichment 的富集结果，见 [string_functional_enrichment](37_functional_enrichment.md#dnakit.databases.string_functional_enrichment)。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.string_get_functional_enrichment(
    identifiers='TP53\nBRCA1\nBRCA2',
    species=9606,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 STRING_functional_enrichment 的富集结果，见 [string_functional_enrichment](37_functional_enrichment.md#dnakit.databases.string_functional_enrichment)。

[完整参数与使用说明](../../faq/api/37_functional_enrichment.md#dnakit.databases.string_get_functional_enrichment-parameters)

<!-- END GENERATED API: enrichment -->
