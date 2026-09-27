# 肿瘤变异与治疗证据查询

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: cancer-variants -->

<span id="cancer-variants"></span>

## 1) 查询 CIViC 分子特征组合详情 {#dnakit.databases.civic_get_molecular_profile}

- **作用：** 查询 CIViC 分子特征组合详情。
- **处理方法：** 通过 `CIViC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.civic_get_molecular_profile(*, molecular_profile_id, evidence_limit=None, verbose=True)`
- **输入：** 必须提供 `molecular_profile_id`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.civic_get_molecular_profile(
    molecular_profile_id=12,
    evidence_limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "molecularProfile": {
      "id": "<integer>",
      "name": "<string>",
      "molecularProfileScore": "<number>",
      "description": "<string>",
      "evidenceItems": "<object>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_get_molecular_profile-parameters)

## 2) 查询 CIViC 肿瘤变异详情 {#dnakit.databases.civic_get_variant}

- **作用：** 查询 CIViC 肿瘤变异详情。
- **处理方法：** 通过 `CIViC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.civic_get_variant(*, variant_id, verbose=True)`
- **输入：** 必须提供 `variant_id`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.civic_get_variant(
    variant_id=2,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "variant": {
      "id": "<integer>",
      "name": "<string>",
      "molecularProfiles": "<object>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_get_variant-parameters)

## 3) 按基因查询 CIViC 肿瘤变异 {#dnakit.databases.civic_get_variants_by_gene}

- **作用：** 按基因查询 CIViC 肿瘤变异。
- **处理方法：** 通过 `CIViC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.civic_get_variants_by_gene(*, gene_id=None, gene_name=None, gene_symbol=None, gene=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.civic_get_variants_by_gene(
    gene_name='EGFR',
    limit=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "gene": {
      "id": "<integer>",
      "name": "<string>",
      "variants": "<object>"
    }
  },
  "metadata": {
    "source": "<string>",
    "format": "<string>",
    "note": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_get_variants_by_gene-parameters)

## 4) 检索 CIViC 临床解释结论 {#dnakit.databases.civic_search_assertions}

- **作用：** 检索 CIViC 临床解释结论。
- **处理方法：** 通过 `CIViC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.civic_search_assertions(*, limit=None, therapy=None, therapy_name=None, disease=None, disease_name=None, variant_name=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.civic_search_assertions(
    limit=5,
    therapy='imatinib',
    disease='leukemia',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "assertions": {
      "totalCount": "<integer>",
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_assertions-parameters)

## 5) 检索 CIViC 变异临床证据 {#dnakit.databases.civic_search_evidence_items}

- **作用：** 检索 CIViC 变异临床证据。
- **处理方法：** 通过 `CIViC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.civic_search_evidence_items(*, limit=None, status=None, therapy=None, therapy_name=None, disease=None, disease_name=None, molecular_profile=None, evidence_type=None, significance=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.civic_search_evidence_items(
    limit=5,
    therapy='imatinib',
    disease='leukemia',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "evidenceItems": {
      "totalCount": "<integer>",
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_evidence_items-parameters)

## 6) 检索 CIViC 收录的肿瘤相关基因 {#dnakit.databases.civic_search_genes}

- **作用：** 检索 CIViC 收录的肿瘤相关基因。
- **处理方法：** 通过 `CIViC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.civic_search_genes(*, limit=None, name=None, query=None, gene_name=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.civic_search_genes(
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "genes": {
      "totalCount": "<integer>",
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  },
  "metadata": {
    "source": "<string>",
    "format": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_genes-parameters)

## 7) 检索 CIViC 收录的药物与治疗方案 {#dnakit.databases.civic_search_therapies}

- **作用：** 检索 CIViC 收录的药物与治疗方案。
- **处理方法：** 通过 `CIViC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.civic_search_therapies(*, limit=None, name=None, query=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.civic_search_therapies(
    name='imatinib',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "browseTherapies": {
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_therapies-parameters)

## 8) 按基因查询 COSMIC 体细胞突变 {#dnakit.databases.cosmic_get_mutations_by_gene}

- **作用：** 按基因查询 COSMIC 体细胞突变。
- **处理方法：** 通过 `COSMIC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cosmic_get_mutations_by_gene(*, operation=None, gene=None, gene_name=None, max_results=None, genome_build=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `gene`, `total_count`, `results`, `genome_build` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cosmic_get_mutations_by_gene(
    operation='get_by_gene',
    gene='BRAF',
    max_results=50,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene": "<string>",
  "total_count": "<integer>",
  "results": [
    {
      "mutation_id": "<string>",
      "GeneName": "<string>",
      "MutationCDS": "<string>",
      "MutationAA": "<string>",
      "MutationGenomePosition": "<string>",
      "PrimarySite": "<string>",
      "PrimaryHistology": "<string>",
      "Fathmm": "<string>"
    }
  ],
  "genome_build": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.cosmic_get_mutations_by_gene-parameters)

## 9) 检索 COSMIC 肿瘤体细胞突变 {#dnakit.databases.cosmic_search_mutations}

- **作用：** 检索 COSMIC 肿瘤体细胞突变。
- **处理方法：** 通过 `COSMIC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cosmic_search_mutations(*, operation=None, terms=None, query=None, max_results=None, genome_build=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total_count`, `results`, `genome_build` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cosmic_search_mutations(
    operation='search',
    terms='BRAF V600E',
    max_results=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total_count": "<integer>",
  "results": [
    {
      "mutation_id": "<string>",
      "display": "<array>",
      "GeneName": "<string>",
      "MutationCDS": "<string>",
      "MutationAA": "<string>",
      "PrimarySite": "<string>",
      "PrimaryHistology": "<string>"
    }
  ],
  "genome_build": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.cosmic_search_mutations-parameters)

## 10) 查询 OncoKB 肿瘤变异作用与治疗证据 {#dnakit.databases.onco_kb_annotate_variant}

- **作用：** 查询 OncoKB 肿瘤变异作用与治疗证据。
- **处理方法：** 通过 `OncoKB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.onco_kb_annotate_variant(*, operation=None, gene=None, gene_symbol=None, variant=None, tumor_type=None, alteration=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `query`, `geneExist`, `variantExist`, `mutationEffect` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.onco_kb_annotate_variant(
    gene='BRAF',
    variant='V600E',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "query": {
    "id": "<null>",
    "referenceGenome": "<string>",
    "hugoSymbol": "<string>",
    "entrezGeneId": "<integer>",
    "alteration": "<string>",
    "alterationType": "<null>",
    "svType": "<null>",
    "tumorType": "<string>"
  },
  "geneExist": "<boolean>",
  "variantExist": "<boolean>",
  "mutationEffect": {
    "knownEffect": "<string>",
    "description": "<string>",
    "citations": {
      "pmids": "<array>",
      "abstracts": "<array>"
    }
  },
  "highestSensitiveLevel": "<string>",
  "highestResistanceLevel": "<null>",
  "highestDiagnosticImplicationLevel": "<string>",
  "highestPrognosticImplicationLevel": "<null>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.onco_kb_annotate_variant-parameters)

## 11) 按区间检索肿瘤拷贝数变异样本 {#dnakit.databases.progenetix_cnv_search}

- **作用：** 按区间检索肿瘤拷贝数变异样本。
- **处理方法：** 通过 `Progenetix` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.progenetix_cnv_search(*, reference_name, start, end, variant_type=None, filters=None, limit=None, verbose=True)`
- **输入：** 必须提供 `reference_name`, `start`, `end`。
- **输出：** 返回字典，包含 `region`, `variant_type`, `filters`, `total_count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.progenetix_cnv_search(
    reference_name='refseq:NC_000007.14',
    start=55019017,
    end=55211628,
    variant_type='DUP',
    filters='NCIT:C4017',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "region": "<string>",
  "variant_type": "<string>",
  "filters": "<string>",
  "total_count": "<integer>",
  "returned_count": "<integer>",
  "biosamples": [
    {
      "id": "<string>",
      "biosample_status": "<string>",
      "histological_diagnosis": "<string>",
      "histological_diagnosis_id": "<string>",
      "external_references": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_cnv_search-parameters)

## 12) 列出 Progenetix 肿瘤拷贝数变异队列 {#dnakit.databases.progenetix_list_cohorts}

- **作用：** 列出 Progenetix 肿瘤拷贝数变异队列。
- **处理方法：** 通过 `Progenetix` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.progenetix_list_cohorts(*, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total_cohorts`, `cohorts` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.progenetix_list_cohorts(
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total_cohorts": "<integer>",
  "cohorts": [
    {
      "id": "<string>",
      "name": "<string>",
      "cohort_size": "<integer>",
      "cohort_type": "<string>",
      "data_types": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_list_cohorts-parameters)

## 13) 列出 Progenetix 可用筛选术语 {#dnakit.databases.progenetix_list_filtering_terms}

- **作用：** 列出 Progenetix 可用筛选术语。
- **处理方法：** 通过 `Progenetix` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.progenetix_list_filtering_terms(*, prefixes=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `prefixes`, `total_terms`, `terms` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.progenetix_list_filtering_terms(
    prefixes='NCIT',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "prefixes": "<string>",
  "total_terms": "<integer>",
  "terms": [
    {
      "id": "<string>",
      "label": "<string>",
      "count": "<integer>",
      "type": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_list_filtering_terms-parameters)

## 14) 按疾病本体检索肿瘤生物样本 {#dnakit.databases.progenetix_search_biosamples}

- **作用：** 按疾病本体检索肿瘤生物样本。
- **处理方法：** 通过 `Progenetix` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.progenetix_search_biosamples(*, filters, limit=None, skip=None, verbose=True)`
- **输入：** 必须提供 `filters`。
- **输出：** 返回字典，包含 `filters`, `total_count`, `returned_count`, `biosamples` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.progenetix_search_biosamples(
    filters='NCIT:C3058',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "filters": "<string>",
  "total_count": "<integer>",
  "returned_count": "<integer>",
  "biosamples": [
    {
      "id": "<string>",
      "biosample_status": "<string>",
      "histological_diagnosis": "<string>",
      "histological_diagnosis_id": "<string>",
      "pathological_stage": "<string>",
      "pathological_tnm": "<string>",
      "sample_origin_type": "<string>",
      "external_references": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_search_biosamples-parameters)

<!-- END GENERATED API: cancer-variants -->
