# GWAS 研究与关联查询

依赖与返回值见[FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics)。

<!-- BEGIN GENERATED API: gwas -->

<span id="gwas"></span>

## 1) 查询指定 SNP 的 GWAS 关联 {#dnakit.databases.gwas_get_associations_for_snp}

- **作用：** 查询指定 SNP 的 GWAS 关联。
- **处理方法：** 通过 `GWASAssociationsForSNP` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_associations_for_snp(*, rs_id, sort=None, direction=None, size=None, page=None, verbose=True)`
- **输入：** 必须提供 `rs_id`。
- **输出：** 返回记录列表；每条记录包含 `association_id`, `p_value`, `beta`, `efo_traits` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_associations_for_snp(
    rs_id='rs7903146',
    sort='p_value',
    direction='asc',
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_associations_for_snp-parameters)

## 2) 查询指定 GWAS 研究的关联结果 {#dnakit.databases.gwas_get_associations_for_study}

- **作用：** 查询指定 GWAS 研究的关联结果。
- **处理方法：** 通过 `GWASAssociationsForStudy` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_associations_for_study(*, accession_id, size=None, page=None, verbose=True)`
- **输入：** 必须提供 `accession_id`。
- **输出：** 返回记录列表；每条记录包含 `association_id`, `p_value`, `beta`, `efo_traits` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_associations_for_study(
    accession_id='GCST000392',
    size=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_associations_for_study-parameters)

## 3) 查询指定性状的 GWAS 关联结果 {#dnakit.databases.gwas_get_associations_for_trait}

- **作用：** 查询指定性状的 GWAS 关联结果。
- **处理方法：** 通过 `GWASAssociationsForTrait` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_associations_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, page=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组恰好选一组：`disease_trait` / `efo_uri` / `efo_id` / `efo_trait`。
- **输出：** 返回记录列表；每条记录包含 `association_id`, `p_value`, `beta`, `efo_traits` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_associations_for_trait(
    disease_trait='type 2 diabetes',
    size=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_associations_for_trait-parameters)

## 4) 按 rsID 查询 GWAS SNP 详情 {#dnakit.databases.gwas_get_snp_by_id}

- **作用：** 按 rsID 查询 GWAS SNP 详情。
- **处理方法：** 通过 `GWASSNPByID` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_snp_by_id(*, rs_id, verbose=True)`
- **输入：** 必须提供 `rs_id`。
- **输出：** 返回字典，包含 `rs_id`, `merged`, `functional_class`, `last_update_date` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_snp_by_id(
    rs_id='rs429358',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "rs_id": "<string>",
  "merged": "<integer>",
  "functional_class": "<string>",
  "last_update_date": "<string>",
  "locations": [
    {
      "chromosome_name": "<string>",
      "chromosome_position": "<integer>",
      "region": "<object>"
    }
  ],
  "alleles": "<string>",
  "most_severe_consequence": "<string>",
  "mapped_genes": [
    "<string>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_snp_by_id-parameters)

## 5) 查询映射到指定基因的 GWAS SNP {#dnakit.databases.gwas_get_snps_for_gene}

- **作用：** 查询映射到指定基因的 GWAS SNP。
- **处理方法：** 通过 `GWASSNPsForGene` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_snps_for_gene(*, size=None, page=None, gene_symbol=None, gene=None, mapped_gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `rs_id`, `merged`, `functional_class`, `last_update_date` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_snps_for_gene(
    mapped_gene='BRCA1',
    size=10,
    page=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "rs_id": "<string>",
    "merged": "<integer>",
    "functional_class": "<string>",
    "last_update_date": "<string>",
    "locations": [
      "<object>"
    ],
    "alleles": "<string>",
    "most_severe_consequence": "<string>",
    "mapped_genes": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_snps_for_gene-parameters)

## 6) 查询指定性状的 GWAS 研究 {#dnakit.databases.gwas_get_studies_for_trait}

- **作用：** 查询指定性状的 GWAS 研究。
- **处理方法：** 通过 `GWASStudiesForTrait` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_studies_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组恰好选一组：`disease_trait` / `efo_uri` / `efo_id` / `efo_trait`。
- **输出：** 返回记录列表；每条记录包含 `accession_id`, `disease_trait`, `efo_traits`, `initial_sample_size` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_studies_for_trait(
    disease_trait='hypertension',
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "accession_id": "<string>",
    "disease_trait": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "initial_sample_size": "<string>",
    "replication_sample_size": "<string>",
    "discovery_ancestry": [
      "<string>"
    ],
    "replication_ancestry": [
      "<string>"
    ],
    "cohort": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_studies_for_trait-parameters)

## 7) 按研究编号查询 GWAS 详情 {#dnakit.databases.gwas_get_study_by_id}

- **作用：** 按研究编号查询 GWAS 详情。
- **处理方法：** 通过 `GWASStudyByID` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_study_by_id(*, study_id=None, accession_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `accession_id`, `disease_trait`, `efo_traits`, `initial_sample_size` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_study_by_id(
    study_id='GCST000392',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "accession_id": "<string>",
  "disease_trait": "<string>",
  "efo_traits": [
    {
      "efo_id": "<string>",
      "efo_trait": "<string>"
    }
  ],
  "initial_sample_size": "<string>",
  "replication_sample_size": "<string>",
  "discovery_ancestry": [
    "<string>"
  ],
  "replication_ancestry": [
    "<string>"
  ],
  "cohort": [
    "<string>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_study_by_id-parameters)

## 8) 检索疾病或性状相关的 GWAS 变异 {#dnakit.databases.gwas_get_variants_for_trait}

- **作用：** 检索疾病或性状相关的 GWAS 变异。
- **处理方法：** 通过 `GWASVariantsForTrait` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_get_variants_for_trait(*, disease_trait=None, trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, limit=None, page=None, sort=None, direction=None, p_value_threshold=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组恰好选一组：`disease_trait` / `trait` / `efo_uri` / `efo_id` / `efo_trait`。
- **输出：** 返回记录列表；每条记录包含 `association_id`, `risk_frequency`, `pvalue_description`, `range` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_get_variants_for_trait(
    disease_trait='breast cancer',
    size=20,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "association_id": "<integer>",
    "risk_frequency": "<string>",
    "pvalue_description": "<string>",
    "range": "<string>",
    "locations": [
      "<object>"
    ],
    "rs_id": "<string>",
    "merged": "<integer>",
    "functional_class": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_variants_for_trait-parameters)

## 9) 按关键词检索 GWAS 关联 {#dnakit.databases.gwas_search_associations}

- **作用：** 按关键词检索 GWAS 关联。
- **处理方法：** 通过 `GWASAssociationSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_search_associations(*, disease_trait=None, query=None, efo_uri=None, efo_id=None, efo_trait=None, rs_id=None, accession_id=None, p_value=None, p_value_threshold=None, sort=None, direction=None, size=None, page=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `association_id`, `p_value`, `beta`, `efo_traits` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_search_associations(
    disease_trait='coronary artery disease',
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_associations-parameters)

## 10) 按基因检索 GWAS 风险等位基因与显著性 {#dnakit.databases.gwas_search_associations_by_gene}

- **作用：** 按基因检索 GWAS 风险等位基因与显著性。
- **处理方法：** 通过 `GWASGeneSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_search_associations_by_gene(*, gene_name, size=None, verbose=True)`
- **输入：** 必须提供 `gene_name`。
- **输出：** 返回字典，包含 `gene_name`, `association_count`, `associations`, `total_found` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_search_associations_by_gene(
    gene_name='TP53',
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_name": "<string>",
  "association_count": "<integer>",
  "associations": [
    "<object>"
  ],
  "total_found": "<integer>",
  "queried_gene_name": "<string>",
  "gene_name_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_associations_by_gene-parameters)

## 11) 按 rsID 或基因检索 GWAS SNP {#dnakit.databases.gwas_search_snps}

- **作用：** 按 rsID 或基因检索 GWAS SNP。
- **处理方法：** 通过 `GWASSNPSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_search_snps(*, rs_id=None, mapped_gene=None, size=None, page=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `rs_id`, `merged`, `functional_class`, `last_update_date` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_search_snps(
    mapped_gene='BRCA1',
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "rs_id": "<string>",
    "merged": "<integer>",
    "functional_class": "<string>",
    "last_update_date": "<string>",
    "locations": [
      "<object>"
    ],
    "alleles": "<string>",
    "most_severe_consequence": "<string>",
    "mapped_genes": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_snps-parameters)

## 12) 按性状与队列条件检索 GWAS 研究 {#dnakit.databases.gwas_search_studies}

- **作用：** 按性状与队列条件检索 GWAS 研究。
- **处理方法：** 通过 `GWASStudySearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gwas_search_studies(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `accession_id`, `disease_trait`, `efo_traits`, `initial_sample_size` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gwas_search_studies(
    disease_trait='diabetes',
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "accession_id": "<string>",
    "disease_trait": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "initial_sample_size": "<string>",
    "replication_sample_size": "<string>",
    "discovery_ancestry": [
      "<string>"
    ],
    "replication_ancestry": [
      "<string>"
    ],
    "cohort": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_studies-parameters)

<!-- END GENERATED API: gwas -->
