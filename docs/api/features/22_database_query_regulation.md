# 表达、调控与肿瘤数据查询

查询基因表达、调控和肿瘤数据。

## 1) 按关键词查询 GEO 表达数据 {#search-expression}

- **作用：** 查找匹配关键词的 GEO 表达数据记录，获取摘要与来源信息。
- **计算方法：** 通过 NCBI Entrez 查询 GEO，返回统一的 QueryResult。
- **API：** `dnakit.search.expression(term[必须], retmax[可选], config[可选], verbose[可选])`
- **输入：** Entrez 查询词；可选返回记录上限 retmax 和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import expression

result = expression("single cell RNA sequencing[All Fields]", retmax=20, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
expression NCBI Entrez <返回记录数>
```

<!-- BEGIN GENERATED API: regulation -->

<span id="regulation"></span>

## 2) 分页查询 cBioPortal 肿瘤研究目录 {#dnakit.databases.c_bio_portal_get_cancer_studies}

- **作用：** 分页查询 cBioPortal 肿瘤研究目录。
- **处理方法：** 通过 `CBioPortal` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.c_bio_portal_get_cancer_studies(*, limit=None, offset=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `studyId`, `name`, `description`, `cancerTypeId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_cancer_studies(
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
    "studyId": "<string>",
    "name": "<string>",
    "description": "<string>",
    "cancerTypeId": "<string>",
    "publicStudy": "<boolean>",
    "pmid": "<string>",
    "citation": "<string>",
    "groups": [
      "<value>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_cancer_studies-parameters)

## 3) 按 Entrez 编号查询肿瘤数据库基因信息 {#dnakit.databases.c_bio_portal_get_gene_info}

- **作用：** 按 Entrez 编号查询肿瘤数据库基因信息。
- **处理方法：** 通过 `CBioPortal` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.c_bio_portal_get_gene_info(*, entrez_gene_id, verbose=True)`
- **输入：** 必须提供 `entrez_gene_id`。
- **输出：** 返回字典，包含 `entrezGeneId`, `hugoGeneSymbol`, `type`, `cytoband` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_gene_info(
    entrez_gene_id=672,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "entrezGeneId": "<integer>",
  "hugoGeneSymbol": "<string>",
  "type": "<string>",
  "cytoband": "<string>",
  "length": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_gene_info-parameters)

## 4) 查询肿瘤研究的分子数据类型 {#dnakit.databases.c_bio_portal_get_molecular_profiles}

- **作用：** 查询肿瘤研究的分子数据类型。
- **处理方法：** 通过 `CBioPortal` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.c_bio_portal_get_molecular_profiles(*, study_id, verbose=True)`
- **输入：** 必须提供 `study_id`。
- **输出：** 返回记录列表；每条记录包含 `molecularProfileId`, `molecularAlterationType`, `datatype`, `name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_molecular_profiles(
    study_id='brca_tcga',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "molecularProfileId": "<string>",
    "molecularAlterationType": "<string>",
    "datatype": "<string>",
    "name": "<string>",
    "description": "<string>",
    "showProfileInAnalysisTab": "<boolean>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_molecular_profiles-parameters)

## 5) 查询肿瘤研究中指定基因的突变 {#dnakit.databases.c_bio_portal_get_mutations}

- **作用：** 查询肿瘤研究中指定基因的突变。
- **处理方法：** 通过 `CBioPortal` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.c_bio_portal_get_mutations(*, study_id, gene_list, sample_list_id=None, verbose=True)`
- **输入：** 必须提供 `study_id`, `gene_list`。
- **输出：** 返回记录列表；每条记录包含 `molecularProfileId`, `sampleId`, `patientId`, `studyId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_mutations(
    study_id='brca_tcga',
    gene_list='BRCA2',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "molecularProfileId": "<string>",
    "sampleId": "<string>",
    "patientId": "<string>",
    "studyId": "<string>",
    "entrezGeneId": "<integer>",
    "gene": {
      "entrezGeneId": "<integer>",
      "hugoGeneSymbol": "<string>",
      "type": "<string>"
    },
    "mutationType": "<string>",
    "proteinChange": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_mutations-parameters)

## 6) 查询肿瘤研究的样本集合 {#dnakit.databases.c_bio_portal_get_samples}

- **作用：** 查询肿瘤研究的样本集合。
- **处理方法：** 通过 `CBioPortal` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.c_bio_portal_get_samples(*, study_id, page_size=None, page_number=None, verbose=True)`
- **输入：** 必须提供 `study_id`。
- **输出：** 返回记录列表；每条记录包含 `sampleId`, `sampleType`, `patientId`, `studyId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_samples(
    study_id='brca_tcga',
    page_size=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "sampleId": "<string>",
    "sampleType": "<string>",
    "patientId": "<string>",
    "studyId": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_samples-parameters)

## 7) 获取 CELLxGENE 原始 H5AD 数据文件 {#dnakit.databases.cellx_gene_download_h5ad}

- **作用：** 获取 CELLxGENE 原始 H5AD 数据文件。
- **处理方法：** 通过 CELLxGENE Census 定位数据集，下载对应的 H5AD 文件。
- **API：** `dnakit.databases.cellx_gene_download_h5ad(*, operation=None, dataset_id, output_path=None, census_version=None, verbose=True)`
- **输入：** 必须提供 `dataset_id`。 需安装 `cellxgene-census`。 运行前准备示例涉及的本地文件与外部程序，并将路径替换为自己的输入或输出位置。
- **输出：** 返回字典，包含 `dataset_id`, `output_path`, `uri`, `message` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cellx_gene_download_h5ad(
    operation='download_h5ad',
    dataset_id='4eb29386-de81-452f-b3c0-e00844e8c7fd',
    census_version='stable',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "dataset_id": "<string>",
  "output_path": "<string>",
  "uri": "<string>",
  "message": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_download_h5ad-parameters)

## 8) 按筛选条件查询单细胞元数据 {#dnakit.databases.cellx_gene_get_cell_metadata}

- **作用：** 按筛选条件查询单细胞元数据。
- **处理方法：** 通过 `CELLxGENECensus` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cellx_gene_get_cell_metadata(*, operation=None, organism=None, obs_value_filter=None, column_names=None, census_version=None, gene=None, tissue=None, max_cells=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需安装 `cellxgene-census`。
- **输出：** 返回记录列表；每条记录包含 `cell_type`, `tissue_general`, `disease`, `assay` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cellx_gene_get_cell_metadata(
    operation='get_obs_metadata',
    organism='Homo sapiens',
    obs_value_filter="tissue_general == 'brain' and disease == 'normal'",
    census_version='stable',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "cell_type": "<string>",
    "tissue_general": "<string>",
    "disease": "<string>",
    "assay": "<string>",
    "sex": "<string>",
    "donor_id": "<string>",
    "dataset_id": "<string>",
    "soma_joinid": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_get_cell_metadata-parameters)

## 9) 获取预计算的单细胞嵌入表示 {#dnakit.databases.cellx_gene_get_embeddings}

- **作用：** 获取预计算的单细胞嵌入表示。
- **处理方法：** 通过 `CELLxGENECensus` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cellx_gene_get_embeddings(*, operation=None, organism=None, embedding_name=None, census_version=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需安装 `cellxgene-census`。
- **输出：** 返回字典，包含 `message`, `organism`, `embedding_name`, `available` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cellx_gene_get_embeddings(
    operation='get_embeddings',
    organism='Homo sapiens',
    census_version='stable',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "message": "<string>",
  "organism": "<string>",
  "embedding_name": "<string>",
  "available": "<boolean>",
  "shape": [
    "<value>"
  ],
  "available_embeddings": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_get_embeddings-parameters)

## 10) 按细胞与基因条件查询表达矩阵 {#dnakit.databases.cellx_gene_get_expression_data}

- **作用：** 按细胞与基因条件查询表达矩阵。
- **处理方法：** 通过 `CELLxGENECensus` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.cellx_gene_get_expression_data(*, operation=None, organism=None, obs_value_filter=None, var_value_filter=None, obs_column_names=None, var_column_names=None, census_version=None, gene=None, tissue=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需安装 `cellxgene-census`。
- **输出：** 返回字典，包含 `n_obs`, `n_vars`, `obs_columns`, `var_columns` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.cellx_gene_get_expression_data(
    operation='get_anndata',
    organism='Homo sapiens',
    obs_value_filter="tissue_general == 'brain' and disease == 'normal'",
    var_value_filter="feature_name in ['SOX2', 'GFAP', 'OLIG2']",
    census_version='stable',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "n_obs": "<integer>",
  "n_vars": "<integer>",
  "obs_columns": [
    "<value>"
  ],
  "var_columns": [
    "<value>"
  ],
  "message": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_get_expression_data-parameters)

## 11) 生成 ChIP-Atlas 富集任务提交说明 {#dnakit.databases.ch_ipatlas_enrichment_analysis}

- **作用：** 根据 BED、motif 或基因列表生成 ChIP-Atlas 富集任务的提交说明；本函数不计算富集统计量。
- **处理方法：** 通过 `ChIPAtlas` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_ipatlas_enrichment_analysis(*, operation=None, bed_data=None, motif=None, gene_list=None, genome=None, antigen_class=None, cell_type_class=None, threshold=None, distance=None, gene=None, cell_type=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 字典包含 `status=requires_submission`、`analysis_completed=False` 和任务提交说明；本次调用尚未完成富集检验。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_enrichment_analysis(
    operation='enrichment_analysis',
    gene_list=['TP53', 'MDM2', 'CDKN1A'],
    genome='hg38',
    distance='5000',
    verbose=False,
)
print(result)
```

- **示例结果：** 字典包含 `status=requires_submission`、`analysis_completed=False` 和任务提交说明；本次调用尚未完成富集检验。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "message": "<string>",
  "genes": [
    "<string>"
  ],
  "distance_from_tss": "<string>",
  "url": "<string>",
  "note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_enrichment_analysis-parameters)

## 12) 检索 ChIP-Atlas 实验元数据 {#dnakit.databases.ch_ipatlas_get_experiments}

- **作用：** 检索 ChIP-Atlas 实验元数据。
- **处理方法：** 通过 `ChIPAtlas` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_ipatlas_get_experiments(*, operation=None, genome=None, antigen=None, cell_type=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `message`, `metadata_file_url`, `web_search_url`, `file_size` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_get_experiments(
    operation='get_experiment_list',
    genome='hg38',
    antigen='CTCF',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "message": "<string>",
  "metadata_file_url": "<string>",
  "web_search_url": "<string>",
  "file_size": "<string>",
  "total_experiments": "<string>",
  "filters_requested": "<object>",
  "recommendation": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_get_experiments-parameters)

## 13) 查询 ChIP-Atlas 峰资源与本地区间重叠 {#dnakit.databases.ch_ipatlas_get_peak_data}

- **作用：** 查询 ChIP-Atlas 实验峰数据入口；gene 模式先定位 GRCh38 基因区域，提供 peak_file 后计算 BED 重叠，否则返回数据获取计划。
- **处理方法：** 通过 `ChIPAtlas` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_ipatlas_get_peak_data(*, operation=None, experiment_id=None, genome=None, format=None, threshold=None, gene=None, experiment_type=None, limit=None, peak_file=None, flank_bp=None, max_file_bytes=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `experiment_id`, `genome`, `format`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_get_peak_data(
    operation='get_peak_data',
    experiment_id='SRX097088',
    genome='hg19',
    format='bigwig',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "experiment_id": "<string>",
  "genome": "<string>",
  "format": "<string>",
  "url": "<string>",
  "message": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_get_peak_data-parameters)

## 14) 按抗原或细胞类型检索 ChIP-Atlas 数据集 {#dnakit.databases.ch_ipatlas_search_datasets}

- **作用：** 按抗原或细胞类型检索 ChIP-Atlas 数据集。
- **处理方法：** 通过 `ChIPAtlas` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_ipatlas_search_datasets(*, operation=None, antigen=None, cell_type=None, genome=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `search_key`, `search_value`, `genome`, `num_results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_search_datasets(
    operation='search_datasets',
    antigen='CTCF',
    genome='hg38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "search_key": "<string>",
  "search_value": "<string>",
  "genome": "<string>",
  "num_results": "<integer>",
  "results": [
    {
      "genome": "<string>",
      "class": "<string>",
      "name": "<string>",
      "num_experiments": "<string>",
      "experiment_ids": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_search_datasets-parameters)

## 15) 查询 CRISPR 筛选的基因依赖信息 {#dnakit.databases.dep_map_get_gene_dependencies}

- **作用：** 查询 CRISPR 筛选的基因依赖信息。
- **处理方法：** 通过 `DepMap` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dep_map_get_gene_dependencies(*, gene_symbol, model_id=None, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** 返回字典，包含 `status`, `data`, `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dep_map_get_gene_dependencies(
    gene_symbol='EGFR',
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
    "genes": [
      "<value>"
    ]
  },
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.dep_map_get_gene_dependencies-parameters)

## 16) 按名称检索肿瘤细胞系 {#dnakit.databases.dep_map_search_cell_lines}

- **作用：** 按名称检索肿瘤细胞系。
- **处理方法：** 通过 `DepMap` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dep_map_search_cell_lines(*, query, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `status`, `data`, `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dep_map_search_cell_lines(
    query='A549',
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
    "query": "<string>",
    "cell_lines": [
      "<value>"
    ],
    "count": "<integer>"
  },
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.dep_map_search_cell_lines-parameters)

## 17) 查询 eQTL 数据集中的变异表达关联 {#dnakit.databases.e_qtl_get_associations}

- **作用：** 查询 eQTL 数据集中的变异表达关联。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.e_qtl_get_associations(*, dataset_id, size=None, gene_id=None, variant=None, molecular_trait_id=None, verbose=True)`
- **输入：** 必须提供 `dataset_id`。
- **输出：** 返回记录列表；每条记录包含 `molecular_trait_id`, `gene_id`, `variant`, `rsid` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.e_qtl_get_associations(
    dataset_id='QTD000001',
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
    "molecular_trait_id": "<string>",
    "gene_id": "<string>",
    "variant": "<string>",
    "rsid": "<string>",
    "chromosome": "<string>",
    "position": "<integer>",
    "ref": "<string>",
    "alt": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.e_qtl_get_associations-parameters)

## 18) 列出可用的 eQTL 研究与组织数据集 {#dnakit.databases.e_qtl_list_datasets}

- **作用：** 列出可用的 eQTL 研究与组织数据集。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.e_qtl_list_datasets(*, size=None, quant_method=None, tissue_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `dataset_id`, `study_id`, `study_label`, `quant_method` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.e_qtl_list_datasets(
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
    "dataset_id": "<string>",
    "study_id": "<string>",
    "study_label": "<string>",
    "quant_method": "<string>",
    "tissue_id": "<string>",
    "tissue_label": "<string>",
    "condition_label": "<string>",
    "sample_size": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.e_qtl_list_datasets-parameters)

## 19) 查询 ENCODE 生物样本完整元数据 {#dnakit.databases.encode_get_biosample}

- **作用：** 查询 ENCODE 生物样本完整元数据。
- **处理方法：** 通过 `ENCODE` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_get_biosample(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字典，包含 `accession`, `url`, `aliases`, `schema_version` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_get_biosample(
    accession='ENCBS000AAA',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "accession": "<string>",
  "url": "<string>",
  "aliases": [
    "<string>"
  ],
  "schema_version": "<string>",
  "status": "<string>",
  "lab": "<string>",
  "award": "<string>",
  "date_created": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_get_biosample-parameters)

## 20) 检索染色质状态分割注释 {#dnakit.databases.encode_get_chromatin_state}

- **作用：** 检索染色质状态分割注释。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_get_chromatin_state(*, biosample_term_name=None, organism=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `annotations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_get_chromatin_state(
    biosample_term_name='K562',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "annotations": [
    {
      "accession": "<string>",
      "annotation_type": "<string>",
      "description": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_get_chromatin_state-parameters)

## 21) 查询 ENCODE 实验元数据与分析信息 {#dnakit.databases.encode_get_experiment}

- **作用：** 查询 ENCODE 实验元数据与分析信息。
- **处理方法：** 通过 `ENCODE` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_get_experiment(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字典，包含 `assay_term_name`, `biosample_ontology`, `documents`, `references` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_get_experiment(
    accession='ENCSR000AKS',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "assay_term_name": "<string>",
  "biosample_ontology": "<string>",
  "documents": [
    "<value>"
  ],
  "references": [
    "<string>"
  ],
  "schema_version": "<string>",
  "accession": "<string>",
  "alternate_accessions": [
    "<value>"
  ],
  "analyses": [
    "<string>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_get_experiment-parameters)

## 22) 按格式与实验类型列出 ENCODE 文件 {#dnakit.databases.encode_list_files}

- **作用：** 按格式与实验类型列出 ENCODE 文件。
- **处理方法：** 通过 `ENCODEFiles` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_list_files(*, file_type=None, assay_title=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `@graph`, `@type`, `all`, `notification` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_list_files(
    file_type='fastq',
    limit=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "@graph": [
    {
      "@id": "<string>",
      "@type": "<array>",
      "accession": "<string>",
      "assay_term_name": "<string>",
      "award": "<object>",
      "biological_replicates": "<array>",
      "biosample_ontology": "<object>",
      "dataset": "<string>"
    }
  ],
  "@type": [
    "<string>"
  ],
  "all": "<string>",
  "notification": "<string>",
  "title": "<string>",
  "total": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_list_files-parameters)

## 23) 检索顺式调控元件与表观组注释 {#dnakit.databases.encode_search_annotations}

- **作用：** 检索顺式调控元件与表观组注释。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_search_annotations(*, annotation_type=None, biosample_term_name=None, organism=None, assembly=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `annotations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_search_annotations(
    annotation_type='candidate Cis-Regulatory Elements',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "annotations": [
    {
      "accession": "<string>",
      "annotation_type": "<string>",
      "description": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_annotations-parameters)

## 24) 按物种与处理条件检索 ENCODE 样本 {#dnakit.databases.encode_search_biosamples}

- **作用：** 按物种与处理条件检索 ENCODE 样本。
- **处理方法：** 通过 `ENCODESearchBiosamples` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_search_biosamples(*, organism=None, biosample_type=None, treatment=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `@graph`, `@type`, `notification`, `title` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_search_biosamples(
    organism='Homo sapiens',
    limit=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "@graph": [
    {
      "@id": "<string>",
      "@type": "<array>",
      "accession": "<string>",
      "age": "<string>",
      "applied_modifications": "<array>",
      "biosample_ontology": "<object>",
      "lab": "<object>",
      "life_stage": "<string>"
    }
  ],
  "@type": [
    "<string>"
  ],
  "notification": "<string>",
  "title": "<string>",
  "total": "<integer>",
  "all": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_biosamples-parameters)

## 25) 检索 ATAC-seq 与 DNase-seq 实验 {#dnakit.databases.encode_search_chromatin_accessibility}

- **作用：** 检索 ATAC-seq 与 DNase-seq 实验。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_search_chromatin_accessibility(*, assay_type=None, biosample_term_name=None, organism=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `experiments` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_search_chromatin_accessibility(
    assay_type='ATAC-seq',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "experiments": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>",
      "lab": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_chromatin_accessibility-parameters)

## 26) 通过数据库模块检索 ENCODE 实验 {#dnakit.databases.encode_search_experiments}

- **作用：** 通过数据库模块检索 ENCODE 实验。
- **处理方法：** 通过 `ENCODESearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_search_experiments(*, assay_title=None, target=None, organism=None, status=None, limit=None, biosample=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回 ENCODE 实验搜索结果；普通查询在 data 中保留响应，biosample 查询返回解析后的 ENCODE 搜索结果。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_search_experiments(
    assay_title='ATAC-seq',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 ENCODE 实验搜索结果；普通查询在 data 中保留响应，biosample 查询返回解析后的 ENCODE 搜索结果。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_experiments-parameters)

## 27) 按组蛋白修饰检索 ChIP-seq 实验 {#dnakit.databases.encode_search_histone_experiments}

- **作用：** 按组蛋白修饰检索 ChIP-seq 实验。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_search_histone_experiments(*, histone_mark=None, target=None, biosample_term_name=None, biosample=None, organism=None, limit=None, biosample_term=None, cell_type=None, tissue=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `experiments` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_search_histone_experiments(
    histone_mark='H3K27ac',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "experiments": [
    {
      "accession": "<string>",
      "histone_mark": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>",
      "lab": "<string>",
      "date_released": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_histone_experiments-parameters)

## 28) 按样本与实验类型检索 ENCODE RNA-seq {#dnakit.databases.encode_search_rnaseq_experiments}

- **作用：** 按样本与实验类型检索 ENCODE RNA-seq。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.encode_search_rnaseq_experiments(*, biosample_term_name=None, biosample=None, cell_type=None, tissue=None, organism=None, assay_type=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `experiments` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.encode_search_rnaseq_experiments(
    biosample_term_name='K562',
    assay_type='total RNA-seq',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "experiments": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>",
      "lab": "<string>",
      "date_released": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_rnaseq_experiments-parameters)

## 29) 检索差异表达实验的项目记录 {#dnakit.databases.expression_atlas_search_differential}

- **作用：** 检索差异表达实验的项目记录。
- **处理方法：** 通过 `ExpressionAtlas` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.expression_atlas_search_differential(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `gene`, `condition`, `species`, `experiments` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.expression_atlas_search_differential(
    gene='TP53',
    species='homo sapiens',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene": "<string>",
  "condition": "<string>",
  "species": "<string>",
  "experiments": [
    {
      "experiment_accession": "<string>",
      "experiment_type": "<string>",
      "experiment_description": "<string>",
      "species": "<string>",
      "num_assays": "<integer>",
      "gene_mentioned": "<boolean>",
      "factors": "<array>"
    }
  ],
  "experiment_count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.expression_atlas_search_differential-parameters)

## 30) 检索基线与差异表达实验 {#dnakit.databases.expression_atlas_search_experiments}

- **作用：** 检索基线与差异表达实验。
- **处理方法：** 通过 `ExpressionAtlas` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.expression_atlas_search_experiments(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `gene`, `condition`, `species`, `experiments` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.expression_atlas_search_experiments(
    gene='TP53',
    species='homo sapiens',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene": "<string>",
  "condition": "<string>",
  "species": "<string>",
  "experiments": [
    {
      "experiment_accession": "<string>",
      "experiment_type": "<string>",
      "experiment_description": "<string>",
      "species": "<string>",
      "num_assays": "<integer>",
      "gene_mentioned": "<boolean>"
    }
  ],
  "total_count": "<integer>",
  "gene_specific_count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.expression_atlas_search_experiments-parameters)

## 31) 检索 Hi-C 与染色质构象数据 {#dnakit.databases.four_dn_search_data}

- **作用：** 检索 Hi-C 与染色质构象数据。
- **处理方法：** 通过 `FourDN` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.four_dn_search_data(*, operation=None, query=None, item_type=None, file_type=None, assay_title=None, biosource_name=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `num_results`, `total`, `results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.four_dn_search_data(
    operation='search',
    query='Hi-C',
    item_type='File',
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
  "num_results": "<integer>",
  "total": "<integer>",
  "results": [
    "<value>"
  ],
  "search_url": "<string>",
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.four_dn_search_data-parameters)

## 32) 查询 GDC 肿瘤病例的临床资料 {#dnakit.databases.gdc_get_clinical_data}

- **作用：** 查询 GDC 肿瘤病例的临床资料。
- **处理方法：** 通过 `GDCClinicalData` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gdc_get_clinical_data(*, project_id=None, primary_site=None, disease_type=None, vital_status=None, gender=None, size=None, offset=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** data.cases 包含病例、诊断和治疗记录，data.pagination 包含分页信息。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gdc_get_clinical_data(
    project_id='TCGA-BRCA',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data.cases 包含病例、诊断和治疗记录，data.pagination 包含分页信息。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_clinical_data-parameters)

## 33) 查询基因在泛癌队列中的突变频率 {#dnakit.databases.gdc_get_mutation_frequency}

- **作用：** 查询基因在泛癌队列中的突变频率。
- **处理方法：** 通过 `GDCMutationFrequency` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gdc_get_mutation_frequency(*, gene_symbol=None, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `source`, `gene`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gdc_get_mutation_frequency(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "success",
  "source": "GDC",
  "gene": "<str>",
  "data": {
    "gene_info": "<dict>",
    "total_ssm_occurrences": "<int>",
    "is_cancer_gene_census": "<bool | None>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_mutation_frequency-parameters)

## 34) 按基因查询 GDC 体细胞突变 {#dnakit.databases.gdc_get_ssm_by_gene}

- **作用：** 按基因查询 GDC 体细胞突变。
- **处理方法：** 通过 `GDCSSM` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gdc_get_ssm_by_gene(*, gene_symbol, project_id=None, size=None, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** data 包含 GDC 体细胞变异记录，gene 标明查询基因。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gdc_get_ssm_by_gene(
    gene_symbol='TP53',
    project_id='TCGA-BRCA',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data 包含 GDC 体细胞变异记录，gene 标明查询基因。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_ssm_by_gene-parameters)

## 35) 查询肿瘤队列的生存与删失数据 {#dnakit.databases.gdc_get_survival}

- **作用：** 查询肿瘤队列的生存与删失数据。
- **处理方法：** 通过 `GDCSurvival` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gdc_get_survival(*, project_id, gene_symbol=None, verbose=True)`
- **输入：** 必须提供 `project_id`。
- **输出：** data 包含供体随访时间、删失状态、生存估计和队列汇总；donors 最多展示 50 条。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gdc_get_survival(
    project_id='TCGA-BRCA',
    verbose=False,
)
print(result)
```

- **示例结果：** data 包含供体随访时间、删失状态、生存估计和队列汇总；donors 最多展示 50 条。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_survival-parameters)

## 36) 按数据类型列出 GDC 文件 {#dnakit.databases.gdc_list_files}

- **作用：** 按数据类型列出 GDC 文件。
- **处理方法：** 通过 `GDCFiles` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gdc_list_files(*, data_type=None, size=None, offset=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** data 中保留 GDC 文件列表及分页响应，记录含文件标识与元数据。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gdc_list_files(
    data_type='Gene Expression Quantification',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data 中保留 GDC 文件列表及分页响应，记录含文件标识与元数据。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_list_files-parameters)

## 37) 列出 GDC 肿瘤项目与汇总信息 {#dnakit.databases.gdc_list_projects}

- **作用：** 列出 GDC 肿瘤项目与汇总信息。
- **处理方法：** 通过 `GDCProjects` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gdc_list_projects(*, program=None, size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** data 中保留 GDC 项目列表及项目统计。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gdc_list_projects(
    program='TCGA',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data 中保留 GDC 项目列表及项目统计。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_list_projects-parameters)

## 38) 按筛选条件检索 GDC 肿瘤病例 {#dnakit.databases.gdc_search_cases}

- **作用：** 按筛选条件检索 GDC 肿瘤病例。
- **处理方法：** 通过 `GDCCases` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gdc_search_cases(*, project_id=None, size=None, offset=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** data 中保留 GDC 病例查询响应及分页信息。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gdc_search_cases(
    project_id='TCGA-BRCA',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data 中保留 GDC 病例查询响应及分页信息。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_search_cases-parameters)

## 39) 检索 GEO 染色质可及性数据集 {#dnakit.databases.geo_search_atacseq_datasets}

- **作用：** 检索 GEO 染色质可及性数据集。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.geo_search_atacseq_datasets(*, query=None, organism=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `datasets` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.geo_search_atacseq_datasets(
    query='T cell differentiation',
    organism='Homo sapiens',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "datasets": [
    {
      "accession": "<string>",
      "title": "<string>",
      "summary": "<string>",
      "organism": "<string>",
      "n_samples": "<integer>",
      "date_published": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.geo_search_atacseq_datasets-parameters)

## 40) 检索 GEO 蛋白质 DNA 结合数据集 {#dnakit.databases.geo_search_chipseq_datasets}

- **作用：** 检索 GEO 蛋白质 DNA 结合数据集。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.geo_search_chipseq_datasets(*, query, organism=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `total`, `datasets` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.geo_search_chipseq_datasets(
    query='H3K27ac liver',
    organism='Homo sapiens',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "datasets": [
    {
      "accession": "<string>",
      "title": "<string>",
      "summary": "<string>",
      "organism": "<string>",
      "n_samples": "<integer>",
      "date_published": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.geo_search_chipseq_datasets-parameters)

## 41) 检索 GEO 转录组表达数据集 {#dnakit.databases.geo_search_rnaseq_datasets}

- **作用：** 检索 GEO 转录组表达数据集。
- **处理方法：** 通过 `Epigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.geo_search_rnaseq_datasets(*, query=None, organism=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `total`, `datasets` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.geo_search_rnaseq_datasets(
    query='breast cancer BRCA1',
    organism='Homo sapiens',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total": "<integer>",
  "datasets": [
    {
      "accession": "<string>",
      "title": "<string>",
      "summary": "<string>",
      "organism": "<string>",
      "n_samples": "<integer>",
      "date_published": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.geo_search_rnaseq_datasets-parameters)

## 42) 计算指定组织中基因与变异的 eQTL 关联 {#dnakit.databases.gtex_calculate_eqtl}

- **作用：** 计算指定组织中基因与变异的 eQTL 关联。
- **处理方法：** 通过 `GTExV2` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtex_calculate_eqtl(*, operation=None, gencode_id, variant_id, tissue_site_detail_id, dataset_id=None, verbose=True)`
- **输入：** 必须提供 `gencode_id`, `variant_id`, `tissue_site_detail_id`。
- **输出：** 返回字典，包含 `data`, `error`, `gencodeId`, `geneSymbol` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtex_calculate_eqtl(
    operation='calculate_eqtl',
    gencode_id='ENSG00000227232.5',
    variant_id='chr1_64764_C_T_b38',
    tissue_site_detail_id='Whole_Blood',
    dataset_id='gtex_v8',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": [
    "<number>"
  ],
  "error": "<number>",
  "gencodeId": "<string>",
  "geneSymbol": "<string>",
  "genotypes": [
    "<integer>"
  ],
  "hetCount": "<integer>",
  "homoAltCount": "<integer>",
  "homoRefCount": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_calculate_eqtl-parameters)

## 43) 查询基因的 GTEx 表达摘要 {#dnakit.databases.gtex_get_expression_summary}

- **作用：** 查询基因的 GTEx 表达摘要。
- **处理方法：** 通过 `GTExExpression` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtex_get_expression_summary(*, gene_symbol=None, ensembl_gene_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回解析到的基因标识及组织表达数据；无表达记录时会附 note。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtex_get_expression_summary(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回解析到的基因标识及组织表达数据；无表达记录时会附 note。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_expression_summary-parameters)

## 44) 查询基因在不同组织中的表达中位数 {#dnakit.databases.gtex_get_median_gene_expression}

- **作用：** 查询基因在不同组织中的表达中位数。
- **处理方法：** 通过 `GTExV2` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtex_get_median_gene_expression(*, operation=None, gencode_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `median`, `tissueSiteDetailId`, `ontologyId`, `datasetId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtex_get_median_gene_expression(
    operation='get_median_gene_expression',
    gencode_id='ENSG00000141510.18',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "median": "<number>",
    "tissueSiteDetailId": "<string>",
    "ontologyId": "<string>",
    "datasetId": "<string>",
    "gencodeId": "<string>",
    "geneSymbol": "<string>",
    "unit": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_median_gene_expression-parameters)

## 45) 查询跨组织 eQTL 荟萃分析结果 {#dnakit.databases.gtex_get_multi_tissue_eqtls}

- **作用：** 查询跨组织 eQTL 荟萃分析结果。
- **处理方法：** 通过 `GTExV2` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtex_get_multi_tissue_eqtls(*, operation=None, gencode_id, variant_id=None, dataset_id=None, page=None, items_per_page=None, verbose=True)`
- **输入：** 必须提供 `gencode_id`。
- **输出：** 返回记录列表；每条记录包含 `gencodeId`, `datasetId`, `metaP`, `variantId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtex_get_multi_tissue_eqtls(
    operation='get_multi_tissue_eqtls',
    gencode_id='ENSG00000141510.16',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "gencodeId": "<string>",
    "datasetId": "<string>",
    "metaP": "<number>",
    "variantId": "<string>",
    "tissues": {
      "Thyroid": "<object>",
      "Testis": "<object>",
      "Small_Intestine_Terminal_Ileum": "<object>",
      "Brain_Frontal_Cortex_BA9": "<object>",
      "Skin_Not_Sun_Exposed_Suprapubic": "<object>",
      "Vagina": "<object>",
      "Whole_Blood": "<object>",
      "Breast_Mammary_Tissue": "<object>"
    }
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_multi_tissue_eqtls-parameters)

## 46) 查询显著的单组织 eQTL 关联 {#dnakit.databases.gtex_get_single_tissue_eqtls}

- **作用：** 查询显著的单组织 eQTL 关联。
- **处理方法：** 通过 `GTExV2` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtex_get_single_tissue_eqtls(*, operation=None, gencode_id=None, variant_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `snpId`, `pos`, `snpIdUpper`, `variantId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtex_get_single_tissue_eqtls(
    operation='get_single_tissue_eqtls',
    gencode_id=['ENSG00000227232.5'],
    tissue_site_detail_id=['Whole_Blood'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "snpId": "<string>",
    "pos": "<integer>",
    "snpIdUpper": "<string>",
    "variantId": "<string>",
    "geneSymbol": "<string>",
    "pValue": "<number>",
    "geneSymbolUpper": "<string>",
    "datasetId": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_single_tissue_eqtls-parameters)

## 47) 列出 GTEx 组织类型与样本统计 {#dnakit.databases.gtex_get_tissue_sites}

- **作用：** 列出 GTEx 组织类型与样本统计。
- **处理方法：** 通过 `GTExV2` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtex_get_tissue_sites(*, operation=None, dataset_id=None, page=None, items_per_page=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `tissueSiteDetailId`, `colorHex`, `colorRgb`, `datasetId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtex_get_tissue_sites(
    operation='get_tissue_sites',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "tissueSiteDetailId": "<string>",
    "colorHex": "<string>",
    "colorRgb": "<string>",
    "datasetId": "<string>",
    "eGeneCount": "<integer>",
    "expressedGeneCount": "<integer>",
    "hasEGenes": "<boolean>",
    "hasSGenes": "<boolean>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_tissue_sites-parameters)

## 48) 按条件查询 GTEx 单组织 eQTL {#dnakit.databases.gtex_query_eqtl}

- **作用：** 按条件查询 GTEx 单组织 eQTL。
- **处理方法：** 通过 `GTExEQTL` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtex_query_eqtl(*, gene_symbol=None, ensembl_gene_id=None, page=None, size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** data.singleTissueEqtl 包含单组织 eQTL 记录，metadata 记录查询来源与条件。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtex_query_eqtl(
    gene_symbol='TP53',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data.singleTissueEqtl 包含单组织 eQTL 记录，metadata 记录查询来源与条件。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_query_eqtl-parameters)

## 49) 汇总多组学数据中的基因信息 {#dnakit.databases.harmonizome_get_gene}

- **作用：** 汇总多组学数据中的基因信息。
- **处理方法：** 通过 `Harmonizome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.harmonizome_get_gene(*, gene_symbol, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** 返回字典，包含 `symbol`, `name`, `ncbi_entrez_gene_id`, `ncbi_entrez_gene_url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.harmonizome_get_gene(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "symbol": "<string>",
  "name": "<string>",
  "ncbi_entrez_gene_id": "<integer>",
  "ncbi_entrez_gene_url": "<string>",
  "description": "<string>",
  "synonyms": [
    "<string>"
  ],
  "proteins": [
    {
      "symbol": "<string>",
      "href": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.harmonizome_get_gene-parameters)

## 50) 查询基因表达与肿瘤生存的关联 {#dnakit.databases.hpa_get_cancer_prognostics_by_gene}

- **作用：** 查询基因表达与肿瘤生存的关联。
- **处理方法：** 通过 `HPAGetCancerPrognostics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpa_get_cancer_prognostics_by_gene(*, ensembl_id=None, gene_name=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `ensembl_id`, `gene`, `gene_synonym`, `cancers_assessed_count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpa_get_cancer_prognostics_by_gene(
    ensembl_id='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "ensembl_id": "<string>",
  "gene": "<string>",
  "gene_synonym": [
    "<string>"
  ],
  "cancers_assessed_count": "<integer>",
  "prognostic_cancers_count": "<integer>",
  "non_prognostic_cancers_count": "<integer>",
  "prognostic_summary": [
    {
      "cancer_type": "<string>",
      "prognostic_type": "<string>",
      "p_value": "<string>",
      "is_prognostic": "<boolean>"
    }
  ],
  "non_prognostic_summary": [
    {
      "cancer_type": "<string>",
      "prognostic_type": "<string>",
      "p_value": "<string>",
      "is_prognostic": "<boolean>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_cancer_prognostics_by_gene-parameters)

## 51) 查看已停用的 HPA 蛋白互作查询 {#dnakit.databases.hpa_get_protein_interactions_by_gene}

- **作用：** 查看已停用的 HPA 蛋白互作查询。
- **处理方法：** 通过 `HPAGetProteinInteractions` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpa_get_protein_interactions_by_gene(*, gene_name, verbose=True)`
- **输入：** 必须提供 `gene_name`。 该数据源入口当前不可用，调用会返回错误；以下代码仅说明调用形式。
- **输出：** 此入口当前报告服务不可用，不返回有效互作结果。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpa_get_protein_interactions_by_gene(
    gene_name='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 此入口当前报告服务不可用，不返回有效互作结果。

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_protein_interactions_by_gene-parameters)

## 52) 按组织或细胞来源查询 RNA 表达 {#dnakit.databases.hpa_get_rna_expression_by_source}

- **作用：** 按组织或细胞来源查询 RNA 表达。
- **处理方法：** 通过 `HPAGetRnaExpressionBySource` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpa_get_rna_expression_by_source(*, gene_name, source_type, source_name, verbose=True)`
- **输入：** 必须提供 `gene_name`, `source_type`, `source_name`。
- **输出：** 返回字典，包含 `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpa_get_rna_expression_by_source(
    gene_name='TP53',
    source_type='tissue',
    source_name='liver',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_rna_expression_by_source-parameters)

## 53) 查询基因在指定组织中的 RNA 表达 {#dnakit.databases.hpa_get_rna_expression_in_specific_tissues}

- **作用：** 查询基因在指定组织中的 RNA 表达。
- **处理方法：** 通过 `HPAGetRnaExpressionByTissue` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpa_get_rna_expression_in_specific_tissues(*, ensembl_id=None, tissue_names=None, tissues=None, gene_name=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpa_get_rna_expression_in_specific_tissues(
    ensembl_id='ENSG00000254647',
    tissue_names=['pancreas'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_rna_expression_in_specific_tissues-parameters)

## 54) 查询蛋白质亚细胞定位注释 {#dnakit.databases.hpa_get_subcellular_location}

- **作用：** 查询蛋白质亚细胞定位注释。
- **处理方法：** 通过 `HPAGetSubcellularLocation` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpa_get_subcellular_location(*, gene_name, verbose=True)`
- **输入：** 必须提供 `gene_name`。
- **输出：** 返回字典，包含 `gene_name`, `gene_synonym`, `main_locations`, `additional_locations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpa_get_subcellular_location(
    gene_name='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_name": "<string>",
  "gene_synonym": [
    "<string>"
  ],
  "main_locations": [
    "<string>"
  ],
  "additional_locations": [
    "<string>"
  ],
  "total_locations": "<integer>",
  "location_summary": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_subcellular_location-parameters)

## 55) 检索 HPA 基因与 Ensembl 标识 {#dnakit.databases.hpa_search_genes_by_query}

- **作用：** 检索 HPA 基因与 Ensembl 标识。
- **处理方法：** 通过 `HPASearchGenes` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hpa_search_genes_by_query(*, search_query, max_results=None, verbose=True)`
- **输入：** 必须提供 `search_query`。
- **输出：** 返回字典，包含 `error` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hpa_search_genes_by_query(
    search_query='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "error": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_search_genes_by_query-parameters)

## 56) 获取转录因子结合矩阵与元数据 {#dnakit.databases.jaspar_get_matrix}

- **作用：** 获取转录因子结合矩阵与元数据。
- **处理方法：** 通过 `JASPAR` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.jaspar_get_matrix(*, matrix_id, verbose=True)`
- **输入：** 必须提供 `matrix_id`。
- **输出：** 返回字典，包含 `matrix_id`, `base_id`, `name`, `collection` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.jaspar_get_matrix(
    matrix_id='MA0002.2',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "matrix_id": "<string>",
  "base_id": "<string>",
  "name": "<string>",
  "collection": "<string>",
  "version": "<string>",
  "pfm": {
    "A": [
      "<number>"
    ],
    "C": [
      "<number>"
    ],
    "G": [
      "<number>"
    ],
    "T": [
      "<number>"
    ]
  },
  "tax_group": "<string>",
  "species": [
    {
      "tax_id": "<integer>",
      "name": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_get_matrix-parameters)

## 57) 列出 JASPAR 结合矩阵的历史版本 {#dnakit.databases.jaspar_get_matrix_versions}

- **作用：** 列出 JASPAR 结合矩阵的历史版本。
- **处理方法：** 通过 `JASPAR` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.jaspar_get_matrix_versions(*, base_id, page=None, page_size=None, verbose=True)`
- **输入：** 必须提供 `base_id`。
- **输出：** 返回字典，包含 `count`, `next`, `previous`, `results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.jaspar_get_matrix_versions(
    base_id='MA0002',
    page=1,
    page_size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "count": "<integer>",
  "next": "<string>",
  "previous": "<string>",
  "results": [
    {
      "matrix_id": "<string>",
      "base_id": "<string>",
      "name": "<string>",
      "collection": "<string>",
      "version": "<string>",
      "url": "<string>",
      "sequence_logo": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_get_matrix_versions-parameters)

## 58) 分页列出 JASPAR 转录因子结合谱 {#dnakit.databases.jaspar_get_transcription_factors}

- **作用：** 分页列出 JASPAR 转录因子结合谱。
- **处理方法：** 通过 `JASPAR` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.jaspar_get_transcription_factors(*, collection=None, page=None, page_size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `count`, `next`, `previous`, `results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.jaspar_get_transcription_factors(
    collection='CORE',
    page=1,
    page_size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "count": "<integer>",
  "next": "<string>",
  "previous": "<string>",
  "results": [
    {
      "matrix_id": "<string>",
      "base_id": "<string>",
      "name": "<string>",
      "collection": "<string>",
      "version": "<string>",
      "url": "<string>",
      "sequence_logo": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_get_transcription_factors-parameters)

## 59) 按名称、物种或分类检索结合矩阵 {#dnakit.databases.jaspar_search_matrices}

- **作用：** 按名称、物种或分类检索结合矩阵。
- **处理方法：** 通过 `JASPAR` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.jaspar_search_matrices(*, search=None, name=None, collection=None, tax_group=None, species=None, page=None, page_size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `count`, `next`, `previous`, `results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.jaspar_search_matrices(
    search='CTCF',
    page=1,
    page_size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "count": "<integer>",
  "next": "<string>",
  "previous": "<string>",
  "results": [
    {
      "matrix_id": "<string>",
      "base_id": "<string>",
      "name": "<string>",
      "collection": "<string>",
      "version": "<string>",
      "url": "<string>",
      "sequence_logo": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_search_matrices-parameters)

## 60) 查询基因与细胞类型对应的转录因子结合位点 {#dnakit.databases.re_map_get_transcription_factor_binding}

- **作用：** 查询基因与细胞类型对应的转录因子结合位点。
- **处理方法：** 通过 `ReMap` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.re_map_get_transcription_factor_binding(*, gene_name, cell_type=None, limit=None, verbose=True)`
- **输入：** 必须提供 `gene_name`。
- **输出：** 返回字典，包含 `experiments`, `gene_name`, `cell_type`, `count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.re_map_get_transcription_factor_binding(
    gene_name='BRCA1',
    cell_type='HepG2',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "experiments": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "target": "<object>",
      "biosample_ontology": "<object>",
      "description": "<string>",
      "status": "<string>"
    }
  ],
  "gene_name": "<string>",
  "cell_type": "<string>",
  "count": "<integer>",
  "url": "<string>",
  "note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.re_map_get_transcription_factor_binding-parameters)

## 61) 按 rsID 查询变异调控注释 {#dnakit.databases.regulome_db_query_variant}

- **作用：** 按 rsID 查询变异调控注释。
- **处理方法：** 通过 `RegulomeDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.regulome_db_query_variant(*, rsid=None, variant=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `rsid`, `assembly`, `query_coordinates`, `regulome_score` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.regulome_db_query_variant(
    rsid='rs4994',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "rsid": "<string>",
  "assembly": "<string>",
  "query_coordinates": [
    "<string>"
  ],
  "regulome_score": {
    "probability": "<string>",
    "ranking": "<string>",
    "tissue_specific_scores": "<object>"
  },
  "variants": [
    {
      "chrom": "<string>",
      "start": "<integer>",
      "end": "<integer>",
      "rsids": "<array>"
    }
  ],
  "features": "<object>",
  "nearby_snps": [
    "<value>"
  ],
  "notifications": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.regulome_db_query_variant-parameters)

## 62) 查询基因组区间的候选顺式调控元件 {#dnakit.databases.screen_get_regulatory_elements}

- **作用：** 查询基因组区间的候选顺式调控元件。
- **处理方法：** 通过 `SCREEN` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.screen_get_regulatory_elements(*, gene_name=None, query=None, element_type=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `gene_name`, `element_type`, `count`, `regulatory_elements` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.screen_get_regulatory_elements(
    gene_name='BRCA1',
    element_type='enhancer',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_name": "<string>",
  "element_type": "<string>",
  "count": "<integer>",
  "regulatory_elements": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "target": "<string>",
      "biosample": "<string>",
      "status": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.screen_get_regulatory_elements-parameters)

## 63) 从 UCSC 获取 ENCODE 候选顺式调控元件 {#dnakit.databases.ucsc_get_encode_c_cres}

- **作用：** 从 UCSC 获取 ENCODE 候选顺式调控元件。
- **处理方法：** 通过 `UCSCEpigenomics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ucsc_get_encode_c_cres(*, genome=None, chrom, start, end, verbose=True)`
- **输入：** 必须提供 `chrom`, `start`, `end`。
- **输出：** 返回字典，包含 `genome`, `region`, `ccre_count`, `ccres` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ucsc_get_encode_c_cres(
    genome='hg38',
    chrom='chr17',
    start=7668421,
    end=7687490,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "genome": "<string>",
  "region": "<string>",
  "ccre_count": "<integer>",
  "ccres": [
    {
      "name": "<string>",
      "chrom": "<string>",
      "chromStart": "<integer>",
      "chromEnd": "<integer>",
      "cCRE_class": "<string>",
      "DNase_maxZ": "<number>",
      "H3K4me3_maxZ": "<number>",
      "H3K27ac_maxZ": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.ucsc_get_encode_c_cres-parameters)

## 64) 查询单个直接转录因子结合数据集 {#dnakit.databases.uni_bind_get_dataset}

- **作用：** 查询单个直接转录因子结合数据集。
- **处理方法：** 通过 `UniBind` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_bind_get_dataset(*, operation=None, dataset_id, verbose=True)`
- **输入：** 必须提供 `dataset_id`。
- **输出：** 返回字典，包含 `dataset_id`, `tf_name`, `cell_line`, `biological_condition` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_bind_get_dataset(
    operation='get_dataset',
    dataset_id='EXP030726.neural_stem_cells.SMAD3',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "dataset_id": "<string>",
  "tf_name": "<string>",
  "cell_line": [
    "<string>"
  ],
  "biological_condition": [
    "<string>"
  ],
  "identifier": [
    "<string>"
  ],
  "jaspar_id": [
    "<string>"
  ],
  "prediction_models": [
    "<string>"
  ],
  "tfbs_models": [
    {
      "prediction_model": "<string>",
      "jaspar_id": "<string>",
      "jaspar_version": "<string>",
      "total_tfbs": "<string>",
      "score_threshold": "<string>",
      "distance_threshold": "<string>",
      "adj_centrimo_pvalue": "<number>",
      "bed_url": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.uni_bind_get_dataset-parameters)

## 65) 列出 UniBind 收录的转录因子 {#dnakit.databases.uni_bind_list_tfs}

- **作用：** 列出 UniBind 收录的转录因子。
- **处理方法：** 通过 `UniBind` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_bind_list_tfs(*, operation=None, search=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回“列出 UniBind 收录的转录因子”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_bind_list_tfs(
    operation='list_tfs',
    search='SMAD',
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

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.uni_bind_list_tfs-parameters)

## 66) 按转录因子与细胞条件检索结合数据集 {#dnakit.databases.uni_bind_search_datasets}

- **作用：** 按转录因子与细胞条件检索结合数据集。
- **处理方法：** 通过 `UniBind` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_bind_search_datasets(*, operation=None, tf_name=None, species=None, cell_line=None, collection=None, order=None, page=None, page_size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `tf_name`, `total_peaks`, `dataset_url`, `dataset_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_bind_search_datasets(
    operation='search_datasets',
    tf_name='SMAD3',
    page_size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "tf_name": "<string>",
    "total_peaks": "<string>",
    "dataset_url": "<string>",
    "dataset_id": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_regulation.md#dnakit.databases.uni_bind_search_datasets-parameters)

<!-- END GENERATED API: regulation -->
