# 通路、互作、蛋白与结构查询

查询通路、分子互作、蛋白和结构数据。

<!-- BEGIN GENERATED API: pathways -->

<span id="pathways"></span>

## 1) 获取已收录的 AlphaFold 蛋白质结构预测 {#dnakit.databases.alphafold_get_prediction}

- **作用：** 获取已收录的 AlphaFold 蛋白质结构预测。
- **处理方法：** 通过 `AlphaFold` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.alphafold_get_prediction(*, qualifier=None, uniprot_id=None, uniprot_accession=None, accession=None, sequence_checksum=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `toolUsed`, `providerId`, `entityType`, `isUniProt` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.alphafold_get_prediction(
    qualifier='P69905',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "toolUsed": "<string>",
    "providerId": "<string>",
    "entityType": "<string>",
    "isUniProt": "<boolean>",
    "modelEntityId": "<string>",
    "modelCreatedDate": "<string>",
    "sequenceVersionDate": "<string>",
    "globalMetricValue": "<number>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.alphafold_get_prediction-parameters)

## 2) 查询实验验证的蛋白质与遗传互作 {#dnakit.databases.bio_grid_get_interactions}

- **作用：** 查询实验验证的蛋白质与遗传互作。
- **处理方法：** 通过 `BioGRID` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.bio_grid_get_interactions(*, gene_names=None, gene_name=None, organism=None, interaction_type=None, evidence_types=None, limit=None, throughput=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 需配置：`BIOGRID_ACCESS_KEY`.
- **输出：** 返回以记录键为索引的字典；每条记录包含 `BIOGRID_INTERACTION_ID`, `ENTREZ_GENE_A`, `ENTREZ_GENE_B`, `OFFICIAL_SYMBOL_A` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.bio_grid_get_interactions(
    gene_names=['TP53'],
    organism='9606',
    interaction_type='physical',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "<record key>": {
    "BIOGRID_INTERACTION_ID": "<integer>",
    "ENTREZ_GENE_A": "<string>",
    "ENTREZ_GENE_B": "<string>",
    "OFFICIAL_SYMBOL_A": "<string>",
    "OFFICIAL_SYMBOL_B": "<string>",
    "SYNONYMS_A": "<string>",
    "SYNONYMS_B": "<string>",
    "EXPERIMENTAL_SYSTEM": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.bio_grid_get_interactions-parameters)

## 3) 检索病原体 B 细胞与 T 细胞表位 {#dnakit.databases.bvbrc_search_epitopes}

- **作用：** 检索病原体 B 细胞与 T 细胞表位。
- **处理方法：** 通过 `BVBRC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.bvbrc_search_epitopes(*, taxon_id=None, protein_name=None, epitope_type=None, organism=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `epitope_id`, `epitope_type`, `epitope_sequence`, `organism` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.bvbrc_search_epitopes(
    taxon_id='2697049',
    protein_name='Spike glycoprotein',
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
    "epitope_id": "<string>",
    "epitope_type": "<string>",
    "epitope_sequence": "<string>",
    "organism": "<string>",
    "protein_name": "<string>",
    "start": "<integer>",
    "end": "<integer>",
    "bcell_assays": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.bvbrc_search_epitopes-parameters)

## 4) 查询基因或蛋白质的 GO 注释 {#dnakit.databases.go_get_annotations_for_gene}

- **作用：** 查询基因或蛋白质的 GO 注释。
- **处理方法：** 通过 `GeneOntology` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.go_get_annotations_for_gene(*, gene_id, rows=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `document_category`, `source`, `bioentity`, `bioentity_internal_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.go_get_annotations_for_gene(
    gene_id='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "document_category": "<string>",
    "source": "<string>",
    "bioentity": "<string>",
    "bioentity_internal_id": "<string>",
    "bioentity_label": "<string>",
    "bioentity_label_searchable": "<string>",
    "annotation_class": "<string>",
    "annotation_class_label": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.go_get_annotations_for_gene-parameters)

## 5) 按关键词检索 GO 术语 {#dnakit.databases.go_search_terms}

- **作用：** 按关键词检索 GO 术语。
- **处理方法：** 通过 `GeneOntology` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.go_search_terms(*, query, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `responseHeader`, `response` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.go_search_terms(
    query='gene:TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "responseHeader": {
    "status": "<integer>",
    "QTime": "<integer>",
    "params": {
      "q": "<string>",
      "fq": "<array>",
      "rows": "<string>",
      "wt": "<string>"
    }
  },
  "response": {
    "numFound": "<integer>",
    "start": "<integer>",
    "docs": [
      "<value>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.go_search_terms-parameters)

## 6) 获取组织特异性蛋白互作网络 {#dnakit.databases.humanbase_ppi_analysis}

- **作用：** 获取组织特异性蛋白互作网络。
- **处理方法：** 通过 `HumanBase` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.humanbase_ppi_analysis(*, gene_list=None, tissue=None, max_node=None, interaction=None, string_mode=None, genes=None, top_n=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字符串。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.humanbase_ppi_analysis(
    gene_list=['TP53', 'BRCA1'],
    tissue='brain',
    max_node=3,
    interaction='co-expression',
    string_mode=True,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
"<string>"
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.humanbase_ppi_analysis-parameters)

## 7) 查询免疫表位关联的抗原 {#dnakit.databases.iedb_get_epitope_antigens}

- **作用：** 查询免疫表位关联的抗原。
- **处理方法：** 通过 `IEDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.iedb_get_epitope_antigens(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **输入：** 必须提供 `structure_id`。
- **输出：** 返回记录列表；每条记录包含 `structure_id`, `parent_source_antigen_iri` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.iedb_get_epitope_antigens(
    structure_id=2,
    limit=2,
    select=['structure_id', 'parent_source_antigen_iri'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "structure_id": "<integer>",
    "parent_source_antigen_iri": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_get_epitope_antigens-parameters)

## 8) 查询免疫表位的 MHC 限制性关联 {#dnakit.databases.iedb_get_epitope_mhc}

- **作用：** 查询免疫表位的 MHC 限制性关联。
- **处理方法：** 通过 `IEDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.iedb_get_epitope_mhc(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **输入：** 必须提供 `structure_id`。
- **输出：** 返回记录列表；每条记录包含 `structure_id`, `elution_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.iedb_get_epitope_mhc(
    structure_id=2,
    limit=2,
    select=['structure_id', 'elution_id'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "structure_id": "<integer>",
    "elution_id": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_get_epitope_mhc-parameters)

## 9) 检索 B 细胞表位实验记录 {#dnakit.databases.iedb_search_bcell}

- **作用：** 检索 B 细胞表位实验记录。
- **处理方法：** 通过 `IEDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.iedb_search_bcell(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `bcell_id`, `structure_id`, `assay_description`, `reference_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.iedb_search_bcell(
    limit=2,
    select=['bcell_id', 'structure_id', 'assay_description', 'reference_id'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "bcell_id": "<integer>",
    "structure_id": "<integer>",
    "assay_description": "<string>",
    "reference_id": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_bcell-parameters)

## 10) 检索免疫表位及其结构编号 {#dnakit.databases.iedb_search_epitopes}

- **作用：** 检索免疫表位及其结构编号。
- **处理方法：** 通过 `IEDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.iedb_search_epitopes(*, sequence_contains=None, structure_type=None, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `structure_id`, `structure_type`, `linear_sequence` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.iedb_search_epitopes(
    sequence_contains='KVF',
    structure_type='Linear peptide',
    limit=2,
    select=['structure_id', 'structure_type', 'linear_sequence'],
    order='structure_id.asc',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "structure_id": "<integer>",
    "structure_type": "<string>",
    "linear_sequence": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_epitopes-parameters)

## 11) 检索 MHC 相关实验记录 {#dnakit.databases.iedb_search_mhc}

- **作用：** 检索 MHC 相关实验记录。
- **处理方法：** 通过 `IEDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.iedb_search_mhc(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `structure_id`, `mhc_restriction`, `mhc_class`, `mhc_allele_resolution` 等字段。
- **示例代码：**

```python
import dnakit

example_select = ['structure_id',
 'mhc_restriction',
 'mhc_class',
 'mhc_allele_resolution',
 'reference_id']

result = dnakit.databases.iedb_search_mhc(
    limit=2,
    select=example_select,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "structure_id": "<integer>",
    "mhc_restriction": "<string>",
    "mhc_class": "<string>",
    "mhc_allele_resolution": "<string>",
    "reference_id": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_mhc-parameters)

## 12) 检索 T 细胞免疫反应实验 {#dnakit.databases.iedb_search_tcell_assays}

- **作用：** 检索 T 细胞免疫反应实验。
- **处理方法：** 通过 `IEDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.iedb_search_tcell_assays(*, sequence=None, sequence_contains=None, mhc_class=None, qualitative_measure=None, limit=None, offset=None, select=None, filters=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `tcell_id`, `linear_sequence`, `structure_type`, `mhc_restriction` 等字段。
- **示例代码：**

```python
import dnakit

example_select = ['tcell_id',
 'linear_sequence',
 'mhc_restriction',
 'qualitative_measure',
 'pubmed_id']

result = dnakit.databases.iedb_search_tcell_assays(
    sequence='GILGFVFTL',
    limit=3,
    select=example_select,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "tcell_id": "<integer>",
    "linear_sequence": "<string>",
    "structure_type": "<string>",
    "mhc_restriction": "<string>",
    "mhc_class": "<string>",
    "qualitative_measure": "<string>",
    "pubmed_id": "<string>",
    "journal_name": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_tcell_assays-parameters)

## 13) 查询免疫遗传学资源与基因命名规则 {#dnakit.databases.imgt_get_gene_info}

- **作用：** 查询免疫遗传学资源与基因命名规则。
- **处理方法：** 通过 `IMGT` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.imgt_get_gene_info(*, operation=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.imgt_get_gene_info(
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "data": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.imgt_get_gene_info-parameters)

## 14) 定位免疫球蛋白与 TCR 胚系基因检索入口 {#dnakit.databases.imgt_search_genes}

- **作用：** 定位免疫球蛋白与 TCR 胚系基因检索入口。
- **处理方法：** 通过 `IMGT` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.imgt_search_genes(*, operation=None, query=None, gene_type=None, species=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.imgt_search_genes(
    gene_type='IGHV',
    species='Homo sapiens',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "data": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.imgt_search_genes-parameters)

## 15) 按疾病与实验方法检索免疫学研究 {#dnakit.databases.imm_port_search_studies}

- **作用：** 按疾病与实验方法检索免疫学研究。
- **处理方法：** 通过 `ImmPort` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.imm_port_search_studies(*, query, condition_or_disease=None, assay_method=None, research_focus=None, species=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `study_accession`, `title`, `brief_description`, `condition_or_disease` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.imm_port_search_studies(
    query='influenza vaccine',
    species='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "study_accession": "<string>",
    "title": "<string>",
    "brief_description": "<string>",
    "condition_or_disease": [
      "<value>"
    ],
    "research_focus": [
      "<value>"
    ],
    "species": [
      "<value>"
    ],
    "assay_methods": [
      "<value>"
    ],
    "assay_method_counts": [
      "<value>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.imm_port_search_studies-parameters)

## 16) 查询蛋白质复合物组成与功能 {#dnakit.databases.intact_get_complex_details}

- **作用：** 查询蛋白质复合物组成与功能。
- **处理方法：** 通过 `IntAct` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.intact_get_complex_details(*, complex_ac, format=None, verbose=True)`
- **输入：** 必须提供 `complex_ac`。
- **输出：** 返回字典，包含 `complexAc`, `name`, `systematicName`, `synonyms` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.intact_get_complex_details(
    complex_ac='CPX-915',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "complexAc": "<string>",
  "name": "<string>",
  "systematicName": "<string>",
  "synonyms": [
    "<string>"
  ],
  "functions": [
    "<string>"
  ],
  "properties": [
    "<string>"
  ],
  "ac": "<string>",
  "species": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_get_complex_details-parameters)

## 17) 查询以指定蛋白质为中心的互作网络 {#dnakit.databases.intact_get_interaction_network}

- **作用：** 查询以指定蛋白质为中心的互作网络。
- **处理方法：** 通过 `IntAct` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.intact_get_interaction_network(*, identifier=None, uniprot_id=None, protein_id=None, gene_symbol=None, gene_name=None, depth=None, limit=None, size=None, format=None, protein_name=None, protein=None, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回“查询以指定蛋白质为中心的互作网络”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.intact_get_interaction_network(
    identifier='P04637',
    depth=1,
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

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_get_interaction_network-parameters)

## 18) 查询指定蛋白质或基因的互作记录 {#dnakit.databases.intact_get_interactions}

- **作用：** 查询指定蛋白质或基因的互作记录。
- **处理方法：** 通过 `IntAct` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.intact_get_interactions(*, identifier, format=None, size=None, verbose=True)`
- **输入：** 必须提供 `identifier`。
- **输出：** 返回“查询指定蛋白质或基因的互作记录”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.intact_get_interactions(
    identifier='P04637',
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

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_get_interactions-parameters)

## 19) 按关键词检索 IntAct 分子互作 {#dnakit.databases.intact_search_interactions}

- **作用：** 按关键词检索 IntAct 分子互作。
- **处理方法：** 通过 `IntAct` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.intact_search_interactions(*, query, first=None, max=None, format=None, species=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回“按关键词检索 IntAct 分子互作”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.intact_search_interactions(
    query='P04637',
    max=10,
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

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_search_interactions-parameters)

## 20) 查询蛋白质关联的结构域与家族条目 {#dnakit.databases.inter_pro_get_entries_for_protein}

- **作用：** 查询蛋白质关联的结构域与家族条目。
- **处理方法：** 通过 `InterProEntry` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.inter_pro_get_entries_for_protein(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字典，包含 `protein_accession`, `total_entries`, `entries` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.inter_pro_get_entries_for_protein(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "protein_accession": "<string>",
  "total_entries": "<integer>",
  "entries": [
    {
      "accession": "<string>",
      "name": "<string>",
      "type": "<string>",
      "source_database": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.inter_pro_get_entries_for_protein-parameters)

## 21) 查询蛋白质结构域与功能注释 {#dnakit.databases.inter_pro_get_protein_domains}

- **作用：** 查询蛋白质结构域与功能注释。
- **处理方法：** 通过 `InterPro` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.inter_pro_get_protein_domains(*, protein_id=None, accession=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `metadata`, `proteins` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.inter_pro_get_protein_domains(
    protein_id='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "metadata": {
      "accession": "<string>",
      "name": "<string>",
      "source_database": "<string>",
      "type": "<string>",
      "integrated": "<null>",
      "member_databases": "<object>",
      "go_terms": "<array>"
    },
    "proteins": [
      "<object>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.inter_pro_get_protein_domains-parameters)

## 22) 按关键词与物种检索 KEGG 基因 {#dnakit.databases.kegg_find_genes}

- **作用：** 按关键词与物种检索 KEGG 基因。
- **处理方法：** 通过 `KEGGFindGenes` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.kegg_find_genes(*, keyword=None, organism=None, query=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `gene_id`, `description` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.kegg_find_genes(
    keyword='insulin',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "gene_id": "<string>",
    "description": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_find_genes-parameters)

## 23) 查询 KEGG 基因序列、功能与通路 {#dnakit.databases.kegg_get_gene_info}

- **作用：** 查询 KEGG 基因序列、功能与通路。
- **处理方法：** 通过 `KEGGGetGeneInfo` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.kegg_get_gene_info(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `status`, `data`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.kegg_get_gene_info(
    gene_id='hsa:348',
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
    "gene_id": "<string>",
    "raw_data": "<string>",
    "line_count": "<integer>"
  },
  "url": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_gene_info-parameters)

## 24) 查询基因参与的 KEGG 通路 {#dnakit.databases.kegg_get_gene_pathways}

- **作用：** 查询基因参与的 KEGG 通路。
- **处理方法：** 通过 `KEGGExt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.kegg_get_gene_pathways(*, gene_id=None, gene_symbol=None, gene=None, organism=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `gene_id`, `pathway_count`, `pathways` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.kegg_get_gene_pathways(
    gene_id='hsa:7157',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_id": "<string>",
  "pathway_count": "<integer>",
  "pathways": [
    {
      "pathway_id": "<string>",
      "pathway_name": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_gene_pathways-parameters)

## 25) 列出 KEGG 通路中的成员基因 {#dnakit.databases.kegg_get_pathway_genes}

- **作用：** 列出 KEGG 通路中的成员基因。
- **处理方法：** 通过 `KEGGExt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.kegg_get_pathway_genes(*, pathway_id, verbose=True)`
- **输入：** 必须提供 `pathway_id`。
- **输出：** 返回字典，包含 `pathway_id`, `pathway_name`, `gene_count`, `genes` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.kegg_get_pathway_genes(
    pathway_id='hsa04115',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "pathway_id": "<string>",
  "pathway_name": "<string>",
  "gene_count": "<integer>",
  "genes": [
    "<string>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_pathway_genes-parameters)

## 26) 查询 KEGG 通路的基因、化合物与反应 {#dnakit.databases.kegg_get_pathway_info}

- **作用：** 查询 KEGG 通路的基因、化合物与反应。
- **处理方法：** 通过 `KEGGGetPathwayInfo` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.kegg_get_pathway_info(*, pathway_id, verbose=True)`
- **输入：** 必须提供 `pathway_id`。
- **输出：** 返回字典，包含 `status`, `data`, `url` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.kegg_get_pathway_info(
    pathway_id='hsa00010',
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
    "pathway_id": "<string>",
    "raw_data": "<string>",
    "line_count": "<integer>"
  },
  "url": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_pathway_info-parameters)

## 27) 按关键词检索 KEGG 通路 {#dnakit.databases.kegg_search_pathway}

- **作用：** 按关键词检索 KEGG 通路。
- **处理方法：** 通过 `KEGGSearchPathway` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.kegg_search_pathway(*, keyword, verbose=True)`
- **输入：** 必须提供 `keyword`。
- **输出：** 返回记录列表；每条记录包含 `pathway_id`, `description` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.kegg_search_pathway(
    keyword='diabetes',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "pathway_id": "<string>",
    "description": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_search_pathway-parameters)

## 28) 按编号查询 NCBI 蛋白质记录摘要 {#dnakit.databases.ncbiprotein_get_summary}

- **作用：** 按编号查询 NCBI 蛋白质记录摘要。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbiprotein_get_summary(*, id, verbose=True)`
- **输入：** 必须提供 `id`。
- **输出：** 返回字典，包含 `result` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbiprotein_get_summary(
    id='119395750',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "result": {
    "uids": [
      "<string>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.ncbiprotein_get_summary-parameters)

## 29) 查询 PANTHER 跨物种直系同源基因 {#dnakit.databases.panther_ortholog}

- **作用：** 查询 PANTHER 跨物种直系同源基因。
- **处理方法：** 通过 `PANTHER` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.panther_ortholog(*, gene_id, organism=None, target_organism=None, ortholog_type=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `gene_id`, `source_organism`, `target_organism`, `ortholog_type` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.panther_ortholog(
    gene_id='P04637',
    organism=9606,
    target_organism=10090,
    ortholog_type='LDO',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_id": "<string>",
  "source_organism": "<integer>",
  "target_organism": "<integer>",
  "ortholog_type": "<string>",
  "mapping": {
    "source_gene": "<string>",
    "target_gene": "<string>",
    "target_gene_symbol": "<string>",
    "ortholog_type": "<string>",
    "persistent_id": "<string>",
    "target_persistent_id": "<string>"
  },
  "mappings": [
    {
      "source_gene": "<string>",
      "target_gene": "<string>",
      "target_gene_symbol": "<string>",
      "ortholog_type": "<string>",
      "persistent_id": "<string>",
      "target_persistent_id": "<string>"
    }
  ],
  "total_mappings": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.panther_ortholog-parameters)

## 30) 读取 PDB 结构坐标文件文本 {#dnakit.databases.pdb_get_structure}

- **作用：** 读取 PDB 结构坐标文件文本。
- **处理方法：** 通过 `API` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pdb_get_structure(*, pdb_id, format='cif', max_bytes=10000000, verbose=True)`
- **输入：** 必须提供 `pdb_id`。
- **输出：** 字典含 pdb_id、format、data（坐标文件文本）、url、sha256；不写入磁盘。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pdb_get_structure(
    pdb_id='1CRN',
    format='cif',
    verbose=False,
)
print(result)
```

- **示例结果：** 字典含 pdb_id、format、data（坐标文件文本）、url、sha256；不写入磁盘。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "pdb_id": "<str>",
  "format": "<str>",
  "data": "<coordinate-file text>",
  "url": "<str>",
  "sha256": "<str>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.pdb_get_structure-parameters)

## 31) 映射 PDB 链与 UniProt 残基区间 {#dnakit.databases.pdbe_get_uniprot_mappings}

- **作用：** 映射 PDB 链与 UniProt 残基区间。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pdbe_get_uniprot_mappings(*, pdb_id=None, uniprot_id=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回以记录键为索引的字典；每条记录包含 `UniProt` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pdbe_get_uniprot_mappings(
    pdb_id='4hhb',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "<record key>": {
    "UniProt": {
      "<record key>": "<object>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.pdbe_get_uniprot_mappings-parameters)

## 32) 按覆盖度与分辨率筛选蛋白质结构 {#dnakit.databases.pdbe_sifts_get_best_structures}

- **作用：** 按覆盖度与分辨率筛选蛋白质结构。
- **处理方法：** 通过 `PDBeSIFTS` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pdbe_sifts_get_best_structures(*, uniprot_accession, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `uniprot_accession`。
- **输出：** 返回字典，包含 `uniprot_accession`, `structures`, `total_structures`, `returned_structures` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pdbe_sifts_get_best_structures(
    uniprot_accession='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "uniprot_accession": "<string>",
  "structures": [
    {
      "pdb_id": "<string>",
      "chain_id": "<string>",
      "uniprot_start": "<integer>",
      "uniprot_end": "<integer>",
      "resolution": "<number>",
      "experimental_method": "<string>",
      "coverage": "<number>",
      "tax_id": "<integer>"
    }
  ],
  "total_structures": "<integer>",
  "returned_structures": "<integer>",
  "offset": "<integer>",
  "note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.pdbe_sifts_get_best_structures-parameters)

## 33) 查询蛋白质 Pfam 结构域与位置 {#dnakit.databases.pfam_get_protein_annotations}

- **作用：** 查询蛋白质 Pfam 结构域与位置。
- **处理方法：** 通过 `Pfam` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pfam_get_protein_annotations(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字典，包含 `accession`, `protein_length`, `domain_count`, `domains` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pfam_get_protein_annotations(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "accession": "<string>",
  "protein_length": "<integer>",
  "domain_count": "<integer>",
  "domains": [
    {
      "pfam_accession": "<string>",
      "name": "<string>",
      "short_name": "<string>",
      "type": "<string>",
      "integrated_interpro": "<string>",
      "start": "<integer>",
      "end": "<integer>",
      "score": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.pfam_get_protein_annotations-parameters)

## 34) 查询植物通路的反应与注释 {#dnakit.databases.plant_reactome_get_pathway}

- **作用：** 查询植物通路的反应与注释。
- **处理方法：** 通过 `PlantReactome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.plant_reactome_get_pathway(*, pathway_id, verbose=True)`
- **输入：** 必须提供 `pathway_id`。
- **输出：** 返回字典，包含 `stId`, `displayName`, `speciesName`, `isInDisease` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.plant_reactome_get_pathway(
    pathway_id='R-OSA-1119312',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "stId": "<string>",
  "displayName": "<string>",
  "speciesName": "<string>",
  "isInDisease": "<boolean>",
  "isInferred": "<boolean>",
  "hasEvent": [
    "<value>"
  ],
  "compartment": [
    "<value>"
  ],
  "goBiologicalProcess": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.plant_reactome_get_pathway-parameters)

## 35) 列出植物通路数据库覆盖的物种 {#dnakit.databases.plant_reactome_list_species}

- **作用：** 列出植物通路数据库覆盖的物种。
- **处理方法：** 通过 `PlantReactome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.plant_reactome_list_species(*, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `dbId`, `displayName`, `taxId`, `abbreviation` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.plant_reactome_list_species(
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "dbId": "<integer>",
    "displayName": "<string>",
    "taxId": "<string>",
    "abbreviation": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.plant_reactome_list_species-parameters)

## 36) 检索植物代谢与调控通路 {#dnakit.databases.plant_reactome_search_pathways}

- **作用：** 检索植物代谢与调控通路。
- **处理方法：** 通过 `PlantReactome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.plant_reactome_search_pathways(*, query, species=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `stId`, `name`, `species`, `exact_type` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.plant_reactome_search_pathways(
    query='photosynthesis',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "stId": "<string>",
    "name": "<string>",
    "species": "<string>",
    "exact_type": "<string>",
    "compartment_names": [
      "<value>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.plant_reactome_search_pathways-parameters)

## 37) 查询蛋白质序列特征与修饰位点 {#dnakit.databases.proteins_api_get_features}

- **作用：** 查询蛋白质序列特征与修饰位点。
- **处理方法：** 通过 `ProteinsAPI` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.proteins_api_get_features(*, accession, format=None, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回“查询蛋白质序列特征与修饰位点”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.proteins_api_get_features(
    accession='P05067',
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

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.proteins_api_get_features-parameters)

## 38) 按基因、蛋白名称或编号检索蛋白质 {#dnakit.databases.proteins_api_search}

- **作用：** 按基因、蛋白名称或编号检索蛋白质。
- **处理方法：** 通过 `ProteinsAPI` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.proteins_api_search(*, query=None, size=None, offset=None, format=None, gene_symbol=None, gene=None, reviewed=None, organism=None, taxid=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回“按基因、蛋白名称或编号检索蛋白质”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.proteins_api_search(
    query='BRCA1',
    size=5,
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

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.proteins_api_search-parameters)

## 39) 查询基因产物的 GO 注释与证据代码 {#dnakit.databases.quick_go_annotations_by_gene}

- **作用：** 查询基因产物的 GO 注释与证据代码。
- **处理方法：** 通过 `QuickGO` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.quick_go_annotations_by_gene(*, gene_product_id, aspect=None, taxon_id=None, limit=None, verbose=True)`
- **输入：** 必须提供 `gene_product_id`。
- **输出：** 返回记录列表；每条记录包含 `go_id`, `go_name`, `go_aspect`, `gene_product_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.quick_go_annotations_by_gene(
    gene_product_id='UniProtKB:P04637',
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "go_id": "<string>",
    "go_name": "<string>",
    "go_aspect": "<string>",
    "gene_product_id": "<string>",
    "symbol": "<string>",
    "qualifier": "<string>",
    "go_evidence": "<string>",
    "assigned_by": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.quick_go_annotations_by_gene-parameters)

## 40) 查询 PDB 结构实验方法与元数据 {#dnakit.databases.rcsbdata_get_entry}

- **作用：** 查询 PDB 结构实验方法与元数据。
- **处理方法：** 通过 `RCSBData` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rcsbdata_get_entry(*, pdb_id=None, entry_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组至少选一组：`pdb_id` / `entry_id`。
- **输出：** 返回字典，包含 `pdb_id`, `title`, `method`, `resolution` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rcsbdata_get_entry(
    pdb_id='4HHB',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "pdb_id": "<string>",
  "title": "<string>",
  "method": "<string>",
  "resolution": "<number>",
  "deposit_date": "<string>",
  "release_date": "<string>",
  "revision_date": "<string>",
  "polymer_entity_count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.rcsbdata_get_entry-parameters)

## 41) 检索 Reactome 通路、反应与生物实体 {#dnakit.databases.reactome_content_search}

- **作用：** 检索 Reactome 通路、反应与生物实体。
- **处理方法：** 通过 `ReactomeContent` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.reactome_content_search(*, query, species=None, types=None, rows=None, start=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `query`, `species`, `count`, `total_results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.reactome_content_search(
    query='apoptosis',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "query": "<string>",
  "species": "<string>",
  "count": "<integer>",
  "total_results": "<integer>",
  "start": "<integer>",
  "has_more": "<boolean>",
  "results": [
    {
      "type": "<string>",
      "stId": "<string>",
      "name": "<string>",
      "species": "<array>",
      "is_disease": "<boolean>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_content_search-parameters)

## 42) 查询生物反应的参与分子与复合物 {#dnakit.databases.reactome_get_participants}

- **作用：** 查询生物反应的参与分子与复合物。
- **处理方法：** 通过 `Reactome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.reactome_get_participants(*, stId, verbose=True)`
- **输入：** 必须提供 `stId`。
- **输出：** 返回记录列表；每条记录包含 `peDbId`, `displayName`, `schemaClass`, `refEntities` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.reactome_get_participants(
    stId='R-HSA-111289',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "peDbId": "<integer>",
    "displayName": "<string>",
    "schemaClass": "<string>",
    "refEntities": [
      "<object>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_get_participants-parameters)

## 43) 查询 Reactome 通路完整信息 {#dnakit.databases.reactome_get_pathway}

- **作用：** 查询 Reactome 通路完整信息。
- **处理方法：** 通过 `Reactome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.reactome_get_pathway(*, stId, verbose=True)`
- **输入：** 必须提供 `stId`。
- **输出：** 返回字典，包含 `dbId`, `displayName`, `stId`, `stIdVersion` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.reactome_get_pathway(
    stId='R-HSA-73817',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "dbId": "<integer>",
  "displayName": "<string>",
  "stId": "<string>",
  "stIdVersion": "<string>",
  "isInDisease": "<boolean>",
  "isInferred": "<boolean>",
  "maxDepth": "<integer>",
  "name": [
    "<string>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_get_pathway-parameters)

## 44) 列出通路包含的反应与子通路 {#dnakit.databases.reactome_get_pathway_reactions}

- **作用：** 列出通路包含的反应与子通路。
- **处理方法：** 通过 `Reactome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.reactome_get_pathway_reactions(*, stId, verbose=True)`
- **输入：** 必须提供 `stId`。
- **输出：** 返回“列出通路包含的反应与子通路”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.reactome_get_pathway_reactions(
    stId='R-HSA-73817',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  "<integer>"
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_get_pathway_reactions-parameters)

## 45) 将 UniProt 蛋白映射到 Reactome 通路 {#dnakit.databases.reactome_map_uniprot_to_pathways}

- **作用：** 将 UniProt 蛋白映射到 Reactome 通路。
- **处理方法：** 通过 `Reactome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.reactome_map_uniprot_to_pathways(*, uniprot_id=None, id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `dbId`, `displayName`, `stId`, `stIdVersion` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.reactome_map_uniprot_to_pathways(
    uniprot_id='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "dbId": "<integer>",
    "displayName": "<string>",
    "stId": "<string>",
    "stIdVersion": "<string>",
    "isInDisease": "<boolean>",
    "isInferred": "<boolean>",
    "maxDepth": "<integer>",
    "name": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_map_uniprot_to_pathways-parameters)

## 46) 查询 STRING 蛋白质综合功能注释 {#dnakit.databases.string_get_functional_annotations}

- **作用：** 查询 STRING 蛋白质综合功能注释。
- **处理方法：** 通过 `STRINGExt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.string_get_functional_annotations(*, identifiers, species=None, category=None, verbose=True)`
- **输入：** 必须提供 `identifiers`。
- **输出：** 返回字典，包含 `identifiers`, `species`, `total_annotations`, `category_summary` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.string_get_functional_annotations(
    identifiers='TP53',
    species=9606,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "identifiers": "<string>",
  "species": "<integer>",
  "total_annotations": "<integer>",
  "category_summary": "<object>",
  "annotations": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_get_functional_annotations-parameters)

## 47) 查询蛋白质互作伙伴与置信度 {#dnakit.databases.string_get_interaction_partners}

- **作用：** 查询蛋白质互作伙伴与置信度。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.string_get_interaction_partners(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `stringId_A`, `stringId_B`, `preferredName_A`, `preferredName_B` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.string_get_interaction_partners(
    identifiers='TP53',
    species=9606,
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "stringId_A": "<string>",
    "stringId_B": "<string>",
    "preferredName_A": "<string>",
    "preferredName_B": "<string>",
    "ncbiTaxonId": "<string>",
    "score": "<number>",
    "nscore": "<number>",
    "fscore": "<number>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_get_interaction_partners-parameters)

## 48) 构建查询蛋白及其邻居的互作网络 {#dnakit.databases.string_get_network}

- **作用：** 构建查询蛋白及其邻居的互作网络。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.string_get_network(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `stringId_A`, `stringId_B`, `preferredName_A`, `preferredName_B` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.string_get_network(
    identifiers='TP53',
    species=9606,
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
    "stringId_A": "<string>",
    "stringId_B": "<string>",
    "preferredName_A": "<string>",
    "preferredName_B": "<string>",
    "ncbiTaxonId": "<string>",
    "score": "<number>",
    "nscore": "<number>",
    "fscore": "<number>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_get_network-parameters)

## 49) 检验蛋白质集合的互作富集程度 {#dnakit.databases.string_ppi_enrichment}

- **作用：** 检验蛋白质集合的互作富集程度。
- **处理方法：** 通过 `STRING` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.string_ppi_enrichment(*, protein_ids, species=None, confidence_score=None, verbose=True)`
- **输入：** 必须提供 `protein_ids`。
- **输出：** 返回字典，包含 `number_of_nodes`, `number_of_edges`, `expected_number_of_edges`, `p_value` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.string_ppi_enrichment(
    protein_ids=['TP53', 'MDM2', 'ATM', 'CHEK2', 'BRCA1', 'BRCA2'],
    species=9606,
    confidence_score=0.4,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "number_of_nodes": "<integer>",
  "number_of_edges": "<integer>",
  "expected_number_of_edges": "<number>",
  "p_value": "<number>",
  "avg_node_degree": "<number>",
  "local_clustering_coefficient": "<number>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_ppi_enrichment-parameters)

## 50) 逐残基标注结合界面、口袋与溶剂暴露 {#dnakit.databases.structure_annotate_per_residue}

- **作用：** 逐残基标注结合界面、口袋与溶剂暴露。
- **处理方法：** 通过 `StructureAnnotation` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.structure_annotate_per_residue(*, operation=None, pdb_id=None, pdb_content=None, target_chain=None, partner_chains=None, ligand_resnames=None, distance_cutoff=None, core_rsa_cutoff=None, include_secondary_structure=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `pdb_id`, `target_chain`, `partner_chains` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.structure_annotate_per_residue(
    operation='annotate_per_residue',
    pdb_id='1CRN',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "status": "<string>",
  "pdb_id": "<string>",
  "target_chain": "<string>",
  "partner_chains": [
    "<string>"
  ],
  "ligand_resnames": [
    "<string>"
  ],
  "distance_cutoff": "<number>",
  "core_rsa_cutoff": "<number>",
  "n_residues": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.structure_annotate_per_residue-parameters)

## 51) 提取 UniProt 蛋白质变异与相关注释 {#dnakit.databases.uni_prot_get_disease_variants_by_accession}

- **作用：** 提取 UniProt 蛋白质变异与相关注释。
- **处理方法：** 通过 `UniProt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_prot_get_disease_variants_by_accession(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回记录列表；每条记录包含 `type`, `location`, `description`, `featureId` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_prot_get_disease_variants_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "type": "<string>",
    "location": {
      "start": "<object>",
      "end": "<object>"
    },
    "description": "<string>",
    "featureId": "<string>",
    "alternativeSequence": {
      "originalSequence": "<string>",
      "alternativeSequences": "<array>"
    },
    "featureCrossReferences": [
      "<object>"
    ],
    "evidences": [
      "<object>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_disease_variants_by_accession-parameters)

## 52) 按编号读取完整 UniProt 条目 {#dnakit.databases.uni_prot_get_entry_by_accession}

- **作用：** 按编号读取完整 UniProt 条目。
- **处理方法：** 通过 `UniProt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_prot_get_entry_by_accession(*, accession, compact=None, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字典，包含 `entryType`, `primaryAccession`, `uniProtkbId`, `protein_name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_prot_get_entry_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "entryType": "<string>",
  "primaryAccession": "<string>",
  "uniProtkbId": "<string>",
  "protein_name": "<string>",
  "gene_names": [
    "<string>"
  ],
  "organism": {
    "scientificName": "<string>",
    "commonName": "<string>",
    "taxonId": "<integer>",
    "lineage": [
      "<string>"
    ]
  },
  "sequence": {
    "value": "<string>",
    "length": "<integer>",
    "molWeight": "<integer>",
    "crc64": "<string>",
    "md5": "<string>"
  },
  "annotationScore": "<number>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_entry_by_accession-parameters)

## 53) 提取蛋白质结构域、活性位点与序列特征 {#dnakit.databases.uni_prot_get_features_by_accession}

- **作用：** 提取蛋白质结构域、活性位点与序列特征。
- **处理方法：** 通过 `UniProt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_prot_get_features_by_accession(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回记录列表；每条记录包含 `type`, `location`, `description` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_prot_get_features_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "type": "<string>",
    "location": "<object>",
    "description": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_features_by_accession-parameters)

## 54) 提取 UniProt 蛋白质功能注释 {#dnakit.databases.uni_prot_get_function_by_accession}

- **作用：** 提取 UniProt 蛋白质功能注释。
- **处理方法：** 通过 `UniProt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_prot_get_function_by_accession(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回“提取 UniProt 蛋白质功能注释”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_prot_get_function_by_accession(
    accession='P04637',
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

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_function_by_accession-parameters)

## 55) 提取 UniProt 标准蛋白质序列 {#dnakit.databases.uni_prot_get_sequence_by_accession}

- **作用：** 提取 UniProt 标准蛋白质序列。
- **处理方法：** 通过 `UniProt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_prot_get_sequence_by_accession(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字符串。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_prot_get_sequence_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
"<string>"
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_sequence_by_accession-parameters)

## 56) 按组合检索式查询 UniProt 蛋白质 {#dnakit.databases.uni_prot_search}

- **作用：** 按组合检索式查询 UniProt 蛋白质。
- **处理方法：** 通过 `UniProt` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.uni_prot_search(*, query, organism=None, limit=None, min_length=None, max_length=None, fields=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `total_results`, `returned`, `results`, `normalization_note` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.uni_prot_search(
    query='gene:TP53',
    organism='homo sapiens',
    limit=5,
    min_length=0,
    max_length=10000,
    fields=[],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total_results": "<integer>",
  "returned": "<integer>",
  "results": [
    {
      "accession": "<string>",
      "id": "<string>",
      "protein_name": "<string>",
      "gene_names": "<array>",
      "organism": "<string>",
      "length": "<integer>"
    }
  ],
  "normalization_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_search-parameters)

## 57) 按通路标题检索 WikiPathways {#dnakit.databases.wiki_pathways_search}

- **作用：** 按通路标题检索 WikiPathways。
- **处理方法：** 通过 `WikiPathwaysSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.wiki_pathways_search(*, query, organism=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `status`, `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.wiki_pathways_search(
    query='apoptosis',
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
  "status": "success",
  "data": [
    {
      "wpid": "<str>",
      "uri": "<str>",
      "title": "<str>",
      "organism": "<str>"
    }
  ],
  "metadata": {
    "query": "<str>",
    "count": "<int>",
    "backend": "WikiPathways SPARQL"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_pathways.md#dnakit.databases.wiki_pathways_search-parameters)

<!-- END GENERATED API: pathways -->
