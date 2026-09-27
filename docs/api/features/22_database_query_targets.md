# 靶点、药物与关联证据查询

查询靶点、药物及关联证据。

<!-- BEGIN GENERATED API: targets -->

<span id="targets"></span>

## 1) 按药物名称或编号查询作用机制 {#dnakit.databases.ch_embl_get_drug_mechanisms}

- **作用：** 按药物名称或编号查询作用机制。
- **处理方法：** 通过 `ChEMBL` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_embl_get_drug_mechanisms(*, limit=None, offset=None, drug_chembl_id=None, drug_name=None, molecule_chembl_id=None, chembl_id=None, drug_chembl_id__exact=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_embl_get_drug_mechanisms(
    drug_chembl_id='CHEMBL1201581',
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
    "mechanisms": [
      "<object>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_get_drug_mechanisms-parameters)

## 2) 查询 ChEMBL 分子结构与性质 {#dnakit.databases.ch_embl_get_molecule}

- **作用：** 查询 ChEMBL 分子结构与性质。
- **处理方法：** 通过 `ChEMBL` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_embl_get_molecule(*, chembl_id=None, format=None, molecule_chembl_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_embl_get_molecule(
    chembl_id='CHEMBL25',
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

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_get_molecule-parameters)

## 3) 查询靶点相关的生物活性测量值 {#dnakit.databases.ch_embl_get_target_activities}

- **作用：** 查询靶点相关的生物活性测量值。
- **处理方法：** 通过 `ChEMBL` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_embl_get_target_activities(*, target_chembl_id__exact=None, limit=None, offset=None, target_chembl_id=None, chembl_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。 参数组恰好选一组：`target_chembl_id__exact` / `target_chembl_id`。
- **输出：** 返回字典，包含 `status`, `data`, `count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_embl_get_target_activities(
    target_chembl_id__exact='CHEMBL2074',
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
    "activities": [
      "<object>"
    ]
  },
  "count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_get_target_activities-parameters)

## 4) 按名称与研发阶段检索药物 {#dnakit.databases.ch_embl_search_drugs}

- **作用：** 按名称与研发阶段检索药物。
- **处理方法：** 通过 `ChEMBL` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_embl_search_drugs(*, query=None, max_phase=None, limit=None, offset=None, molecule_chembl_id=None, q=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data`, `count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_embl_search_drugs(
    max_phase=4,
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
    "drugs": [
      "<object>"
    ]
  },
  "count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_search_drugs-parameters)

## 5) 按名称、物种或类型检索药物靶点 {#dnakit.databases.ch_embl_search_targets}

- **作用：** 按名称、物种或类型检索药物靶点。
- **处理方法：** 通过 `ChEMBL` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ch_embl_search_targets(*, target_chembl_id=None, pref_name__contains=None, organism=None, target_type=None, fields=None, limit=None, offset=None, q=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `status`, `data`, `count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ch_embl_search_targets(
    pref_name__contains='kinase',
    organism='Homo sapiens',
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
    "targets": [
      "<object>"
    ]
  },
  "count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_search_targets-parameters)

## 6) 查看已不可用的 CTD 化学物基因互作查询 {#dnakit.databases.ctd_get_chemical_gene_interactions}

- **作用：** 查看已不可用的 CTD 化学物基因互作查询。
- **处理方法：** 通过 `CTD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ctd_get_chemical_gene_interactions(*, input_terms, verbose=True)`
- **输入：** 必须提供 `input_terms`。 该数据源入口当前不可用，调用会返回错误；以下代码仅说明调用形式。
- **输出：** 此入口当前报告服务不可用，不返回有效互作结果。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ctd_get_chemical_gene_interactions(
    input_terms='bisphenol A',
    verbose=False,
)
print(result)
```

- **示例结果：** 此入口当前报告服务不可用，不返回有效互作结果。

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.ctd_get_chemical_gene_interactions-parameters)

## 7) 批量查询药物与基因互作证据 {#dnakit.databases.dgidb_get_drug_gene_interactions}

- **作用：** 批量查询药物与基因互作证据。
- **处理方法：** 通过 `DGIdb` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dgidb_get_drug_gene_interactions(*, genes=None, gene_name=None, gene=None, interaction_sources=None, interaction_types=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dgidb_get_drug_gene_interactions(
    genes=['EGFR'],
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
      "nodes": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.dgidb_get_drug_gene_interactions-parameters)

## 8) 查询基因的可成药类别 {#dnakit.databases.dgidb_get_gene_druggability}

- **作用：** 查询基因的可成药类别。
- **处理方法：** 通过 `DGIdb` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.dgidb_get_gene_druggability(*, genes=None, gene_name=None, gene=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.dgidb_get_gene_druggability(
    genes=['EGFR', 'TP53', 'BRCA1'],
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
      "nodes": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.dgidb_get_gene_druggability-parameters)

## 9) 查询药物基本信息与批准状态 {#dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id}

- **作用：** 查询药物基本信息与批准状态。
- **处理方法：** 通过 `XML` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `error`, `query`, `total_matches`, `total_returned_results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id(
    query='Aspirin',
    case_sensitive=False,
    exact_match=False,
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "error": "<string>",
  "query": "<string>",
  "total_matches": "<integer>",
  "total_returned_results": "<integer>",
  "results": [
    {
      "_tag": "<string>",
      "_text": "<string>",
      "_attributes": "<object>",
      "drug_name": "<string>",
      "drugbank_id": "<string>",
      "secondary_ids": "<string>",
      "description": "<string>",
      "cas_number": "<string>"
    }
  ],
  "search_parameters": {
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id-parameters)

## 10) 查询药效学、作用机制与药代动力学 {#dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id}

- **作用：** 查询药效学、作用机制与药代动力学。
- **处理方法：** 通过 `XML` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `error`, `query`, `total_matches`, `total_returned_results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id(
    query='Aspirin',
    case_sensitive=False,
    exact_match=False,
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "error": "<string>",
  "query": "<string>",
  "total_matches": "<integer>",
  "total_returned_results": "<integer>",
  "results": [
    {
      "_tag": "<string>",
      "_text": "<string>",
      "_attributes": "<object>",
      "drug_name": "<string>",
      "drugbank_id": "<string>",
      "synonyms": "<string>",
      "pharmacodynamics": "<string>",
      "mechanism_of_action": "<string>"
    }
  ],
  "search_parameters": {
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id-parameters)

## 11) 查询药物靶点、酶、载体与转运体 {#dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id}

- **作用：** 查询药物靶点、酶、载体与转运体。
- **处理方法：** 通过 `XML` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `error`, `query`, `total_matches`, `total_returned_results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id(
    query='Aspirin',
    case_sensitive=False,
    exact_match=False,
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "error": "<string>",
  "query": "<string>",
  "total_matches": "<integer>",
  "total_returned_results": "<integer>",
  "results": [
    {
      "_tag": "<string>",
      "_text": "<string>",
      "_attributes": "<object>",
      "drug_name": "<string>",
      "drugbank_id": "<string>",
      "synonyms": "<string>",
      "targets": "<array>",
      "enzymes": "<array>"
    }
  ],
  "search_parameters": {
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id-parameters)

## 12) 按名称与同义词检索 DrugBank 词表 {#dnakit.databases.drugbank_vocab_search}

- **作用：** 按名称与同义词检索 DrugBank 词表。
- **处理方法：** 通过 `Dataset` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.drugbank_vocab_search(*, query, search_fields=None, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `query`, `total_results`, `results`, `search_parameters` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.drugbank_vocab_search(
    query='aspirin',
    search_fields=['Common name', 'Synonyms'],
    case_sensitive=False,
    exact_match=False,
    limit=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "query": "<string>",
  "total_results": "<integer>",
  "results": [
    {
      "DrugBank ID": "<string>",
      "Accession Numbers": "<string>",
      "Common name": "<string>",
      "CAS": "<string>",
      "UNII": "<string>",
      "Synonyms": "<string>",
      "Standard InChI Key": "<string>",
      "matched_fields": "<array>"
    }
  ],
  "search_parameters": {
    "search_fields": [
      "<string>"
    ],
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_vocab_search-parameters)

## 13) 查询文献支持的疾病靶点关联评分 {#dnakit.databases.europepmc_disease_target_score}

- **作用：** 查询文献支持的疾病靶点关联评分。
- **处理方法：** 通过 `DiseaseTargetScore` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.europepmc_disease_target_score(*, efoId, pageSize=None, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `disease_info`, `datasource`, `total_targets_with_scores`, `target_scores` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.europepmc_disease_target_score(
    efoId='MONDO_0011996',
    pageSize=100,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "disease_info": {
    "disease_id": "<string>",
    "disease_name": "<string>"
  },
  "datasource": "<string>",
  "total_targets_with_scores": "<integer>",
  "target_scores": [
    {
      "target_symbol": "<string>",
      "target_id": "<string>",
      "datasource": "<string>",
      "score": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.europepmc_disease_target_score-parameters)

## 14) 查询表达数据支持的疾病靶点关联评分 {#dnakit.databases.expression_atlas_disease_target_score}

- **作用：** 查询表达数据支持的疾病靶点关联评分。
- **处理方法：** 通过 `DiseaseTargetScore` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.expression_atlas_disease_target_score(*, efoId, pageSize=None, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `disease_info`, `datasource`, `total_targets_with_scores`, `target_scores` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.expression_atlas_disease_target_score(
    efoId='MONDO_0011996',
    pageSize=100,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "disease_info": {
    "disease_id": "<string>",
    "disease_name": "<string>"
  },
  "datasource": "<string>",
  "total_targets_with_scores": "<integer>",
  "target_scores": [
    {
      "target_symbol": "<string>",
      "target_id": "<string>",
      "datasource": "<string>",
      "score": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.expression_atlas_disease_target_score-parameters)

## 15) 查询药品有效成分 {#dnakit.databases.fda_get_active_ingredient_info_by_drug_name}

- **作用：** 查询药品有效成分。
- **处理方法：** 通过 `FDADrugLabel` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fda_get_active_ingredient_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **输入：** 必须提供 `drug_name`。
- **输出：** 返回字典，包含 `meta`, `results`, `result_count`, `duplicates_removed` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fda_get_active_ingredient_info_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "active_ingredient": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_active_ingredient_info_by_drug_name-parameters)

## 16) 查询药品标签中的不良反应信息 {#dnakit.databases.fda_get_adverse_reactions_by_drug_name}

- **作用：** 查询药品标签中的不良反应信息。
- **处理方法：** 通过 `FDADrugLabel` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fda_get_adverse_reactions_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **输入：** 必须提供 `drug_name`。
- **输出：** 返回字典，包含 `meta`, `results`, `result_count`, `duplicates_removed` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fda_get_adverse_reactions_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "adverse_reactions": "<array>",
      "warnings_and_cautions": "<array>",
      "spl_product_data_elements": "<array>",
      "warnings": "<array>",
      "boxed_warning": "<array>",
      "precautions": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>",
  "section_note": "<string>",
  "note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_adverse_reactions_by_drug_name-parameters)

## 17) 查询药品黑框警告与相关安全信息 {#dnakit.databases.fda_get_boxed_warning_info_by_drug_name}

- **作用：** 查询药品黑框警告与相关安全信息。
- **处理方法：** 通过 `FDADrugLabel` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fda_get_boxed_warning_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **输入：** 必须提供 `drug_name`。
- **输出：** 返回字典，包含 `meta`, `results`, `result_count`, `duplicates_removed` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fda_get_boxed_warning_info_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "boxed_warning": "<array>",
      "spl_product_data_elements": "<array>",
      "warnings_and_cautions": "<array>",
      "warnings": "<array>",
      "related_sections_present": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>",
  "section_note": "<string>",
  "note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_boxed_warning_info_by_drug_name-parameters)

## 18) 查询药品标签中的临床研究信息 {#dnakit.databases.fda_get_clinical_studies_info_by_drug_name}

- **作用：** 查询药品标签中的临床研究信息。
- **处理方法：** 通过 `FDADrugLabel` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fda_get_clinical_studies_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **输入：** 必须提供 `drug_name`。
- **输出：** 返回字典，包含 `meta`, `results`, `result_count`, `duplicates_removed` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fda_get_clinical_studies_info_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "clinical_studies": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_clinical_studies_info_by_drug_name-parameters)

## 19) 查询药品适应证与使用说明 {#dnakit.databases.fda_get_indications_by_drug_name}

- **作用：** 查询药品适应证与使用说明。
- **处理方法：** 通过 `FDADrugLabel` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fda_get_indications_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **输入：** 必须提供 `drug_name`。
- **输出：** 返回字典，包含 `meta`, `results`, `result_count`, `duplicates_removed` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fda_get_indications_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "indications_and_usage": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_indications_by_drug_name-parameters)

## 20) 查询药品标签中的作用机制 {#dnakit.databases.fda_get_mechanism_of_action_by_drug_name}

- **作用：** 查询药品标签中的作用机制。
- **处理方法：** 通过 `FDADrugLabel` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fda_get_mechanism_of_action_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **输入：** 必须提供 `drug_name`。
- **输出：** 返回字典，包含 `meta`, `results`, `result_count`, `duplicates_removed` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fda_get_mechanism_of_action_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "mechanism_of_action": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_mechanism_of_action_by_drug_name-parameters)

## 21) 查询代谢物的疾病关联 {#dnakit.databases.hmdb_get_diseases}

- **作用：** 查询代谢物的疾病关联。
- **处理方法：** 通过 `Metabolite` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hmdb_get_diseases(*, operation=None, hmdb_id=None, compound_name=None, pubchem_cid=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `identifier`, `compound_name`, `pubchem_cid`, `ctd_query_term` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hmdb_get_diseases(
    hmdb_id='HMDB0000122',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "identifier": "<string>",
  "compound_name": "<string>",
  "pubchem_cid": "<integer>",
  "ctd_query_term": "<string>",
  "disease_count": "<integer>",
  "diseases": [
    {
      "disease_name": "<string>",
      "disease_id": "<string>",
      "disease_categories": "<string>",
      "direct_evidence": "<string>",
      "pubmed_ids": "<array>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.hmdb_get_diseases-parameters)

## 22) 按编号或名称查询代谢物信息 {#dnakit.databases.hmdb_get_metabolite}

- **作用：** 按编号或名称查询代谢物信息。
- **处理方法：** 通过 `Metabolite` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hmdb_get_metabolite(*, operation=None, hmdb_id=None, compound_name=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `pubchem_cid`, `name`, `iupac_name`, `formula` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hmdb_get_metabolite(
    hmdb_id='HMDB0000122',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "pubchem_cid": "<integer>",
  "name": "<string>",
  "iupac_name": "<string>",
  "formula": "<string>",
  "molecular_weight": "<string>",
  "smiles": "<string>",
  "inchikey": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.hmdb_get_metabolite-parameters)

## 23) 按名称或分子式检索代谢物 {#dnakit.databases.hmdb_search}

- **作用：** 按名称或分子式检索代谢物。
- **处理方法：** 通过 `Metabolite` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hmdb_search(*, operation=None, query, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `query`, `results`, `count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hmdb_search(
    query='glucose',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "query": "<string>",
  "results": [
    {
      "pubchem_cid": "<integer>",
      "name": "<string>",
      "formula": "<string>",
      "molecular_weight": "<string>",
      "smiles": "<string>"
    }
  ],
  "count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.hmdb_search-parameters)

## 24) 按疾病查询已知相关药物 {#dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id}

- **作用：** 按疾病查询已知相关药物。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id(*, efoId, size=None, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id(
    efoId='MONDO_0004975',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "drugAndClinicalCandidates": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id-parameters)

## 25) 按靶点查询相关药物与研发阶段 {#dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id}

- **作用：** 按靶点查询相关药物与研发阶段。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id(*, ensemblId, size=None, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id(
    ensemblId='ENSG00000121879',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "drugAndClinicalCandidates": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id-parameters)

## 26) 查询疾病相关的 HPO 表型 {#dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id}

- **作用：** 查询疾病相关的 HPO 表型。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id(*, efoId, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "phenotypes": "<object>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id-parameters)

## 27) 分页查询疾病关联靶点 {#dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id}

- **作用：** 分页查询疾病关联靶点。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id(*, efoId, size=None, index=None, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "associatedTargets": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id-parameters)

## 28) 查询 GWAS 可信集合与位点基因映射 {#dnakit.databases.open_targets_get_credible_set_detail}

- **作用：** 查询 GWAS 可信集合与位点基因映射。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_credible_set_detail(*, studyLocusId, verbose=True)`
- **输入：** 必须提供 `studyLocusId`。
- **输出：** 返回字典，包含 `credibleSet` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_credible_set_detail(
    studyLocusId='b758d8fb10924f5338cbad8d27c7dee8',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "credibleSet": {
    "studyLocusId": "<string>",
    "studyId": "<string>",
    "region": "<string>",
    "chromosome": "<string>",
    "position": "<integer>",
    "finemappingMethod": "<string>",
    "pValueMantissa": "<number>",
    "pValueExponent": "<integer>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_credible_set_detail-parameters)

## 29) 查询疾病本体的父级与祖先概念 {#dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id}

- **作用：** 查询疾病本体的父级与祖先概念。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id(*, efoId, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "ancestors": "<array>",
      "parents": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id-parameters)

## 30) 查询疾病本体的子级与后代概念 {#dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id}

- **作用：** 查询疾病本体的子级与后代概念。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id(*, efoId, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "descendants": "<array>",
      "children": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id-parameters)

## 31) 查询疾病名称、描述与跨库引用 {#dnakit.databases.open_targets_get_disease_description_by_efo_id}

- **作用：** 查询疾病名称、描述与跨库引用。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_disease_description_by_efo_id(*, efoId, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_description_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "description": "<string>",
      "dbXRefs": "<array>",
      "isTherapeuticArea": "<boolean>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_description_by_efo_id-parameters)

## 32) 按疾病名称解析本体编号与描述 {#dnakit.databases.open_targets_get_disease_id_description_by_name}

- **作用：** 按疾病名称解析本体编号与描述。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_disease_id_description_by_name(*, diseaseName=None, query=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_id_description_by_name(
    diseaseName='Diabetes Mellitus',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "search": {
      "hits": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_id_description_by_name-parameters)

## 33) 查询疾病同义名称 {#dnakit.databases.open_targets_get_disease_synonyms_by_efo_id}

- **作用：** 查询疾病同义名称。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_disease_synonyms_by_efo_id(*, efoId, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_synonyms_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "synonyms": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_synonyms_by_efo_id-parameters)

## 34) 查询疾病所属的治疗领域 {#dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id}

- **作用：** 查询疾病所属的治疗领域。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id(*, efoId, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "therapeuticAreas": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id-parameters)

## 35) 分页查询靶点关联疾病与表型 {#dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl}

- **作用：** 分页查询靶点关联疾病与表型。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl(*, ensemblId, page=None, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "associatedDiseases": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl-parameters)

## 36) 查询药物显著不良事件记录 {#dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id}

- **作用：** 查询药物显著不良事件记录。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id(*, chemblId, page=None, verbose=True)`
- **输入：** 必须提供 `chemblId`。
- **输出：** 返回字典，包含 `drug`, `metadata`, `error`, `error_details` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id(
    chemblId='CHEMBL25',
    page={'index': 0, 'size': 10},
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "drug": {
    "id": "<string>",
    "name": "<string>",
    "adverseEvents": {
      "count": "<integer>",
      "returned": "<integer>",
      "truncated": "<boolean>",
      "criticalValue": "<number>",
      "rows": "<array>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  },
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id-parameters)

## 37) 按通用名称解析药物 ChEMBL 编号 {#dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name}

- **作用：** 按通用名称解析药物 ChEMBL 编号。
- **处理方法：** 通过 `OpentargetToolDrugNameMatch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name(*, drugName=None, query=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name(
    drugName='Aspirin',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "search": {
      "hits": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name-parameters)

## 38) 按 ChEMBL 编号查询药物适应证 {#dnakit.databases.open_targets_get_drug_indications_by_chembl_id}

- **作用：** 按 ChEMBL 编号查询药物适应证。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_drug_indications_by_chembl_id(*, chemblId, verbose=True)`
- **输入：** 必须提供 `chemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_indications_by_chembl_id(
    chemblId='CHEMBL25',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "drug": {
      "id": "<string>",
      "name": "<string>",
      "indications": "<object>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_indications_by_chembl_id-parameters)

## 39) 按 ChEMBL 编号查询药物作用机制 {#dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id}

- **作用：** 按 ChEMBL 编号查询药物作用机制。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id(*, chemblId, verbose=True)`
- **输入：** 必须提供 `chemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id(
    chemblId='CHEMBL25',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "drug": {
      "id": "<string>",
      "name": "<string>",
      "mechanismsOfAction": "<object>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id-parameters)

## 40) 按 ChEMBL 编号查询药物警告 {#dnakit.databases.open_targets_get_drug_warnings_by_chembl_id}

- **作用：** 按 ChEMBL 编号查询药物警告。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_drug_warnings_by_chembl_id(*, chemblId, verbose=True)`
- **输入：** 必须提供 `chemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_warnings_by_chembl_id(
    chemblId='CHEMBL25',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "drug": {
      "id": "<string>",
      "name": "<string>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_warnings_by_chembl_id-parameters)

## 41) 按证据来源筛选靶点疾病关联 {#dnakit.databases.open_targets_get_evidence_by_datasource}

- **作用：** 按证据来源筛选靶点疾病关联。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_evidence_by_datasource(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, datasourceIds=None, size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_evidence_by_datasource(
    efoId='MONDO_0005011',
    ensemblId='ENSG00000141510',
    datasourceIds=['europepmc'],
    size=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "evidences": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>",
    "datasource_rename_note": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_evidence_by_datasource-parameters)

## 42) 查询 GWAS 研究队列与统计元数据 {#dnakit.databases.open_targets_get_gwas_study}

- **作用：** 查询 GWAS 研究队列与统计元数据。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_gwas_study(*, studyId, verbose=True)`
- **输入：** 必须提供 `studyId`。
- **输出：** 返回字典，包含 `study` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_gwas_study(
    studyId='GCST000392',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "study": {
    "id": "<string>",
    "studyType": "<string>",
    "traitFromSource": "<string>",
    "publicationFirstAuthor": "<string>",
    "publicationTitle": "<string>",
    "publicationDate": "<string>",
    "publicationJournal": "<string>",
    "pubmedId": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_gwas_study-parameters)

## 43) 查询靶点相关文献记录 {#dnakit.databases.open_targets_get_publications_by_target_ensembl_id}

- **作用：** 查询靶点相关文献记录。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_publications_by_target_ensembl_id(*, entityId, additionalIds=None, startYear=None, startMonth=None, endYear=None, endMonth=None, verbose=True)`
- **输入：** 必须提供 `entityId`。
- **输出：** 返回字典，包含 `data`, `metadata`, `error`, `error_details` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_publications_by_target_ensembl_id(
    entityId='ENSG00000141510',
    startYear=2023,
    startMonth=0,
    endYear=2023,
    endMonth=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "literatureOcurrences": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  },
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_publications_by_target_ensembl_id-parameters)

## 44) 查询与疾病相似的生物医学实体 {#dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id}

- **作用：** 查询与疾病相似的生物医学实体。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id(*, efoId, threshold, size, verbose=True)`
- **输入：** 必须提供 `efoId`, `threshold`, `size`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id(
    efoId='MONDO_0005011',
    threshold=1.0,
    size=1,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "similarEntities": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id-parameters)

## 45) 查询指定 GWAS 研究的可信集合 {#dnakit.databases.open_targets_get_study_credible_sets}

- **作用：** 查询指定 GWAS 研究的可信集合。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_study_credible_sets(*, studyIds, size=None, index=None, verbose=True)`
- **输入：** 必须提供 `studyIds`。
- **输出：** 返回字典，包含 `credibleSets` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_study_credible_sets(
    studyIds=['GCST000392'],
    size=10,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "credibleSets": {
    "count": "<integer>",
    "rows": [
      "<object>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_study_credible_sets-parameters)

## 46) 查询靶点的功能类别 {#dnakit.databases.open_targets_get_target_classes_by_ensembl_id}

- **作用：** 查询靶点的功能类别。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_classes_by_ensembl_id(*, ensemblId, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_classes_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "targetClass": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_classes_by_ensembl_id-parameters)

## 47) 查询靶点的 GO 功能注释 {#dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id}

- **作用：** 查询靶点的 GO 功能注释。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id(*, ensemblId, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "geneOntology": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id-parameters)

## 48) 查询靶点的基因组位置与方向 {#dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id}

- **作用：** 查询靶点的基因组位置与方向。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id(*, ensemblId, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "genomicLocation": "<object>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id-parameters)

## 49) 查询靶点的同源基因 {#dnakit.databases.open_targets_get_target_homologues_by_ensembl_id}

- **作用：** 查询靶点的同源基因。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_homologues_by_ensembl_id(*, ensemblId, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_homologues_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "homologues": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_homologues_by_ensembl_id-parameters)

## 50) 按名称解析靶点标识与描述 {#dnakit.databases.open_targets_get_target_id_description_by_name}

- **作用：** 按名称解析靶点标识与描述。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_id_description_by_name(*, targetName, verbose=True)`
- **输入：** 必须提供 `targetName`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_id_description_by_name(
    targetName='BRAF',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "search": {
      "hits": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_id_description_by_name-parameters)

## 51) 分页查询靶点互作伙伴与证据 {#dnakit.databases.open_targets_get_target_interactions_by_ensembl_id}

- **作用：** 分页查询靶点互作伙伴与证据。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_interactions_by_ensembl_id(*, ensemblId, page=None, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data`, `metadata`, `error`, `error_details` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_interactions_by_ensembl_id(
    ensemblId='ENSG00000141510',
    page={'index': 0, 'size': 10},
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "interactions": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  },
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_interactions_by_ensembl_id-parameters)

## 52) 查询靶点已知安全性风险记录 {#dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id}

- **作用：** 查询靶点已知安全性风险记录。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id(*, ensemblId, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "safetyLiabilities": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id-parameters)

## 53) 查询靶点对不同治疗方式的可开发性 {#dnakit.databases.open_targets_get_target_tractability_by_ensembl_id}

- **作用：** 查询靶点对不同治疗方式的可开发性。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_target_tractability_by_ensembl_id(*, ensemblId, verbose=True)`
- **输入：** 必须提供 `ensemblId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_tractability_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "tractability": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_tractability_by_ensembl_id-parameters)

## 54) 查询包含指定变异的 GWAS 可信集合 {#dnakit.databases.open_targets_get_variant_credible_sets}

- **作用：** 查询包含指定变异的 GWAS 可信集合。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_variant_credible_sets(*, variantId=None, size=None, variant_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `variant` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_variant_credible_sets(
    variantId='10_112998590_C_T',
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": {
    "id": "<string>",
    "rsIds": [
      "<string>"
    ],
    "credibleSets": {
      "count": "<integer>",
      "rows": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_variant_credible_sets-parameters)

## 55) 查询变异坐标、频率与功能后果 {#dnakit.databases.open_targets_get_variant_info}

- **作用：** 查询变异坐标、频率与功能后果。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_get_variant_info(*, variantId, verbose=True)`
- **输入：** 必须提供 `variantId`。
- **输出：** 返回字典，包含 `variant` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_get_variant_info(
    variantId='10_112998590_C_T',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "variant": {
    "id": "<string>",
    "rsIds": [
      "<string>"
    ],
    "chromosome": "<string>",
    "position": "<integer>",
    "referenceAllele": "<string>",
    "alternateAllele": "<string>",
    "hgvsId": "<string>",
    "variantDescription": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_variant_info-parameters)

## 56) 转换疾病与表型的跨数据库标识 {#dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids}

- **作用：** 转换疾病与表型的跨数据库标识。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids(*, inputId, verbose=True)`
- **输入：** 必须提供 `inputId`。
- **输出：** 返回字典，包含 `data` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids(
    inputId='1',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "mapIds": {
      "mappings": "<array>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids-parameters)

## 57) 联合检索疾病、靶点与药物实体 {#dnakit.databases.open_targets_multi_entity_search_by_query_string}

- **作用：** 联合检索疾病、靶点与药物实体。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_multi_entity_search_by_query_string(*, queryString, entityNames=None, page=None, verbose=True)`
- **输入：** 必须提供 `queryString`。
- **输出：** 返回字典，包含 `error`, `error_details` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_multi_entity_search_by_query_string(
    queryString='gene:TP53',
    entityNames=[],
    page={'index': 0, 'size': 10},
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_multi_entity_search_by_query_string-parameters)

## 58) 按疾病或性状查找 GWAS 研究 {#dnakit.databases.open_targets_search_gwas_studies_by_disease}

- **作用：** 按疾病或性状查找 GWAS 研究。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_search_gwas_studies_by_disease(*, diseaseIds=None, enableIndirect=None, disease_name=None, disease=None, trait=None, size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `studies` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_search_gwas_studies_by_disease(
    diseaseIds=['MONDO_0005148'],
    size=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "studies": {
    "count": "<integer>",
    "rows": [
      "<object>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_search_gwas_studies_by_disease-parameters)

## 59) 查询 IntOGen 肿瘤驱动基因关联证据 {#dnakit.databases.open_targets_target_disease_evidence}

- **作用：** 查询 IntOGen 肿瘤驱动基因关联证据。
- **处理方法：** 通过 `OpenTarget` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.open_targets_target_disease_evidence(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.open_targets_target_disease_evidence(
    efoId='MONDO_0005011',
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "evidences": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_target_disease_evidence-parameters)

## 60) 将化合物名称解析为 PubChem CID {#dnakit.databases.pub_chem_get_cid_by_compound_name}

- **作用：** 将化合物名称解析为 PubChem CID。
- **处理方法：** 通过 `PubChem` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pub_chem_get_cid_by_compound_name(*, name=None, compound_name=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `IdentifierList` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pub_chem_get_cid_by_compound_name(
    name='aspirin',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "IdentifierList": {
    "CID": [
      "<integer>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.pub_chem_get_cid_by_compound_name-parameters)

## 61) 按 CID 查询化合物分子性质 {#dnakit.databases.pub_chem_get_compound_properties_by_cid}

- **作用：** 按 CID 查询化合物分子性质。
- **处理方法：** 通过 `PubChem` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.pub_chem_get_compound_properties_by_cid(*, cid, properties=None, verbose=True)`
- **输入：** 必须提供 `cid`。
- **输出：** 返回字典，包含 `PropertyTable` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.pub_chem_get_compound_properties_by_cid(
    cid=2244,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "PropertyTable": {
    "Properties": [
      "<object>"
    ]
  }
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.pub_chem_get_compound_properties_by_cid-parameters)

## 62) 查询通路支持的疾病靶点关联评分 {#dnakit.databases.reactome_disease_target_score}

- **作用：** 查询通路支持的疾病靶点关联评分。
- **处理方法：** 通过 `DiseaseTargetScore` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.reactome_disease_target_score(*, efoId, pageSize=None, verbose=True)`
- **输入：** 必须提供 `efoId`。
- **输出：** 返回字典，包含 `disease_info`, `datasource`, `total_targets_with_scores`, `target_scores` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.reactome_disease_target_score(
    efoId='MONDO_0011996',
    pageSize=100,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "disease_info": {
    "disease_id": "<string>",
    "disease_name": "<string>"
  },
  "datasource": "<string>",
  "total_targets_with_scores": "<integer>",
  "target_scores": [
    {
      "target_symbol": "<string>",
      "target_id": "<string>",
      "datasource": "<string>",
      "score": "<number>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_targets.md#dnakit.databases.reactome_disease_target_score-parameters)

<!-- END GENERATED API: targets -->
