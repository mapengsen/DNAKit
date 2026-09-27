# 模式生物、物种与微生物组查询

查询模式生物、物种分类和微生物组数据。

## 1) 查询物种分类与 TaxID {#search-taxonomy}

- **作用：** 按物种名称、TaxID 或别名获取分类记录，用于确定查询和下载目标。
- **计算方法：** 查询 NCBI Taxonomy，将分类报告整理为 QueryResult。
- **API：** `dnakit.search.taxonomy(query[必须], report[可选], children[可选], include_lineage[可选], ranks[可选], page_size[可选], page_token[可选], config[可选], verbose[可选])`
- **输入：** 物种名称或 TaxID；可选报告类型、子节点、谱系、等级、分页和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import SearchConfig, taxonomy

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)
result = taxonomy("human", page_size=5, config=config, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
taxonomy NCBI Datasets <返回记录数>
```

<!-- BEGIN GENERATED API: organisms -->

<span id="organisms"></span>

## 2) 查询果蝇基因详情 {#dnakit.databases.fly_base_get_gene}

- **作用：** 查询果蝇基因详情。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fly_base_get_gene(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `id`, `symbol`, `name`, `species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene(
    gene_id='FB:FBgn0003996',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "synonyms": [
    "<string>"
  ],
  "so_term": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene-parameters)

## 3) 查询果蝇基因等位变异与表型 {#dnakit.databases.fly_base_get_gene_alleles}

- **作用：** 查询果蝇基因等位变异与表型。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fly_base_get_gene_alleles(*, gene_id, limit=None, page=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `id`, `symbol`, `symbol_text`, `category` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_alleles(
    gene_id='FB:FBgn0000490',
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
    "symbol": "<string>",
    "symbol_text": "<string>",
    "category": "<string>",
    "has_disease": "<boolean>",
    "has_phenotype": "<boolean>",
    "variants": [
      "<object>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_alleles-parameters)

## 4) 查询果蝇基因相关疾病模型 {#dnakit.databases.fly_base_get_gene_disease_models}

- **作用：** 查询果蝇基因相关疾病模型。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fly_base_get_gene_disease_models(*, gene_id, limit=None, page=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `model_id`, `model_name`, `gene_symbol`, `gene_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_disease_models(
    gene_id='FB:FBgn0000490',
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
    "model_id": "<string>",
    "model_name": "<string>",
    "gene_symbol": "<string>",
    "gene_id": "<string>",
    "data_provider": "<string>",
    "diseases": [
      "<object>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_disease_models-parameters)

## 5) 查询果蝇基因的发育阶段与组织表达 {#dnakit.databases.fly_base_get_gene_expression}

- **作用：** 查询果蝇基因的发育阶段与组织表达。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fly_base_get_gene_expression(*, gene_id, limit=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `total_annotations`, `returned`, `annotations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_expression(
    gene_id='FB:FBgn0000490',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total_annotations": "<integer>",
  "returned": "<integer>",
  "annotations": [
    {
      "stage": "<string>",
      "location": "<string>",
      "data_provider": "<string>",
      "reference": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_expression-parameters)

## 6) 查询果蝇基因分子互作与遗传互作 {#dnakit.databases.fly_base_get_gene_interactions}

- **作用：** 查询果蝇基因分子互作与遗传互作。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fly_base_get_gene_interactions(*, gene_id, interaction_type=None, limit=None, page=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `subject_gene_id`, `subject_symbol`, `interactor_gene_id`, `interactor_symbol` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_interactions(
    gene_id='FB:FBgn0000490',
    interaction_type='genetic',
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
    "subject_gene_id": "<string>",
    "subject_symbol": "<string>",
    "interactor_gene_id": "<string>",
    "interactor_symbol": "<string>",
    "interactor_species": "<string>",
    "interaction_type": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_interactions-parameters)

## 7) 查询果蝇基因的跨物种直系同源关系 {#dnakit.databases.fly_base_get_gene_orthologs}

- **作用：** 查询果蝇基因的跨物种直系同源关系。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fly_base_get_gene_orthologs(*, gene_id, limit=None, page=None, stringency=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `subject_gene_id`, `subject_symbol`, `subject_species`, `ortholog_gene_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fly_base_get_gene_orthologs(
    gene_id='FB:FBgn0004647',
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
    "subject_gene_id": "<string>",
    "subject_symbol": "<string>",
    "subject_species": "<string>",
    "ortholog_gene_id": "<string>",
    "ortholog_symbol": "<string>",
    "ortholog_species": "<string>",
    "is_best_score": "<string>",
    "is_best_score_reverse": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_base_get_gene_orthologs-parameters)

## 8) 检索果蝇综合基因组学记录 {#dnakit.databases.fly_mine_search}

- **作用：** 检索果蝇综合基因组学记录。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.fly_mine_search(*, q, size=None, format=None, facet_Category=None, verbose=True)`
- **输入：** 必须提供 `q`。
- **输出：** 返回字典，包含 `totalHits`, `facets`, `results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.fly_mine_search(
    q='Notch',
    size=5,
    format='json',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "totalHits": "<integer>",
  "facets": "<object>",
  "results": [
    {
      "id": "<integer>",
      "type": "<string>",
      "relevance": "<number>",
      "fields": "<object>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.fly_mine_search-parameters)

## 9) 按分类编号查询 GBIF 物种详情 {#dnakit.databases.gbif_get_species}

- **作用：** 按分类编号查询 GBIF 物种详情。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gbif_get_species(*, speciesKey, verbose=True)`
- **输入：** 必须提供 `speciesKey`。
- **输出：** 返回字典，包含 `key`, `scientificName`, `canonicalName`, `vernacularName` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gbif_get_species(
    speciesKey=2435099,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "key": "<integer>",
  "scientificName": "<string>",
  "canonicalName": "<string>",
  "vernacularName": "<string>",
  "authorship": "<string>",
  "rank": "<string>",
  "taxonomicStatus": "<string>",
  "kingdom": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.gbif_get_species-parameters)

## 10) 检索 GBIF 物种名称与分类 {#dnakit.databases.gbif_search_species}

- **作用：** 检索 GBIF 物种名称与分类。
- **处理方法：** 通过 `GBIF` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gbif_search_species(*, query, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 字典的 data 包含 GBIF 物种搜索响应及分页信息。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gbif_search_species(
    query='Arabidopsis thaliana',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 字典的 data 包含 GBIF 物种搜索响应及分页信息。

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.gbif_search_species-parameters)

## 11) 查询 GTDB 物种簇及其基因组 {#dnakit.databases.gtdb_get_species}

- **作用：** 查询 GTDB 物种簇及其基因组。
- **处理方法：** 通过 `GTDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtdb_get_species(*, operation=None, species, max_genomes=None, verbose=True)`
- **输入：** 必须提供 `species`。
- **输出：** 返回字典，包含 `species_name`, `total_genomes`, `genomes_shown`, `genomes` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtdb_get_species(
    operation='get_species',
    species='Escherichia coli',
    max_genomes=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "species_name": "<string>",
  "total_genomes": "<integer>",
  "genomes_shown": "<integer>",
  "genomes": [
    {
      "accession": "<string>",
      "ncbi_org_name": "<string>",
      "ncbi_tax": "<string>",
      "gtdb_species_rep": "<boolean>",
      "ncbi_type_material": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_get_species-parameters)

## 12) 查询 GTDB 分类层级与版本历史 {#dnakit.databases.gtdb_get_taxon_info}

- **作用：** 查询 GTDB 分类层级与版本历史。
- **处理方法：** 通过 `GTDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtdb_get_taxon_info(*, operation=None, taxon, verbose=True)`
- **输入：** 必须提供 `taxon`。
- **输出：** 返回字典，包含 `taxon`, `rank`, `n_genomes`, `higher_ranks` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtdb_get_taxon_info(
    operation='get_taxon_info',
    taxon='f__Lachnospiraceae',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "taxon": "<string>",
  "rank": "<string>",
  "n_genomes": "<integer>",
  "higher_ranks": [
    "<string>"
  ],
  "in_releases": [
    "<string>"
  ],
  "lineage": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_get_taxon_info-parameters)

## 13) 按物种名称检索原核基因组 {#dnakit.databases.gtdb_search_genomes}

- **作用：** 按物种名称检索原核基因组。
- **处理方法：** 通过 `GTDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtdb_search_genomes(*, operation=None, query, page=None, items_per_page=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `query`, `match_type`, `total_results`, `total_results_scope` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtdb_search_genomes(
    operation='search_genomes',
    query='Escherichia',
    items_per_page=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "query": "<string>",
  "match_type": "<string>",
  "total_results": "<integer>",
  "total_results_scope": "<string>",
  "note": "<string>",
  "page": "<integer>",
  "results": [
    {
      "gid": "<string>",
      "accession": "<string>",
      "ncbiOrgName": "<string>",
      "ncbiTaxonomy": "<string>",
      "gtdbTaxonomy": "<string>",
      "isGtdbSpeciesRep": "<boolean>",
      "isNcbiTypeMaterial": "<boolean>"
    }
  ],
  "count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_search_genomes-parameters)

## 14) 按部分名称检索原核分类单元 {#dnakit.databases.gtdb_search_taxon}

- **作用：** 按部分名称检索原核分类单元。
- **处理方法：** 通过 `GTDB` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.gtdb_search_taxon(*, operation=None, query, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回字典，包含 `query`, `matches`, `count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.gtdb_search_taxon(
    operation='search_taxon',
    query='Lachnospiraceae',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "query": "<string>",
  "matches": [
    "<string>"
  ],
  "count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.gtdb_search_taxon-parameters)

## 15) 按 MGI 编号查询小鼠基因详情 {#dnakit.databases.mgi_get_gene}

- **作用：** 按 MGI 编号查询小鼠基因详情。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgi_get_gene(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `id`, `symbol`, `name`, `species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgi_get_gene(
    gene_id='MGI:98834',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "synonyms": [
    "<string>"
  ],
  "so_term": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgi_get_gene-parameters)

## 16) 查询小鼠基因的突变表型注释 {#dnakit.databases.mgi_get_phenotypes}

- **作用：** 查询小鼠基因的突变表型注释。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgi_get_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `gene_symbol`, `gene_id`, `phenotype_statement` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgi_get_phenotypes(
    gene_id='MGI:98834',
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
    "gene_symbol": "<string>",
    "gene_id": "<string>",
    "phenotype_statement": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgi_get_phenotypes-parameters)

## 17) 按名称或关键词检索小鼠基因 {#dnakit.databases.mgi_search_genes}

- **作用：** 按名称或关键词检索小鼠基因。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgi_search_genes(*, query, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `symbol`, `name`, `gene_id`, `category` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgi_search_genes(
    query='Trp53',
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
    "symbol": "<string>",
    "name": "<string>",
    "gene_id": "<string>",
    "category": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgi_search_genes-parameters)

## 18) 查询微生物基因组质量与分类信息 {#dnakit.databases.mgnify_get_genome}

- **作用：** 查询微生物基因组质量与分类信息。
- **处理方法：** 通过 `MGnifyExpanded` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_get_genome(*, genome_id, verbose=True)`
- **输入：** 必须提供 `genome_id`。
- **输出：** 返回字典，包含 `genome_id`, `type`, `taxonomy`, `completeness` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_get_genome(
    genome_id='MGYG000000001',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "genome_id": "<string>",
  "type": "<string>",
  "taxonomy": "<string>",
  "completeness": "<number>",
  "contamination": "<number>",
  "length": "<integer>",
  "num_proteins": "<integer>",
  "gc_content": "<number>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_genome-parameters)

## 19) 查询宏基因组分析的 GO 功能注释 {#dnakit.databases.mgnify_get_go_terms}

- **作用：** 查询宏基因组分析的 GO 功能注释。
- **处理方法：** 通过 `MGnifyExpanded` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_get_go_terms(*, analysis_id, page_size=None, verbose=True)`
- **输入：** 必须提供 `analysis_id`。
- **输出：** 返回记录列表；每条记录包含 `go_id`, `description`, `count`, `category` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_get_go_terms(
    analysis_id='MGYA00585482',
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
    "go_id": "<string>",
    "description": "<string>",
    "count": "<integer>",
    "category": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_go_terms-parameters)

## 20) 查询宏基因组分析的蛋白结构域组成 {#dnakit.databases.mgnify_get_interpro}

- **作用：** 查询宏基因组分析的蛋白结构域组成。
- **处理方法：** 通过 `MGnifyExpanded` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_get_interpro(*, analysis_id, page_size=None, verbose=True)`
- **输入：** 必须提供 `analysis_id`。
- **输出：** 返回记录列表；每条记录包含 `interpro_id`, `description`, `count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_get_interpro(
    analysis_id='MGYA00585482',
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
    "interpro_id": "<string>",
    "description": "<string>",
    "count": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_interpro-parameters)

## 21) 查询 MGnify 宏基因组研究详情 {#dnakit.databases.mgnify_get_study_detail}

- **作用：** 查询 MGnify 宏基因组研究详情。
- **处理方法：** 通过 `MGnifyExpanded` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_get_study_detail(*, study_accession, verbose=True)`
- **输入：** 必须提供 `study_accession`。
- **输出：** 返回字典，包含 `study_id`, `study_name`, `study_abstract`, `bioproject` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_get_study_detail(
    study_accession='MGYS00002008',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "study_id": "<string>",
  "study_name": "<string>",
  "study_abstract": "<string>",
  "bioproject": "<string>",
  "centre_name": "<string>",
  "biomes": [
    "<string>"
  ],
  "is_public": "<boolean>",
  "samples_count": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_study_detail-parameters)

## 22) 查询宏基因组分类组成与丰度 {#dnakit.databases.mgnify_get_taxonomy}

- **作用：** 查询宏基因组分类组成与丰度。
- **处理方法：** 通过 `MGnifyExpanded` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_get_taxonomy(*, analysis_id, rna_type=None, page_size=None, verbose=True)`
- **输入：** 必须提供 `analysis_id`。
- **输出：** 返回记录列表；每条记录包含 `organism_id`, `name`, `rank`, `domain` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_get_taxonomy(
    analysis_id='MGYA00585482',
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
    "organism_id": "<string>",
    "name": "<string>",
    "rank": "<string>",
    "domain": "<string>",
    "count": "<integer>",
    "lineage": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_get_taxonomy-parameters)

## 23) 列出指定宏基因组研究的分析记录 {#dnakit.databases.mgnify_list_analyses}

- **作用：** 列出指定宏基因组研究的分析记录。
- **处理方法：** 通过 `MGnifyAnalyses` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_list_analyses(*, study_accession, size=None, verbose=True)`
- **输入：** 必须提供 `study_accession`。
- **输出：** data 中包含指定研究的 MGnify 分析记录及服务返回的分页信息。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_list_analyses(
    study_accession='MGYS00000596',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data 中包含指定研究的 MGnify 分析记录及服务返回的分页信息。

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_list_analyses-parameters)

## 24) 浏览宏基因组生境分类 {#dnakit.databases.mgnify_list_biomes}

- **作用：** 浏览宏基因组生境分类。
- **处理方法：** 通过 `MGnifyExpanded` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_list_biomes(*, depth=None, page=None, page_size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `biome_id`, `biome_name`, `samples_count` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_list_biomes(
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
    "biome_id": "<string>",
    "biome_name": "<string>",
    "samples_count": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_list_biomes-parameters)

## 25) 检索宏基因组组装与分离株基因组 {#dnakit.databases.mgnify_search_genomes}

- **作用：** 检索宏基因组组装与分离株基因组。
- **处理方法：** 通过 `MGnifyExpanded` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_search_genomes(*, taxonomy=None, page=None, page_size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `genome_id`, `type`, `taxonomy`, `completeness` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_search_genomes(
    page=1,
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
    "genome_id": "<string>",
    "type": "<string>",
    "taxonomy": "<string>",
    "completeness": "<number>",
    "contamination": "<number>",
    "length": "<integer>",
    "num_proteins": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_search_genomes-parameters)

## 26) 按生境与关键词检索宏基因组研究 {#dnakit.databases.mgnify_search_studies}

- **作用：** 按生境与关键词检索宏基因组研究。
- **处理方法：** 通过 `MGnifyStudies` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.mgnify_search_studies(*, biome=None, search=None, size=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** data 中包含 MGnify 研究检索结果及服务返回的分页信息。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.mgnify_search_studies(
    search='human gut',
    size=3,
    verbose=False,
)
print(result)
```

- **示例结果：** data 中包含 MGnify 研究检索结果及服务返回的分页信息。

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.mgnify_search_studies-parameters)

## 27) 检索全球植物名称与分类记录 {#dnakit.databases.powo_search_plants}

- **作用：** 检索全球植物名称与分类记录。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.powo_search_plants(*, q, f=None, perPage=None, cursor=None, verbose=True)`
- **输入：** 必须提供 `q`。
- **输出：** 返回字典，包含 `totalResults`, `page`, `totalPages`, `perPage` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.powo_search_plants(
    q='Rosa canina',
    f='species_f:true',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "totalResults": "<integer>",
  "page": "<integer>",
  "totalPages": "<integer>",
  "perPage": "<integer>",
  "cursor": "<string>",
  "results": [
    {
      "accepted": "<boolean>",
      "author": "<string>",
      "kingdom": "<string>",
      "family": "<string>",
      "name": "<string>",
      "rank": "<string>",
      "snippet": "<string>",
      "url": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.powo_search_plants-parameters)

## 28) 查询大鼠基因的疾病、表型与功能注释 {#dnakit.databases.rgd_get_annotations}

- **作用：** 查询大鼠基因的疾病、表型与功能注释。
- **处理方法：** 通过 `RGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_get_annotations(*, rgd_id, aspect=None, verbose=True)`
- **输入：** 必须提供 `rgd_id`。
- **输出：** 返回“查询大鼠基因的疾病、表型与功能注释”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_get_annotations(
    rgd_id='2218',
    aspect='D',
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

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_annotations-parameters)

## 29) 按 RGD 编号查询大鼠基因详情 {#dnakit.databases.rgd_get_gene}

- **作用：** 按 RGD 编号查询大鼠基因详情。
- **处理方法：** 通过 `RGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_get_gene(*, rgd_id, verbose=True)`
- **输入：** 必须提供 `rgd_id`。
- **输出：** 返回大鼠基因名称、功能描述、基因类型及 Ensembl 交叉引用。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_get_gene(
    rgd_id='2218',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回大鼠基因名称、功能描述、基因类型及 Ensembl 交叉引用。

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_gene-parameters)

## 30) 查询大鼠基因的跨物种直系同源关系 {#dnakit.databases.rgd_get_orthologs}

- **作用：** 查询大鼠基因的跨物种直系同源关系。
- **处理方法：** 通过 `RGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_get_orthologs(*, rgd_id, verbose=True)`
- **输入：** 必须提供 `rgd_id`。
- **输出：** 返回“查询大鼠基因的跨物种直系同源关系”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_get_orthologs(
    rgd_id='2218',
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

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_orthologs-parameters)

## 31) 查询大鼠基因组区间的数量性状位点 {#dnakit.databases.rgd_get_qtls_in_region}

- **作用：** 查询大鼠基因组区间的数量性状位点。
- **处理方法：** 通过 `RGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_get_qtls_in_region(*, chromosome, start, stop, map_key, verbose=True)`
- **输入：** 必须提供 `chromosome`, `start`, `stop`, `map_key`。
- **输出：** 返回记录列表；每条记录包含 `rgd_id`, `symbol`, `name`, `chromosome` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_get_qtls_in_region(
    chromosome='10',
    start=1,
    stop=50000000,
    map_key=360,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "rgd_id": "<integer>",
    "symbol": "<string>",
    "name": "<string>",
    "chromosome": "<string>",
    "lod": "<number>",
    "p_value": "<number>",
    "variance": "<number>",
    "inheritance_type": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_qtls_in_region-parameters)

## 32) 按品系名称查询大鼠品系详情 {#dnakit.databases.rgd_get_strain}

- **作用：** 按品系名称查询大鼠品系详情。
- **处理方法：** 通过 `RGDStrain` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_get_strain(*, symbol, verbose=True)`
- **输入：** 必须提供 `symbol`。
- **输出：** 返回字典，包含 `rgd_id`, `symbol`, `name`, `full_strain_name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_get_strain(
    symbol='SHR',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "rgd_id": "<integer>",
  "symbol": "<string>",
  "name": "<string>",
  "full_strain_name": "<string>",
  "substrain": "<string>",
  "strain_type": "<string>",
  "genetics": "<string>",
  "origin": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_strain-parameters)

## 33) 查询大鼠品系的疾病与表型注释 {#dnakit.databases.rgd_get_strain_annotations}

- **作用：** 查询大鼠品系的疾病与表型注释。
- **处理方法：** 通过 `RGDStrain` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_get_strain_annotations(*, symbol, category=None, verbose=True)`
- **输入：** 必须提供 `symbol`。
- **输出：** 返回记录列表；每条记录包含 `term`, `term_acc`, `category`, `aspect` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_get_strain_annotations(
    symbol='SHR',
    category='disease',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "term": "<string>",
    "term_acc": "<string>",
    "category": "<string>",
    "aspect": "<string>",
    "qualifier": "<string>",
    "evidence": "<string>",
    "data_source": "<string>",
    "ref_rgd_id": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_get_strain_annotations-parameters)

## 34) 按名称或关键词检索 RGD 基因 {#dnakit.databases.rgd_search_genes}

- **作用：** 按名称或关键词检索 RGD 基因。
- **处理方法：** 通过 `RGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_search_genes(*, query, species=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回“按名称或关键词检索 RGD 基因”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_search_genes(
    query='Brca1',
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

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_search_genes-parameters)

## 35) 按关键词与类型检索大鼠品系 {#dnakit.databases.rgd_search_strains}

- **作用：** 按关键词与类型检索大鼠品系。
- **处理方法：** 通过 `RGDStrain` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.rgd_search_strains(*, query=None, strain_type=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `rgd_id`, `symbol`, `name`, `full_strain_name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.rgd_search_strains(
    query='hypertensive',
    strain_type='inbred',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "rgd_id": "<integer>",
    "symbol": "<string>",
    "name": "<string>",
    "full_strain_name": "<string>",
    "substrain": "<string>",
    "strain_type": "<string>",
    "genetics": "<string>",
    "origin": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.rgd_search_strains-parameters)

## 36) 查询酵母基因与基因座概况 {#dnakit.databases.sgd_get_gene}

- **作用：** 查询酵母基因与基因座概况。
- **处理方法：** 通过 `SGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.sgd_get_gene(*, sgd_id, verbose=True)`
- **输入：** 必须提供 `sgd_id`。
- **输出：** 返回字典，包含 `sgd_id`, `display_name`, `gene_name`, `systematic_name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sgd_get_gene(
    sgd_id='S000000259',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "sgd_id": "<string>",
  "display_name": "<string>",
  "gene_name": "<string>",
  "systematic_name": "<string>",
  "locus_type": "<string>",
  "description": "<string>",
  "uniprot_id": "<string>",
  "aliases": [
    "<string>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_gene-parameters)

## 37) 查询酵母基因的 GO 注释与证据 {#dnakit.databases.sgd_get_go_annotations}

- **作用：** 查询酵母基因的 GO 注释与证据。
- **处理方法：** 通过 `SGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.sgd_get_go_annotations(*, sgd_id, verbose=True)`
- **输入：** 必须提供 `sgd_id`。
- **输出：** 返回记录列表；每条记录包含 `go_id`, `go_term`, `go_aspect`, `evidence_code` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sgd_get_go_annotations(
    sgd_id='S000003219',
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
    "go_term": "<string>",
    "go_aspect": "<string>",
    "evidence_code": "<string>",
    "qualifier": "<string>",
    "annotation_type": "<string>",
    "reference": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_go_annotations-parameters)

## 38) 查询酵母基因的遗传与物理互作 {#dnakit.databases.sgd_get_interactions}

- **作用：** 查询酵母基因的遗传与物理互作。
- **处理方法：** 通过 `SGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.sgd_get_interactions(*, sgd_id, verbose=True)`
- **输入：** 必须提供 `sgd_id`。
- **输出：** 返回记录列表；每条记录包含 `interaction_type`, `experiment_type`, `bait_gene`, `hit_gene` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sgd_get_interactions(
    sgd_id='S000003219',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "interaction_type": "<string>",
    "experiment_type": "<string>",
    "bait_gene": "<string>",
    "hit_gene": "<string>",
    "annotation_type": "<string>",
    "reference": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_interactions-parameters)

## 39) 查询酵母突变表型与实验条件 {#dnakit.databases.sgd_get_phenotypes}

- **作用：** 查询酵母突变表型与实验条件。
- **处理方法：** 通过 `SGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.sgd_get_phenotypes(*, sgd_id, verbose=True)`
- **输入：** 必须提供 `sgd_id`。
- **输出：** 返回记录列表；每条记录包含 `phenotype`, `mutant_type`, `experiment_type`, `strain_name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sgd_get_phenotypes(
    sgd_id='S000003219',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "phenotype": "<string>",
    "mutant_type": "<string>",
    "experiment_type": "<string>",
    "strain_name": "<string>",
    "allele": "<string>",
    "chemical": "<string>",
    "reference": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_get_phenotypes-parameters)

## 40) 跨数据类型检索酵母数据库 {#dnakit.databases.sgd_search}

- **作用：** 跨数据类型检索酵母数据库。
- **处理方法：** 通过 `SGD` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.sgd_search(*, query, category=None, limit=None, offset=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `name`, `category`, `description`, `href` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sgd_search(
    query='GAL4',
    category='locus',
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
    "name": "<string>",
    "category": "<string>",
    "description": "<string>",
    "href": "<string>",
    "aliases": [
      "<value>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.sgd_search-parameters)

## 41) 查询植物形态、生长与繁殖特征 {#dnakit.databases.usda_plants_get_characteristics}

- **作用：** 查询植物形态、生长与繁殖特征。
- **处理方法：** 通过 `USDAPlantsCharacteristics` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.usda_plants_get_characteristics(*, symbol, verbose=True)`
- **输入：** 必须提供 `symbol`。
- **输出：** 返回“查询植物形态、生长与繁殖特征”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.usda_plants_get_characteristics(
    symbol='ABBA',
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

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.usda_plants_get_characteristics-parameters)

## 42) 查询植物名称、分类与生活型 {#dnakit.databases.usda_plants_get_profile}

- **作用：** 查询植物名称、分类与生活型。
- **处理方法：** 通过 `USDAPlantsProfile` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.usda_plants_get_profile(*, symbol, verbose=True)`
- **输入：** 必须提供 `symbol`。
- **输出：** 返回植物学名、俗名、分类、生长习性和区域原生状态；精细分布并非所有记录均有。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.usda_plants_get_profile(
    symbol='ABBA',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回植物学名、俗名、分类、生长习性和区域原生状态；精细分布并非所有记录均有。

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.usda_plants_get_profile-parameters)

## 43) 查询线虫基因的组织与发育表达 {#dnakit.databases.worm_base_get_expression}

- **作用：** 查询线虫基因的组织与发育表达。
- **处理方法：** 通过 `WormBase` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.worm_base_get_expression(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `wormbase_id`, `gene_name`, `expressed_in_count`, `expressed_in` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.worm_base_get_expression(
    gene_id='WBGene00006763',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "wormbase_id": "<string>",
  "gene_name": "<string>",
  "expressed_in_count": "<integer>",
  "expressed_in": [
    {
      "term_id": "<string>",
      "term_name": "<string>"
    }
  ],
  "expressed_during": [
    {
      "term_id": "<string>",
      "term_name": "<string>"
    }
  ],
  "subcellular_localization": [
    {
      "term_id": "<string>",
      "term_name": "<string>"
    }
  ],
  "expression_clusters_count": "<integer>",
  "expression_clusters": [
    {
      "cluster_id": "<string>",
      "cluster_label": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.worm_base_get_expression-parameters)

## 44) 查询秀丽隐杆线虫基因详情 {#dnakit.databases.worm_base_get_gene}

- **作用：** 查询秀丽隐杆线虫基因详情。
- **处理方法：** 通过 `WormBase` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.worm_base_get_gene(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `wormbase_id`, `gene_name`, `sequence_name`, `species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.worm_base_get_gene(
    gene_id='WBGene00006763',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "wormbase_id": "<string>",
  "gene_name": "<string>",
  "sequence_name": "<string>",
  "species": "<string>",
  "description": "<string>",
  "gene_type": "<string>",
  "status": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.worm_base_get_gene-parameters)

## 45) 查询线虫基因的实验表型 {#dnakit.databases.worm_base_get_phenotypes}

- **作用：** 查询线虫基因的实验表型。
- **处理方法：** 通过 `WormBase` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.worm_base_get_phenotypes(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `wormbase_id`, `gene_name`, `phenotype_count`, `phenotypes` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.worm_base_get_phenotypes(
    gene_id='WBGene00006763',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "wormbase_id": "<string>",
  "gene_name": "<string>",
  "phenotype_count": "<integer>",
  "phenotypes": [
    {
      "phenotype_id": "<string>",
      "phenotype_name": "<string>",
      "evidence_type": "<string>"
    }
  ],
  "not_observed_count": "<integer>",
  "phenotypes_not_observed": [
    {
      "phenotype_id": "<string>",
      "phenotype_name": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.worm_base_get_phenotypes-parameters)

## 46) 查询爪蟾基因详情与跨库标识 {#dnakit.databases.xenbase_get_gene}

- **作用：** 查询爪蟾基因详情与跨库标识。
- **处理方法：** 通过 `Xenbase` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.xenbase_get_gene(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `gene_id`, `symbol`, `name`, `species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.xenbase_get_gene(
    gene_id='XB-GENE-484286',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "gene_id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "synonyms": [
    "<string>"
  ],
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "so_term": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.xenbase_get_gene-parameters)

## 47) 检索热带爪蟾与非洲爪蟾基因 {#dnakit.databases.xenbase_search_genes}

- **作用：** 检索热带爪蟾与非洲爪蟾基因。
- **处理方法：** 通过 `Xenbase` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.xenbase_search_genes(*, query, species=None, limit=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `gene_id`, `symbol`, `name`, `species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.xenbase_search_genes(
    query='tp53',
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
    "symbol": "<string>",
    "name": "<string>",
    "species": "<string>",
    "synonyms": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.xenbase_search_genes-parameters)

## 48) 查询斑马鱼基因详情 {#dnakit.databases.zfin_get_gene}

- **作用：** 查询斑马鱼基因详情。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.zfin_get_gene(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `id`, `symbol`, `name`, `species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.zfin_get_gene(
    gene_id='ZFIN:ZDB-GENE-990415-8',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "species": {
    "name": "<string>",
    "short_name": "<string>",
    "taxon_id": "<string>",
    "data_provider": "<string>"
  },
  "gene_synopsis": "<string>",
  "automated_gene_synopsis": "<string>",
  "synonyms": [
    "<string>"
  ],
  "so_term": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.zfin_get_gene-parameters)

## 49) 查询斑马鱼基因的发育与组织表达 {#dnakit.databases.zfin_get_gene_expression}

- **作用：** 查询斑马鱼基因的发育与组织表达。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.zfin_get_gene_expression(*, gene_id, limit=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `total_annotations`, `returned`, `annotations` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.zfin_get_gene_expression(
    gene_id='ZFIN:ZDB-GENE-980526-166',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "total_annotations": "<integer>",
  "returned": "<integer>",
  "annotations": [
    {
      "stage": "<string>",
      "location": "<string>",
      "data_provider": "<string>",
      "reference": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.zfin_get_gene_expression-parameters)

## 50) 查询斑马鱼基因突变表型 {#dnakit.databases.zfin_get_gene_phenotypes}

- **作用：** 查询斑马鱼基因突变表型。
- **处理方法：** 通过 `AllianceGenome` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.zfin_get_gene_phenotypes(*, gene_id, limit=None, page=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `gene_symbol`, `gene_id`, `phenotype_statement` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.zfin_get_gene_phenotypes(
    gene_id='ZFIN:ZDB-GENE-990415-8',
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
    "gene_symbol": "<string>",
    "gene_id": "<string>",
    "phenotype_statement": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_organisms.md#dnakit.databases.zfin_get_gene_phenotypes-parameters)

<!-- END GENERATED API: organisms -->
