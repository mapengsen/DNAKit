# 原始 DNA 的基因名称与区间注释

将基因组 DNA 映射到物种参考序列，返回重叠的基因；安装与详细约定见[FAQ：详细说明](../../faq/tutorial_integrations.md#gene-annotation)。

## 1) 基因名称与区间注释 {#annotate-genes}

- **作用：** 返回基因名称、ID、方向及其在输入序列中的区间。
- **处理方法：** 用 minimap2 比对基因组 DNA，再将参考 GFF3 的已有基因注释映射到输入序列。
- **API：** `dnakit.annotation.annotate_genes(sequence, *, species, assembly=None, sequence_id="query", config=None, ...)`
- **输入：** 一条基因组 DNA 和物种名或 TaxID；需要外部程序 `minimap2`，首次调用可能下载参考。适用于已有参考注释映射，不用于从头预测或拼接后的转录本。
- **输出：** 返回 `GeneAnnotationResult`，示例展示 `reference`、`gene_records`、`segments`。
- **示例代码：**

```python
from dnakit.annotation import GeneAnnotationConfig, annotate_genes
from dnakit.io import read

with read("query.fasta", format="fasta") as source:
    record = next(iter(source))

result = annotate_genes(
    record.sequence.symbols,
    species="Homo sapiens",  # 也接受 NCBI TaxID，例如 9606
    assembly=None,
    sequence_id=record.id,
    config=GeneAnnotationConfig(
        cache_dir=".dnakit/gene_annotation",
        threads=4,
    ),
)
print(result.reference["assembly"])
for gene in result.gene_records:
    print(gene.to_dict())
for segment in result.segments:
    print(segment.to_dict())
```

- **示例结果：** 实际参考组装、基因记录及 `gene` / `intergenic` / `unknown` 分段；内容由 `query.fasta` 和参考决定，缺失名称为 `null`。Python 区间为 0-based 半开，`to_dict()` 与 CLI JSON 为 1-based 闭区间。

## 2) 按基因符号或 GeneID 查询基因记录 {#search-gene}

- **作用：** 获取 NCBI Gene 的基因标识、摘要或产物报告，解析基因查询目标。
- **计算方法：** 使用 NCBI Datasets 查询基因，并按报告类型整理为 QueryResult。
- **API：** `dnakit.search.gene(query[必须], taxon[可选], by[可选], report[可选], gene_types[可选], page_size[可选], page_token[可选], config[可选], verbose[可选])`
- **输入：** GeneID 或基因符号；基因符号查询必须提供 taxon。可选匹配方式、报告类型、基因类型、分页和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import SearchConfig, gene

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)
result = gene("BRCA2", taxon="human", report="product", page_size=5, config=config, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
gene NCBI Datasets <返回记录数>
```

## 3) 按基因查询转录本模型 {#search-transcripts}

- **作用：** 按基因名称或标识获取 Ensembl 转录本模型及关联注释。
- **计算方法：** 解析 Ensembl 查询标识，并根据 expand 参数获取展开的转录本信息。
- **API：** `dnakit.search.transcripts(query[必须], species[可选], by[可选], expand[可选], config[可选], verbose[可选])`
- **输入：** 基因名称或标识；可选物种、匹配方式、是否展开和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import transcripts

result = transcripts("BRCA2", species="human", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
transcripts Ensembl <返回记录数>
```

## 4) 查询基因组区间的注释特征 {#search-annotation}

- **作用：** 获取与指定区域重叠的基因等注释特征，用于区域注释。
- **计算方法：** 转换 DNAKit 区间坐标，查询 Ensembl overlap 接口并保留坐标信息。
- **API：** `dnakit.search.annotation(species[必须], region[必须], features[可选], strand[可选], biotype[可选], config[可选], verbose[可选])`
- **输入：** 物种和 **0-based 半开** 区间；可选特征类型、链方向、biotype 和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import annotation

result = annotation("human", "13:32315000-32320000", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
annotation Ensembl <返回记录数>
```

## 5) 查询 UCSC 轨道的区间命中 {#search-ucsc-track-data}

- **作用：** 获取指定 UCSC 轨道与基因组区间重叠的记录，例如 CpG 岛。
- **计算方法：** 按 genome、track 和区域调用 UCSC API，将轨道命中整理为 QueryResult。
- **API：** `dnakit.search.ucsc_track_data(genome[必须], track[必须], region[必须], max_items[可选], config[可选], verbose[可选])`
- **输入：** UCSC genome 名称、轨道名称和 **0-based 半开** 区间；可选返回上限 max_items 和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import ucsc_track_data

result = ucsc_track_data("hg38", "cpgIslandExt", "chr13:32315000-32320000", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
annotation UCSC <返回记录数>
```

## 6) 在参考组装之间映射区间坐标 {#search-map-coordinates}

- **作用：** 将一个参考组装上的基因组区间映射到另一个组装。
- **计算方法：** 调用 Ensembl 坐标映射接口，并在返回结果中保留源组装、目标组装与坐标约定。
- **API：** `dnakit.search.map_coordinates(species[必须], region[必须], source_assembly[必须], target_assembly[必须], strand[可选], config[可选], verbose[可选])`
- **输入：** 物种、**0-based 半开** 区间、源组装和目标组装；可选链方向和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import map_coordinates

result = map_coordinates(
    "human",
    "1:100000-101000",
    source_assembly="GRCh38",
    target_assembly="GRCh37",
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
coordinate_mapping Ensembl <返回记录数>
```

<!-- BEGIN GENERATED API: gene-annotation -->

<span id="gene-annotation"></span>

## 7) 查询基因家族的跨物种系统发育树 {#dnakit.databases.ensembl_compara_get_gene_tree}

- **作用：** 查询基因家族的跨物种系统发育树。
- **处理方法：** 通过 `EnsemblCompara` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_compara_get_gene_tree(*, gene, species=None, verbose=True)`
- **输入：** 必须提供 `gene`。
- **输出：** 返回字典，包含 `tree_id`, `newick`, `rooted`, `members` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_compara_get_gene_tree(
    gene='BRCA1',
    species='human',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "tree_id": "<string>",
  "newick": "<string>",
  "rooted": "<boolean>",
  "members": [
    "<object>"
  ],
  "total_members": "<integer>"
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_compara_get_gene_tree-parameters)

## 8) 查询跨物种直系同源基因 {#dnakit.databases.ensembl_compara_get_orthologues}

- **作用：** 查询跨物种直系同源基因。
- **处理方法：** 通过 `EnsemblCompara` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_compara_get_orthologues(*, gene, species=None, target_species=None, target_taxon=None, verbose=True)`
- **输入：** 必须提供 `gene`。
- **输出：** 返回记录列表；每条记录包含 `source_gene`, `target_gene`, `target_protein`, `target_species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_compara_get_orthologues(
    gene='BRCA1',
    target_species='mouse',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "source_gene": "<string>",
    "target_gene": "<string>",
    "target_protein": "<string>",
    "target_species": "<string>",
    "homology_type": "<string>",
    "taxonomy_level": "<string>",
    "method": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_compara_get_orthologues-parameters)

## 9) 查询基因复制产生的旁系同源基因 {#dnakit.databases.ensembl_compara_get_paralogues}

- **作用：** 查询基因复制产生的旁系同源基因。
- **处理方法：** 通过 `EnsemblCompara` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_compara_get_paralogues(*, gene, species=None, verbose=True)`
- **输入：** 必须提供 `gene`。
- **输出：** 返回记录列表；每条记录包含 `source_gene`, `paralogue_gene`, `paralogue_protein`, `species` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_compara_get_paralogues(
    gene='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "source_gene": "<string>",
    "paralogue_gene": "<string>",
    "paralogue_protein": "<string>",
    "species": "<string>",
    "paralogy_type": "<string>",
    "taxonomy_level": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_compara_get_paralogues-parameters)

## 10) 按基因名称查询同源关系 {#dnakit.databases.ensembl_get_homology}

- **作用：** 按基因名称查询同源关系。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_get_homology(*, species, symbol, target_species=None, target_taxon=None, type=None, sequence=None, aligned=None, verbose=True)`
- **输入：** 必须提供 `species`, `symbol`。
- **输出：** 返回记录列表；每条记录包含 `homologies` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_get_homology(
    species='human',
    symbol='BRCA2',
    target_species='mouse',
    type='orthologues',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "homologies": [
      "<object>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_get_homology-parameters)

## 11) 查询基因组区间重叠的注释特征 {#dnakit.databases.ensembl_get_overlap_features}

- **作用：** 查询基因组区间重叠的注释特征。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_get_overlap_features(*, species, region, feature=None, biotype=None, verbose=True)`
- **输入：** 必须提供 `species`, `region`。
- **输出：** 返回记录列表；每条记录包含 `id`, `external_name`, `biotype`, `start` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_get_overlap_features(
    species='human',
    region='7:140424943-140624564',
    feature='gene',
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
    "external_name": "<string>",
    "biotype": "<string>",
    "start": "<integer>",
    "end": "<integer>"
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_get_overlap_features-parameters)

## 12) 查询基因组区间内的调控元件 {#dnakit.databases.ensembl_get_regulatory_features}

- **作用：** 查询基因组区间内的调控元件。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_get_regulatory_features(*, region, species=None, feature, verbose=True)`
- **输入：** 必须提供 `region`, `feature`。
- **输出：** 返回记录列表；每条记录包含 `id`, `feature_type`, `start`, `end` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_get_regulatory_features(
    region='1:1000000-2000000',
    species='human',
    feature='regulatory',
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
    "feature_type": "<string>",
    "start": "<integer>",
    "end": "<integer>",
    "strand": "<integer>",
    "seq_region_name": "<string>",
    "description": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_get_regulatory_features-parameters)

## 13) 按编号或名称查询 Ensembl 基因信息 {#dnakit.databases.ensembl_lookup_gene}

- **作用：** 按编号或名称查询 Ensembl 基因信息。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_lookup_gene(*, gene_id=None, species=None, id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `id`, `display_name`, `description`, `biotype` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_lookup_gene(
    gene_id='ENSG00000012048',
    species='homo_sapiens',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "id": "<string>",
  "display_name": "<string>",
  "description": "<string>",
  "biotype": "<string>",
  "start": "<integer>",
  "end": "<integer>",
  "strand": "<integer>",
  "seq_region_name": "<string>"
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ensembl_lookup_gene-parameters)

## 14) 按官方基因名称查询 HGNC 注释 {#dnakit.databases.hgnc_fetch_gene_by_symbol}

- **作用：** 按官方基因名称查询 HGNC 注释。
- **处理方法：** 通过 `HGNC` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.hgnc_fetch_gene_by_symbol(*, symbol, verbose=True)`
- **输入：** 必须提供 `symbol`。
- **输出：** 返回字典，包含 `hgnc_id`, `symbol`, `name`, `locus_type` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.hgnc_fetch_gene_by_symbol(
    symbol='TP53',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "hgnc_id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "locus_type": "<string>",
  "locus_group": "<string>",
  "status": "<string>",
  "location": "<string>",
  "entrez_id": "<string>"
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.hgnc_fetch_gene_by_symbol-parameters)

## 15) 批量查询基因注释与标识映射 {#dnakit.databases.my_gene_batch_query}

- **作用：** 批量查询基因注释与标识映射。
- **处理方法：** 通过 `MyGene` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.my_gene_batch_query(*, gene_ids, species=None, fields=None, verbose=True)`
- **输入：** 必须提供 `gene_ids`。
- **输出：** 返回字典，包含 `results` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.my_gene_batch_query(
    gene_ids=['CDK2', 'TP53', 'BRCA1'],
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "results": [
    {
      "_id": "<string>",
      "query": "<string>",
      "symbol": "<string>",
      "name": "<string>",
      "entrezgene": "<integer>",
      "notfound": "<boolean>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.my_gene_batch_query-parameters)

## 16) 按基因编号获取综合注释 {#dnakit.databases.my_gene_get_gene_annotation}

- **作用：** 按基因编号获取综合注释。
- **处理方法：** 通过 `MyGene` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.my_gene_get_gene_annotation(*, gene_id, fields=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `_id`, `symbol`, `name`, `entrezgene` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.my_gene_get_gene_annotation(
    gene_id='1017',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "_id": "<string>",
  "symbol": "<string>",
  "name": "<string>",
  "entrezgene": "<integer>",
  "ensembl": {
    "gene": "<string>",
    "transcript": [
      "<string>"
    ],
    "protein": [
      "<string>"
    ]
  },
  "summary": "<string>",
  "go": {
    "BP": [
      "<value>"
    ],
    "CC": [
      "<value>"
    ],
    "MF": [
      "<value>"
    ]
  },
  "pathway": "<object>"
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.my_gene_get_gene_annotation-parameters)

## 17) 按关键词或标识检索基因 {#dnakit.databases.my_gene_query_genes}

- **作用：** 按关键词或标识检索基因。
- **处理方法：** 通过 `MyGene` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.my_gene_query_genes(*, query=None, species=None, fields=None, size=None, q=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `took`, `total`, `max_score`, `hits` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.my_gene_query_genes(
    query='CDK2',
    species='human',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "took": "<integer>",
  "total": "<integer>",
  "max_score": "<number>",
  "hits": [
    {
      "_id": "<string>",
      "_score": "<number>",
      "symbol": "<string>",
      "name": "<string>",
      "entrezgene": "<integer>",
      "ensembl": "<object>",
      "summary": "<string>"
    }
  ]
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.my_gene_query_genes-parameters)

## 18) 按 NCBI Gene ID 获取基因详情 {#dnakit.databases.ncbidatasets_get_gene}

- **作用：** 按 NCBI Gene ID 获取基因详情。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_get_gene(*, gene_id, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回字典，包含 `gene_id`, `symbol`, `description`, `tax_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_gene(
    gene_id='7157',
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
  "description": "<string>",
  "tax_id": "<string>",
  "taxname": "<string>",
  "common_name": "<string>",
  "type": "<string>",
  "chromosomes": [
    "<string>"
  ]
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbidatasets_get_gene-parameters)

## 19) 按基因名称与物种查询 NCBI 基因 {#dnakit.databases.ncbidatasets_get_gene_by_symbol}

- **作用：** 按基因名称与物种查询 NCBI 基因。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_get_gene_by_symbol(*, symbol, taxon=None, organism=None, species=None, verbose=True)`
- **输入：** 必须提供 `symbol`。
- **输出：** 返回记录列表；每条记录包含 `gene_id`, `symbol`, `description`, `tax_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_gene_by_symbol(
    symbol='BRCA1',
    taxon='human',
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
    "description": "<string>",
    "tax_id": "<string>",
    "taxname": "<string>",
    "type": "<string>",
    "chromosomes": [
      "<string>"
    ],
    "swiss_prot_accessions": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbidatasets_get_gene_by_symbol-parameters)

## 20) 按 NCBI 基因编号查询直系同源基因 {#dnakit.databases.ncbidatasets_get_orthologs}

- **作用：** 按 NCBI 基因编号查询直系同源基因。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_get_orthologs(*, gene_id, page_size=None, verbose=True)`
- **输入：** 必须提供 `gene_id`。
- **输出：** 返回记录列表；每条记录包含 `gene_id`, `symbol`, `description`, `tax_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_orthologs(
    gene_id='7157',
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
    "gene_id": "<string>",
    "symbol": "<string>",
    "description": "<string>",
    "tax_id": "<string>",
    "taxname": "<string>",
    "common_name": "<string>",
    "type": "<string>",
    "chromosomes": [
      "<string>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbidatasets_get_orthologs-parameters)

## 21) 获取 NCBI 基因摘要与功能信息 {#dnakit.databases.ncbigene_get_summary}

- **作用：** 获取 NCBI 基因摘要与功能信息。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbigene_get_summary(*, id, verbose=True)`
- **输入：** 必须提供 `id`。
- **输出：** 返回字典，包含 `result` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbigene_get_summary(
    id='7157',
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

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbigene_get_summary-parameters)

## 22) 按物种和关键词检索 NCBI 基因编号 {#dnakit.databases.ncbigene_search}

- **作用：** 按物种和关键词检索 NCBI 基因编号。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbigene_search(*, term, retmax=None, verbose=True)`
- **输入：** 必须提供 `term`。
- **输出：** 返回字典，包含 `header`, `esearchresult` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbigene_search(
    term='TP53[Symbol] AND Homo sapiens[Organism]',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "header": {
    "type": "<string>",
    "version": "<string>"
  },
  "esearchresult": {
    "count": "<string>",
    "retmax": "<string>",
    "retstart": "<string>",
    "idlist": [
      "<string>"
    ],
    "querytranslation": "<string>"
  }
}
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.ncbigene_search-parameters)

## 23) 查询 MANE 转录本与 RefSeq 对应关系 {#dnakit.databases.tark_get_mane_transcripts}

- **作用：** 查询 MANE 转录本与 RefSeq 对应关系。
- **处理方法：** 通过 `TarkManeTranscripts` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.tark_get_mane_transcripts(*, gene=None, ensembl_id=None, refseq_id=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回记录列表；每条记录包含 `gene`, `mane_type`, `ensembl_transcript`, `refseq_transcript` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.tark_get_mane_transcripts(
    gene='BRCA2',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "gene": "<string>",
    "mane_type": "<string>",
    "ensembl_transcript": "<string>",
    "refseq_transcript": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.tark_get_mane_transcripts-parameters)

## 24) 查询 Ensembl 转录本历史记录 {#dnakit.databases.tark_get_transcript}

- **作用：** 查询 Ensembl 转录本历史记录。
- **处理方法：** 通过 `TarkTranscript` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.tark_get_transcript(*, stable_id, verbose=True)`
- **输入：** 必须提供 `stable_id`。
- **输出：** 返回“查询 Ensembl 转录本历史记录”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.tark_get_transcript(
    stable_id='ENST00000380152',
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

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.tark_get_transcript-parameters)

## 25) 查询基因的转录本与 MANE 注释 {#dnakit.databases.variant_validator_gene2transcripts}

- **作用：** 查询基因的转录本与 MANE 注释。
- **处理方法：** 通过 `Base` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.variant_validator_gene2transcripts(*, gene_symbol, gene=None, gene_name=None, transcript_set=None, genome_build=None, verbose=True)`
- **输入：** 必须提供 `gene_symbol`。
- **输出：** 返回记录列表；每条记录包含 `current_name`, `current_symbol`, `hgnc`, `previous_symbol` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.variant_validator_gene2transcripts(
    gene_symbol='TP53',
    transcript_set='mane',
    genome_build='GRCh38',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "current_name": "<string>",
    "current_symbol": "<string>",
    "hgnc": "<string>",
    "previous_symbol": "<string>",
    "requested_symbol": "<string>",
    "transcripts": [
      "<object>"
    ]
  }
]
```

[完整参数与使用说明](../../faq/api/26_gene_annotation.md#dnakit.databases.variant_validator_gene2transcripts-parameters)

<!-- END GENERATED API: gene-annotation -->
