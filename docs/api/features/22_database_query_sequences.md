# 序列、组装与测序数据查询

查询序列、参考组装和测序项目。

## 1) 查询参考组装与组装 accession {#search-assembly}

- **作用：** 按物种或组装 accession 获取组装记录，用于选择参考基因组。
- **计算方法：** 查询 NCBI Genome，并按参考组装、组装来源和当前版本条件筛选。
- **API：** `dnakit.search.assembly(query[必须], by[可选], reference_only[可选], assembly_source[可选], current_only[可选], page_size[可选], page_token[可选], config[可选], verbose[可选])`
- **输入：** 物种名称、TaxID 或组装 accession；可选匹配方式、参考组装筛选、分页和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import SearchConfig, assembly

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)
result = assembly("human", reference_only=True, page_size=5, config=config, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
assembly NCBI Datasets <返回记录数>
```

## 2) 按基因组区间查询序列 {#search-sequence}

- **作用：** 获取指定基因组区间的 DNA 序列，可指定链方向和上下游扩展。
- **计算方法：** 将 0-based 半开区间转换为 Ensembl 坐标，查询序列并保留坐标信息。
- **API：** `dnakit.search.sequence(species[必须], region[必须], strand[可选], upstream[可选], downstream[可选], mask[可选], progress[可选], config[可选], verbose[可选])`
- **输入：** 物种和一个或多个 chromosome:start-end 区间，坐标为 **0-based 半开区间**；可选链方向、上下游扩展、掩码、进度回调和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import sequence

result = sequence(
    "human",
    ("1:0-100", "X:1000-1100"),
    strand=-1,
    upstream=20,
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
sequence Ensembl <返回记录数>
```

## 3) 查询 ENA 测序运行与文件信息 {#search-reads}

- **作用：** 按测序策略、平台等条件查找 ENA 测序运行，获取元数据及文件信息。
- **计算方法：** 向 ENA Portal 提交查询表达式，按指定字段、偏移量和记录上限返回 QueryResult。
- **API：** `dnakit.search.reads(query[必须], fields[可选], result[可选], limit[可选], offset[可选], config[可选], verbose[可选])`
- **输入：** ENA 查询表达式；可选字段、结果类型、记录上限、偏移量和 SearchConfig。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import reads

result = reads(
    'library_strategy="RNA-Seq" AND instrument_platform="ILLUMINA"',
    limit=20,
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
reads ENA <返回记录数>
```

## 4) 查询跨物种基因组比对区块 {#search-comparative-alignment}

- **作用：** 获取指定基因组区域的跨物种比对区块，用于比较基因组序列。
- **计算方法：** 通过 Ensembl Compara 获取区域比对，按比对方法和物种集合筛选。
- **API：** `dnakit.search.comparative_alignment(species[必须], region[必须], strand[可选], method[可选], species_set_group[可选], display_species[可选], aligned[可选], mask[可选], config[可选], verbose[可选])`
- **输入：** 物种和 **0-based 半开** 基因组区间；可选链方向、方法、物种集合、展示物种、是否保留比对和掩码。需要联网。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
- **示例代码：**

```python
from dnakit.search import comparative_alignment

result = comparative_alignment("human", "1:100000-101000", verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下为返回格式，未执行联网查询；实际记录数由在线数据库决定。

```text
comparative_alignment Ensembl <返回记录数>
```

<!-- BEGIN GENERATED API: sequences -->

<span id="sequences"></span>

## 5) 按分类与关键词检索 ENA 生物样本 {#dnakit.databases.enaportal_search_samples}

- **作用：** 按分类与关键词检索 ENA 生物样本。
- **处理方法：** 通过 `ENAPortal` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.enaportal_search_samples(*, query, limit=None, fields=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `sample_accession`, `sample_alias`, `description`, `tax_id` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.enaportal_search_samples(
    query='tax_tree(562)',
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
    "sample_accession": "<string>",
    "sample_alias": "<string>",
    "description": "<string>",
    "tax_id": "<string>",
    "scientific_name": "<string>",
    "first_public": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.enaportal_search_samples-parameters)

## 6) 检索 ENA 测序研究项目 {#dnakit.databases.enaportal_search_studies}

- **作用：** 检索 ENA 测序研究项目。
- **处理方法：** 通过 `ENAPortal` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.enaportal_search_studies(*, query, limit=None, fields=None, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `study_accession`, `study_title`, `center_name`, `first_public` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.enaportal_search_studies(
    query='description="SARS-CoV-2"',
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
    "study_title": "<string>",
    "center_name": "<string>",
    "first_public": "<string>",
    "description": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.enaportal_search_studies-parameters)

## 7) 获取基因、转录本或区间序列 {#dnakit.databases.ensembl_get_sequence}

- **作用：** 获取基因、转录本或区间序列。
- **处理方法：** 通过 `Ensembl` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_get_sequence(*, id, type=None, species=None, multiple_sequences=None, verbose=True)`
- **输入：** 必须提供 `id`。
- **输出：** 返回记录列表；每条记录包含 `id`, `seq`, `desc`, `molecule` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_get_sequence(
    id='ENSG00000012048',
    type='genomic',
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
    "seq": "<string>",
    "desc": "<string>",
    "molecule": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ensembl_get_sequence-parameters)

## 8) 按 Ensembl 编号获取核酸或蛋白质序列 {#dnakit.databases.ensembl_seq_get_id_sequence}

- **作用：** 按 Ensembl 编号获取核酸或蛋白质序列。
- **处理方法：** 通过 `EnsemblSequence` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_seq_get_id_sequence(*, ensembl_id, type=None, verbose=True)`
- **输入：** 必须提供 `ensembl_id`。
- **输出：** 返回字典，包含 `ensembl_id`, `molecule`, `description`, `sequence_length` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_seq_get_id_sequence(
    ensembl_id='ENSP00000269305',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "ensembl_id": "<string>",
  "molecule": "<string>",
  "description": "<string>",
  "sequence_length": "<integer>",
  "sequence": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ensembl_seq_get_id_sequence-parameters)

## 9) 按基因组坐标获取 Ensembl 参考序列 {#dnakit.databases.ensembl_seq_get_region_sequence}

- **作用：** 按基因组坐标获取 Ensembl 参考序列。
- **处理方法：** 通过 `EnsemblSequence` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ensembl_seq_get_region_sequence(*, region, species=None, verbose=True)`
- **输入：** 必须提供 `region`。
- **输出：** 返回字典，包含 `species`, `region`, `molecule`, `description` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ensembl_seq_get_region_sequence(
    region='17:7668421-7668520',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "species": "<string>",
  "region": "<string>",
  "molecule": "<string>",
  "description": "<string>",
  "sequence_length": "<integer>",
  "sequence": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ensembl_seq_get_region_sequence-parameters)

## 10) 将 GenBank UID 转换为序列编号 {#dnakit.databases.ncbi_fetch_accessions}

- **作用：** 将 GenBank UID 转换为序列编号。
- **处理方法：** 通过 `NCBINucleotideSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbi_fetch_accessions(*, operation=None, uids, verbose=True)`
- **输入：** 必须提供 `uids`。
- **输出：** 返回“将 GenBank UID 转换为序列编号”的列表；元素结构见下方示例。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbi_fetch_accessions(
    uids=['545778205'],
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

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_fetch_accessions-parameters)

## 11) 按 NCBI 序列编号获取核酸记录 {#dnakit.databases.ncbi_get_sequence}

- **作用：** 按 NCBI 序列编号获取核酸记录。
- **处理方法：** 通过 `NCBINucleotideSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbi_get_sequence(*, operation=None, accession, format=None, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字符串。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbi_get_sequence(
    accession='U00096',
    format='fasta',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
"<string>"
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_get_sequence-parameters)

## 12) 按物种、基因或菌株检索核酸记录 {#dnakit.databases.ncbi_search_nucleotide}

- **作用：** 按物种、基因或菌株检索核酸记录。
- **处理方法：** 通过 `NCBINucleotideSearch` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbi_search_nucleotide(*, operation=None, organism=None, gene=None, strain=None, keywords=None, seq_type=None, query=None, limit=None, sort=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `uids`, `count`, `returned`, `search_term` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbi_search_nucleotide(
    organism='Escherichia coli',
    strain='K-12',
    keywords='complete genome',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "uids": [
    "<string>"
  ],
  "count": "<integer>",
  "returned": "<integer>",
  "search_term": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_search_nucleotide-parameters)

## 13) 按研究与实验条件检索 SRA 测序运行 {#dnakit.databases.ncbi_sra_search_runs}

- **作用：** 按研究与实验条件检索 SRA 测序运行。
- **处理方法：** 通过 `NCBISRA` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbi_sra_search_runs(*, operation=None, study=None, organism=None, strategy=None, platform=None, source=None, query=None, limit=None, sort=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `uids`, `count`, `returned`, `search_term` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbi_sra_search_runs(
    operation='search',
    organism='Homo sapiens',
    strategy='RNA-Seq',
    limit=5,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "uids": [
    "<string>"
  ],
  "count": "<integer>",
  "returned": "<integer>",
  "search_term": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbi_sra_search_runs-parameters)

## 14) 查询基因组组装的版本与元数据 {#dnakit.databases.ncbidatasets_get_genome_assembly}

- **作用：** 查询基因组组装的版本与元数据。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_get_genome_assembly(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字典，包含 `accession`, `paired_accession`, `source_database`, `organism_name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_genome_assembly(
    accession='GCF_000005845.2',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "accession": "<string>",
  "paired_accession": "<string>",
  "source_database": "<string>",
  "organism_name": "<string>",
  "tax_id": "<integer>",
  "strain": "<string>",
  "assembly_name": "<string>",
  "assembly_level": "<string>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_get_genome_assembly-parameters)

## 15) 获取组装中染色体与支架的序列报告 {#dnakit.databases.ncbidatasets_get_sequence_reports}

- **作用：** 获取组装中染色体与支架的序列报告。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_get_sequence_reports(*, accession, page_size=None, page_token=None, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回记录列表；每条记录包含 `chr_name`, `sequence_name`, `role`, `location_type` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_sequence_reports(
    accession='GCF_000005845.2',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "chr_name": "<string>",
    "sequence_name": "<string>",
    "role": "<string>",
    "location_type": "<string>",
    "refseq_accession": "<string>",
    "genbank_accession": "<string>",
    "length": "<integer>",
    "gc_percent": "<number>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_get_sequence_reports-parameters)

## 16) 按分类编号查询物种谱系与统计 {#dnakit.databases.ncbidatasets_get_taxonomy}

- **作用：** 按分类编号查询物种谱系与统计。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_get_taxonomy(*, tax_id, verbose=True)`
- **输入：** 必须提供 `tax_id`。
- **输出：** 返回字典，包含 `tax_id`, `organism_name`, `genbank_common_name`, `rank` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_get_taxonomy(
    tax_id='9606',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "tax_id": "<integer>",
  "organism_name": "<string>",
  "genbank_common_name": "<string>",
  "rank": "<string>",
  "blast_name": "<string>",
  "lineage": [
    "<integer>"
  ],
  "lineage_names": [
    {
      "tax_id": "<integer>",
      "organism_name": "<string>",
      "rank": "<string>"
    }
  ],
  "children": [
    "<integer>"
  ]
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_get_taxonomy-parameters)

## 17) 列出指定分类单元的基因组组装 {#dnakit.databases.ncbidatasets_list_genomes_by_taxon}

- **作用：** 列出指定分类单元的基因组组装。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_list_genomes_by_taxon(*, taxon, limit=None, reference_only=None, verbose=True)`
- **输入：** 必须提供 `taxon`。
- **输出：** 返回记录列表；每条记录包含 `accession`, `paired_accession`, `source_database`, `organism_name` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_list_genomes_by_taxon(
    taxon='Saccharomyces cerevisiae',
    limit=5,
    reference_only=True,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "accession": "<string>",
    "paired_accession": "<string>",
    "source_database": "<string>",
    "organism_name": "<string>",
    "tax_id": "<integer>",
    "strain": "<string>",
    "assembly_name": "<string>",
    "assembly_level": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_list_genomes_by_taxon-parameters)

## 18) 按部分名称查找物种分类候选 {#dnakit.databases.ncbidatasets_suggest_taxonomy}

- **作用：** 按部分名称查找物种分类候选。
- **处理方法：** 通过 `NCBIDatasets` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.ncbidatasets_suggest_taxonomy(*, query, verbose=True)`
- **输入：** 必须提供 `query`。
- **输出：** 返回记录列表；每条记录包含 `scientific_name`, `tax_id`, `common_name`, `rank` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.ncbidatasets_suggest_taxonomy(
    query='drosophila',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
[
  {
    "scientific_name": "<string>",
    "tax_id": "<string>",
    "common_name": "<string>",
    "rank": "<string>",
    "group_name": "<string>",
    "matched_term": "<string>"
  }
]
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.ncbidatasets_suggest_taxonomy-parameters)

## 19) 按实验编号查询 SRA 测序元数据 {#dnakit.databases.sra_get_experiment}

- **作用：** 按实验编号查询 SRA 测序元数据。
- **处理方法：** 通过 `SRA` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.sra_get_experiment(*, accession, verbose=True)`
- **输入：** 必须提供 `accession`。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sra_get_experiment(
    accession='SRX515811',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "uid": "<string>",
    "title": "<string>",
    "organism": "<string>",
    "taxid": "<string>",
    "platform": "<string>",
    "library_strategy": "<string>",
    "library_source": "<string>",
    "library_selection": "<string>"
  },
  "metadata": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.sra_get_experiment-parameters)

## 20) 检索 SRA 测序实验 {#dnakit.databases.sra_search_experiments}

- **作用：** 检索 SRA 测序实验。
- **处理方法：** 通过 `SRA` 按输入条件检索或提交请求，并解析响应。
- **API：** `dnakit.databases.sra_search_experiments(*, query=None, organism=None, library_strategy=None, platform=None, limit=None, verbose=True)`
- **输入：** 按示例选择查询条件；可选参数使用接口默认值。
- **输出：** 返回字典，包含 `data`, `metadata` 等字段。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.sra_search_experiments(
    query='BRCA1',
    organism='Homo sapiens',
    library_strategy='RNA-Seq',
    limit=3,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回结构如下，实际记录随输入和数据源变化。

返回结构节选；`<…>` 表示类型占位符，非本次查询的实测数据。

```json
{
  "data": {
    "total": "<integer>",
    "returned": "<integer>",
    "query_used": "<string>",
    "experiments": [
      "<object>"
    ]
  },
  "metadata": "<object>"
}
```

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.databases.sra_search_experiments-parameters)

## 21) 批量将序列 UID 解析为带版本的编号 {#dnakit.search.ncbi.accessions}

- **作用：** 把 NCBI 序列 UID 转为版本化 accession，保留输入顺序和重复项。
- **处理方法：** 用 NCBI ESummary 将 UID 转换成带版本的 accession。
- **API：** `dnakit.search.ncbi.accessions(uids, *, database='nuccore', config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)`
- **输入：** 必须提供 `uids`。
- **输出：** tuple[str, ...]；空输入返回空元组，缺失 UID 或版本号抛出 QueryError。
- **示例代码：**

```python
import dnakit

result = dnakit.search.ncbi.accessions(
    uids=['556503834'],
    verbose=False,
)
print(result)
```

- **示例结果：** tuple[str, ...]；空输入返回空元组，缺失 UID 或版本号抛出 QueryError。

[完整参数与使用说明](../../faq/api/22_database_query_sequences.md#dnakit.search.ncbi.accessions-parameters)

<!-- END GENERATED API: sequences -->
