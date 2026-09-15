# 数据查询

使用统一接口查询 NCBI、Ensembl、ENA、ENCODE 和 UCSC 等公共数据库，从而获取 DNA 序列及相关元数据。

## 新增数据库接口入口

新增的 `dnakit.databases` 提供 455 个公开函数；`dnakit.encode` 提供 20 个独立的实验查询、下载和管理接口。它们的完整参数与返回值已逐项收录：

- [基因、序列与测序数据](../native/sequences.md) · [ENCODE](../native/encode.md)
- [遗传变异与 GWAS](../native/variants.md) · [疾病、表型与本体](../native/phenotypes.md)
- [表达与调控](../native/regulation.md) · [物种与微生物组](../native/organisms.md)
- [通路、互作与结构](../native/pathways.md) · [靶点与药物证据](../native/targets.md)
- [文献与研究检索](../native/literature.md) · [接口搜索与分发](../native/execution.md)

完整检索入口见[新增 484 个函数索引](../native/index.md)，从搜索函数到实际调用的示例见[调用教程](../native/usage.md)。下方继续介绍原有 `dnakit.search` 的统一查询结果与扩充后的序列检索。

<span id="1"></span>**提供方与入口**

- **作用：** 汇总 DNAKit 可访问的公共数据库，并统一超时、重试、缓存、限速及返回对象，使不同来源的查询结果具有一致的 provenance 和错误处理。
- **API：** `dnakit.search.SearchConfig(timeout[可选], max_response_bytes[可选], max_records[可选], api_key[可选], email[可选], tool[可选])`、`dnakit.search.QueryResult`。
- **输入：** 可选超时时间、最大响应字节数、最大记录数以及 NCBI 联系信息；具体查询条件由后续查询 API 提供。

| 数据来源            | 主要入口                                                                           | 用途                                                         |
| ------------------- | ---------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| NCBI Datasets       | `taxonomy()`、`assembly()`、`gene()`、`virus()`                            | 分类、组装、基因和病毒查询                                   |
| NCBI Entrez         | `entrez()` 及 project/sample/variant/expression/literature 包装                  | accession、BioProject、BioSample、ClinVar、GEO、PubMed、版本 |
| NCBI BLAST          | `submit_blast()`、`identify()`、`novelty()`                                  | 异步相似性、来源线索和远程新颖性检查                         |
| Ensembl REST        | `sequence()`、`transcripts()`、`annotation()`、`variant()`、`homology()` | 坐标序列、转录本、注释、VEP、同源和比较基因组                |
| ENA Portal          | `reads()`、`ena_project()`、`ena_sample()`                                   | Study/Sample/Experiment/Run/Analysis 和公开文件记录          |
| ENCODE Portal       | `encode_search()`                                                                | 表观基因组实验和文件记录                                     |
| UCSC REST/Downloads | `ucsc_*()`                                                                       | 染色体、序列、CpG/重复/调控/保守性轨道及文件目录             |

这些入口依据 [NCBI Datasets REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)、[NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/)、[Ensembl REST](https://rest.ensembl.org/)、[ENA Portal API](https://www.ebi.ac.uk/ena/portal/api)、[ENCODE REST API](https://www.encodeproject.org/help/rest-api/) 和 [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html) 实现。

## 按物种、基因和菌株检索序列

`dnakit.search.ncbi.entrez(database, term=None, *, organism=None, gene=None, strain=None,
keywords=None, seq_type=None, source="all", retmax=20, retstart=0, config=None)`
保留原来的 `entrez("nuccore", "查询表达式")` 用法，也支持结构化序列条件。
可选条件与 `term` 通过 AND 组合；`organism/gene/strain` 为字面值，`keywords/term` 支持 Entrez 表达式。
结构化过滤适用于 `nuccore` 和 `protein`；菌株通过 All Fields 查询，需核对记录中的实际菌株。

```python
import dnakit

result = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(result.metadata["uids"])
print(result.total_count, accessions)
for summary in result.records:
    print(summary)
```

`seq_type` 支持 `complete_genome/genomic/mrna/ncrna/refseq/protein`；蛋白序列使用
`database="protein"`。`source="refseq"` 可与分子类型组合；`source="genbank"` 排除 RefSeq，
不代表实验验证等级。空结果不会自动扩大查询。

`accessions(uids, *, database="nuccore", config=None)` 每批最多解析 200 个 UID，
按输入顺序返回版本化 accession 元组，保留重复项；空输入返回空元组，不发请求；缺失 UID 报错。
`accession(identifier)` 仍用于查询一个 accession 的摘要。
`QueryResult.records` 是摘要，`metadata["uids"]` 是 UID，`total_count` 是总命中数；
`metadata["warnings"]` 保留 NCBI 的查询提示。这里不返回原始序列，文件下载见
[按 accession 下载序列](15_download.md#accession-sequences)。

MCP 中的对应名称为 `dnakit.search.entrez` 和 `dnakit.search.accessions`，
新增参数会出现在工具描述的输入 Schema 中。

## 1) 基础查询 <span id="2"></span>

- **作用：** 按物种、TaxID、组装 accession、GeneID 或基因符号查询 NCBI Taxonomy、Genome 和 Gene，返回标准化标识、摘要及来源链接，用于解析后续下载目标。
- **API：** `dnakit.search.taxonomy(query[必须], report[可选], children[可选], include_lineage[可选], ranks[可选], page_size[可选], page_token[可选], config[可选])`、`dnakit.search.assembly(query[必须], by[可选], reference_only[可选], assembly_source[可选], current_only[可选], page_size[可选], page_token[可选], config[可选])`、`dnakit.search.gene(query[必须], taxon[可选], by[可选], report[可选], gene_types[可选], page_size[可选], page_token[可选], config[可选])`。
- **输入：** 必填物种名称、TaxID、组装 accession、GeneID 或基因名称；可选报告类型、物种范围、分页参数和统一 `SearchConfig`。基因符号查询必须同时提供 `taxon`。
- **示例代码：**

```python
from dnakit.search import SearchConfig, assembly, gene, taxonomy

config = SearchConfig(timeout=30, max_response_bytes=20_000_000, max_records=100)

taxon = taxonomy("human", page_size=5, config=config)
genomes = assembly("human", reference_only=True, page_size=5, config=config)
brca2 = gene("BRCA2", taxon="human", report="product", page_size=5, config=config)

for result in (taxon, genomes, brca2):
    print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下仅展示稳定的返回结构；实时记录数由查询时的数据库决定。

```text
taxonomy NCBI Datasets <返回记录数>
assembly NCBI Datasets <返回记录数>
gene NCBI Datasets <返回记录数>
```

## 2) 坐标、转录本与区域注释 <span id="3"></span>

- **作用：** 按组装和基因组坐标查询区域序列、转录本、基因注释及 UCSC 轨道命中，返回绑定坐标系的结构化记录，用于局部区域注释。
- **API：** `dnakit.search.sequence(species[必须], region[必须], strand[可选], upstream[可选], downstream[可选], mask[可选], progress[可选], config[可选])`、`dnakit.search.transcripts(query[必须], species[可选], by[可选], expand[可选], config[可选])`、`dnakit.search.annotation(species[必须], region[必须], features[可选], strand[可选], biotype[可选], config[可选])`、`dnakit.search.ucsc_track_data(genome[必须], track[必须], region[必须], max_items[可选], config[可选])`。
- **输入：** 必填物种或 UCSC genome 名称、`chromosome:start-end` 区域以及查询目标；可选链方向、上下游扩展、注释类型、轨道名称、进度回调和资源配置。DNAKit 输入坐标统一为 **0-based 半开区间**。
- **示例代码：**

```python
from dnakit.search import annotation, sequence, transcripts, ucsc_track_data

region = sequence(
    "human",
    ("1:0-100", "X:1000-1100"),
    strand=-1,
    upstream=20,
)

gene_model = transcripts("BRCA2", species="human")
features = annotation("human", "13:32315000-32320000")
cpg = ucsc_track_data("hg38", "cpgIslandExt", "chr13:32315000-32320000")

for result in (region, gene_model, features, cpg):
    print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下仅展示稳定的返回结构；实时记录数由查询时的数据库决定。

```text
sequence Ensembl <返回记录数>
transcripts Ensembl <返回记录数>
annotation Ensembl <返回记录数>
annotation UCSC <返回记录数>
```

批量 Ensembl 序列查询支持进度回调。下面使用项目已有的 Rich 依赖显示进度条：

```python
from rich.progress import Progress

from dnakit.search import sequence

regions = ("1:0-100", "2:0-100", "X:0-100")
with Progress() as bar:
    task = bar.add_task("查询 Ensembl", total=len(regions))

    def update(event):
        bar.update(task, completed=event.completed, description=event.item)

    result = sequence("human", regions, progress=update)
```

## 3) 测序、表达、调控与比较基因组 <span id="4"></span>

- **作用：** 查询 ENA 测序运行、GEO 表达数据、ENCODE 调控记录、跨组装坐标映射及比较基因组比对，返回可筛选 metadata 和文件链接。
- **API：** `dnakit.search.reads(query[必须], fields[可选], result[可选], limit[可选], offset[可选], config[可选])`、`dnakit.search.expression(term[必须], retmax[可选], config[可选])`、`dnakit.search.encode_search(object_type[可选], search_term[可选], filters[可选], limit[可选], config[可选])`、`dnakit.search.map_coordinates(species[必须], region[必须], source_assembly[必须], target_assembly[必须], strand[可选], config[可选])`、`dnakit.search.comparative_alignment(species[必须], region[必须], strand[可选], method[可选], species_set_group[可选], display_species[可选], aligned[可选], mask[可选], config[可选])`。
- **输入：** 必填对应提供方的查询表达式或基因组区域；坐标映射还必须提供源 assembly 和目标 assembly。可选返回字段、过滤条件、记录上限、比较方法和统一 `SearchConfig`。
- **示例代码：**

```python
from dnakit.search import (
    comparative_alignment,
    encode_search,
    expression,
    map_coordinates,
    reads,
)

runs = reads(
    'library_strategy="RNA-Seq" AND instrument_platform="ILLUMINA"',
    limit=20,
)
geo = expression("single cell RNA sequencing[All Fields]", retmax=20)
encode_files = encode_search(
    object_type="File",
    filters={"file_format": "bigWig", "status": "released"},
    limit=20,
)
mapped = map_coordinates(
    "human",
    "1:100000-101000",
    source_assembly="GRCh38",
    target_assembly="GRCh37",
)
alignment = comparative_alignment("human", "1:100000-101000")

for result in (runs, geo, encode_files, mapped, alignment):
    print(result.query_type, result.provider, len(result.records))
```

- **示例结果：** 以下仅展示稳定的返回结构；实时记录数由查询时的数据库决定。

```text
reads ENA <返回记录数>
expression NCBI Entrez <返回记录数>
regulation ENCODE <返回记录数>
coordinate_mapping Ensembl <返回记录数>
comparative_alignment Ensembl <返回记录数>
```
