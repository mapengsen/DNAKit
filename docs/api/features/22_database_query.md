# 数据查询

使用统一接口查询 NCBI、Ensembl、ENA、ENCODE 和 UCSC 等公共数据库，从而获取 DNA 序列及相关元数据。

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

- **限制：** 返回结果为不可变 `QueryResult`，包含脱敏后的请求 URL、记录、总数、分页 token、metadata 和 provenance。记录内容及数量取决于查询时的 NCBI 数据；响应字节数、记录数和分页大小均受硬上限约束。

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

- **限制：** Ensembl 请求由 adapter 转为 1-based 闭区间，并在结果中记录两套坐标口径；UCSC 原生使用 0-based start 和半开 end。区域跨度、批量数量、返回记录数和响应大小均有硬上限；实时注释和轨道内容会随提供方版本变化。

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

- **限制：** ENA 查询表达式、NCBI Entrez 语法和 ENCODE filters 均由对应提供方定义。坐标映射和比较比对只返回提供方已有的数据，不保证每个区域都有结果；所有调用均受网络、记录数、区域跨度和响应大小上限约束。

## 4) BLAST、来源识别与远程新颖性 <span id="5-blast"></span>

- **作用：** 提交并轮询 NCBI BLAST 异步任务，返回数据库命中、Identity、Coverage 和可能来源；也可相对所选远程数据库计算基于最佳命中的新颖性指标。
- **API：** `dnakit.search.identify(sequence[必须], wait[可选], config[可选], database[可选], program[可选], hitlist_size[可选], expect[可选], megablast[可选], poll_interval[可选], timeout[可选], progress[可选])`、`dnakit.search.novelty(sequence[必须], wait[可选], config[可选], database[可选], program[可选], hitlist_size[可选], expect[可选], megablast[可选], poll_interval[可选], timeout[可选], progress[可选])`、`dnakit.search.wait_for_blast(job[必须], poll_interval[可选], timeout[可选], progress[可选], config[可选])`。
- **输入：** 必填 DNA 字符串或无 Gap 的 `DNASequence`；必须在 `SearchConfig` 中提供联系邮箱。可选 BLAST 数据库、程序、命中上限、E-value、轮询间隔和总超时。
- **示例代码：**

```python
from dnakit.search import SearchConfig, identify, wait_for_blast

config = SearchConfig(email="researcher@example.org", max_records=20)
job = identify("ACGTACGTACGT", config=config, hitlist_size=20)
result = wait_for_blast(
    job,
    config=config,
    poll_interval=60,
    timeout=900,
)

print(result.query_type, result.provider, len(result.records))
if result.records:
    print(result.records[0]["identity"], result.records[0]["query_coverage"])
print(result.metadata["novelty_score"])
```

- **示例结果：** 以下仅展示稳定的返回结构；实时命中数量和数值由 NCBI BLAST 决定。

```text
sequence_similarity NCBI BLAST <命中记录数>
<identity: 0～1> <query_coverage: 0～1>
<novelty_score: 0～1 或 None>
```

- **限制：** NCBI BLAST 使用异步作业。默认 `identify()`/`novelty()` 只返回 `BlastJob`，不会隐式等待；`wait_for_blast()` 的轮询间隔不得小于 60 秒。远程 `novelty_score` 定义为 `1 - 最大 identity`，必须同时检查 coverage；它不替代 `dnakit.evaluation.evaluate_novelty()` 相对于版本化本地参考库的定义。BLAST 调用遵循 [NCBI BLAST URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html)。

<span id="6"></span>**明确边界**

- 只查询公开数据；dbGaP 等受控数据仍需用户正式授权，DNAKit 不接收或绕过凭据。
- 远程 adapter 的可用性受网络、提供方限流、维护和数据更新影响。单元测试验证请求、解析和上限契约，不把替身响应写成实时数据库结果。
- `QueryResult` 保留提供方原始字段，只做必要的坐标和 BLAST 命中标准化；不自动生成临床结论，也不把 annotation payload 改写成未经证据支持的自然语言解释。
- BLAST 来源识别返回最相似记录、物种和命中统计；只有提供方命中记录包含对应字段时，结果才会附带组装与坐标。
- 需要把查询结果保存为 JSON、JSONL、CSV、TSV 或 XML 时，使用[下载](15_download.md#query-result-export)页中的 `dnakit.download.metadata()`。
