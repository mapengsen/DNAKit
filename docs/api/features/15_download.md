# 下载

根据链接、物种名称、数据库编号或查询结果下载文件、参考基因组、序列和注释。

本轮新增接口的完整参数见[基因、序列与下载](../native/sequences.md)，包括 `dnakit.download.ncbi.sequence`、`dnakit.download.ena.sequence` 和数据库文件获取入口。ENCODE 的 `download_files`、`batch_download` 及文件追踪见[ENCODE 数据管理](../native/encode.md)。[484 个新增函数完整索引](../native/index.md)可按函数名查找。

## 按 accession 下载序列 <span id="accession-sequences"></span>

`dnakit.download.ncbi.sequence(accession, output_path, *, database="nuccore", format="fasta",
seq_start=None, seq_stop=None, strand=1, config=None, search_config=None)`
直接访问 NCBI EFetch；`dnakit.download.ena.sequence(accession, output_path, *, format="fasta",
config=None, search_config=None)` 访问 ENA Browser。这两个接口无需 ToolUniverse。

```python
import dnakit

result = dnakit.download.ncbi.sequence(
    "NC_000913.3", "ecoli_k12.fasta", format="fasta",
)
print(result.files[0].path)
print(result.metadata["resolved_accession"], result.metadata["sequence_length"])
print(result.manifest_path)
```

| 提供方 | 基础包直接获取 | 需要 `dnakit[bio-formats]` 的格式转换 |
| --- | --- | --- |
| NCBI `nuccore`（DNA/RNA） | `fasta`、`genbank` | `embl`：通过 Biopython 转换 GenBank |
| NCBI `protein` | `fasta`、`genbank`（GenPept） | 不支持 EMBL 蛋白输出 |
| ENA INSDC 核酸记录 | `fasta`、`embl` | `genbank`：通过 Biopython 转换 EMBL |

NCBI 支持 RefSeq（包括 `NZ_`）和 INSDC 序列 accession；ENA 需要 INSDC accession，例如
`U00096.3`。格式转换不代表 accession 已转成另一个数据库的编号。`GCA_/GCF_` 是组装编号，
请用已有的 `genome_package()`。

返回 `DatasetDownloadResult`，保存序列和相邻的 `.manifest.json`，记录请求及解析出的
accession、序列长度、来源 URL、下载/输出格式、转换标记和校验值。正式写入前校验记录的
accession/version、格式和序列长度；异常响应不会作为序列文件保存。这里的校验值用于本地完整性
记录，并非提供方独立签名。已有文件默认受保护；需要覆盖时传入 `DownloadConfig(overwrite=True)`。

可选 `seq_start/seq_stop` 是 **1-based 闭区间**，必须同时指定；`strand=2` 获取反向互补序列，
只适用于 NCBI 核酸序列。蛋白下载目前获取完整记录。此处与 Ensembl 坐标下载的 0-based 半开区间不同。

`DownloadConfig` 管理保存、覆盖和文件大小上限；`SearchConfig` 管理请求超时、响应大小
（默认 20 MB、最高 100 MB）及 NCBI `api_key/email`。显式提供 `search_config` 时以它的
凭据为准，否则使用 `DownloadConfig.api_key`；持久化 URL 会隐藏密钥和邮箱。
较大的组装请使用分包下载。`format="embl"` 的可选依赖可通过
`python -m pip install 'dnakit[bio-formats]'` 一次安装，数据本身仍从上游数据库获取。

MCP 对应名称为 `dnakit.download.ncbi_sequence` 和 `dnakit.download.ena_sequence`；
Python 也支持这两个扁平别名。原有 `dnakit.download.sequence` 继续用于 Ensembl 坐标序列。

### 下载进度

下载默认向终端标准错误流显示开始提示、进度条和完成路径。参考基因组还显示解压进度；
服务器未提供文件大小时显示已下载字节数，不估算百分比。只有下载、校验和文件写入成功后才提示完成。
`download_genome(..., show_progress=False)` 可关闭显示；`dnakit.download` 中接受 `config` 的网络下载接口
使用 `DownloadConfig(show_progress=False)`。显式传入 `progress`（序列下载为 `query_progress`）回调时，
由回调负责显示，避免重复进度条。NCBI 数据包、ENA、ENCODE、UCSC 和通用文件下载共享这一行为。

### 下载失败与重试

基因组包下载、通用 HTTPS 文件下载及共享远程 GET 查询遇到超时、连接中断或 HTTP
408/429/500/502/503/504 时，最多尝试 3 次，重试前分别等待 0.25、0.5 秒。
404、参数错误、校验和不匹配和无效 ZIP 不自动重试；POST 查询也不自动重发。
每次失败的下载从暂存文件重新开始，尚不支持程序退出后的 HTTP Range 断点续传。

服务器提供 `Content-Length` 时会核对实际字节数；基因组包还核对 FASTA 校验和及
报告中的 assembly accession，不能把其他版本误报为当前请求成功。
只有完整下载和校验通过后才安装正式文件。覆盖冲突保留 `FileExistsError`，目录权限错误保留
相应文件系统异常；损坏输入不会替换已有完整结果。瞬时故障重试耗尽后，异常的
`context["attempts"]` 为 3，并提示检查网络后重新运行。

<span id="common-genome-names"></span>**常见物种名称**

本节列出可以直接复制到 `download_genome(query, ...)` 或
`dnakit download-genome` 命令中的常见 NCBI 物种名称。

- **API**：`dnakit.references.download_genome(query[必须], output_dir[必须], overwrite[可选], keep_package[可选], timeout[可选], api_key[可选], api_base_url[可选], chunk_size[可选], max_download_bytes[可选], progress[可选])`、`dnakit download-genome query[必须] output_dir[必须] --overwrite[可选] --keep-package[可选] --api-key[可选] --progress/--no-progress[可选]`。

<span id="1"></span>**人和常用模式动物**

| 中文名称 | 可直接复制的`query`      | 当前参考组装 accession |
| -------- | -------------------------- | ---------------------- |
| 人类     | `Homo sapiens`           | `GCF_000001405.40`   |
| 小鼠     | `Mus musculus`           | `GCF_000001635.27`   |
| 大鼠     | `Rattus norvegicus`      | `GCF_036323735.1`    |
| 猪       | `Sus scrofa`             | `GCF_054392235.1`    |
| 牛       | `Bos taurus`             | `GCF_002263795.3`    |
| 羊       | `Ovis aries`             | `GCF_016772045.2`    |
| 山羊     | `Capra hircus`           | `GCF_001704415.2`    |
| 马       | `Equus caballus`         | `GCF_041296265.1`    |
| 狗       | `Canis lupus familiaris` | `GCF_011100685.1`    |
| 猕猴     | `Macaca mulatta`         | `GCF_049350105.2`    |
| 鸡       | `Gallus gallus`          | `GCF_016699485.2`    |
| 斑马鱼   | `Danio rerio`            | `GCF_049306965.2`    |
| 爪蟾     | `Xenopus tropicalis`     | `GCF_000004195.4`    |

<span id="2"></span>**常见农作物和植物**

| 中文名称     | 可直接复制的`query`           | 当前参考组装 accession |
| ------------ | ------------------------------- | ---------------------- |
| 拟南芥       | `Arabidopsis thaliana`        | `GCF_000001735.4`    |
| 水稻         | `Oryza sativa Japonica Group` | `GCF_034140825.1`    |
| 玉米         | `Zea mays`                    | `GCF_902167145.1`    |
| 小麦         | `Triticum aestivum`           | `GCF_018294505.1`    |
| 大豆         | `Glycine max`                 | `GCF_000004515.6`    |
| 番茄         | `Solanum lycopersicum`        | `GCF_036512215.1`    |
| 马铃薯       | `Solanum tuberosum`           | `GCF_000226075.1`    |
| 陆地棉       | `Gossypium hirsutum`          | `GCF_007990345.1`    |
| 烟草         | `Nicotiana tabacum`           | `GCF_000715075.1`    |
| 毛果杨       | `Populus trichocarpa`         | `GCF_000002775.5`    |
| 阿拉比卡咖啡 | `Coffea arabica`              | `GCF_036785885.1`    |

<span id="3"></span>**其他常用模式生物**

| 中文名称     | 可直接复制的`query`              | 当前参考组装 accession |
| ------------ | ---------------------------------- | ---------------------- |
| 果蝇         | `Drosophila melanogaster`        | `GCF_000001215.4`    |
| 秀丽隐杆线虫 | `Caenorhabditis elegans`         | `GCF_000002985.6`    |
| 家蚕         | `Bombyx mori`                    | `GCF_030269925.1`    |
| 蜜蜂         | `Apis mellifera`                 | `GCF_003254395.2`    |
| 酿酒酵母     | `Saccharomyces cerevisiae S288C` | `GCF_000146045.2`    |

## 1) `DBD-001` 完整基因组

- **作用：** 按组装 accession 或常见物种名称下载完整参考基因组 FASTA，同时保存组装报告、来源 URL 和校验信息，用于建立可追溯的本地参考序列。
- **API：** `dnakit.references.resolve_genome_assembly(query[必须], timeout[可选], api_key[可选], api_base_url[可选])`、`dnakit.references.download_genome(query[必须], output_dir[必须], overwrite[可选], keep_package[可选], timeout[可选], api_key[可选], api_base_url[可选], chunk_size[可选], max_download_bytes[可选], progress[可选])`。
- **输入：** NCBI 版本化 accession、固定别名（如 `hg38`）或物种名；需要本地输出目录。
- **Agent 参数：** `query` 是组装标识或物种名，不是工具搜索词。请求下载人类 hg38 时使用 `query="hg38"`，不要传入 `"human hg38 reference genome FASTA"`。工具 Schema 包含此约束；下载接口的 `api_key` 指 NCBI 密钥，不是 DeepSeek 密钥。
- **示例代码：**

```python
from pathlib import Path

from dnakit.references import download_genome

result = download_genome("hg38", Path("references/hg38"))
print(result.fasta_path)
print(result.fasta_md5, result.fasta_sha256)
```

- **示例结果：**

```text
references/hg38/<genomic FASTA>
<NCBI MD5> <local SHA-256>
```

## 2) `DBD-002` 部分基因组数据包

- **作用：** 从指定组装中只下载需要的染色体、序列类别或文件类型，并生成清单，用于减少大型基因组数据包的传输和存储量。
- **API：** `dnakit.download.genome_package(accessions[必须], output_dir[必须], include[可选], chromosomes[可选], hydrated[可选], config[可选], progress[可选])`。
- **输入：** 一个或多个带版本号的 `GCA_`/`GCF_` accession；可选染色体名称和 `include` 文件类型。
- **示例代码：**

```python
from dnakit.download import genome_package

result = genome_package(
    "GCF_000001405.40",
    "downloads/human_chrM",
    chromosomes=("MT",),
    include=("GENOME_FASTA", "GENOME_GFF", "SEQUENCE_REPORT"),
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/human_chrM/ncbi_genome_package_manifest.json
```

## 3) `DBD-003` 区域序列

- **作用：** 按组装、染色体和起止坐标下载目标区域序列，返回 FASTA 及坐标来源记录，用于提取基因座或验证局部分析结果。
- **API：** `dnakit.download.sequence(species[必须], region[必须], output_path[必须], strand[可选], upstream[可选], downstream[可选], mask[可选], config[可选], query_progress[可选])`。
- **输入：** Ensembl 物种名和 0-based 半开区域；支持链方向、上下游扩展和多个区域。
- **示例代码：**

```python
from dnakit.download import sequence

result = sequence(
    "homo_sapiens",
    "7:140424943-140424950",
    "downloads/region.fa",
)
print(result.files[0].path)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/region.fa
downloads/region.fa.manifest.json
```

## 4) `DBD-004` 基因组注释

- **作用：** 下载与指定组装匹配的 GFF、GTF、GBFF 等基因组注释文件，并记录版本和校验值，用于基因、转录本及功能区域分析。
- **API：** `dnakit.download.annotation(accession[必须], output_dir[必须], formats[可选], chromosomes[可选], config[可选], progress[可选])`。
- **输入：** 带版本号的组装 accession；可选注释格式和染色体筛选。
- **示例代码：**

```python
from dnakit.download import annotation

result = annotation(
    "GCF_000001405.40",
    "downloads/hg38_annotation",
    formats=("GENOME_GFF", "GENOME_GTF", "SEQUENCE_REPORT"),
)
print(result.files)
print(result.manifest_path)
```

- **示例结果：**

```text
(DownloadedFile(...),)
downloads/hg38_annotation/ncbi_genome_package_manifest.json
```

## 5) `DBD-005` 基因序列

- **作用：** 按 GeneID 选择性下载基因组区段、RNA、CDS、UTR 和蛋白 FASTA，保留来源 accession，用于构建指定基因的本地序列集合。
- **API：** `dnakit.download.gene(gene_ids[必须], output_dir[必须], include[可选], accession_filter[可选], include_product_report[可选], config[可选], progress[可选])`。
- **输入：** 一个或多个数字型 NCBI GeneID；可选 FASTA 类型和 accession 过滤器。
- **示例代码：**

```python
from dnakit.download import gene

result = gene(
    672,
    "downloads/gene_672",
    include=("FASTA_GENE", "FASTA_RNA", "FASTA_CDS"),
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/gene_672/ncbi_gene_package_manifest.json
```

## 6) `DBD-006` 蛋白序列

- **作用：** 下载基因或病毒记录关联的公开蛋白 FASTA，保留蛋白 accession 和来源信息，用于翻译结果核对或蛋白参考库准备。
- **API：** `dnakit.download.gene(..., include=("FASTA_PROTEIN",))`；病毒蛋白使用 `dnakit.download.virus_package(..., include=("PROTEIN",))`。
- **输入：** 基因下载使用数字型 GeneID；病毒下载使用病毒 accession 或单个 taxon。
- **示例代码：**

```python
from dnakit.download import gene

result = gene(
    672,
    "downloads/gene_672_protein",
    include=("FASTA_PROTEIN",),
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/gene_672_protein/ncbi_gene_package_manifest.json
```

## 7) `DBD-007` 组装信息

- **作用：** 下载组装报告、序列报告和注释统计等元数据文件，用于核对染色体名称、组装层级、文件版本及注释完整性。
- **API：** `dnakit.download.genome_package(accessions[必须], output_dir[必须], include=("SEQUENCE_REPORT",), config[可选], progress[可选])`。
- **输入：** 带版本号的组装 accession；通常只选择 `SEQUENCE_REPORT` 或其他报告文件。
- **示例代码：**

```python
from dnakit.download import genome_package

result = genome_package(
    "GCF_000001405.40",
    "downloads/hg38_assembly_report",
    include=("SEQUENCE_REPORT",),
)
print(result.metadata)
print(result.manifest_path)
```

- **示例结果：**

```text
{'accessions': ('GCF_000001405.40',), 'include': ('SEQUENCE_REPORT',), ...}
downloads/hg38_assembly_report/ncbi_genome_package_manifest.json
```

## 8) `DBD-008` 分类信息

- **作用：** 下载指定 TaxID 对应的科学名称、别名和完整分类谱系，保存结构化摘要，用于统一物种标识及按分类层级整理数据。
- **API：** `dnakit.download.taxonomy(tax_ids[必须], output_dir[必须], include_names[可选], include_summary[可选], config[可选], progress[可选])`。
- **输入：** 一个或多个数字型 NCBI TaxID；可选择名称报告和摘要报告。
- **示例代码：**

```python
from dnakit.download import taxonomy

result = taxonomy(
    (9606, 10090),
    "downloads/taxonomy",
    include_names=True,
    include_summary=True,
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/taxonomy/ncbi_taxonomy_package_manifest.json
```

## 9) `DBD-009` 基因元数据 {#query-result-export}

- **作用：** 把公共数据库查询得到的基因、样本或运行元数据导出为 CSV、TSV 或 JSON，便于离线筛选、审计和与下载文件关联。
- **API：** `dnakit.download.metadata(query[必须], output_path[必须], format[可选], config[可选])`。
- **输入：** `dnakit.search` 返回的 `QueryResult`；格式支持 JSON、JSONL、CSV、TSV 和 XML。
- **示例代码：**

```python
from dnakit.download import metadata
from dnakit.search import sample

samples = sample("Homo sapiens[Organism]", retmax=20)
result = metadata(samples, "downloads/samples.csv", format="csv")
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/samples.csv
downloads/samples.csv.manifest.json
```

## 10) `DBD-010` 病毒数据包

- **作用：** 按病毒名称、TaxID 或 accession 下载基因组、CDS、蛋白和注释组成的数据包，并保存来源清单，用于建立版本明确的病毒参考集。
- **API：** `dnakit.download.virus_package(query[必须], output_dir[必须], by[可选], include[可选], aux_reports[可选], refseq_only[可选], annotated_only[可选], complete_only[可选], host[可选], geo_location[可选], released_since[可选], updated_since[可选], config[可选], progress[可选])`。
- **输入：** 病毒 accession 或单个 taxon；可选 RefSeq、已注释、完整基因组、宿主、地域和时间筛选。
- **示例代码：**

```python
from dnakit.download import virus_package

result = virus_package(
    "NC_045512.2",
    "downloads/sars_cov_2",
    by="accession",
    include=("GENOME", "CDS", "PROTEIN"),
    aux_reports=("ANNOTATION", "BIOSAMPLE_REPORT"),
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/sars_cov_2/ncbi_virus_NC_045512.2_manifest.json
```

## 11) `DBD-011` 原始测序数据

- **作用：** 按 ENA Run accession 下载单端或双端 FASTQ，依据 ENA 提供的 MD5 校验完整性，并记录样本与文件对应关系。
- **API：** `dnakit.download.reads(query[必须], output_dir[必须], file_kind="fastq", config[可选], progress[可选])`。
- **输入：** ENA Run accession、ENA 查询表达式或 `dnakit.search.reads()` 返回的 `QueryResult`；`file_kind` 为 `fastq`。
- **示例代码：**

```python
from dnakit.download import reads

result = reads(
    "SRR390728",
    "downloads/SRR390728_fastq",
    file_kind="fastq",
)
print(result.files)
print(result.manifest_path)
```

- **示例结果：**

```text
(DownloadedFile(...),)
downloads/SRR390728_fastq/ena_fastq_manifest.json
```

## 12) `DBD-012` 比对数据

- **作用：** 下载 ENA 公开的 BAM、CRAM 或 SAM 文件及可用校验信息，用于复用已有 reads-to-reference 比对结果而无需重新比对。
- **API：** `dnakit.download.reads(query[必须], output_dir[必须], file_kind="submitted", config[可选], progress[可选])`。
- **输入：** ENA Run accession、查询表达式或 ENA `QueryResult`；`file_kind` 为 `submitted`。
- **示例代码：**

```python
from dnakit.download import reads

result = reads(
    "SRR390728",
    "downloads/SRR390728_submitted",
    file_kind="submitted",
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/SRR390728_submitted/ena_submitted_manifest.json
```

## 13) `DBD-013` 测序元数据

- **作用：** 查询并保存 ENA 的 Study、Sample、Experiment、Run 和文件链接等结构化信息，用于筛选测序数据及追踪样本来源。
- **API：** `dnakit.search.reads(query[必须], fields[可选], result[可选], limit[可选], offset[可选], config[可选])` 加 `dnakit.download.metadata(query[必须], output_path[必须], format[可选], config[可选])`。
- **输入：** ENA 查询表达式和可选字段；查询结果可导出为 JSON、JSONL、CSV、TSV 或 XML。
- **示例代码：**

```python
from dnakit.download import metadata
from dnakit.search import reads as search_reads

runs = search_reads('run_accession="SRR390728"', limit=20)
result = metadata(runs, "downloads/ena_runs.json", format="json")
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/ena_runs.json
downloads/ena_runs.json.manifest.json
```

## 14) `DBD-014` 变异数据

- **作用：** 下载指定参考版本的 dbSNP VCF、索引和校验文件，用于本地查询已收录变异及对序列变异结果进行注释。
- **API：** `dnakit.download.variants(output_dir[必须], source="dbsnp", format="vcf", assembly[可选], include_index[可选], include_checksum_files[可选], config[可选], progress[可选])`。
- **输入：** `source="dbsnp"`、`format="vcf"` 和 `GRCh37`/`GRCh38` 组装版本。
- **示例代码：**

```python
from dnakit.download import variants

result = variants(
    "downloads/dbsnp_grch38",
    source="dbsnp",
    format="vcf",
    assembly="GRCh38",
)
print(result.files)
print(result.manifest_path)
```

- **示例结果：**

```text
(dbsnp_GRCh38.vcf.gz, dbsnp_GRCh38.vcf.gz.tbi, ...)
downloads/dbsnp_grch38/dbsnp_vcf_manifest.json
```

## 15) `DBD-015` ClinVar 数据

- **作用：** 下载 ClinVar 的 VCF、XML 或 TSV 发布文件及版本信息，用于本地关联变异与公开临床意义注释。
- **API：** `dnakit.download.variants(output_dir[必须], source="clinvar", format[可选], assembly[可选], include_index[可选], config[可选], progress[可选])`。
- **输入：** `source="clinvar"`；格式为 `vcf`、`xml` 或 `tsv`，VCF 可选择 `GRCh37`/`GRCh38` 和索引。
- **示例代码：**

```python
from dnakit.download import variants

result = variants(
    "downloads/clinvar",
    source="clinvar",
    format="vcf",
    assembly="GRCh38",
    include_index=True,
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/clinvar/clinvar_vcf_manifest.json
```

## 16) `DBD-016` 表达数据

- **作用：** 按 GEO accession 下载公开的表达矩阵或补充文件，并保存样本与来源信息，用于离线表达数据分析。
- **API：** `dnakit.download.expression(accession[必须], output_dir[必须], format[可选], platform[可选], config[可选], progress[可选])`。
- **输入：** `GSE`、`GDS`、`GPL` 或 `GSM` accession；Series Matrix 和 common RAW 通常需要 `GSE`，多平台矩阵可指定 `GPL`。
- **示例代码：**

```python
from dnakit.download import expression

result = expression(
    "GSE100",
    "downloads/GSE100",
    format="matrix",
    platform="GPL96",
)
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/GSE100/geo_matrix_manifest.json
```

## 17) `DBD-017` 调控数据

- **作用：** 按物种、组装、实验或轨道条件下载 ENCODE/UCSC 调控数据，保留文件 metadata，用于注释开放染色质、结合位点等区域。
- **API：** `dnakit.search.encode_search(...)` 加 `dnakit.download.encode_files(query[必须], output_dir[必须], config[可选], progress[可选])`；显式资源也可使用 `dnakit.download.tracks(resources[必须], output_dir[必须], source[可选], config[可选], progress[可选])`。
- **输入：** ENCODE `object_type="File"` 的查询结果，或调用方已确认的 HTTPS `RemoteFile` 列表。
- **示例代码：**

```python
from dnakit.download import encode_files
from dnakit.search import encode_search

found = encode_search(
    object_type="File",
    filters={"file_format": "bigWig", "status": "released"},
    limit=10,
)
result = encode_files(found, "downloads/encode_regulation")
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/encode_regulation/encode_files_manifest.json
```

## 18) `DBD-018` 重复序列数据

- **作用：** 下载与指定组装匹配的 UCSC RepeatMasker 等重复序列轨道，用于标记基因组中的重复类别和坐标区间。
- **API：** `dnakit.search.ucsc_files(genome[必须], pattern[可选], limit[可选], config[可选])` 加 `dnakit.download.ucsc_files(query[必须], output_dir[必须], config[可选], progress[可选])`。
- **输入：** UCSC assembly 和文件 glob；下载接口接收对应的 `QueryResult`。
- **示例代码：**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*repeat*", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_repeats")
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/hg38_repeats/ucsc_files_manifest.json
```

## 19) `DBD-019` 保守性数据

- **作用：** 下载指定组装的 phyloP、phastCons 等保守性轨道及说明信息，用于为基因组位置关联跨物种保守性分数。
- **API：** `dnakit.search.ucsc_files(...)` 加 `dnakit.download.ucsc_files(...)`；已明确 URL 的单个轨道也可用 `dnakit.download.tracks(...)`。
- **输入：** UCSC assembly 和 conservation 文件 glob，或调用方确认的显式 HTTPS 文件。
- **示例代码：**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*phastCons*", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_conservation")
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/hg38_conservation/ucsc_files_manifest.json
```

## 20) `DBD-020` 多物种比对

- **作用：** 下载指定组装和物种集合的 UCSC 多物种 MAF 比对文件，用于离线查看同源区域及开展比较基因组分析。
- **API：** `dnakit.search.ucsc_files(...)` 加 `dnakit.download.ucsc_files(...)`；查询结果元数据可用 `dnakit.download.metadata(...)` 导出。
- **输入：** UCSC assembly 和 MAF/多序列文件 glob，或显式 HTTPS 文件。
- **示例代码：**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*.maf.gz", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_maf")
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/hg38_maf/ucsc_files_manifest.json
```

## 21) `DBD-021` 坐标转换文件

- **作用：** 下载两个基因组组装版本之间的 UCSC chain 文件并记录方向，用于将区域、变异和注释坐标从一个组装映射到另一个组装。
- **API：** `dnakit.search.ucsc_files(...)` 加 `dnakit.download.ucsc_files(...)`。
- **输入：** UCSC assembly 和 `liftOver`/`chain` 文件 glob。
- **示例代码：**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*liftOver*chain.gz", limit=1000)
result = download_ucsc_files(catalog, "downloads/ucsc_chain")
print(result.manifest_path)
```

- **示例结果：**

```text
downloads/ucsc_chain/ucsc_files_manifest.json
```

## 22) `DBD-022` 通用链接下载

`download_url(url, output_dir)` 只需提供 HTTPS 文件直链和目标文件夹，目录不存在时会自动创建。
默认使用链接路径末尾的文件名，自动解码 URL 转义字符，不把 `?` 后的查询参数放入文件名。

```python
from dnakit.download import download_url

result = download_url(
    "https://example.org/files/data.fa.gz",  # 替换为需要下载的文件直链
    "downloads/my_data",
)
print(result.path)       # 下载文件的绝对路径
print(result.byte_size)  # 文件字节数
```

可选参数 `filename="reference.fa.gz"` 可以指定保存名称；链接没有文件名时必须提供该参数。
默认同名文件已存在就抛出 `FileExistsError`；需要覆盖时使用
`config=DownloadConfig(overwrite=True)`（从 `dnakit.download` 导入 `DownloadConfig`）。
超时、单文件大小上限和进度显示也通过 `DownloadConfig` 设置。
可传入 `expected_md5` / `expected_sha256` 验证下载内容，或用 `headers` 提供请求头。

此接口不限文件类型，按响应原样保存；例如 `.gz` 下载后仍是压缩文件。
它接收文件直链，不解析网页或网盘分享页面，不自动解压。
Agent/MCP 可通过 `dnakit.download.download_url` 调用，必填参数同样只有 `url` 和 `output_dir`。
需要批量下载并生成来源清单时，可以使用 `dataset([RemoteFile(url1), RemoteFile(url2)], output_dir)`。
