# 下载

下载序列、参考基因组和注释；格式、重试与物种名称见 [FAQ](../../faq/tutorial_integrations.md#downloads)。

<span id="accession-sequences"></span>按 accession 下载：[NCBI](#dnakit.download.ncbi.sequence) · [ENA](#dnakit.download.ena.sequence)；格式与坐标约定见 [FAQ](../../faq/tutorial_integrations.md#downloads)。

## 1) 完整基因组 {#1-dbd-001}

- **作用：** 按组装 accession 或常见物种名称下载完整参考基因组 FASTA，同时保存组装报告、来源 URL 和校验信息，用于建立可追溯的本地参考序列。
- **处理方法：** 解析组装编号，从 NCBI 获取基因组包并校验 FASTA 和来源报告。
- **API：** `dnakit.references.resolve_genome_assembly(query[必须], timeout[可选], api_key[可选], api_base_url[可选])`、`dnakit.references.download_genome(query[必须], output_dir[必须], overwrite[可选], keep_package[可选], timeout[可选], api_key[可选], api_base_url[可选], chunk_size[可选], max_download_bytes[可选], progress[可选])`。
- **输入：** NCBI 版本化 accession、固定别名（如 `hg38`）或物种名；需要本地输出目录。
- **输出：** `resolve_genome_assembly` 返回组装信息；`download_genome` 返回 `GenomeAssembly`，包含 FASTA 路径和校验值。
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

## 2) 部分基因组数据包 {#2-dbd-002}

- **作用：** 从指定组装中只下载需要的染色体、序列类别或文件类型，并生成清单，用于减少大型基因组数据包的传输和存储量。
- **处理方法：** 按染色体和文件类别筛选 NCBI 组装包内容。
- **API：** `dnakit.download.genome_package(accessions[必须], output_dir[必须], include[可选], chromosomes[可选], hydrated[可选], config[可选], progress[可选])`。
- **输入：** 一个或多个带版本号的 `GCA_`/`GCF_` accession；可选染色体名称和 `include` 文件类型。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 3) 区域序列 {#3-dbd-003}

- **作用：** 按组装、染色体和起止坐标下载目标区域序列，返回 FASTA 及坐标来源记录，用于提取基因座或验证局部分析结果。
- **处理方法：** 将区间和方向参数提交 Ensembl，保存返回序列。
- **API：** `dnakit.download.sequence(species[必须], region[必须], output_path[必须], strand[可选], upstream[可选], downstream[可选], mask[可选], config[可选], query_progress[可选])`。
- **输入：** Ensembl 物种名和 0-based 半开区域；支持链方向、上下游扩展和多个区域。
- **输出：** 返回 `DatasetDownloadResult`，包含区域 FASTA 文件和来源清单路径。
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

## 4) 基因组注释 {#4-dbd-004}

- **作用：** 下载与指定组装匹配的 GFF、GTF、GBFF 等基因组注释文件，并记录版本和校验值，用于基因、转录本及功能区域分析。
- **处理方法：** 从 NCBI 组装包提取所选注释格式并保存来源清单。
- **API：** `dnakit.download.annotation(accession[必须], output_dir[必须], formats[可选], chromosomes[可选], config[可选], progress[可选])`。
- **输入：** 带版本号的组装 accession；可选注释格式和染色体筛选。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 5) 基因序列 {#5-dbd-005}

- **作用：** 按 GeneID 选择性下载基因组区段、RNA、CDS、UTR 和蛋白 FASTA，保留来源 accession，用于构建指定基因的本地序列集合。
- **处理方法：** 按 GeneID 和 include 选项请求 NCBI 基因数据包。
- **API：** `dnakit.download.gene(gene_ids[必须], output_dir[必须], include[可选], accession_filter[可选], include_product_report[可选], config[可选], progress[可选])`。
- **输入：** 一个或多个数字型 NCBI GeneID；可选 FASTA 类型和 accession 过滤器。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 6) 蛋白序列 {#6-dbd-006}

- **作用：** 下载基因或病毒记录关联的公开蛋白 FASTA，保留蛋白 accession 和来源信息，用于翻译结果核对或蛋白参考库准备。
- **处理方法：** 按蛋白 accession 获取序列并保存文件和来源信息。
- **API：** `dnakit.download.gene(..., include=("FASTA_PROTEIN",))`；病毒蛋白使用 `dnakit.download.virus_package(..., include=("PROTEIN",))`。
- **输入：** 基因下载使用数字型 GeneID；病毒下载使用病毒 accession 或单个 taxon。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 7) 组装信息 {#7-dbd-007}

- **作用：** 下载组装报告、序列报告和注释统计等元数据文件，用于核对染色体名称、组装层级、文件版本及注释完整性。
- **处理方法：** 获取组装的结构化元数据，再按所选格式导出。
- **API：** `dnakit.download.genome_package(accessions[必须], output_dir[必须], include=("SEQUENCE_REPORT",), config[可选], progress[可选])`。
- **输入：** 带版本号的组装 accession；通常只选择 `SEQUENCE_REPORT` 或其他报告文件。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 8) 分类信息 {#8-dbd-008}

- **作用：** 下载指定 TaxID 对应的科学名称、别名和完整分类谱系，保存结构化摘要，用于统一物种标识及按分类层级整理数据。
- **处理方法：** 查询分类记录及谱系，保存结构化结果。
- **API：** `dnakit.download.taxonomy(tax_ids[必须], output_dir[必须], include_names[可选], include_summary[可选], config[可选], progress[可选])`。
- **输入：** 一个或多个数字型 NCBI TaxID；可选择名称报告和摘要报告。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 9) 基因元数据 {#9-dbd-009}

- **作用：** 把公共数据库查询得到的基因、样本或运行元数据导出为 CSV、TSV 或 JSON，便于离线筛选、审计和与下载文件关联。
- **处理方法：** 查询基因元数据，将结果表导出到本地文件。
- **API：** `dnakit.download.metadata(query[必须], output_path[必须], format[可选], config[可选])`。
- **输入：** `dnakit.search` 返回的 `QueryResult`；格式支持 JSON、JSONL、CSV、TSV 和 XML。
- **输出：** 返回 `DatasetDownloadResult`，包含导出的元数据文件和来源清单路径。
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

## 10) 病毒数据包 {#10-dbd-010}

- **作用：** 按病毒名称、TaxID 或 accession 下载基因组、CDS、蛋白和注释组成的数据包，并保存来源清单，用于建立版本明确的病毒参考集。
- **处理方法：** 从 NCBI 病毒数据包获取指定基因组及可选注释。
- **API：** `dnakit.download.virus_package(query[必须], output_dir[必须], by[可选], include[可选], aux_reports[可选], refseq_only[可选], annotated_only[可选], complete_only[可选], host[可选], geo_location[可选], released_since[可选], updated_since[可选], config[可选], progress[可选])`。
- **输入：** 病毒 accession 或单个 taxon；可选 RefSeq、已注释、完整基因组、宿主、地域和时间筛选。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 11) 原始测序数据 {#11-dbd-011}

- **作用：** 按 ENA Run accession 下载单端或双端 FASTQ，依据 ENA 提供的 MD5 校验完整性，并记录样本与文件对应关系。
- **处理方法：** 解析 ENA Run 的 FASTQ 下载链接，下载后按提供方 MD5 校验。
- **API：** `dnakit.download.reads(query[必须], output_dir[必须], file_kind="fastq", config[可选], progress[可选])`。
- **输入：** ENA Run accession、ENA 查询表达式或 `dnakit.search.reads()` 返回的 `QueryResult`；`file_kind` 为 `fastq`。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 12) 比对数据 {#12-dbd-012}

- **作用：** 下载 ENA 公开的 BAM、CRAM 或 SAM 文件及可用校验信息，用于复用已有 reads-to-reference 比对结果而无需重新比对。
- **处理方法：** 读取 ENA submitted 文件列表，下载已公开的比对文件。
- **API：** `dnakit.download.reads(query[必须], output_dir[必须], file_kind="submitted", config[可选], progress[可选])`。
- **输入：** ENA Run accession、查询表达式或 ENA `QueryResult`；`file_kind` 为 `submitted`。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 13) 测序元数据 {#13-dbd-013}

- **作用：** 查询并保存 ENA 的 Study、Sample、Experiment、Run 和文件链接等结构化信息，用于筛选测序数据及追踪样本来源。
- **处理方法：** 查询 ENA 元数据，再序列化为所选表格或结构化格式。
- **API：** `dnakit.search.reads(query[必须], fields[可选], result[可选], limit[可选], offset[可选], config[可选])` 加 `dnakit.download.metadata(query[必须], output_path[必须], format[可选], config[可选])`。
- **输入：** ENA 查询表达式和可选字段；查询结果可导出为 JSON、JSONL、CSV、TSV 或 XML。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
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

## 14) 变异数据 {#14-dbd-014}

- **作用：** 下载指定参考版本的 dbSNP VCF、索引和校验文件，用于本地查询已收录变异及对序列变异结果进行注释。
- **处理方法：** 按参考版本选择 dbSNP 发布文件并下载 VCF、索引和校验文件。
- **API：** `dnakit.download.variants(output_dir[必须], source="dbsnp", format="vcf", assembly[可选], include_index[可选], include_checksum_files[可选], config[可选], progress[可选])`。
- **输入：** `source="dbsnp"`、`format="vcf"` 和 `GRCh37`/`GRCh38` 组装版本。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 15) ClinVar 数据 {#15-dbd-015-clinvar}

- **作用：** 下载 ClinVar 的 VCF、XML 或 TSV 发布文件及版本信息，用于本地关联变异与公开临床意义注释。
- **处理方法：** 选择 ClinVar 的发布格式与参考版本，再下载对应文件。
- **API：** `dnakit.download.variants(output_dir[必须], source="clinvar", format[可选], assembly[可选], include_index[可选], config[可选], progress[可选])`。
- **输入：** `source="clinvar"`；格式为 `vcf`、`xml` 或 `tsv`，VCF 可选择 `GRCh37`/`GRCh38` 和索引。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 16) 表达数据 {#16-dbd-016}

- **作用：** 按 GEO accession 下载公开的表达矩阵或补充文件，并保存样本与来源信息，用于离线表达数据分析。
- **处理方法：** 按 GEO accession 定位并下载表达矩阵或补充文件。
- **API：** `dnakit.download.expression(accession[必须], output_dir[必须], format[可选], platform[可选], config[可选], progress[可选])`。
- **输入：** `GSE`、`GDS`、`GPL` 或 `GSM` accession；Series Matrix 和 common RAW 通常需要 `GSE`，多平台矩阵可指定 `GPL`。
- **输出：** 返回 `DatasetDownloadResult`，包含下载文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
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

## 17) 调控数据 {#17-dbd-017}

- **作用：** 按物种、组装、实验或轨道条件下载 ENCODE/UCSC 调控数据，保留文件 metadata，用于注释开放染色质、结合位点等区域。
- **处理方法：** 筛选 ENCODE 文件或显式轨道资源，按查询结果下载。
- **API：** `dnakit.search.encode_search(...)` 加 `dnakit.download.encode_files(query[必须], output_dir[必须], config[可选], progress[可选])`；显式资源也可使用 `dnakit.download.tracks(resources[必须], output_dir[必须], source[可选], config[可选], progress[可选])`。
- **输入：** ENCODE `object_type="File"` 的查询结果，或调用方已确认的 HTTPS `RemoteFile` 列表。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
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

## 18) 重复序列数据 {#18-dbd-018}

- **作用：** 下载与指定组装匹配的 UCSC RepeatMasker 等重复序列轨道，用于标记基因组中的重复类别和坐标区间。
- **处理方法：** 按文件模式筛选 UCSC RepeatMasker 资源，再下载匹配文件。
- **API：** `dnakit.search.ucsc_files(genome[必须], pattern[可选], limit[可选], config[可选])` 加 `dnakit.download.ucsc_files(query[必须], output_dir[必须], config[可选], progress[可选])`。
- **输入：** UCSC assembly 和文件 glob；下载接口接收对应的 `QueryResult`。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
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

## 19) 保守性数据 {#19-dbd-019}

- **作用：** 下载指定组装的 phyloP、phastCons 等保守性轨道及说明信息，用于为基因组位置关联跨物种保守性分数。
- **处理方法：** 筛选 UCSC phyloP 或 phastCons 文件，下载已有保守性轨道。
- **API：** `dnakit.search.ucsc_files(...)` 加 `dnakit.download.ucsc_files(...)`；已明确 URL 的单个轨道也可用 `dnakit.download.tracks(...)`。
- **输入：** UCSC assembly 和 conservation 文件 glob，或调用方确认的显式 HTTPS 文件。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
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

## 20) 多物种比对 {#20-dbd-020}

- **作用：** 下载指定组装和物种集合的 UCSC 多物种 MAF 比对文件，用于离线查看同源区域及开展比较基因组分析。
- **处理方法：** 查询 UCSC 多物种比对资源并下载所选文件。
- **API：** `dnakit.search.ucsc_files(...)` 加 `dnakit.download.ucsc_files(...)`；查询结果元数据可用 `dnakit.download.metadata(...)` 导出。
- **输入：** UCSC assembly 和 MAF/多序列文件 glob，或显式 HTTPS 文件。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
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

## 21) 坐标转换文件 {#21-dbd-021}

- **作用：** 下载两个基因组组装版本之间的 UCSC chain 文件并记录方向，用于将区域、变异和注释坐标从一个组装映射到另一个组装。
- **处理方法：** 根据来源与目标组装选择 UCSC chain 文件并下载。
- **API：** `dnakit.search.ucsc_files(...)` 加 `dnakit.download.ucsc_files(...)`。
- **输入：** UCSC assembly 和 `liftOver`/`chain` 文件 glob。
- **输出：** 返回 `QueryResult`，包含查询记录 `records`、数据源 `provider` 和分页信息。
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

## 22) 通用链接下载 {#22-dbd-022}

- **作用：** 从文件直链下载任意类型的文件。
- **处理方法：** 原样保存响应内容，可校验 MD5 或 SHA-256。
- **API：** `dnakit.download.download_url(url, output_dir, *, filename=None, config=None, ...)`
- **输入：** HTTPS 文件直链和输出目录；替换示例链接。默认保护已有文件，压缩包不会自动解压。
- **输出：** 返回 `DownloadedFile`，包含本地路径、大小和校验值。
- **示例代码：**

```python
from dnakit.download import download_url

result = download_url(
    "https://example.org/files/data.fa.gz",  # 替换为需要下载的文件直链
    "downloads/my_data",
)
print(result.path)       # 下载文件的绝对路径
print(result.byte_size)  # 文件字节数
```

- **示例结果：** 下载文件的绝对路径和字节数；具体结果由实际链接决定。

## 23) 从检索结果下载序列 {#从检索结果下载序列}

- **作用：** 按物种检索序列，再下载第一条匹配记录。
- **处理方法：** 用 NCBI Entrez 检索，将记录 ID 转为 accession，再调用序列下载接口。
- **API：** `dnakit.search.ncbi.entrez(...)`；`dnakit.search.ncbi.accessions(...)`；`dnakit.download.ncbi.sequence(...)`
- **输入：** 物种及筛选条件、输出路径；需要联网，会写入本地文件。
- **输出：** 检索返回 `QueryResult`；下载返回 `DatasetDownloadResult`，包含序列文件、序列长度和来源清单。
- **示例代码：**

```python
import dnakit

found = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(found.metadata["uids"])
if accessions:
    downloaded = dnakit.download.ncbi.sequence(
        accessions[0], "ecoli_k12.fasta", format="fasta",
    )
    print(downloaded.files[0].path)
    print(downloaded.metadata["sequence_length"])
    print(downloaded.manifest_path)
else:
    print("没有匹配序列，请核对查询条件。")
```

- **示例结果：** 下载路径、序列长度和来源清单；没有命中时输出提示。具体命中随数据库变化。

<span id="common-genome-names"></span>常见物种名称、ENA 示例及更多参数见 [FAQ](../../faq/tutorial_integrations.md#downloads)。

<span id="encode-download-example"></span>[ENCODE 批量下载与预览](#dnakit.encode.batch_download)。

<!-- BEGIN GENERATED API: downloads -->

<span id="downloads"></span>

## 24) 通过 HTTP 或 HTTPS 下载文件 {#dnakit.databases.download_file}

- **作用：** 通过 HTTP 或 HTTPS 下载文件。
- **处理方法：** 通过 HTTP(S) 分块写入指定路径，或在 return_content=True 时返回文本。
- **API：** `dnakit.databases.download_file(*, url, output_path=None, timeout=None, return_content=None, chunk_size=None, follow_redirects=None, verbose=True)`
- **输入：** 必须提供 `url`。 运行前准备示例涉及的本地文件与外部程序，并将路径替换为自己的输入或输出位置。
- **输出：** 保存模式返回文件路径；`return_content=True` 时返回文件的文本内容。
- **示例代码：**

```python
import dnakit

result = dnakit.databases.download_file(
    url='https://www.wikipedia.org',
    output_path='./downloads/output.txt',
    timeout=30,
    return_content=False,
    chunk_size=8192,
    follow_redirects=True,
    verbose=False,
)
print(result)
```

- **示例结果：** 保存模式返回文件路径；`return_content=True` 时返回文件的文本内容。

[完整参数与使用说明](../../faq/api/15_download.md#dnakit.databases.download_file-parameters)

## 25) 下载 ENA 核酸序列并保存校验清单 {#dnakit.download.ena.sequence}

- **作用：** 按 INSDC accession 下载 ENA 序列文件。
- **处理方法：** 解析 accession，向对应序列服务请求数据，写出文件与来源清单。
- **API：** `dnakit.download.ena.sequence(accession, output_path, *, format='fasta', config=None, search_config=None, api_base_url='https://www.ebi.ac.uk/ena/browser/api', verbose=True)`
- **输入：** 必须提供 `accession`, `output_path`。
- **输出：** 返回 `DatasetDownloadResult`，包含文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
- **示例代码：**

```python
import dnakit

result = dnakit.download.ena.sequence(
    accession='U00096.3',
    output_path='ecoli_ena.fasta',
    format='fasta',
    verbose=False,
)
print(result)
```

- **示例结果：** DatasetDownloadResult：files 为文件记录，metadata 含来源、解析出的 accession 和序列长度，manifest_path 指向来源清单。详细格式支持与坐标约定见[序列下载教程](15_download.md#accession-sequences)。

[完整参数与使用说明](../../faq/api/15_download.md#dnakit.download.ena.sequence-parameters)

## 26) 下载 NCBI 核酸或蛋白质序列 {#dnakit.download.ncbi.sequence}

- **作用：** 按 accession 下载 NCBI 序列，可选择格式、核酸子区间和链方向。
- **处理方法：** 解析 accession，向对应序列服务请求数据，写出文件与来源清单。
- **API：** `dnakit.download.ncbi.sequence(accession, output_path, *, database='nuccore', format='fasta', seq_start=None, seq_stop=None, strand=1, config=None, search_config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)`
- **输入：** 必须提供 `accession`, `output_path`。
- **输出：** 返回 `DatasetDownloadResult`，包含文件 `files`、来源信息 `metadata` 和清单路径 `manifest_path`。
- **示例代码：**

```python
import dnakit

result = dnakit.download.ncbi.sequence(
    accession='NC_000913.3',
    output_path='ecoli_ncbi.fasta',
    format='fasta',
    verbose=False,
)
print(result)
```

- **示例结果：** DatasetDownloadResult：files 为文件记录，metadata 含来源、解析出的 accession 和序列长度，manifest_path 指向来源清单。详细格式支持与坐标约定见[序列下载教程](15_download.md#accession-sequences)。

[完整参数与使用说明](../../faq/api/15_download.md#dnakit.download.ncbi.sequence-parameters)

<!-- END GENERATED API: downloads -->

<!-- BEGIN GENERATED API: encode-downloads -->

<span id="encode-downloads"></span>

## 27) 批量检索并下载 ENCODE 文件 {#dnakit.encode.batch_download}

- **作用：** 批量下载匹配实验的文件。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.batch_download(download_dir, file_format=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, target=None, preferred_default=None, organize_by='experiment', verify_md5=True, limit=100, dry_run=True, *, verbose=True)`
- **输入：** 必须提供 `download_dir`。
- **输出：** 返回 JSON 文本，包含下载预览（dry_run=True）或逐文件下载结果。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.batch_download(
    download_dir='downloads/encode',
    assay_title='ATAC-seq',
    assembly='GRCh38',
    file_format='bed',
    limit=3,
    dry_run=True,
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含下载预览（dry_run=True）或逐文件下载结果。

[完整参数与使用说明](../../faq/api/15_download.md#dnakit.encode.batch_download-parameters)

## 28) 按文件编号下载 ENCODE 数据 {#dnakit.encode.download_files}

- **作用：** 下载指定文件并记录元数据。
- **处理方法：** 按指定条件访问 ENCODE，解析实验、文件或下载结果。
- **API：** `dnakit.encode.download_files(file_accessions, download_dir, organize_by='flat', verify_md5=True, *, verbose=True)`
- **输入：** 必须提供 `file_accessions`, `download_dir`。 先用文件查询功能取得实际编号，替换 `ENCFF_REPLACE_ME`。
- **输出：** 返回 JSON 文本，包含逐文件路径、大小、成功状态与 MD5 校验结果。
- **示例代码：**

```python
import dnakit

result = dnakit.encode.download_files(
    file_accessions=['ENCFF_REPLACE_ME'],
    download_dir='downloads/encode',
    verbose=False,
)
print(result)
```

- **示例结果：** 返回 JSON 文本，包含逐文件路径、大小、成功状态与 MD5 校验结果。

[完整参数与使用说明](../../faq/api/15_download.md#dnakit.encode.download_files-parameters)

<!-- END GENERATED API: encode-downloads -->
