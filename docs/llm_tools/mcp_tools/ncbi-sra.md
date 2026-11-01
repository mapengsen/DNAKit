# NCBI SRA MCP Server (`ncbi-sra`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **4** 个 DNA/基因相关工具（服务器共列出 6 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [CSI-Genomics-and-Data-Analytics-Core/ncbi-sra-mcp-server @ `5569215`](https://github.com/CSI-Genomics-and-Data-Analytics-Core/ncbi-sra-mcp-server/tree/5569215722ba8554c80aa66f8d971b5bfe09a984) |
| 许可证 | MIT |
| 启动命令 | `python -c from sra_mcp_server import mcp; mcp.run()` |
| 传输方式 | stdio |
| 安装 | Not on PyPI: install separately from the pinned source checkout in a Python &gt;=3.11 environment (pip install -r requirements.txt, which provides fastmcp), set workdir to the checkout and executable_path to that Python; optionally set NCBI_EMAIL and NCBI_API_KEY, and install sra-tools (prefetch, fasterq-dump) for toolkit downloads. |
| 副作用 | Queries NCBI Entrez E-utilities and the ENA portal API; sra_download starts background jobs that create the destination directory (default ./downloads) and write .sra or FASTQ files via urllib or the sra-tools prefetch/fasterq-dump subprocesses. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`sra_download`](#sra_download) | `sequencing_data` | Launches a download/conversion job for an SRA run accession to SRA or FASTQ using SRA Toolkit or direct HTTP/S3. |
| [`sra_files`](#sra_files) | `sequencing_data` | Retrieves direct FASTQ/SRA file URLs, MD5 checksums and sizes for run accessions via the ENA mirror. |
| [`sra_metadata`](#sra_metadata) | `sequencing_data` | Retrieves SRA run/study metadata (platform, library strategy, spots, bases, organism) via Entrez esummary. |
| [`sra_search`](#sra_search) | `sequencing_data` | Searches SRA runs and studies with Entrez query syntax (organism, platform, strategy, dates). |

## `sra_download` {#sra_download}

Launches a download/conversion job for an SRA run accession to SRA or FASTQ using SRA Toolkit or direct HTTP/S3.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`id`
- **关键词：**fastq, sra toolkit, fasterq-dump, srr, paired-end

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `id` | string | 是 | SRA run accession (e.g., "SRR1234567") |
| `format` | string | 否 | Output format - "sra" or "fastq" (default: "fastq") 默认 `fastq` |
| `split` | string | 否 | Split mode - "none", "spots", or "files" (default: "files") 默认 `files` |
| `threads` | integer | 否 | Number of threads (default: 4) 默认 `4` |
| `dest` | string | 否 | Destination directory (default: "./downloads") 默认 `./downloads` |
| `prefer_direct_http` | boolean | 否 | Use direct HTTP download when possible (default: True) 默认 `True` |

??? note "上游说明（英文原文）"

    Launch download/conversion job using SRA Toolkit or direct HTTP/S3.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_download',
    {'id': '...'},
)
```

## `sra_files` {#sra_files}

Retrieves direct FASTQ/SRA file URLs, MD5 checksums and sizes for run accessions via the ENA mirror.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`ids`
- **关键词：**fastq, ena, file urls, md5, srr

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ids` | array of string | 是 | List of run accessions (e.g., \["SRR1234567", "SRR7890123"\]) Note: Must be run accessions (SRR/ERR/DRR), not study accessions |

??? note "上游说明（英文原文）"

    Retrieve direct file endpoints via ENA mirror.

    RECOMMENDED LIMITS:
    - For download planning: ≤ 50 IDs per call
    - For bulk file listing: ≤ 100 IDs per call

    FILE FORMATS PROVIDED:
    - FASTQ: Compressed FASTQ files (.fastq.gz)
    - SRA: Original SRA format files (.sra)
    - MD5: Checksums for verification
    - File sizes: In bytes for download planning

    DOWNLOAD OPTIONS:
    - Direct HTTP: Use provided URLs with wget/curl
    - ENA mirror: Fast, reliable European server
    - Multiple FASTQ files: For paired-end reads (fastq_1, fastq_2)

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_files',
    {'ids': []},
)
```

## `sra_metadata` {#sra_metadata}

Retrieves SRA run/study metadata (platform, library strategy, spots, bases, organism) via Entrez esummary.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`ids`
- **关键词：**esummary, library strategy, platform, spots, bases

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ids` | array of string | 是 | List of SRA IDs (UIDs or accessions) Examples: \["123456789", "987654321"\] or \["SRR1234567", "SRR7890123"\] |

??? note "上游说明（英文原文）"

    Retrieve metadata for SRA runs/studies via Entrez esummary.

    RECOMMENDED LIMITS:
    - For detailed analysis: ≤ 50 IDs per call
    - For bulk metadata: ≤ 200 IDs per call
    - Use pagination if you have more IDs

    METADATA FIELDS RETURNED:
    - accession: SRA run/study accession (e.g., SRR1234567)
    - title: Experiment title
    - organism: Scientific name of organism
    - platform: Sequencing platform (e.g., ILLUMINA, PACBIO)
    - instrument: Specific instrument model
    - library_strategy: Library strategy (RNA-Seq, WGS, etc.)
    - library_source: Library source (GENOMIC, TRANSCRIPTOMIC, etc.)
    - spots: Total number of spots/reads
    - bases: Total number of bases
    - study_accession: Associated study accession

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_metadata',
    {'ids': []},
)
```

## `sra_search` {#sra_search}

Searches SRA runs and studies with Entrez query syntax (organism, platform, strategy, dates).

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`term`
- **关键词：**esearch, entrez, rna-seq, wgs, illumina

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `term` | string | 是 | Search term using Entrez query syntax Examples: - "illumina\[platform\] AND human\[orgn\]" - "SARS-CoV-2\[orgn\] AND RNA-Seq\[strategy\]" - "breast cancer\[all\] AND illumina\[platform\]" |
| `retmax` | integer | 否 | Maximum results to return (default: 200, recommended ≤ 200) 默认 `200` |
| `retstart` | integer | 否 | Starting position for pagination (default: 0) 默认 `0` |
| `sort` | string | 否 | Sort order - "relevance" (default) or "date" 默认 `relevance` |
| `datetype` | string \| null | 否 | Date type filter - "pdat" (publication), "mdat" (modification), "edat" (entry) 默认 `None` |
| `mindate` | string \| null | 否 | Minimum date in YYYY/MM/DD format (e.g., "2020/01/01") 默认 `None` |
| `maxdate` | string \| null | 否 | Maximum date in YYYY/MM/DD format (e.g., "2023/12/31") 默认 `None` |

??? note "上游说明（英文原文）"

    Search SRA runs/studies using Entrez esearch.

    RECOMMENDED LIMITS:
    - For exploration: retmax ≤ 50
    - For bulk retrieval: Use pagination with retmax ≤ 200
    - Maximum supported: 100,000 records per query

    SEARCH SYNTAX:
    - Organism: "Homo sapiens\[orgn\]" or "human\[orgn\]"
    - Platform: "illumina\[platform\]" or "pacbio\[platform\]"
    - Strategy: "RNA-Seq\[strategy\]" or "WGS\[strategy\]"
    - Combine: Use AND, OR, NOT operators
    - Date ranges: Use datetype, mindate, maxdate parameters

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_search',
    {'term': '...'},
)
```
