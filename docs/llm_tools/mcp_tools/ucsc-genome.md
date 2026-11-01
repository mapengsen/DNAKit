# UCSC Genome MCP (`ucsc-genome`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **11** 个 DNA/基因相关工具（服务器共列出 12 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [hlydecker/ucsc-genome-mcp @ `7f0f83a`](https://github.com/hlydecker/ucsc-genome-mcp/tree/7f0f83a2f1935916f8cd35ddfb8dd3c36dac2e84) |
| 许可证 | MIT |
| 启动命令 | `python ucsc-genome-mcp.py` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned source checkout and its requirements in Python &gt;=3.10, set workdir to the checkout and executable_path to that Python. Use the existing hyphenated ucsc-genome-mcp.py script; upstream README/packaging also mention nonexistent underscore filenames. Pin 'mcp&lt;2': under mcp 2.x the script fails to import. |
| 副作用 | Sends genome, track, region and sequence queries to UCSC and caller-selected public assembly hubs. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`find_genome`](#find_genome) | `genome_annotation` | Searches UCSC Genome Browser genome assemblies by name or accession with filters for assembly level, category and year. |
| [`get_sequence`](#get_sequence) | `sequence` | Retrieves DNA sequence for a chromosome or coordinate range from a UCSC or hub genome assembly, optionally reverse complemented. |
| [`get_track_data`](#get_track_data) | `genome_annotation` | Retrieves records from a UCSC or hub data track, optionally restricted to a chromosome region. |
| [`list_chromosomes`](#list_chromosomes) | `genome_annotation` | Lists chromosomes (and sizes) of a UCSC database genome, assembly hub or track. |
| [`list_files`](#list_files) | `genome_annotation` | Lists download files available for a UCSC genome assembly such as hg38 or mm10. |
| [`list_genark_genomes`](#list_genark_genomes) | `genome_annotation` | Lists UCSC GenArk assembly hub genomes or tests whether a specific genome exists. |
| [`list_hub_genomes`](#list_hub_genomes) | `genome_annotation` | Lists all genomes available in a given track hub or assembly hub URL. |
| [`list_public_hubs`](#list_public_hubs) | `genome_annotation` | Lists all public track hubs registered with the UCSC Genome Browser. |
| [`list_tracks`](#list_tracks) | `genome_annotation` | Lists all data tracks available for a UCSC database genome or hub genome. |
| [`list_ucsc_genomes`](#list_ucsc_genomes) | `genome_annotation` | Lists all UCSC Genome Browser database genome assemblies. |
| [`search_genome`](#search_genome) | `genome_annotation` | Searches within a UCSC genome assembly across track data, help docs and public hubs for a term. |

## `find_genome` {#find_genome}

Searches UCSC Genome Browser genome assemblies by name or accession with filters for assembly level, category and year.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**genome assembly, gcf accession, species, assembly level, grch38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search string to find genomes (e.g., 'dog', 'GRCh38', 'GCF_028858775.2') |
| `browser` | string | 否 | Filter by browser availability (default: mustExist) 可选值：`mustExist`, `mayExist`, `notExist` |
| `stats_only` | boolean | 否 | Only show statistics about search results |
| `year` | integer | 否 | Filter results by year |
| `category` | string | 否 | Filter by NCBI category 可选值：`reference`, `representative` |
| `status` | string | 否 | Filter by NCBI status 可选值：`reference`, `representative` |
| `level` | string | 否 | Filter by NCBI assembly level 可选值：`complete`, `chromosome`, `scaffold`, `contig` |
| `max_items` | integer | 否 | Maximum number of items to return (default: 1000000, use -1 for max) |

??? note "上游说明（英文原文）"

    Search for a genome in the UCSC browser using a search string. Supports advanced search with +word (force inclusion), -word (exclusion), and word* (wildcard).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'find_genome',
    {'query': '...'},
)
```

## `get_sequence` {#get_sequence}

Retrieves DNA sequence for a chromosome or coordinate range from a UCSC or hub genome assembly, optionally reverse complemented.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`genome`, `chrom`
- **关键词：**dna sequence, chromosome, coordinates, reverse complement, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome` | string | 是 | Genome assembly name (e.g., 'hg38') |
| `chrom` | string | 是 | Chromosome name (e.g., 'chr1', 'chrM') |
| `start` | integer | 否 | Start coordinate (0-based, optional, requires end) |
| `end` | integer | 否 | End coordinate (1-based, optional, requires start) |
| `hub_url` | string | 否 | URL of assembly hub (optional) |
| `reverse_complement` | boolean | 否 | Return reverse complement of sequence |

??? note "上游说明（英文原文）"

    Retrieve DNA sequence from a specified genome assembly. Can retrieve entire chromosome or specific coordinates.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'get_sequence',
    {'genome': '...', 'chrom': '...'},
)
```

## `get_track_data` {#get_track_data}

Retrieves records from a UCSC or hub data track, optionally restricted to a chromosome region.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`genome`, `track`
- **关键词：**annotation track, genomic region, bed, coordinates, ucsc table

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome` | string | 是 | Genome assembly name |
| `track` | string | 是 | Track name |
| `chrom` | string | 否 | Chromosome name (optional) |
| `start` | integer | 否 | Start coordinate (0-based, optional, requires end) |
| `end` | integer | 否 | End coordinate (1-based, optional, requires start) |
| `hub_url` | string | 否 | URL of track/assembly hub (optional) |
| `max_items` | integer | 否 | Maximum number of items to return (default: 1000000) |
| `json_output_arrays` | boolean | 否 | Return data as JSON arrays instead of objects |

??? note "上游说明（英文原文）"

    Retrieve data from a specified track in a hub or UCSC database genome. Can be filtered by chromosome and coordinates.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'get_track_data',
    {'genome': '...', 'track': '...'},
)
```

## `list_chromosomes` {#list_chromosomes}

Lists chromosomes (and sizes) of a UCSC database genome, assembly hub or track.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`genome`
- **关键词：**chromosome sizes, assembly, contigs, hub, chrom

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome` | string | 是 | Genome assembly name |
| `hub_url` | string | 否 | URL of track/assembly hub (optional) |
| `track` | string | 否 | Specific track name to list chromosomes from (optional) |

??? note "上游说明（英文原文）"

    List chromosomes in an assembly hub, track hub, or UCSC database genome. Optionally filter by specific track.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_chromosomes',
    {'genome': '...'},
)
```

## `list_files` {#list_files}

Lists download files available for a UCSC genome assembly such as hg38 or mm10.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`genome`
- **关键词：**downloads, assembly files, 2bit, annotation files, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome` | string | 是 | Genome assembly name (e.g., 'hg38', 'mm10') |
| `format` | string | 否 | Output format (default: json) 可选值：`json`, `text` |
| `max_items` | integer | 否 | Maximum number of items to return |

??? note "上游说明（英文原文）"

    List download files available for a specified UCSC genome assembly.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_files',
    {'genome': '...'},
)
```

## `list_genark_genomes` {#list_genark_genomes}

Lists UCSC GenArk assembly hub genomes or tests whether a specific genome exists.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**genark, assembly hub, genome assembly, gca, species

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome` | string | 否 | Specific genome to test for existence (optional) |
| `max_items` | integer | 否 | Maximum number of items to return (default: 1000000) |

??? note "上游说明（英文原文）"

    List UCSC Genome Browser database genomes from assembly hub host (GenArk). Can also test for existence of a specific genome.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_genark_genomes',
    {},
)
```

## `list_hub_genomes` {#list_hub_genomes}

Lists all genomes available in a given track hub or assembly hub URL.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`hub_url`
- **关键词：**track hub, assembly hub, genome assembly, hub url, genomes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `hub_url` | string | 是 | URL of the track hub or assembly hub |

??? note "上游说明（英文原文）"

    List all genomes available in a specified track or assembly hub.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_hub_genomes',
    {'hub_url': '...'},
)
```

## `list_public_hubs` {#list_public_hubs}

Lists all public track hubs registered with the UCSC Genome Browser.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**track hubs, public hubs, annotation collections, genome browser, datasets

??? note "上游说明（英文原文）"

    List all available public track hubs in the UCSC Genome Browser.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_public_hubs',
    {},
)
```

## `list_tracks` {#list_tracks}

Lists all data tracks available for a UCSC database genome or hub genome.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`genome`
- **关键词：**annotation tracks, track list, genome browser, hub, composite

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genome` | string | 是 | Genome assembly name |
| `hub_url` | string | 否 | URL of track/assembly hub (optional, required with genome for hub tracks) |
| `track_leaves_only` | boolean | 否 | Only show tracks without composite container information |

??? note "上游说明（英文原文）"

    List all data tracks available in a specified hub or UCSC database genome.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_tracks',
    {'genome': '...'},
)
```

## `list_ucsc_genomes` {#list_ucsc_genomes}

Lists all UCSC Genome Browser database genome assemblies.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**genome assembly, hg38, mm10, database genomes, species

??? note "上游说明（英文原文）"

    List all UCSC Genome Browser database genomes available on the database host.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_ucsc_genomes',
    {},
)
```

## `search_genome` {#search_genome}

Searches within a UCSC genome assembly across track data, help docs and public hubs for a term.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`search`, `genome`
- **关键词：**gene search, position, track search, hg38, genome browser

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `search` | string | 是 | Search term |
| `genome` | string | 是 | Genome assembly to search in |
| `categories` | string | 否 | Restrict search to specific categories 可选值：`helpDocs`, `publicHubs`, `trackDb` |

??? note "上游说明（英文原文）"

    Search for matches within a UCSC Genome Browser genome assembly across tracks, help docs, and public hubs.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'search_genome',
    {'search': '...', 'genome': '...'},
)
```
