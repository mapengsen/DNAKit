# AWS Kiro Life Sciences: agriculture (`kiro-agriculture`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **4** 个 DNA/基因相关工具（服务器共列出 4 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| 许可证 | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| 启动命令 | `life-sciences-agriculture` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-agriculture component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-agriculture console script. No uvx or automatic install is used. |
| 副作用 | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`gramene_search`](#gramene_search) | `genome_annotation` | Searches Gramene comparative plant genomics data (genes, orthologs, genomes) by keyword. |
| [`phytozome_search`](#phytozome_search) | `genome_annotation` | Searches the JGI Phytozome plant genome portal by keyword. |
| [`plantgdb_search`](#plantgdb_search) | `genome_annotation` | Searches PlantGDB plant genome and sequence data by keyword. |
| [`tair_search`](#tair_search) | `genome_annotation` | Searches TAIR (The Arabidopsis Information Resource) data by keyword. |

## `gramene_search` {#gramene_search}

Searches Gramene comparative plant genomics data (genes, orthologs, genomes) by keyword.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**plant genomics, comparative genomics, crop, orthologs, plant genes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query. |

??? note "上游说明（英文原文）"

    Search Gramene comparative genomics data by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-agriculture',
    'gramene_search',
    {'query': '...'},
)
```

## `phytozome_search` {#phytozome_search}

Searches the JGI Phytozome plant genome portal by keyword.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**plant genome, jgi, crop, gene families, assembly

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query. |

??? note "上游说明（英文原文）"

    Search Phytozome plant genomes by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-agriculture',
    'phytozome_search',
    {'query': '...'},
)
```

## `plantgdb_search` {#plantgdb_search}

Searches PlantGDB plant genome and sequence data by keyword.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**plant genome, crop, est, gene models, sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query. |

??? note "上游说明（英文原文）"

    Search PlantGDB plant genome data by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-agriculture',
    'plantgdb_search',
    {'query': '...'},
)
```

## `tair_search` {#tair_search}

Searches TAIR (The Arabidopsis Information Resource) data by keyword.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**arabidopsis, thaliana, plant, agi locus, gene

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search query. |

??? note "上游说明（英文原文）"

    Search TAIR Arabidopsis data by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-agriculture',
    'tair_search',
    {'query': '...'},
)
```
