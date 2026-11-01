# Addgene MCP (`addgene`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **4** 个 DNA/基因相关工具（服务器共列出 4 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [longevity-genie/addgene-mcp @ `b402762`](https://github.com/longevity-genie/addgene-mcp/tree/b402762fdf10de6c5a179700f4dd7912a31682f2) |
| 许可证 | MIT |
| 启动命令 | `addgene-mcp stdio` |
| 传输方式 | stdio |
| 安装 | pip install addgene-mcp==0.1.3 (PyPI, Python &gt;=3.12) or install the pinned repository separately; the installed console script is addgene-mcp and 'addgene-mcp stdio' starts the stdio server. No credentials are required; MCP_HOST, MCP_PORT and MCP_TRANSPORT only affect the HTTP/SSE modes. |
| 副作用 | Scrapes https://www.addgene.org through a Scrapy subprocess for plasmid search and sequence information (Scrapy may keep a local HTTP cache); download_sequence fetches a GenBank or SnapGene file from Addgene and writes it to the current working directory or a caller-chosen directory, creating it if needed. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`addgene_download_sequence`](#addgene_download_sequence) | `molecular_tools` | Downloads an Addgene plasmid sequence file in GenBank (.gb) or SnapGene (.dna) format to a local directory and returns the file path. |
| [`addgene_get_popular_plasmids`](#addgene_get_popular_plasmids) | `molecular_tools` | Lists up to 50 of the most-requested Addgene plasmids (100+ requests each). |
| [`addgene_get_sequence_info`](#addgene_get_sequence_info) | `molecular_tools` | Returns the download URL and availability status of the GenBank or SnapGene sequence file for an Addgene plasmid ID. |
| [`addgene_search_plasmids`](#addgene_search_plasmids) | `molecular_tools` | Searches the Addgene plasmid repository by text query with filters for species, expression system, vector type (e.g. CRISPR), plasmid type and resistance marke… |

## `addgene_download_sequence` {#addgene_download_sequence}

Downloads an Addgene plasmid sequence file in GenBank (.gb) or SnapGene (.dna) format to a local directory and returns the file path.

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`plasmid_id`
- **关键词：**plasmid, genbank, snapgene, vector sequence, download

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plasmid_id` | integer | 是 |  |
| `format` | string | 否 | 可选值：`genbank`, `snapgene` 默认 `genbank` |
| `download_directory` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Download a plasmid sequence file to the local filesystem. This actually downloads the file from Addgene and saves it locally. Use EXACTLY format='genbank' for .gb files (default) or format='snapgene' for .dna files. These are the only two valid format values. Returns the local file path and download status.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'addgene',
    'addgene_download_sequence',
    {'plasmid_id': 0},
)
```

## `addgene_get_popular_plasmids` {#addgene_get_popular_plasmids}

Lists up to 50 of the most-requested Addgene plasmids (100+ requests each).

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**无
- **关键词：**plasmid, popular, vector, repository, cloning

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `page_size` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Get the most popular plasmids from Addgene (100+ requests each). This is a convenience method for finding widely-used, well-established plasmids. Returns up to 50 results.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'addgene',
    'addgene_get_popular_plasmids',
    {},
)
```

## `addgene_get_sequence_info` {#addgene_get_sequence_info}

Returns the download URL and availability status of the GenBank or SnapGene sequence file for an Addgene plasmid ID.

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**`plasmid_id`
- **关键词：**plasmid, sequence file, genbank, snapgene, download url

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plasmid_id` | integer | 是 |  |
| `format` | string | 否 | 可选值：`genbank`, `snapgene` 默认 `snapgene` |

??? note "上游说明（英文原文）"

    Get information about downloading a plasmid sequence file. Requires a valid Addgene plasmid ID. Returns download URL and availability status. Use EXACTLY format='snapgene' for SnapGene .dna files or format='genbank' for GenBank .gb files. These are the only two valid format values.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'addgene',
    'addgene_get_sequence_info',
    {'plasmid_id': 0},
)
```

## `addgene_search_plasmids` {#addgene_search_plasmids}

Searches the Addgene plasmid repository by text query with filters for species, expression system, vector type (e.g. CRISPR), plasmid type and resistance marker.

- **类别：**`molecular_tools`（Plasmids, primers, CRISPR guides and cloning）
- **必填参数：**无
- **关键词：**plasmid, crispr, vector, expression system, cloning backbone, resistance marker

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string \| null | 否 | 默认 `None` |
| `page_size` | integer | 否 | 默认 `50` |
| `page_number` | integer | 否 | 默认 `1` |
| `expression` | string \| null | 否 | 默认 `None` |
| `vector_types` | string \| null | 否 | 默认 `None` |
| `species` | string \| null | 否 | 默认 `None` |
| `plasmid_type` | string \| null | 否 | 默认 `None` |
| `resistance_marker` | string \| null | 否 | 默认 `None` |
| `bacterial_resistance` | string \| null | 否 | 默认 `None` |
| `popularity` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Search for plasmids in the Addgene repository using text queries and filters. You MUST use the exact string values specified in the parameter documentation. Common use cases: find plasmids for gene expression (specify species and expression system), gene editing (vector_types='crispr'), cloning vectors (plasmid_type='empty_backbone'), or protein expression (expression='bacterial' or 'mammalian').

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'addgene',
    'addgene_search_plasmids',
    {},
)
```
