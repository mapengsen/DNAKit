# Juicebox MCP (`juicebox`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **6** 个 DNA/基因相关工具（服务器共列出 17 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [aidenlab/juicebox-mcp @ `0149fd1`](https://github.com/aidenlab/juicebox-mcp/tree/0149fd1410ae37448dd899eb33b30082ae2748d5) |
| 许可证 | MIT |
| 启动命令 | `node server.js` |
| 传输方式 | stdio |
| 安装 | Prepare the official source checkout and npm dependencies separately, then set workdir to the checkout and executable_path to Node if necessary. |
| 副作用 | Starts the MCP server and local browser communication service; visualization tools may require an attached Juicebox frontend and network datasets. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`get_data_source_statistics`](#get_data_source_statistics) | `regulation_epigenome` | Returns unfiltered statistics for a Hi-C map source (4DN or ENCODE): total maps and breakdowns by assembly, biosource and other metadata. |
| [`get_map_details`](#get_map_details) | `regulation_epigenome` | Returns detailed metadata for one Hi-C contact map from 4DN or ENCODE search results, selected by result index or URL. |
| [`goto_locus`](#goto_locus) | `regulation_epigenome` | Navigates the loaded Hi-C contact map to a genomic locus given as coordinates, gene name, natural language or a chr/start/end object. |
| [`load_control_map`](#load_control_map) | `regulation_epigenome` | Loads a control Hi-C contact map (.hic file) from a URL with an optional normalization for comparison with the main map. |
| [`load_map`](#load_map) | `regulation_epigenome` | Loads a Hi-C contact map (.hic file) from a URL into Juicebox with optional name, normalization method and starting locus. |
| [`search_maps`](#search_maps) | `regulation_epigenome` | Searches Hi-C contact maps in 4DN and/or ENCODE by natural-language query over metadata such as assembly, biosource and biosample. |

## `get_data_source_statistics` {#get_data_source_statistics}

Returns unfiltered statistics for a Hi-C map source (4DN or ENCODE): total maps and breakdowns by assembly, biosource and other metadata.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`source`
- **关键词：**hi-c, 4dn, encode, genome assembly, cell type, 3d genome, contact map catalog

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 是 | Data source ID ('4dn' or 'encode') |

??? note "上游说明（英文原文）"

    Get statistical overview of a data source including total maps, assemblies covered, and breakdowns by metadata fields. Use this when users ask "what assemblies are available", "how many maps are there", "what cell types are covered", etc. This returns unfiltered statistics without search limits.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'juicebox',
    'get_data_source_statistics',
    {'source': '...'},
)
```

## `get_map_details` {#get_map_details}

Returns detailed metadata for one Hi-C contact map from 4DN or ENCODE search results, selected by result index or URL.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`source`
- **关键词：**hi-c, contact map metadata, 4dn, encode, 3d genome, chromatin conformation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 是 | Data source ID ('4dn' or 'encode') |
| `index` | integer | 否 | Index from search results (0-based). Required if url is not provided. |
| `url` | string | 否 | Direct URL to the map. Required if index is not provided. |

??? note "上游说明（英文原文）"

    Get detailed information about a specific Hi-C contact map. Use this when users want more information about a specific map from search results.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'juicebox',
    'get_map_details',
    {'source': '...'},
)
```

## `goto_locus` {#goto_locus}

Navigates the loaded Hi-C contact map to a genomic locus given as coordinates, gene name, natural language or a chr/start/end object.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`locus`
- **关键词：**hi-c, genomic coordinates, gene name, chromosome region, 3d genome, navigate

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `locus` | string \| object | 是 | Locus to navigate to. Can be a string (natural language or standard format) or an object with chr, start, and end properties. |

??? note "上游说明（英文原文）"

    Navigate to a specific genomic locus in the currently loaded map. Supports natural language, gene names, standard format, and structured objects. Examples: "chr1:1000-2000", "BRCA1", "chromosome 1 from 1000 to 2000", or \{chr: "chr1", start: 1000, end: 2000\}. When a single chromosome is specified, it applies to both axes of the Hi-C contact map.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'juicebox',
    'goto_locus',
    {'locus': '...'},
)
```

## `load_control_map` {#load_control_map}

Loads a control Hi-C contact map (.hic file) from a URL with an optional normalization for comparison with the main map.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`url`
- **关键词：**hi-c, .hic file, control map, normalization, kr, 3d genome

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `url` | string | 是 | URL to the control .hic file |
| `name` | string | 否 | Optional name for the control map |
| `normalization` | string | 否 | Normalization method (e.g., "VC", "VC_SQRT", "KR", "NONE") |

??? note "上游说明（英文原文）"

    Load a control map (.hic file) for comparison

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'juicebox',
    'load_control_map',
    {'url': '...'},
)
```

## `load_map` {#load_map}

Loads a Hi-C contact map (.hic file) from a URL into Juicebox with optional name, normalization method and starting locus.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`url`
- **关键词：**hi-c, .hic file, contact map, normalization, 3d genome, chromatin conformation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `url` | string | 是 | URL to the .hic file |
| `name` | string | 否 | Optional name for the map |
| `normalization` | string | 否 | Normalization method (e.g., "VC", "VC_SQRT", "KR", "NONE") |
| `locus` | string | 否 | Optional genomic locus (e.g., "1:1000000-2000000 1:1000000-2000000") |

??? note "上游说明（英文原文）"

    Load a Hi-C contact map (.hic file) into Juicebox

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'juicebox',
    'load_map',
    {'url': '...'},
)
```

## `search_maps` {#search_maps}

Searches Hi-C contact maps in 4DN and/or ENCODE by natural-language query over metadata such as assembly, biosource and biosample.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`query`
- **关键词：**hi-c, 4dn, encode, cell line, genome assembly, 3d genome, dataset search

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 否 | Data source ID ('4dn', 'encode') or 'all' to search all sources. Default: 'all' |
| `query` | string | 是 | Natural language search query (e.g., "human hg38", "mouse cells", "K562") |
| `limit` | integer | 否 | Maximum number of results to return (default: 50) |

??? note "上游说明（英文原文）"

    Search for Hi-C contact maps using natural language queries. Searches across all metadata fields (Assembly, Biosource, Biosample, Description, etc.). Use this when users want to find specific maps, e.g., "human hg38 maps", "mouse cell lines", "K562 cells", etc. NOTE: Results are limited to 50 by default. For statistical questions like "what assemblies are covered" or "how many maps are there", use get_data_source_statistics instead.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'juicebox',
    'search_maps',
    {'query': '...'},
)
```
