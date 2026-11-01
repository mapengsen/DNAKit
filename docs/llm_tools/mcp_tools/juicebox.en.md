# Juicebox MCP (`juicebox`)

DNAKit's MCP tool index keeps **6** DNA/gene tools of this server (17 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [aidenlab/juicebox-mcp @ `0149fd1`](https://github.com/aidenlab/juicebox-mcp/tree/0149fd1410ae37448dd899eb33b30082ae2748d5) |
| License | MIT |
| Command | `node server.js` |
| Transports | stdio |
| Installation | Prepare the official source checkout and npm dependencies separately, then set workdir to the checkout and executable_path to Node if necessary. |
| Effects | Starts the MCP server and local browser communication service; visualization tools may require an attached Juicebox frontend and network datasets. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`get_data_source_statistics`](#get_data_source_statistics) | `regulation_epigenome` | Returns unfiltered statistics for a Hi-C map source (4DN or ENCODE): total maps and breakdowns by assembly, biosource and other metadata. |
| [`get_map_details`](#get_map_details) | `regulation_epigenome` | Returns detailed metadata for one Hi-C contact map from 4DN or ENCODE search results, selected by result index or URL. |
| [`goto_locus`](#goto_locus) | `regulation_epigenome` | Navigates the loaded Hi-C contact map to a genomic locus given as coordinates, gene name, natural language or a chr/start/end object. |
| [`load_control_map`](#load_control_map) | `regulation_epigenome` | Loads a control Hi-C contact map (.hic file) from a URL with an optional normalization for comparison with the main map. |
| [`load_map`](#load_map) | `regulation_epigenome` | Loads a Hi-C contact map (.hic file) from a URL into Juicebox with optional name, normalization method and starting locus. |
| [`search_maps`](#search_maps) | `regulation_epigenome` | Searches Hi-C contact maps in 4DN and/or ENCODE by natural-language query over metadata such as assembly, biosource and biosample. |

## `get_data_source_statistics` {#get_data_source_statistics}

Returns unfiltered statistics for a Hi-C map source (4DN or ENCODE): total maps and breakdowns by assembly, biosource and other metadata.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `source`
- **Keywords:** hi-c, 4dn, encode, genome assembly, cell type, 3d genome, contact map catalog

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes | Data source ID ('4dn' or 'encode') |

??? note "Upstream description"

    Get statistical overview of a data source including total maps, assemblies covered, and breakdowns by metadata fields. Use this when users ask "what assemblies are available", "how many maps are there", "what cell types are covered", etc. This returns unfiltered statistics without search limits.

Call it through DNAKit once the server is installed:

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

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `source`
- **Keywords:** hi-c, contact map metadata, 4dn, encode, 3d genome, chromatin conformation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes | Data source ID ('4dn' or 'encode') |
| `index` | integer | no | Index from search results (0-based). Required if url is not provided. |
| `url` | string | no | Direct URL to the map. Required if index is not provided. |

??? note "Upstream description"

    Get detailed information about a specific Hi-C contact map. Use this when users want more information about a specific map from search results.

Call it through DNAKit once the server is installed:

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

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `locus`
- **Keywords:** hi-c, genomic coordinates, gene name, chromosome region, 3d genome, navigate

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `locus` | string \| object | yes | Locus to navigate to. Can be a string (natural language or standard format) or an object with chr, start, and end properties. |

??? note "Upstream description"

    Navigate to a specific genomic locus in the currently loaded map. Supports natural language, gene names, standard format, and structured objects. Examples: "chr1:1000-2000", "BRCA1", "chromosome 1 from 1000 to 2000", or \{chr: "chr1", start: 1000, end: 2000\}. When a single chromosome is specified, it applies to both axes of the Hi-C contact map.

Call it through DNAKit once the server is installed:

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

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `url`
- **Keywords:** hi-c, .hic file, control map, normalization, kr, 3d genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `url` | string | yes | URL to the control .hic file |
| `name` | string | no | Optional name for the control map |
| `normalization` | string | no | Normalization method (e.g., "VC", "VC_SQRT", "KR", "NONE") |

??? note "Upstream description"

    Load a control map (.hic file) for comparison

Call it through DNAKit once the server is installed:

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

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `url`
- **Keywords:** hi-c, .hic file, contact map, normalization, 3d genome, chromatin conformation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `url` | string | yes | URL to the .hic file |
| `name` | string | no | Optional name for the map |
| `normalization` | string | no | Normalization method (e.g., "VC", "VC_SQRT", "KR", "NONE") |
| `locus` | string | no | Optional genomic locus (e.g., "1:1000000-2000000 1:1000000-2000000") |

??? note "Upstream description"

    Load a Hi-C contact map (.hic file) into Juicebox

Call it through DNAKit once the server is installed:

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

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `query`
- **Keywords:** hi-c, 4dn, encode, cell line, genome assembly, 3d genome, dataset search

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | no | Data source ID ('4dn', 'encode') or 'all' to search all sources. Default: 'all' |
| `query` | string | yes | Natural language search query (e.g., "human hg38", "mouse cells", "K562") |
| `limit` | integer | no | Maximum number of results to return (default: 50) |

??? note "Upstream description"

    Search for Hi-C contact maps using natural language queries. Searches across all metadata fields (Assembly, Biosource, Biosample, Description, etc.). Use this when users want to find specific maps, e.g., "human hg38 maps", "mouse cell lines", "K562 cells", etc. NOTE: Results are limited to 50 by default. For statistical questions like "what assemblies are covered" or "how many maps are there", use get_data_source_statistics instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'juicebox',
    'search_maps',
    {'query': '...'},
)
```
