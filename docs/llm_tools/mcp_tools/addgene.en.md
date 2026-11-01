# Addgene MCP (`addgene`)

DNAKit's MCP tool index keeps **4** DNA/gene tools of this server (4 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [longevity-genie/addgene-mcp @ `b402762`](https://github.com/longevity-genie/addgene-mcp/tree/b402762fdf10de6c5a179700f4dd7912a31682f2) |
| License | MIT |
| Command | `addgene-mcp stdio` |
| Transports | stdio |
| Installation | pip install addgene-mcp==0.1.3 (PyPI, Python &gt;=3.12) or install the pinned repository separately; the installed console script is addgene-mcp and 'addgene-mcp stdio' starts the stdio server. No credentials are required; MCP_HOST, MCP_PORT and MCP_TRANSPORT only affect the HTTP/SSE modes. |
| Effects | Scrapes https://www.addgene.org through a Scrapy subprocess for plasmid search and sequence information (Scrapy may keep a local HTTP cache); download_sequence fetches a GenBank or SnapGene file from Addgene and writes it to the current working directory or a caller-chosen directory, creating it if needed. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`addgene_download_sequence`](#addgene_download_sequence) | `molecular_tools` | Downloads an Addgene plasmid sequence file in GenBank (.gb) or SnapGene (.dna) format to a local directory and returns the file path. |
| [`addgene_get_popular_plasmids`](#addgene_get_popular_plasmids) | `molecular_tools` | Lists up to 50 of the most-requested Addgene plasmids (100+ requests each). |
| [`addgene_get_sequence_info`](#addgene_get_sequence_info) | `molecular_tools` | Returns the download URL and availability status of the GenBank or SnapGene sequence file for an Addgene plasmid ID. |
| [`addgene_search_plasmids`](#addgene_search_plasmids) | `molecular_tools` | Searches the Addgene plasmid repository by text query with filters for species, expression system, vector type (e.g. CRISPR), plasmid type and resistance marke… |

## `addgene_download_sequence` {#addgene_download_sequence}

Downloads an Addgene plasmid sequence file in GenBank (.gb) or SnapGene (.dna) format to a local directory and returns the file path.

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `plasmid_id`
- **Keywords:** plasmid, genbank, snapgene, vector sequence, download

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `plasmid_id` | integer | yes |  |
| `format` | string | no | One of `genbank`, `snapgene` Default `genbank` |
| `download_directory` | string \| null | no | Default `None` |

??? note "Upstream description"

    Download a plasmid sequence file to the local filesystem. This actually downloads the file from Addgene and saves it locally. Use EXACTLY format='genbank' for .gb files (default) or format='snapgene' for .dna files. These are the only two valid format values. Returns the local file path and download status.

Call it through DNAKit once the server is installed:

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

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** none
- **Keywords:** plasmid, popular, vector, repository, cloning

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `page_size` | integer | no | Default `20` |

??? note "Upstream description"

    Get the most popular plasmids from Addgene (100+ requests each). This is a convenience method for finding widely-used, well-established plasmids. Returns up to 50 results.

Call it through DNAKit once the server is installed:

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

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** `plasmid_id`
- **Keywords:** plasmid, sequence file, genbank, snapgene, download url

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `plasmid_id` | integer | yes |  |
| `format` | string | no | One of `genbank`, `snapgene` Default `snapgene` |

??? note "Upstream description"

    Get information about downloading a plasmid sequence file. Requires a valid Addgene plasmid ID. Returns download URL and availability status. Use EXACTLY format='snapgene' for SnapGene .dna files or format='genbank' for GenBank .gb files. These are the only two valid format values.

Call it through DNAKit once the server is installed:

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

- **Category:** `molecular_tools` (Plasmids, primers, CRISPR guides and cloning)
- **Required arguments:** none
- **Keywords:** plasmid, crispr, vector, expression system, cloning backbone, resistance marker

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string \| null | no | Default `None` |
| `page_size` | integer | no | Default `50` |
| `page_number` | integer | no | Default `1` |
| `expression` | string \| null | no | Default `None` |
| `vector_types` | string \| null | no | Default `None` |
| `species` | string \| null | no | Default `None` |
| `plasmid_type` | string \| null | no | Default `None` |
| `resistance_marker` | string \| null | no | Default `None` |
| `bacterial_resistance` | string \| null | no | Default `None` |
| `popularity` | string \| null | no | Default `None` |

??? note "Upstream description"

    Search for plasmids in the Addgene repository using text queries and filters. You MUST use the exact string values specified in the parameter documentation. Common use cases: find plasmids for gene expression (specify species and expression system), gene editing (vector_types='crispr'), cloning vectors (plasmid_type='empty_backbone'), or protein expression (expression='bacterial' or 'mammalian').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'addgene',
    'addgene_search_plasmids',
    {},
)
```
