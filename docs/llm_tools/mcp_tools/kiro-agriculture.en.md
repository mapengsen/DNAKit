# AWS Kiro Life Sciences: agriculture (`kiro-agriculture`)

DNAKit's MCP tool index keeps **4** DNA/gene tools of this server (4 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| License | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| Command | `life-sciences-agriculture` |
| Transports | stdio |
| Installation | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-agriculture component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-agriculture console script. No uvx or automatic install is used. |
| Effects | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`gramene_search`](#gramene_search) | `genome_annotation` | Searches Gramene comparative plant genomics data (genes, orthologs, genomes) by keyword. |
| [`phytozome_search`](#phytozome_search) | `genome_annotation` | Searches the JGI Phytozome plant genome portal by keyword. |
| [`plantgdb_search`](#plantgdb_search) | `genome_annotation` | Searches PlantGDB plant genome and sequence data by keyword. |
| [`tair_search`](#tair_search) | `genome_annotation` | Searches TAIR (The Arabidopsis Information Resource) data by keyword. |

## `gramene_search` {#gramene_search}

Searches Gramene comparative plant genomics data (genes, orthologs, genomes) by keyword.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** plant genomics, comparative genomics, crop, orthologs, plant genes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query. |

??? note "Upstream description"

    Search Gramene comparative genomics data by keyword.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** plant genome, jgi, crop, gene families, assembly

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query. |

??? note "Upstream description"

    Search Phytozome plant genomes by keyword.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** plant genome, crop, est, gene models, sequence

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query. |

??? note "Upstream description"

    Search PlantGDB plant genome data by keyword.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** arabidopsis, thaliana, plant, agi locus, gene

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search query. |

??? note "Upstream description"

    Search TAIR Arabidopsis data by keyword.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-agriculture',
    'tair_search',
    {'query': '...'},
)
```
