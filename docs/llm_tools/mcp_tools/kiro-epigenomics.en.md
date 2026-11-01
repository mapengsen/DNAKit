# AWS Kiro Life Sciences: epigenomics (`kiro-epigenomics`)

DNAKit's MCP tool index keeps **3** DNA/gene tools of this server (3 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| License | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| Command | `life-sciences-epigenomics` |
| Transports | stdio |
| Installation | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-epigenomics component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-epigenomics console script. No uvx or automatic install is used. |
| Effects | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`ihec_search`](#ihec_search) | `regulation_epigenome` | Searches International Human Epigenome Consortium (IHEC) epigenome datasets by tissue type. |
| [`methbase_search`](#methbase_search) | `regulation_epigenome` | Searches MethBase DNA methylomes by species and optional tissue. |
| [`roadmap_search`](#roadmap_search) | `regulation_epigenome` | Searches Roadmap Epigenomics data files by tissue and optional histone mark (e.g. H3K4me3). |

## `ihec_search` {#ihec_search}

Searches International Human Epigenome Consortium (IHEC) epigenome datasets by tissue type.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `tissue`
- **Keywords:** reference epigenome, histone, methylation, chip-seq, tissue

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `tissue` | string | yes | Tissue type to search for. |

??? note "Upstream description"

    Search IHEC datasets by tissue type.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-epigenomics',
    'ihec_search',
    {'tissue': '...'},
)
```

## `methbase_search` {#methbase_search}

Searches MethBase DNA methylomes by species and optional tissue.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `species`
- **Keywords:** dna methylation, methylome, bisulfite, cpg, wgbs

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | yes | Species name. |
| `tissue` | string | no | Tissue type. Default `` |

??? note "Upstream description"

    Search MethBase methylomes by species and optional tissue.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-epigenomics',
    'methbase_search',
    {'species': '...'},
)
```

## `roadmap_search` {#roadmap_search}

Searches Roadmap Epigenomics data files by tissue and optional histone mark (e.g. H3K4me3).

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `tissue`
- **Keywords:** histone mark, h3k4me3, chromatin, epigenome, tissue

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `tissue` | string | yes | Tissue type to search for. |
| `mark` | string | no | Histone mark (e.g. 'H3K4me3'). Default `` |

??? note "Upstream description"

    Search Roadmap Epigenomics data by tissue and optional histone mark.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-epigenomics',
    'roadmap_search',
    {'tissue': '...'},
)
```
