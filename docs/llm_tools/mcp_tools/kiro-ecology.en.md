# AWS Kiro Life Sciences: ecology (`kiro-ecology`)

DNAKit's MCP tool index keeps **3** DNA/gene tools of this server (7 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| License | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| Command | `life-sciences-ecology` |
| Transports | stdio |
| Installation | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-ecology component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-ecology console script. No uvx or automatic install is used. |
| Effects | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`bold_search`](#bold_search) | `sequence` | Searches BOLD Systems DNA barcode records (e.g. COI sequences) by taxon name. |
| [`genbank_env_search`](#genbank_env_search) | `sequence` | Searches GenBank environmental sample sequences by keyword via NCBI E-utilities. |
| [`mgnify_search`](#mgnify_search) | `sequencing_data` | Searches EBI MGnify metagenomics studies by keyword. |

## `bold_search` {#bold_search}

Searches BOLD Systems DNA barcode records (e.g. COI sequences) by taxon name.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `taxon`
- **Keywords:** dna barcode, coi, barcoding, taxon, species identification

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `taxon` | string | yes | Taxon name. |

??? note "Upstream description"

    Search BOLD barcode records by taxon name.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-ecology',
    'bold_search',
    {'taxon': '...'},
)
```

## `genbank_env_search` {#genbank_env_search}

Searches GenBank environmental sample sequences by keyword via NCBI E-utilities.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `query`
- **Keywords:** environmental sequences, genbank, nucleotide, env_sample, metabarcoding

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search keyword. |

??? note "Upstream description"

    Search GenBank environmental sequences by keyword.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-ecology',
    'genbank_env_search',
    {'query': '...'},
)
```

## `mgnify_search` {#mgnify_search}

Searches EBI MGnify metagenomics studies by keyword.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `keyword`
- **Keywords:** metagenomics, microbiome, amplicon, 16s, ebi

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `keyword` | string | yes | Search keyword. |

??? note "Upstream description"

    Search MGnify metagenomics studies by keyword.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-ecology',
    'mgnify_search',
    {'keyword': '...'},
)
```
