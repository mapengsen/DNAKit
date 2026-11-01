# AWS Kiro Life Sciences: model-organisms (`kiro-model-organisms`)

DNAKit's MCP tool index keeps **5** DNA/gene tools of this server (5 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| License | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| Command | `life-sciences-model-organisms` |
| Transports | stdio |
| Installation | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-model-organisms component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-model-organisms console script. No uvx or automatic install is used. |
| Effects | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`flybase_gene`](#flybase_gene) | `genome_annotation` | Looks up a Drosophila melanogaster gene by symbol in FlyBase. |
| [`mgi_gene`](#mgi_gene) | `genome_annotation` | Looks up a mouse gene by symbol in Mouse Genome Informatics (MGI). |
| [`sgd_gene`](#sgd_gene) | `genome_annotation` | Looks up a Saccharomyces cerevisiae gene by name in the Saccharomyces Genome Database (SGD). |
| [`wormbase_gene`](#wormbase_gene) | `genome_annotation` | Looks up a Caenorhabditis elegans gene by name in WormBase. |
| [`zfin_gene`](#zfin_gene) | `genome_annotation` | Looks up a zebrafish (Danio rerio) gene by name in ZFIN. |

## `flybase_gene` {#flybase_gene}

Looks up a Drosophila melanogaster gene by symbol in FlyBase.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `symbol`
- **Keywords:** drosophila, fruit fly, model organism, fbgn, gene symbol

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Gene symbol (e.g. 'dpp', 'Notch'). |

??? note "Upstream description"

    Lookup a gene by symbol in FlyBase (Drosophila).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'flybase_gene',
    {'symbol': '...'},
)
```

## `mgi_gene` {#mgi_gene}

Looks up a mouse gene by symbol in Mouse Genome Informatics (MGI).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `symbol`
- **Keywords:** mouse, mus musculus, mgi id, model organism, gene symbol

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Gene symbol (e.g. 'Trp53', 'Brca1'). |

??? note "Upstream description"

    Lookup a gene by symbol in MGI (mouse).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'mgi_gene',
    {'symbol': '...'},
)
```

## `sgd_gene` {#sgd_gene}

Looks up a Saccharomyces cerevisiae gene by name in the Saccharomyces Genome Database (SGD).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `name`
- **Keywords:** yeast, saccharomyces cerevisiae, model organism, orf, gene name

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | yes | Gene name (e.g. 'CDC28', 'RAD51'). |

??? note "Upstream description"

    Lookup a gene by name in SGD (yeast).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'sgd_gene',
    {'name': '...'},
)
```

## `wormbase_gene` {#wormbase_gene}

Looks up a Caenorhabditis elegans gene by name in WormBase.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `name`
- **Keywords:** c. elegans, nematode, worm, model organism, gene name

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | yes | Gene name (e.g. 'unc-86', 'lin-12'). |

??? note "Upstream description"

    Lookup a gene by name in WormBase (C. elegans).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'wormbase_gene',
    {'name': '...'},
)
```

## `zfin_gene` {#zfin_gene}

Looks up a zebrafish (Danio rerio) gene by name in ZFIN.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `name`
- **Keywords:** zebrafish, danio rerio, model organism, zdb id, gene name

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | yes | Gene name (e.g. 'shha', 'tp53'). |

??? note "Upstream description"

    Lookup a gene by name in ZFIN (zebrafish).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-model-organisms',
    'zfin_gene',
    {'name': '...'},
)
```
