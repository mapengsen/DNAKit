# SmartAPI MCP (BioThings core) (`smartapi-biothings`)

DNAKit's MCP tool index keeps **3** DNA/gene tools of this server (5 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [biothings/smartapi-mcp @ `2cd3f4d`](https://github.com/biothings/smartapi-mcp/tree/2cd3f4dc4d1db186f8106b257a8aa40d2f994673) |
| License | Apache-2.0 |
| Command | `smartapi-mcp --api_set biothings_core --facade on` |
| Transports | stdio |
| Installation | pip install smartapi-mcp (PyPI, version 0.6.0 at verification) into a Python &gt;=3.10 environment; the catalog arguments select the core BioThings API set with the generic five-tool facade forced on (MyGene.info, MyVariant.info, MyChem.info, MyDisease.info, MyGeneSet.info, MyTaxon.info). No credentials are needed; optionally set SMARTAPI_LOG_LEVEL or LOG_LEVEL. |
| Effects | Fetches API specs from the SmartAPI registry (smart-api.info) at startup and sends query and lookup requests (GET /query, GET and POST by id) to the public BioThings services; no files are written and no credentials are sent. The cataloged --api_set biothings_core selects only BioThings APIs. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`biothings_get`](#biothings_get) | `genome_annotation` | Retrieves one annotation by ID from a chosen BioThings API (mygene, myvariant, mygeneset, mychem, mydisease, mytaxon). |
| [`biothings_getbatch`](#biothings_getbatch) | `genome_annotation` | Retrieves annotations for multiple IDs from a chosen BioThings API, with field and scope selection. |
| [`biothings_query`](#biothings_query) | `genome_annotation` | Searches a chosen BioThings API (e.g. mygene with symbol:CDK2) by query string with fields and paging. |

## `biothings_get` {#biothings_get}

Retrieves one annotation by ID from a chosen BioThings API (mygene, myvariant, mygeneset, mychem, mydisease, mytaxon).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `api`, `id`
- **Keywords:** mygene, myvariant, biothings, annotation, entrez, hgvs

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `api` | string | yes |  |
| `id` | string | yes |  |
| `fields` | string | no | Default `all` |

??? note "Upstream description"

    Retrieve a single annotation by id (GET /\{type\}/\{id\}).

    BioThings API to target, given as its short name. Call list_biothings_apis to discover APIs. Available: mychem, mydisease, mygene, mygeneset, mytaxon, myvariant

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'smartapi-biothings',
    'biothings_get',
    {'api': '...', 'id': '...'},
)
```

## `biothings_getbatch` {#biothings_getbatch}

Retrieves annotations for multiple IDs from a chosen BioThings API, with field and scope selection.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `api`, `ids`
- **Keywords:** mygene, myvariant, batch, biothings, id mapping

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `api` | string | yes |  |
| `ids` | array of string | yes |  |
| `fields` | string | no | Default `all` |
| `scopes` | string \| null | no | Default `None` |

??? note "Upstream description"

    Retrieve annotations for multiple ids (POST /\{type\}).

    BioThings API to target, given as its short name. Call list_biothings_apis to discover APIs. Available: mychem, mydisease, mygene, mygeneset, mytaxon, myvariant

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'smartapi-biothings',
    'biothings_getbatch',
    {'api': '...', 'ids': []},
)
```

## `biothings_query` {#biothings_query}

Searches a chosen BioThings API (e.g. mygene with symbol:CDK2) by query string with fields and paging.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `api`, `q`
- **Keywords:** mygene, myvariant, gene query, biothings, gene symbol

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `api` | string | yes |  |
| `q` | string | yes |  |
| `fields` | string | no | Default `all` |
| `size` | integer | no | Default `10` |
| `from_` | integer | no | Default `0` |

??? note "Upstream description"

    Search a BioThings API. `api` is the API name (see list_biothings_apis); `q` is a query string (e.g. 'symbol:CDK2' or 'diabetes').

    BioThings API to target, given as its short name. Call list_biothings_apis to discover APIs. Available: mychem, mydisease, mygene, mygeneset, mytaxon, myvariant

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'smartapi-biothings',
    'biothings_query',
    {'api': '...', 'q': '...'},
)
```
