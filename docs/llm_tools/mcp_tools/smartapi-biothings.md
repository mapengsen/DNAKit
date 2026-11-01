# SmartAPI MCP (BioThings core) (`smartapi-biothings`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **3** 个 DNA/基因相关工具（服务器共列出 5 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [biothings/smartapi-mcp @ `2cd3f4d`](https://github.com/biothings/smartapi-mcp/tree/2cd3f4dc4d1db186f8106b257a8aa40d2f994673) |
| 许可证 | Apache-2.0 |
| 启动命令 | `smartapi-mcp --api_set biothings_core --facade on` |
| 传输方式 | stdio |
| 安装 | pip install smartapi-mcp (PyPI, version 0.6.0 at verification) into a Python &gt;=3.10 environment; the catalog arguments select the core BioThings API set with the generic five-tool facade forced on (MyGene.info, MyVariant.info, MyChem.info, MyDisease.info, MyGeneSet.info, MyTaxon.info). No credentials are needed; optionally set SMARTAPI_LOG_LEVEL or LOG_LEVEL. |
| 副作用 | Fetches API specs from the SmartAPI registry (smart-api.info) at startup and sends query and lookup requests (GET /query, GET and POST by id) to the public BioThings services; no files are written and no credentials are sent. The cataloged --api_set biothings_core selects only BioThings APIs. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`biothings_get`](#biothings_get) | `genome_annotation` | Retrieves one annotation by ID from a chosen BioThings API (mygene, myvariant, mygeneset, mychem, mydisease, mytaxon). |
| [`biothings_getbatch`](#biothings_getbatch) | `genome_annotation` | Retrieves annotations for multiple IDs from a chosen BioThings API, with field and scope selection. |
| [`biothings_query`](#biothings_query) | `genome_annotation` | Searches a chosen BioThings API (e.g. mygene with symbol:CDK2) by query string with fields and paging. |

## `biothings_get` {#biothings_get}

Retrieves one annotation by ID from a chosen BioThings API (mygene, myvariant, mygeneset, mychem, mydisease, mytaxon).

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`api`, `id`
- **关键词：**mygene, myvariant, biothings, annotation, entrez, hgvs

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `api` | string | 是 |  |
| `id` | string | 是 |  |
| `fields` | string | 否 | 默认 `all` |

??? note "上游说明（英文原文）"

    Retrieve a single annotation by id (GET /\{type\}/\{id\}).

    BioThings API to target, given as its short name. Call list_biothings_apis to discover APIs. Available: mychem, mydisease, mygene, mygeneset, mytaxon, myvariant

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`api`, `ids`
- **关键词：**mygene, myvariant, batch, biothings, id mapping

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `api` | string | 是 |  |
| `ids` | array of string | 是 |  |
| `fields` | string | 否 | 默认 `all` |
| `scopes` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve annotations for multiple ids (POST /\{type\}).

    BioThings API to target, given as its short name. Call list_biothings_apis to discover APIs. Available: mychem, mydisease, mygene, mygeneset, mytaxon, myvariant

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`api`, `q`
- **关键词：**mygene, myvariant, gene query, biothings, gene symbol

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `api` | string | 是 |  |
| `q` | string | 是 |  |
| `fields` | string | 否 | 默认 `all` |
| `size` | integer | 否 | 默认 `10` |
| `from_` | integer | 否 | 默认 `0` |

??? note "上游说明（英文原文）"

    Search a BioThings API. `api` is the API name (see list_biothings_apis); `q` is a query string (e.g. 'symbol:CDK2' or 'diabetes').

    BioThings API to target, given as its short name. Call list_biothings_apis to discover APIs. Available: mychem, mydisease, mygene, mygeneset, mytaxon, myvariant

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'smartapi-biothings',
    'biothings_query',
    {'api': '...', 'q': '...'},
)
```
