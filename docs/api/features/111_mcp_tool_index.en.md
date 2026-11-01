# Search the DNA and gene tools inside third-party MCP servers

Written for: DNAKit users and Agent builders who want to find functions that only exist in third-party MCP servers.

`list_mcp_servers()` lists 49 third-party MCP servers, but until now an Agent could only see each server's name and capability tags. The tools inside a server became visible only after installing and starting it. DNAKit now ships an offline index of those tools, reduced to the ones that work on DNA or genes, so `search_mcp_tools()` and the Agent server's `search_dnakit_tools` can find them without starting anything.

Every indexed tool, with its arguments and the call through DNAKit, is listed in the [third-party MCP tool catalogue](../../llm_tools/mcp_tools/index.md).

## How the index was built

- **Listing.** Each of the 49 servers was installed from its pinned commit and asked for `tools/list` on 2026-10-16. 46 were listed live (stdio or the hosted HTTP endpoint). `ncbi-database`, `seqeron-sequence` and `grebi` were read from the pinned source: the first fails the MCP handshake at that commit, .NET is not installed here, and grebi's hosted endpoint returned 404.
- **Hidden tools.** Two servers expose only a dispatcher in `tools/list`. ToolUniverse starts in compact mode (5 tools) and runs its 2,922 underlying tools through `execute_tool`; genomi runs 75 capability operations through `genomi.invoke`. Both registries were read from the pinned package and indexed with the dispatcher call they need.
- **Selection.** Every one of the 4,005 listed tools was read (description and input schema) and kept only when its main input or output is DNA or genes. 1,402 tools were kept. The other 2,603 are recorded with a reason: server housekeeping, literature search, drugs and chemistry, protein-only, non-genetic clinical data, generic ontology lookups, cell-level single-cell steps, generic platform management, or unrelated topics.
- Kept tools carry one of ten categories: `sequence`, `genome_annotation`, `variant`, `association`, `expression`, `regulation_epigenome`, `sequencing_data`, `prediction_model`, `molecular_tools` and `breeding_population`.

The index is a snapshot. Calling a tool still needs the server installed; `call_mcp_tool` validates arguments against the live schema of that installation.

## 1) `dnakit.integrations.mcp_tool_index_summary()`

- **Purpose：** Count the listed, kept and excluded MCP tools, per category and per server.
- **Method：** Reads the bundled index file; it starts no server and uses no network.
- **API：** `dnakit.integrations.mcp_tool_index_summary(verbose[optional])`.
- **Input：** None besides `verbose`.
- **Output：** A dict with `catalog_servers`, `listed_servers`, `servers_with_indexed_tools`, `listed_tools`, `indexed_tools`, `excluded_tools`, `indexed_by_category`, `excluded_by_reason`, `categories` and `servers`. Each `servers` entry has `server_id`, `listing`, `listed_date`, `source_commit`, `listed_tools`, `indexed_tools` and `note`.
- **Example code：**

```python
from dnakit.integrations import mcp_tool_index_summary

summary = mcp_tool_index_summary(verbose=False)
print(summary["catalog_servers"], summary["listed_tools"], summary["indexed_tools"], summary["excluded_tools"])
for category, count in summary["indexed_by_category"].items():
    print(category, count)
for server in summary["servers"]:
    if server["server_id"] in {"addgene", "genomi", "galaxy-mcp", "tooluniverse"}:
        print(server["server_id"], server["listing"], server["listed_tools"], server["indexed_tools"])
```

- **Example output：**

```text
49 4005 1402 2603
sequence 116
genome_annotation 360
variant 305
association 209
expression 77
regulation_epigenome 72
sequencing_data 122
prediction_model 66
molecular_tools 29
breeding_population 46
tooluniverse live-stdio+source-catalog 2927 845
genomi live-stdio+source-catalog 112 57
galaxy-mcp live-stdio 48 0
addgene live-stdio 4 4
```

## 2) `dnakit.integrations.search_mcp_tools()`

- **Purpose：** Find DNA and gene tools of third-party MCP servers by what they do.
- **Method：** Ranks the indexed tools with local Okapi BM25 over the tool name, the summary, the keywords, the category, the upstream description and the server name. With an empty query and a `server_id` or `category`, it lists those tools in index order.
- **API：** `dnakit.integrations.search_mcp_tools(query, server_id[optional], category[optional], limit[optional], verbose[optional])`.
- **Input：** `query` is English text (at most 1,000 characters). `server_id` is a catalog id, `category` one of the ten categories, `limit` an integer from 1 to 100 (default 10).
- **Output：** A dict with `query`, `count`, `total_matches`, `indexed_tools` and `tools`. Each tool has `server_id`, `tool_name`, `category`, `summary`, `required_arguments`, `listing`, `server_requires` (environment variables and installation note) and `call`, which says how to invoke it with `call_mcp_tool`. ToolUniverse tools that need an API key also carry `tool_requires`.
- **Example code：**

```python
from dnakit.integrations import search_mcp_tools

result = search_mcp_tools("polygenic risk score", limit=3, verbose=False)
print(result["count"], result["total_matches"], result["indexed_tools"])
for tool in result["tools"]:
    print(tool["server_id"], tool["tool_name"], tool["category"], tool["required_arguments"])
```

- **Example output：**

```text
3 170 1402
locus polygenic_risk variant []
genechat calculate_prs variant []
gwas-mcp calculate_prs variant ['sumstats_path', 'target_plink_prefix']
```

List one server's kept tools:

- **Example code：**

```python
from dnakit.integrations import search_mcp_tools

listed = search_mcp_tools("", server_id="gnomad", limit=3, verbose=False)
print(listed["total_matches"])
for tool in listed["tools"]:
    print(tool["tool_name"], "|", tool["summary"])
```

- **Example output：**

```text
11
get_clinvar_variant_info | Retrieves ClinVar information for a variant ID from a gnomAD dataset (v2/v3/v4) on GRCh37 or GRCh38.
get_copy_number_variant_info | Retrieves copy number variant (CNV) information by CNV ID from gnomAD v4 (GRCh38).
get_gene_info | Retrieves gnomAD v4 gene information by Ensembl gene ID or gene symbol.
```

## 3) `dnakit.integrations.describe_mcp_tool()`

- **Purpose：** Show one indexed tool's full description, input schema and server facts before calling it.
- **Method：** Looks the tool up in the index. A tool that was listed but left out raises an error that names the exclusion category and reason.
- **API：** `dnakit.integrations.describe_mcp_tool(server_id, tool_name, verbose[optional])`.
- **Input：** `server_id` and `tool_name` as returned by `search_mcp_tools`.
- **Output：** The search fields plus `description`, `keywords`, `input_schema`, `annotations` and `server` (name, repository, licence, pinned commit, command, arguments, transports, installation, environment variables and effects).
- **Example code：**

```python
from dnakit.integrations import describe_mcp_tool

tool = describe_mcp_tool("tooluniverse", "DNA_reverse_complement", verbose=False)
print(tool["category"], tool["listing"])
print(tool["call"])
print(sorted(tool["input_schema"]["properties"]))
print(tool["server"]["license"], tool["server"]["command"])
```

- **Example output：**

```text
sequence source-catalog
{'function': 'dnakit.integrations.call_mcp_tool', 'server_id': 'tooluniverse', 'tool_name': 'execute_tool', 'arguments': {'tool_name': 'DNA_reverse_complement', 'arguments': "the tool's own arguments (see input_schema)"}}
['operation', 'sequence']
Apache-2.0 tooluniverse
```

An excluded tool explains why it is not indexed:

- **Example code：**

```python
from dnakit.exceptions import ConfigurationError
from dnakit.integrations import describe_mcp_tool

try:
    describe_mcp_tool("galaxy-mcp", "run_tool", verbose=False)
except ConfigurationError as error:
    print(error.code, error.context["excluded_as"])
    print(error.context["reason"])
```

- **Example output：**

```text
MCP_INDEX_TOOL_NOT_INDEXED generic_platform
Galaxy tool discovery/execution is generic platform management, not a DNA function itself.
```

## 4) Calling a tool that was found

`call` in each result gives the arguments for `call_mcp_tool`. For most servers the tool is called by its own name. For ToolUniverse and genomi the call goes to the dispatcher (`execute_tool` or `genomi.invoke`) with the tool name and its arguments nested inside. After installing the server, for example `bioinformatics-mcp`:

```python
from dnakit.integrations import call_mcp_tool, search_mcp_tools

hit = search_mcp_tools("reverse complement", server_id="bioinformatics-mcp", limit=1, verbose=False)["tools"][0]
result = call_mcp_tool(
    hit["server_id"],
    hit["call"]["tool_name"],
    {"sequence": "ATGGCC"},
    executable_path="/path/to/venv/bin/bioinformatics-mcp",
    verbose=False,
)
print(result["status"], result["content"][0]["text"])
```

This path was run on 2026-10-16 for a direct tool (`bioinformatics-mcp` `reverse_complement`), a genomi operation (`sequence.translate` through `genomi.invoke`) and a ToolUniverse tool (`DNA_reverse_complement` through `execute_tool`); all three returned the reverse complement or translation of `ATGGCC`.

## 5) In the Agent server (`dnakit-mcp`)

`search_dnakit_tools` has a new argument `mcp_limit` (default 5, at most 50). When `category` is not set and the query is not empty, the response gains `mcp_tools`: up to `mcp_limit` results of `search_mcp_tools` for the same query. The `tools` list of DNAKit functions is unchanged. An Agent calls a found tool through `call_dnakit_tool` with `dnakit.integrations.call_mcp_tool`, and can read its schema first with `dnakit.integrations.describe_mcp_tool`. Set `mcp_limit=0` to leave MCP tools out. MCP results always use BM25, whatever `mode` is.

`integration_inventory()` now also reports `mcp_listed_servers`, `mcp_listed_tools` and `mcp_indexed_dna_gene_tools`.

## 6) Per-server counts

"Listed" is the full tool list read from each server; "kept" is the DNA/gene subset in the index.

| Server | Listing | Listed | Kept |
| --- | --- | ---: | ---: |
| `tooluniverse` | live + registry | 2,927 | 845 |
| `genomi` | live + registry | 112 | 57 |
| `plant-genomics` | live (stdio) | 50 | 40 |
| `bioinformatics-mcp` | live (stdio) | 70 | 33 |
| `gwas-mcp` | live (stdio) | 44 | 33 |
| `opentargets` | live (stdio) | 68 | 30 |
| `seqeron-sequence` | source | 35 | 30 |
| `gigwa-mcp` | live (stdio) | 35 | 27 |
| `alphagenome-jolab` | live (stdio) | 24 | 23 |
| `locus` | live (stdio) | 22 | 21 |
| `rewire-genomics` | live (stdio) | 24 | 20 |
| `chorus` | live (stdio) | 24 | 17 |
| `kiro-genomics` | live (stdio) | 18 | 17 |
| `just-prs` | live (stdio) | 20 | 16 |
| `gget` | live (stdio) | 20 | 15 |
| `grebi` | source | 33 | 15 |
| `mavedb-link` | live (HTTP) | 15 | 13 |
| `dnallm` | live (stdio) | 16 | 11 |
| `gnomad` | live (stdio) | 12 | 11 |
| `ucsc-genome` | live (stdio) | 12 | 11 |
| `biothings` | live (stdio) | 19 | 10 |
| `genechat` | live (stdio) | 10 | 10 |
| `encode-toolkit` | live (stdio) | 20 | 9 |
| `biocontext-kb` | live (HTTP) | 49 | 7 |
| `ensembl` | live (HTTP) | 7 | 7 |
| `uniprot-provenance` | live (stdio) | 41 | 7 |
| `brapi` | live (HTTP) | 22 | 6 |
| `juicebox` | live (stdio) | 17 | 6 |
| `biomcp` | live (stdio) | 7 | 5 |
| `kiro-model-organisms` | live (stdio) | 5 | 5 |
| `sequence-simulation` | live (stdio) | 6 | 5 |
| `addgene` | live (stdio) | 4 | 4 |
| `bioinformatics-ncbi-uniprot` | live (stdio) | 10 | 4 |
| `kiro-agriculture` | live (stdio) | 4 | 4 |
| `ncbi-database` | source | 4 | 4 |
| `ncbi-sra` | live (stdio) | 6 | 4 |
| `10x-genomics` | live (stdio) | 27 | 3 |
| `kiro-clinical` | live (stdio) | 10 | 3 |
| `kiro-ecology` | live (stdio) | 7 | 3 |
| `kiro-epigenomics` | live (stdio) | 3 | 3 |
| `smartapi-biothings` | live (stdio) | 5 | 3 |
| `just-dna-agents` | live (stdio) | 6 | 2 |
| `opentargets-platform` | live (HTTP) | 5 | 2 |
| `ebi-tools` | live (HTTP) | 4 | 1 |
| `galaxy-mcp` | live (stdio) | 48 | 0 |
| `genepattern` | live (stdio) | 52 | 0 |
| `kiro-pipelines` | live (stdio) | 8 | 0 |
| `pubmed` | live (HTTP) | 11 | 0 |
| `pubmed-anthropic` | live (HTTP) | 7 | 0 |

## Scope

- The selection follows written rules, applied tool by tool. Borderline tools have reasons that start with "borderline:"; `describe_mcp_tool` on an excluded tool shows its category and reason.
- The schemas are snapshots of the pinned commits. Upstream servers can add, rename or remove tools; `list_mcp_tools(server_id)` always shows what an installed server offers now.
- ToolUniverse and genomi entries came from their registries, not from `tools/list`; 208 ToolUniverse tools appear only when their API key, or an LLM key, is set (`tool_requires`). grebi's query-template tools depend on the deployment.
- Some servers need `mcp<2` to start (`genechat`, `gwas-mcp`, `encode-toolkit`, `ucsc-genome`, `bioinformatics-ncbi-uniprot`); their catalog `installation` notes say so.
- `others/tests/unit/test_mcp_tool_index.py` checks that every catalog server has a listing record, that kept and excluded tools add up to each listing, the search filters, the dispatcher call shape and the `mcp_tools` field of the Agent server.
