# BrAPI MCP (`brapi`)

DNAKit's MCP tool index keeps **6** DNA/gene tools of this server (22 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [cyanheads/brapi-mcp-server @ `b7492bd`](https://github.com/cyanheads/brapi-mcp-server/tree/b7492bdf24dcc0abc36626f234eae03ec5e1331b) |
| License | Apache-2.0 |
| Command | `node dist/index.js` |
| Transports | stdio, http |
| Installation | Separately prepare the pinned @cyanheads/brapi-mcp-server release/source with its dist build and dependencies; set workdir to the package root and executable_path to Node &gt;=24. The bridge never runs npm, npx, bunx or an installer. |
| Effects | Queries BrAPI servers; opt-in upstream tools can submit observations or export local data. Connection aliases and dataframes require call_mcp_workflow to preserve a session. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`brapi_export_genotype_matrix`](#brapi_export_genotype_matrix) | `breeding_population` | Pulls genotype calls for germplasm and a variant set and pivots them into a matrix, VCF-subset or PLINK ped/map output. |
| [`brapi_find_genotype_calls`](#brapi_find_genotype_calls) | `breeding_population` | Retrieves genotype calls filtered by variant set, germplasm, call set or variant, spilling large results to a dataframe. |
| [`brapi_find_germplasm`](#brapi_find_germplasm) | `breeding_population` | Finds germplasm by name, synonym, accession number, PUI, crop or collection. |
| [`brapi_find_variants`](#brapi_find_variants) | `variant` | Finds variant records by variant set, reference sequence or genomic region (start/end). |
| [`brapi_get_germplasm`](#brapi_get_germplasm) | `breeding_population` | Fetches one germplasm by DbId with attributes, direct parents and study/descendant counts. |
| [`brapi_walk_pedigree`](#brapi_walk_pedigree) | `breeding_population` | Walks germplasm ancestors and/or descendants as a deduplicated DAG with depth limits and cycle detection. |

## `brapi_export_genotype_matrix` {#brapi_export_genotype_matrix}

Pulls genotype calls for germplasm and a variant set and pivots them into a matrix, VCF-subset or PLINK ped/map output.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `variantSetDbId`, `format`
- **Keywords:** genotype matrix, vcf, plink, genotype calls, variant set, germplasm
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `alias` | string | no | Connection alias registered via brapi_connect. Omit to read the connection registered under alias `default` — i.e. a prior brapi_connect call that did not specify an alias. Calls that used a non-default alias must pass that same alias here. |
| `germplasmDbIds` | array of string | no | Restrict to these germplasm. Omit to pull all germplasm in the variant set (use with caution on large sets). |
| `variantSetDbId` | string | yes | Variant set to pull calls for. Required. |
| `format` | string | yes | Output format. `matrix-json` registers a wide canvas dataframe only. `vcf-lite` returns VCF-subset text and registers the dataframe. `plink` returns .ped/.map text and registers the dataframe. One of `plink`, `vcf-lite`, `matrix-json` |
| `maxCalls` | integer | no | Lower the pull cap for this call. Omit to use the deployment ceiling (BRAPI_GENOTYPE_CALLS_MAX_PULL). Cannot raise it: a value above the deployment ceiling is clamped down to it and the effective cap is reported in `warnings`. |
| `maxColumns` | integer | no | Lower the distinct-variant column cap for this call. Omit to use the deployment ceiling (BRAPI_GENOTYPE_MATRIX_MAX_COLUMNS). Cannot raise it: a value above the deployment ceiling is clamped down to it. When the variant set resolves more distinct variants than the effective cap, the matrix is capped at that many variant columns, `truncated` is set, and the effective cap is reported in `warnings`.… |

??? note "Upstream description"

    Pull genotype calls for a germplasm × variant set and pivot them into a matrix. `format` controls the output: `matrix-json` registers a wide germplasm × variant canvas dataframe for SQL analysis; `vcf-lite` returns VCF-subset text (in the `vcf` field) and also registers the dataframe; `plink` returns .ped/.map text (in the `ped`/`map` fields) and also registers the dataframe. vcf-lite/plink pull /variants metadata for CHROM/POS/REF/ALT (`.`/`0` when the server lacks them). Column names are SQL-safe identifiers; `variantColumnLegend` maps them back to original variant IDs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'brapi',
    'brapi_export_genotype_matrix',
    {'variantSetDbId': '...', 'format': 'plink'},
)
```

## `brapi_find_genotype_calls` {#brapi_find_genotype_calls}

Retrieves genotype calls filtered by variant set, germplasm, call set or variant, spilling large results to a dataframe.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** none
- **Keywords:** genotype calls, variant set, call set, germplasm, genotyping
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `alias` | string | no | Connection alias registered via brapi_connect. Omit to read the connection registered under alias `default` — i.e. a prior brapi_connect call that did not specify an alias. Calls that used a non-default alias must pass that same alias here. |
| `variantSetDbId` | string | no | Scope calls to a single variant set. Strongly recommended. |
| `variantSetDbIds` | array of string | no | Alternative: multiple variant sets at once. |
| `germplasmDbIds` | array of string | no | Restrict to these germplasm (call sets). |
| `callSetDbIds` | array of string | no | Restrict to these call sets directly. |
| `variantDbIds` | array of string | no | Restrict to specific variants. |
| `callFormat` | string | no | Requested call-encoding format, when the server honors it. One of `VCF`, `FLAPJACK`, `DARTSEQ`, `JSON` |
| `loadLimit` | integer | no | Cap on rows returned inline. Omit for the deployment default. When the collected set exceeds this, the full result lands in a dataframe and only the first `loadLimit` rows return inline — query the dataframe with brapi_dataframe_query (SQL) for the rest. Upstream pageSize is fixed for genotype calls, so this knob only affects the inline preview here (no spillover capacity tradeoff). |

??? note "Upstream description"

    Pull genotype calls for a germplasm × variant set. Filter to bound cost — at minimum, set `variantSetDbId` or `germplasmDbIds`. The upstream pull is capped by deployment policy; when the pull is truncated, narrow the filters or query the spilled dataframe. `loadLimit` bounds the rows returned inline; the full collected set is materialized as a dataframe — query it with brapi_dataframe_query (SQL) instead of paging row-by-row.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'brapi',
    'brapi_find_genotype_calls',
    {},
)
```

## `brapi_find_germplasm` {#brapi_find_germplasm}

Finds germplasm by name, synonym, accession number, PUI, crop or collection.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** none
- **Keywords:** germplasm, accession, crop, genebank, plant breeding
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `alias` | string | no | Connection alias registered via brapi_connect. Omit to read the connection registered under alias `default` — i.e. a prior brapi_connect call that did not specify an alias. Calls that used a non-default alias must pass that same alias here. |
| `names` | array of string | no | Filter by germplasm display names. |
| `germplasmDbIds` | array of string | no | Filter by DbIds. |
| `germplasmPUIs` | array of string | no | Persistent unique identifiers. |
| `accessionNumbers` | array of string | no | Filter by accession numbers (gene-bank catalog codes). |
| `crops` | array of string | no | Filter by common crop names. |
| `synonyms` | array of string | no | Match registered synonyms. |
| `collections` | array of string | no | Filter by germplasm collection names. |
| `genus` | string | no | Botanical genus. |
| `species` | string | no | Botanical species. |
| `text` | string | no | Free-text query. Applied client-side as a substring match on returned rows (germplasmName, accessionNumber, defaultDisplayName, registered synonyms) — no BrAPI server reliably supports a server-side free-text filter, so combine with `crops` / `genus` / etc. to narrow the upstream pull first. |
| `loadLimit` | integer | no | Cap on rows returned inline. Omit for the deployment default. Rows beyond the cap land in a dataframe; query with brapi_dataframe_query (SQL) instead of paging row-by-row. |
| `extraFilters` | object | no | Extra BrAPI filters forwarded verbatim. Valid keys vary by endpoint; brapi_describe_filters enumerates them. Named params on this tool take precedence on conflict. |

??? note "Upstream description"

    Find germplasm by name, synonym, accession number, PUI, crop, or free-text query. Matches across registered synonyms. When the upstream total exceeds loadLimit, the full result set is materialized as a dataframe — query it with brapi_dataframe_query (SQL) instead of paging row-by-row.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'brapi',
    'brapi_find_germplasm',
    {},
)
```

## `brapi_find_variants` {#brapi_find_variants}

Finds variant records by variant set, reference sequence or genomic region (start/end).

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** variant records, variant set, genomic region, reference, markers
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `alias` | string | no | Connection alias registered via brapi_connect. Omit to read the connection registered under alias `default` — i.e. a prior brapi_connect call that did not specify an alias. Calls that used a non-default alias must pass that same alias here. |
| `variantSets` | array of string | no | Filter by variantSetDbIds. |
| `variants` | array of string | no | Filter by variantDbIds. |
| `references` | array of string | no | Filter by referenceDbIds. |
| `referenceName` | string | no | Reference display name (e.g. "chr01", "chr1"). |
| `start` | integer | no | Inclusive 1-based start. |
| `end` | integer | no | Exclusive 1-based end. |
| `loadLimit` | integer | no | Cap on rows returned inline. Omit for the deployment default. Rows beyond the cap land in a dataframe; query with brapi_dataframe_query (SQL) instead of paging row-by-row. |
| `extraFilters` | object | no | Extra BrAPI filters forwarded verbatim. Valid keys vary by endpoint; brapi_describe_filters enumerates them. Named params on this tool take precedence on conflict. |

??? note "Upstream description"

    Find variant records by variant set, reference sequence, or genomic region (start/end, 1-based inclusive / exclusive). When the upstream total exceeds loadLimit, the full result set is materialized as a dataframe — query it with brapi_dataframe_query (SQL).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'brapi',
    'brapi_find_variants',
    {},
)
```

## `brapi_get_germplasm` {#brapi_get_germplasm}

Fetches one germplasm by DbId with attributes, direct parents and study/descendant counts.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `germplasmDbId`
- **Keywords:** germplasm, pedigree, parents, accession
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `germplasmDbId` | string | yes | Germplasm identifier. |
| `alias` | string | no | Connection alias registered via brapi_connect. Omit to read the connection registered under alias `default` — i.e. a prior brapi_connect call that did not specify an alias. Calls that used a non-default alias must pass that same alias here. |

??? note "Upstream description"

    Fetch a single germplasm by DbId with attributes and direct parents. Response companions report study count, direct parent count, and direct descendant count — signals for pedigree depth and observation coverage.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'brapi',
    'brapi_get_germplasm',
    {'germplasmDbId': '...'},
)
```

## `brapi_walk_pedigree` {#brapi_walk_pedigree}

Walks germplasm ancestors and/or descendants as a deduplicated DAG with depth limits and cycle detection.

- **Category:** `breeding_population` (Breeding genotyping, markers, population genetics and strain resources)
- **Required arguments:** `germplasmDbIds`
- **Keywords:** pedigree, ancestry, progeny, germplasm, breeding
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `germplasmDbIds` | array of string | yes | Starting germplasm (1–20 roots). All roots are walked concurrently. |
| `direction` | string | no | Which direction to walk: ancestors (parents), descendants (progeny), or both. One of `ancestors`, `descendants`, `both` Default `ancestors` |
| `maxDepth` | integer | no | Max generations to walk per direction (default 3, cap 10). Default `3` |
| `loadLimit` | integer | no | Cap on rows returned inline. Omit for the deployment default. Rows beyond the cap land in a dataframe; query with brapi_dataframe_query (SQL) instead of paging row-by-row. |
| `alias` | string | no | Connection alias registered via brapi_connect. Omit to read the connection registered under alias `default` — i.e. a prior brapi_connect call that did not specify an alias. Calls that used a non-default alias must pass that same alias here. |

??? note "Upstream description"

    Walk germplasm ancestry or descendancy as a deduplicated DAG, with multi-generation traversal, cycle detection, and depth limits. Returns nodes + edges plus traversal stats (depthReached, rootCount, leafCount, cycleCount, deadEndCount).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'brapi',
    'brapi_walk_pedigree',
    {'germplasmDbIds': []},
)
```
