# AWS Kiro Life Sciences: clinical (`kiro-clinical`)

DNAKit's MCP tool index keeps **3** DNA/gene tools of this server (10 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| License | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| Command | `life-sciences-clinical` |
| Transports | stdio |
| Installation | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-clinical component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-clinical console script. No uvx or automatic install is used. |
| Effects | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`omim_entry`](#omim_entry) | `association` | Retrieves a full OMIM entry by MIM number, including gene-phenotype and allelic variant data (requires OMIM_API_KEY). |
| [`omim_search`](#omim_search) | `association` | Searches OMIM gene and genetic disorder entries by keyword (requires OMIM_API_KEY). |
| [`pharmgkb_search`](#pharmgkb_search) | `variant` | Searches PharmGKB pharmacogenomics data (genes, variants, drugs) by keyword. |

## `omim_entry` {#omim_entry}

Retrieves a full OMIM entry by MIM number, including gene-phenotype and allelic variant data (requires OMIM_API_KEY).

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `mim_number`
- **Keywords:** mendelian, mim number, genetic disorder, allelic variants, gene-phenotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `mim_number` | string | yes | MIM number (e.g. '154700'). |

??? note "Upstream description"

    Get OMIM entry by MIM number. Requires OMIM_API_KEY.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-clinical',
    'omim_entry',
    {'mim_number': '...'},
)
```

## `omim_search` {#omim_search}

Searches OMIM gene and genetic disorder entries by keyword (requires OMIM_API_KEY).

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `query`
- **Keywords:** mendelian, genetic disorder, inherited disease, mim, gene-phenotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search keyword (e.g. 'Marfan syndrome'). |
| `max_results` | integer | no | Maximum results to return. Default `10` |

??? note "Upstream description"

    Search OMIM genetic disorder entries. Requires OMIM_API_KEY.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-clinical',
    'omim_search',
    {'query': '...'},
)
```

## `pharmgkb_search` {#pharmgkb_search}

Searches PharmGKB pharmacogenomics data (genes, variants, drugs) by keyword.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query`
- **Keywords:** pharmacogenomics, pgx, cyp2d6, drug response, haplotype

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search keyword (e.g. 'warfarin', 'CYP2D6'). |

??? note "Upstream description"

    Search PharmGKB pharmacogenomics data by keyword.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-clinical',
    'pharmgkb_search',
    {'query': '...'},
)
```
