# BioMCP (`biomcp`)

DNAKit's MCP tool index keeps **5** DNA/gene tools of this server (7 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [genomoncology/biomcp @ `3b1d645`](https://github.com/genomoncology/biomcp/tree/3b1d6452ab63741562d2ef97495a97281b304729) |
| License | MIT |
| Command | `biomcp serve` |
| Transports | stdio |
| Installation | Install genomoncology/biomcp separately; its PyPI distribution is biomcp-cli, not the unrelated biomcp distribution. The installed executable is biomcp. |
| Effects | Contacts upstream biomedical APIs according to the selected tool; some sources require credentials. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`gene_cspec`](#gene_cspec) | `variant` | Retrieves ClinGen CSpec criteria-specification manifests, pages or attachment metadata for a given gene. |
| [`get`](#get) | `variant` | Fetches one biomedical record by ID for an entity type (gene, variant, drug, disease, trial, article, pathway, protein, etc.) with selectable sections. |
| [`search`](#search) | `variant` | Searches a biomedical entity type (gene, variant, gwas, pgx, article, trial, protein, author) with typed, bounded filters. |
| [`variant_erepo`](#variant_erepo) | `variant` | Retrieves ClinGen ERepo expert-panel variant classification assertions by ClinGen allele ID(s) or by gene. |
| [`variant_normalize_car`](#variant_normalize_car) | `variant` | Normalizes 1-50 RefSeq HGVS expressions through the ClinGen Allele Registry to canonical allele identifiers. |

## `gene_cspec` {#gene_cspec}

Retrieves ClinGen CSpec criteria-specification manifests, pages or attachment metadata for a given gene.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene`
- **Keywords:** clingen, acmg criteria, variant curation, vcep, classification rules

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `capture_id` | string \| null | no | Default `None` |
| `files` | boolean | no | Default `False` |
| `gene` | string | yes |  |
| `limit` | integer | no | Default `25` |
| `offset` | integer | no | Default `0` |
| `version_iri` | string \| null | no | Default `None` |

??? note "Upstream description"

    Get ClinGen CSpec manifests, pages, or attachment metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biomcp',
    'gene_cspec',
    {'gene': '...'},
)
```

## `get` {#get}

Fetches one biomedical record by ID for an entity type (gene, variant, drug, disease, trial, article, pathway, protein, etc.) with selectable sections.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clinvar, population frequency, gene record, variant record, civic, cosmic, gencc

??? note "Upstream description"

    Get one biomedical record with typed inputs. Binary assets are CLI-only.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biomcp',
    'get',
    {},
)
```

## `search` {#search}

Searches a biomedical entity type (gene, variant, gwas, pgx, article, trial, protein, author) with typed, bounded filters.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** gene search, variant search, gwas, pharmacogenomics, chromosome region, gene type

??? note "Upstream description"

    Search a biomedical entity with typed, bounded inputs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biomcp',
    'search',
    {},
)
```

## `variant_erepo` {#variant_erepo}

Retrieves ClinGen ERepo expert-panel variant classification assertions by ClinGen allele ID(s) or by gene.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clingen, erepo, caid, pathogenicity, expert panel, acmg

??? note "Upstream description"

    Get ClinGen ERepo expert assertions by CAid or gene.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biomcp',
    'variant_erepo',
    {},
)
```

## `variant_normalize_car` {#variant_normalize_car}

Normalizes 1-50 RefSeq HGVS expressions through the ClinGen Allele Registry to canonical allele identifiers.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `inputs`
- **Keywords:** hgvs, clingen allele registry, caid, normalization, refseq

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `inputs` | array of string | yes |  |

??? note "Upstream description"

    Normalize 1-50 RefSeq HGVS values with ClinGen Allele Registry.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biomcp',
    'variant_normalize_car',
    {'inputs': []},
)
```
