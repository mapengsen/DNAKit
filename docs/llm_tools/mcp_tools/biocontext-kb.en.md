# BioContextAI Knowledgebase MCP (`biocontext-kb`)

DNAKit's MCP tool index keeps **7** DNA/gene tools of this server (49 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [biocontext-ai/knowledgebase-mcp @ `7e99be7`](https://github.com/biocontext-ai/knowledgebase-mcp/tree/7e99be76f69ea89c713c4ddcc13842d70f8da6b2) |
| License | Apache-2.0 |
| Command | `biocontext_kb` |
| Transports | stdio, http |
| Installation | Install biocontext_kb separately in Python &gt;=3.12. DNAKit fixes MCP_ENVIRONMENT=DEVELOPMENT for the verified stdio entry point. |
| Effects | Queries the selected biomedical knowledge services. The public hosted endpoint is described by its authors as testing-only and subject to fair use. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`get_ensembl_id_from_gene_symbol`](#get_ensembl_id_from_gene_symbol) | `genome_annotation` | Returns the stable Ensembl gene ID (ENSG) for a gene symbol and taxonomy ID. |
| [`get_go_terms_by_gene`](#get_go_terms_by_gene) | `genome_annotation` | Searches OLS for Gene Ontology terms related to a gene name or symbol. |
| [`get_human_protein_atlas_info`](#get_human_protein_atlas_info) | `expression` | Retrieves Human Protein Atlas data for a gene (Ensembl ID or symbol): tissue expression, RNA/protein levels, subcellular location and pathology. |
| [`get_panglaodb_marker_genes`](#get_panglaodb_marker_genes) | `expression` | Retrieves PanglaoDB cell-type marker genes filtered by species, organ, cell type, gene symbol and sensitivity/specificity scores. |
| [`get_reactome_info_by_identifier`](#get_reactome_info_by_identifier) | `genome_annotation` | Queries Reactome for pathways associated with a molecule identifier, with species, p-value, sorting and disease-pathway filters. |
| [`get_uniprot_id_by_protein_symbol`](#get_uniprot_id_by_protein_symbol) | `genome_annotation` | Resolves a gene symbol, synonym or protein name to a UniProt accession for a species, preferring reviewed entries. |
| [`query_open_targets_graphql`](#query_open_targets_graphql) | `association` | Executes a GraphQL query against the Open Targets Platform API covering targets, diseases, drugs, variants and studies. |

## `get_ensembl_id_from_gene_symbol` {#get_ensembl_id_from_gene_symbol}

Returns the stable Ensembl gene ID (ENSG) for a gene symbol and taxonomy ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_symbol`
- **Keywords:** ensembl id, gene symbol, id mapping, ensg, taxonomy id

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | Gene name (e.g., 'TP53') |
| `species` | string | no | Taxonomy ID (e.g., 9606 for human, 10090 for mouse) Default `9606` |

??? note "Upstream description"

    Get Ensembl gene ID from gene symbol. Returns the stable Ensembl ID (ENSG*) for the given gene symbol and species.

    Returns:
    dict: Ensembl gene ID in format \{'ensembl_id': 'ENSG...'\} or error message.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'get_ensembl_id_from_gene_symbol',
    {'gene_symbol': '...'},
)
```

## `get_go_terms_by_gene` {#get_go_terms_by_gene}

Searches OLS for Gene Ontology terms related to a gene name or symbol.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_name`
- **Keywords:** gene ontology, go terms, ols, gene function

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string | yes | Gene name or symbol to search for (e.g., 'TP53', 'BRCA1') |
| `size` | integer | no | Maximum number of results to return Default `10` |
| `exact_match` | boolean | no | Whether to perform exact match search Default `False` |

??? note "Upstream description"

    Search OLS for Gene Ontology (GO) terms related to a gene name using structured vocabularies.

    Returns:
    dict: GO terms with go_terms array containing id, label, description, type or error message.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'get_go_terms_by_gene',
    {'gene_name': '...'},
)
```

## `get_human_protein_atlas_info` {#get_human_protein_atlas_info}

Retrieves Human Protein Atlas data for a gene (Ensembl ID or symbol): tissue expression, RNA/protein levels, subcellular location and pathology.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `gene_id`, `gene_symbol`
- **Keywords:** human protein atlas, tissue expression, rna expression, subcellular location, pathology

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string \| null | yes | Ensembl gene ID (e.g., 'ENSG00000141510') |
| `gene_symbol` | string \| null | yes | Gene symbol (e.g., 'TP53') |

??? note "Upstream description"

    Retrieve Human Protein Atlas information including expression, localization, and pathology data. Provide either gene_id or gene_symbol.

    Returns:
    dict: Protein atlas data with tissue_expression, subcellular_location, pathology, antibodies, RNA/protein levels or error message.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'get_human_protein_atlas_info',
    {'gene_id': '...', 'gene_symbol': '...'},
)
```

## `get_panglaodb_marker_genes` {#get_panglaodb_marker_genes}

Retrieves PanglaoDB cell-type marker genes filtered by species, organ, cell type, gene symbol and sensitivity/specificity scores.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `species`
- **Keywords:** panglaodb, marker genes, cell type, single cell, sensitivity, specificity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | yes | Species: 'Hs' for Human or 'Mm' for Mouse |
| `min_sensitivity` | number \| null | no | Minimum sensitivity score (0-1), applied to species-specific column Default `None` |
| `min_specificity` | number \| null | no | Minimum specificity score (0-1), applied to species-specific column Default `None` |
| `organ` | string \| null | no | Organ filter (e.g., 'Brain', 'Lung'), case-insensitive Default `None` |
| `cell_type` | string \| null | no | Cell type filter (e.g., 'Smooth muscle cells', 'T cells'), case-insensitive Default `None` |
| `gene_symbol` | string \| null | no | Gene symbol filter (e.g., 'MAFB', 'SYNPO'), case-insensitive Default `None` |

??? note "Upstream description"

    Retrieve marker genes from PanglaoDB dataset with optional filters. Supports filtering by species, scores, organ, cell type, gene symbol.

    Returns:
    dict: Markers array with gene symbols, cell types, organs, sensitivity/specificity scores or error message.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'get_panglaodb_marker_genes',
    {'species': '...'},
)
```

## `get_reactome_info_by_identifier` {#get_reactome_info_by_identifier}

Queries Reactome for pathways associated with a molecule identifier, with species, p-value, sorting and disease-pathway filters.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `identifier`
- **Keywords:** reactome, pathway, gene to pathway, identifier, species

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `identifier` | string | yes | The identifier of the element to be retrieved |
| `base_url` | string | no | Base URL for the Reactome API Default `https://reactome.org/AnalysisService` |
| `interactors` | boolean | no | Include interactors Default `False` |
| `species` | string \| array \| null | no | List of species to filter the result (accepts taxonomy ids, species names and dbId) Default `None` |
| `page_size` | integer | no | Pathways per page Default `20` |
| `page` | integer | no | Page number Default `1` |
| `sort_by` | string | no | Field to sort results by (e.g., 'ENTITIES_PVALUE', 'ENTITIES_FDR') Default `ENTITIES_PVALUE` |
| `order` | string | no | Sort order ('ASC' or 'DESC') Default `ASC` |
| `resource` | string | no | Resource to filter by (TOTAL includes all molecule types) Default `TOTAL` |
| `p_value` | number | no | P-value threshold (only pathways with p-value &lt;= threshold will be returned) Default `1.0` |
| `include_disease` | boolean | no | Set to False to exclude disease pathways Default `True` |
| `min_entities` | integer \| null | no | Minimum number of contained entities per pathway Default `None` |
| `max_entities` | integer \| null | no | Maximum number of contained entities per pathway Default `None` |
| `importable_only` | boolean | no | Filter to only include importable resources Default `False` |
| `timeout` | integer | no | Request timeout in seconds Default `30` |

??? note "Upstream description"

    Query the Reactome API identifier endpoint.

    Use this endpoint to retrieve pathways associated with a given identifier.
    Always provide the species parameter to ensure the correct protein is returned.

    Args:
    identifier (str): The identifier of the element to be retrieved
    base_url (str): Base URL for the Reactome API
    interactors (bool): Include interactors
    species (str or list): List of species to filter the result (accepts taxonomy ids, species names and dbId)
    page_size (int): Pathways per page
    page (int): Page number
    sort_by (str): Field to sort results by (e.g., "ENTITIES_PVALUE", "ENTITIES_FDR")
    order (str): Sort order ("ASC" or "DESC")
    resource (str): Resource to filter by (TOTAL includes all molecule types)
    p_value (float): P-value threshold (only pathways with p-value &lt;= threshold will be returned)
    include_disease (bool): Set to False to exclude disease pathways
    min_entities (int): Minimum number of contained entities per pathway
    max_entities (int): Maximum number of contained entities per pathway
    importable_only (bool): Filter to only include importable resources
    timeout (int): Request timeout in seconds

    Returns:
    dict: API response data or error information

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'get_reactome_info_by_identifier',
    {'identifier': '...'},
)
```

## `get_uniprot_id_by_protein_symbol` {#get_uniprot_id_by_protein_symbol}

Resolves a gene symbol, synonym or protein name to a UniProt accession for a species, preferring reviewed entries.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `protein_symbol`
- **Keywords:** uniprot accession, gene symbol, id mapping, swiss-prot

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `protein_symbol` | string | yes | Gene symbol (e.g., 'TP53', 'SYNPO'), gene synonym, or exact protein name (e.g., 'Podocin') |
| `species` | string | no | Organism taxonomy ID (e.g., '9606' for human) Default `9606` |

??? note "Upstream description"

    Retrieve the UniProt accession for a gene symbol and species. Prefers reviewed (Swiss-Prot) entries.

    Returns:
    dict: accession, gene, protein_name, reviewed, match_type ('gene', 'synonym' or 'protein_name') and other candidates, or error message (with candidates if ambiguous or only partially matching).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'get_uniprot_id_by_protein_symbol',
    {'protein_symbol': '...'},
)
```

## `query_open_targets_graphql` {#query_open_targets_graphql}

Executes a GraphQL query against the Open Targets Platform API covering targets, diseases, drugs, variants and studies.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `query_string`
- **Keywords:** open targets, target-disease, genetic evidence, graphql, drug target

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query_string` | string | yes | GraphQL query string starting with 'query' keyword |
| `variables` | object \| null | no | Optional variables for the GraphQL query Default `None` |

??? note "Upstream description"

    Execute GraphQL queries against the Open Targets API. Use get_open_targets_query_examples() or get_open_targets_graphql_schema() first.

    Returns:
    dict: GraphQL response with data field containing targets, diseases, drugs, variants, studies or error message.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'query_open_targets_graphql',
    {'query_string': '...'},
)
```
