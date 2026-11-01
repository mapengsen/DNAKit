# BioContextAI Knowledgebase MCP (`biocontext-kb`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **7** 个 DNA/基因相关工具（服务器共列出 49 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [biocontext-ai/knowledgebase-mcp @ `7e99be7`](https://github.com/biocontext-ai/knowledgebase-mcp/tree/7e99be76f69ea89c713c4ddcc13842d70f8da6b2) |
| 许可证 | Apache-2.0 |
| 启动命令 | `biocontext_kb` |
| 传输方式 | stdio, http |
| 安装 | Install biocontext_kb separately in Python &gt;=3.12. DNAKit fixes MCP_ENVIRONMENT=DEVELOPMENT for the verified stdio entry point. |
| 副作用 | Queries the selected biomedical knowledge services. The public hosted endpoint is described by its authors as testing-only and subject to fair use. |

| 工具 | 类别 | 简介 |
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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_symbol`
- **关键词：**ensembl id, gene symbol, id mapping, ensg, taxonomy id

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 | Gene name (e.g., 'TP53') |
| `species` | string | 否 | Taxonomy ID (e.g., 9606 for human, 10090 for mouse) 默认 `9606` |

??? note "上游说明（英文原文）"

    Get Ensembl gene ID from gene symbol. Returns the stable Ensembl ID (ENSG*) for the given gene symbol and species.

    Returns:
    dict: Ensembl gene ID in format \{'ensembl_id': 'ENSG...'\} or error message.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_name`
- **关键词：**gene ontology, go terms, ols, gene function

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene name or symbol to search for (e.g., 'TP53', 'BRCA1') |
| `size` | integer | 否 | Maximum number of results to return 默认 `10` |
| `exact_match` | boolean | 否 | Whether to perform exact match search 默认 `False` |

??? note "上游说明（英文原文）"

    Search OLS for Gene Ontology (GO) terms related to a gene name using structured vocabularies.

    Returns:
    dict: GO terms with go_terms array containing id, label, description, type or error message.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_id`, `gene_symbol`
- **关键词：**human protein atlas, tissue expression, rna expression, subcellular location, pathology

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string \| null | 是 | Ensembl gene ID (e.g., 'ENSG00000141510') |
| `gene_symbol` | string \| null | 是 | Gene symbol (e.g., 'TP53') |

??? note "上游说明（英文原文）"

    Retrieve Human Protein Atlas information including expression, localization, and pathology data. Provide either gene_id or gene_symbol.

    Returns:
    dict: Protein atlas data with tissue_expression, subcellular_location, pathology, antibodies, RNA/protein levels or error message.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`species`
- **关键词：**panglaodb, marker genes, cell type, single cell, sensitivity, specificity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 是 | Species: 'Hs' for Human or 'Mm' for Mouse |
| `min_sensitivity` | number \| null | 否 | Minimum sensitivity score (0-1), applied to species-specific column 默认 `None` |
| `min_specificity` | number \| null | 否 | Minimum specificity score (0-1), applied to species-specific column 默认 `None` |
| `organ` | string \| null | 否 | Organ filter (e.g., 'Brain', 'Lung'), case-insensitive 默认 `None` |
| `cell_type` | string \| null | 否 | Cell type filter (e.g., 'Smooth muscle cells', 'T cells'), case-insensitive 默认 `None` |
| `gene_symbol` | string \| null | 否 | Gene symbol filter (e.g., 'MAFB', 'SYNPO'), case-insensitive 默认 `None` |

??? note "上游说明（英文原文）"

    Retrieve marker genes from PanglaoDB dataset with optional filters. Supports filtering by species, scores, organ, cell type, gene symbol.

    Returns:
    dict: Markers array with gene symbols, cell types, organs, sensitivity/specificity scores or error message.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`identifier`
- **关键词：**reactome, pathway, gene to pathway, identifier, species

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `identifier` | string | 是 | The identifier of the element to be retrieved |
| `base_url` | string | 否 | Base URL for the Reactome API 默认 `https://reactome.org/AnalysisService` |
| `interactors` | boolean | 否 | Include interactors 默认 `False` |
| `species` | string \| array \| null | 否 | List of species to filter the result (accepts taxonomy ids, species names and dbId) 默认 `None` |
| `page_size` | integer | 否 | Pathways per page 默认 `20` |
| `page` | integer | 否 | Page number 默认 `1` |
| `sort_by` | string | 否 | Field to sort results by (e.g., 'ENTITIES_PVALUE', 'ENTITIES_FDR') 默认 `ENTITIES_PVALUE` |
| `order` | string | 否 | Sort order ('ASC' or 'DESC') 默认 `ASC` |
| `resource` | string | 否 | Resource to filter by (TOTAL includes all molecule types) 默认 `TOTAL` |
| `p_value` | number | 否 | P-value threshold (only pathways with p-value &lt;= threshold will be returned) 默认 `1.0` |
| `include_disease` | boolean | 否 | Set to False to exclude disease pathways 默认 `True` |
| `min_entities` | integer \| null | 否 | Minimum number of contained entities per pathway 默认 `None` |
| `max_entities` | integer \| null | 否 | Maximum number of contained entities per pathway 默认 `None` |
| `importable_only` | boolean | 否 | Filter to only include importable resources 默认 `False` |
| `timeout` | integer | 否 | Request timeout in seconds 默认 `30` |

??? note "上游说明（英文原文）"

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

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`protein_symbol`
- **关键词：**uniprot accession, gene symbol, id mapping, swiss-prot

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `protein_symbol` | string | 是 | Gene symbol (e.g., 'TP53', 'SYNPO'), gene synonym, or exact protein name (e.g., 'Podocin') |
| `species` | string | 否 | Organism taxonomy ID (e.g., '9606' for human) 默认 `9606` |

??? note "上游说明（英文原文）"

    Retrieve the UniProt accession for a gene symbol and species. Prefers reviewed (Swiss-Prot) entries.

    Returns:
    dict: accession, gene, protein_name, reviewed, match_type ('gene', 'synonym' or 'protein_name') and other candidates, or error message (with candidates if ambiguous or only partially matching).

安装服务器后通过 DNAKit 调用：

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

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query_string`
- **关键词：**open targets, target-disease, genetic evidence, graphql, drug target

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query_string` | string | 是 | GraphQL query string starting with 'query' keyword |
| `variables` | object \| null | 否 | Optional variables for the GraphQL query 默认 `None` |

??? note "上游说明（英文原文）"

    Execute GraphQL queries against the Open Targets API. Use get_open_targets_query_examples() or get_open_targets_graphql_schema() first.

    Returns:
    dict: GraphQL response with data field containing targets, diseases, drugs, variants, studies or error message.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biocontext-kb',
    'query_open_targets_graphql',
    {'query_string': '...'},
)
```
