# Open Targets Platform official MCP (`opentargets-platform`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **2** 个 DNA/基因相关工具（服务器共列出 5 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [anthropics/life-sciences @ `e96556b`](https://github.com/anthropics/life-sciences/tree/e96556b637b56d6cc3a5ad33987009be9e60aa5c) |
| 许可证 | Remote provider terms; the marketplace does not grant an independent server-code license |
| 启动命令 | — |
| 传输方式 | http |
| 安装 | Select transport=http. The provider-published endpoint supports anonymous MCP initialization and tools/list as verified on 2026-09-20. No local installer is run. |
| 副作用 | Queries are sent to the named public biological/literature service. Returned articles and data retain provider-specific terms. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`batch_query_open_targets_graphql`](#batch_query_open_targets_graphql) | `association` | Executes the same Open Targets Platform GraphQL query repeatedly with different variable sets, keyed by a field. |
| [`query_open_targets_graphql`](#query_open_targets_graphql) | `association` | Executes a GraphQL query against the Open Targets Platform API for targets, diseases, drugs, variants and studies. |

## `batch_query_open_targets_graphql` {#batch_query_open_targets_graphql}

Executes the same Open Targets Platform GraphQL query repeatedly with different variable sets, keyed by a field.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query_string`, `variables_list`, `key_field`
- **关键词：**open targets, graphql, target-disease, batch query, ensembl id, genetic evidence
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query_string` | string | 是 | The GraphQL query string to execute for all variable sets. |
| `variables_list` | array of object | 是 | List of variable dictionaries, one per query execution. |
| `key_field` | string | 是 | Variable field name to use as key in results mapping. |

??? note "上游说明（英文原文）"

    Execute the same GraphQL query multiple times with different variable sets.

    Use this tool instead of the regular query tool when you need to run the same query
    repeatedly with different arguments (e.g., querying multiple drugs, targets, or diseases).

    WORKFLOW - Follow these steps in order:

    Step 1: RESOLVE IDENTIFIERS
    If user provides common names (gene symbols, disease names, drug names),
    use `search_entity` tool FIRST to convert them to standardized IDs:

    - Targets/Genes: "BRCA1", "BRCA2" -&gt; ENSEMBL IDs "ENSG00000012048", "ENSG00000139618"
    - Diseases: "breast cancer" -&gt; EFO/MONDO ID "MONDO_0007254"
    - Drugs: "aspirin", "ibuprofen" -&gt; ChEMBL IDs "CHEMBL1201583", "CHEMBL521"
    - Variants: Use "chr_pos_ref_alt" format or rsIDs

    Example: search_entity(query_string="BRCA1 BRCA2", entity_names=\["target"\])

    Step 2: LEARN QUERY STRUCTURE
    Call `get_open_targets_graphql_schema` with relevant categories to retrieve
    the schema subset needed for your query. Select categories that cover the
    data domains you need - BE INCLUSIVE (it's better to include extra categories
    than to miss required types).

    Example: For a query about drug mechanisms and safety:
    get_open_targets_graphql_schema(categories=\["drug-mechanisms", "drug-safety"\])

    Study the returned schema to understand available types, fields, and their
    relationships, then construct a GraphQL query that fetches the information
    the user needs.

    FALLBACK: If you encounter errors or need detailed information about specific
    types, use `get_type_dependencies` sparingly to explore type relationships.
    This tool provides exhaustive type dependency information but should only be
    used when category-based retrieval is insufficient.

    Step 3: CONSTRUCT AND EXECUTE BATCH QUERY
    Build GraphQL query and variables_list using:
    - Standardized IDs from Step 1 (REQUIRED)
    - Query patterns from Step 2
    - Follow the "COMMON MISTAKES TO AVOID" guidance in the schema output

    Call this tool with query_string, variables_list, and key_field.

    REQUIRED IDENTIFIER FORMATS:
    - Targets/Genes: ENSEMBL IDs (e.g., "ENSG00000139618")
    - Diseases: EFO IDs (e.g., "EFO_0000305") or MONDO IDs (e.g., "MONDO_0007254")
    - Drugs: ChEMBL IDs (e.g., "CHEMBL1201583")
    - Variants: "chr_pos_ref_alt" format (e.g., "19_44908822_C_T") or rsIDs (e.g., "rs7412")
    - Studies: Study IDs (e.g., "GCST90002357")
    - Credible Sets: Study Locus IDs (e.g., "7d68cc9c70351c9dbd2a2c0c145e555d")

    Args:
    query_string (str): The GraphQL query string to execute for all variable sets.
    variables_list (list\[dict\[str, Any\]\]): List of variable dictionaries, one per query execution.
    key_field (str): Variable field name to use as key in results mapping.

    Returns:
    (BatchQueryResult): Results keyed by the specified field value, with execution summary.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets-platform',
    'batch_query_open_targets_graphql',
    {'query_string': '...', 'variables_list': [], 'key_field': '...'},
)
```

## `query_open_targets_graphql` {#query_open_targets_graphql}

Executes a GraphQL query against the Open Targets Platform API for targets, diseases, drugs, variants and studies.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query_string`
- **关键词：**open targets, graphql, target-disease, genetic evidence, ensembl id, l2g
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query_string` | string | 是 | GraphQL query string starting with 'query' keyword. |
| `variables` | object \| null | 否 | Optional dict or JSON string with query variables. 默认 `None` |

??? note "上游说明（英文原文）"

    Execute GraphQL queries against the Open Targets Platform API.

    WORKFLOW - Follow these steps in order:

    Step 1: RESOLVE IDENTIFIERS
    If user provides common names (gene symbols, disease names, drug names),
    use `search_entity` tool FIRST to convert them to standardized IDs:

    - Targets/Genes: "BRCA2" -&gt; ENSEMBL ID "ENSG00000139618"
    - Diseases: "breast cancer" -&gt; EFO/MONDO ID "MONDO_0007254"
    - Drugs: "aspirin" -&gt; ChEMBL ID "CHEMBL1201583"
    - Variants: Use "chr_pos_ref_alt" format or rsIDs

    Example: search_entity(query_string="BRCA2", entity_names=\["target"\])

    Step 2: LEARN QUERY STRUCTURE
    Call `get_open_targets_graphql_schema` with relevant categories to retrieve
    the schema subset needed for your query. Select categories that cover the
    data domains you need - BE INCLUSIVE (it's better to include extra categories
    than to miss required types).

    Example: For a query about drug mechanisms and safety:
    get_open_targets_graphql_schema(categories=\["drug-mechanisms", "drug-safety"\])

    Study the returned schema to understand available types, fields, and their
    relationships, then construct a GraphQL query that fetches the information
    the user needs.

    FALLBACK: If you encounter errors or need detailed information about specific
    types, use `get_type_dependencies` sparingly to explore type relationships.
    This tool provides exhaustive type dependency information but should only be
    used when category-based retrieval is insufficient.

    Step 3: CONSTRUCT AND EXECUTE QUERY
    Build GraphQL query using:
    - Standardized IDs from Step 1 (REQUIRED)
    - Query structure from Step 2
    - Follow the "COMMON MISTAKES TO AVOID" guidance in the schema output

    Call this tool with query_string and optional variables.

    REQUIRED IDENTIFIER FORMATS:
    - Targets/Genes: ENSEMBL IDs (e.g., "ENSG00000139618")
    - Diseases: EFO IDs (e.g., "EFO_0000305") or MONDO IDs (e.g., "MONDO_0007254")
    - Drugs: ChEMBL IDs (e.g., "CHEMBL1201583")
    - Variants: "chr_pos_ref_alt" format (e.g., "19_44908822_C_T") or rsIDs (e.g., "rs7412")
    - Studies: Study IDs (e.g., "GCST90002357")
    - Credible Sets: Study Locus IDs (e.g., "7d68cc9c70351c9dbd2a2c0c145e555d")

    Args:
    query_string (str): GraphQL query string starting with 'query' keyword.
    variables (UnionType\[dict\[str, Any\], None\]): Optional dict or JSON string with query variables.

    Returns:
    (QueryResult): GraphQL response with data field containing targets, diseases, drugs, variants, studies or error message.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'opentargets-platform',
    'query_open_targets_graphql',
    {'query_string': '...'},
)
```
