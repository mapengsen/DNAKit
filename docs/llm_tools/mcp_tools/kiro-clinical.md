# AWS Kiro Life Sciences: clinical (`kiro-clinical`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **3** 个 DNA/基因相关工具（服务器共列出 10 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| 许可证 | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| 启动命令 | `life-sciences-clinical` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-clinical component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-clinical console script. No uvx or automatic install is used. |
| 副作用 | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`omim_entry`](#omim_entry) | `association` | Retrieves a full OMIM entry by MIM number, including gene-phenotype and allelic variant data (requires OMIM_API_KEY). |
| [`omim_search`](#omim_search) | `association` | Searches OMIM gene and genetic disorder entries by keyword (requires OMIM_API_KEY). |
| [`pharmgkb_search`](#pharmgkb_search) | `variant` | Searches PharmGKB pharmacogenomics data (genes, variants, drugs) by keyword. |

## `omim_entry` {#omim_entry}

Retrieves a full OMIM entry by MIM number, including gene-phenotype and allelic variant data (requires OMIM_API_KEY).

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`mim_number`
- **关键词：**mendelian, mim number, genetic disorder, allelic variants, gene-phenotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `mim_number` | string | 是 | MIM number (e.g. '154700'). |

??? note "上游说明（英文原文）"

    Get OMIM entry by MIM number. Requires OMIM_API_KEY.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`query`
- **关键词：**mendelian, genetic disorder, inherited disease, mim, gene-phenotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search keyword (e.g. 'Marfan syndrome'). |
| `max_results` | integer | 否 | Maximum results to return. 默认 `10` |

??? note "上游说明（英文原文）"

    Search OMIM genetic disorder entries. Requires OMIM_API_KEY.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**pharmacogenomics, pgx, cyp2d6, drug response, haplotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search keyword (e.g. 'warfarin', 'CYP2D6'). |

??? note "上游说明（英文原文）"

    Search PharmGKB pharmacogenomics data by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-clinical',
    'pharmgkb_search',
    {'query': '...'},
)
```
