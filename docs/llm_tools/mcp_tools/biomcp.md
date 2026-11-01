# BioMCP (`biomcp`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **5** 个 DNA/基因相关工具（服务器共列出 7 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [genomoncology/biomcp @ `3b1d645`](https://github.com/genomoncology/biomcp/tree/3b1d6452ab63741562d2ef97495a97281b304729) |
| 许可证 | MIT |
| 启动命令 | `biomcp serve` |
| 传输方式 | stdio |
| 安装 | Install genomoncology/biomcp separately; its PyPI distribution is biomcp-cli, not the unrelated biomcp distribution. The installed executable is biomcp. |
| 副作用 | Contacts upstream biomedical APIs according to the selected tool; some sources require credentials. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`gene_cspec`](#gene_cspec) | `variant` | Retrieves ClinGen CSpec criteria-specification manifests, pages or attachment metadata for a given gene. |
| [`get`](#get) | `variant` | Fetches one biomedical record by ID for an entity type (gene, variant, drug, disease, trial, article, pathway, protein, etc.) with selectable sections. |
| [`search`](#search) | `variant` | Searches a biomedical entity type (gene, variant, gwas, pgx, article, trial, protein, author) with typed, bounded filters. |
| [`variant_erepo`](#variant_erepo) | `variant` | Retrieves ClinGen ERepo expert-panel variant classification assertions by ClinGen allele ID(s) or by gene. |
| [`variant_normalize_car`](#variant_normalize_car) | `variant` | Normalizes 1-50 RefSeq HGVS expressions through the ClinGen Allele Registry to canonical allele identifiers. |

## `gene_cspec` {#gene_cspec}

Retrieves ClinGen CSpec criteria-specification manifests, pages or attachment metadata for a given gene.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**clingen, acmg criteria, variant curation, vcep, classification rules

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `capture_id` | string \| null | 否 | 默认 `None` |
| `files` | boolean | 否 | 默认 `False` |
| `gene` | string | 是 |  |
| `limit` | integer | 否 | 默认 `25` |
| `offset` | integer | 否 | 默认 `0` |
| `version_iri` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Get ClinGen CSpec manifests, pages, or attachment metadata.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**clinvar, population frequency, gene record, variant record, civic, cosmic, gencc

??? note "上游说明（英文原文）"

    Get one biomedical record with typed inputs. Binary assets are CLI-only.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**gene search, variant search, gwas, pharmacogenomics, chromosome region, gene type

??? note "上游说明（英文原文）"

    Search a biomedical entity with typed, bounded inputs.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**clingen, erepo, caid, pathogenicity, expert panel, acmg

??? note "上游说明（英文原文）"

    Get ClinGen ERepo expert assertions by CAid or gene.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`inputs`
- **关键词：**hgvs, clingen allele registry, caid, normalization, refseq

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `inputs` | array of string | 是 |  |

??? note "上游说明（英文原文）"

    Normalize 1-50 RefSeq HGVS values with ClinGen Allele Registry.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biomcp',
    'variant_normalize_car',
    {'inputs': []},
)
```
