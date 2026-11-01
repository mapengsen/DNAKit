# gnomAD MCP (`gnomad`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **11** 个 DNA/基因相关工具（服务器共列出 12 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [koido/gnomad-mcp @ `63f09d2`](https://github.com/koido/gnomad-mcp/tree/63f09d21385751430c968ed770bb214f3067d786) |
| 许可证 | Apache-2.0 |
| 启动命令 | `python server.py` |
| 传输方式 | stdio |
| 安装 | Clone the pinned koido/gnomad-mcp checkout (not on PyPI, no console script) and install its dependencies (fastmcp, gql, aiohttp, httpx) in Python &gt;=3.13; set workdir to the checkout and executable_path to that Python. No environment variables or credentials are needed. |
| 副作用 | Sends GraphQL queries over HTTPS to the public gnomAD browser API at https://gnomad.broadinstitute.org/api (gene, variant, region and metadata lookups); the server tools do not write local files or modify remote data. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`get_clinvar_variant_info`](#get_clinvar_variant_info) | `variant` | Retrieves ClinVar information for a variant ID from a gnomAD dataset (v2/v3/v4) on GRCh37 or GRCh38. |
| [`get_copy_number_variant_info`](#get_copy_number_variant_info) | `variant` | Retrieves copy number variant (CNV) information by CNV ID from gnomAD v4 (GRCh38). |
| [`get_gene_info`](#get_gene_info) | `genome_annotation` | Retrieves gnomAD v4 gene information by Ensembl gene ID or gene symbol. |
| [`get_mitochondrial_variant_info`](#get_mitochondrial_variant_info) | `variant` | Retrieves mitochondrial variant information (e.g. M-8602-T-C) from gnomAD v4 on GRCh38. |
| [`get_region_info`](#get_region_info) | `variant` | Retrieves gnomAD v4 information for a GRCh38 chromosomal region given chromosome, start and stop. |
| [`get_str_info`](#get_str_info) | `variant` | Retrieves short tandem repeat (STR) information, e.g. ATXN1, from gnomAD v4. |
| [`get_structural_variant_info`](#get_structural_variant_info) | `variant` | Retrieves structural variant information by SV ID from gnomAD SV datasets (v2/v4). |
| [`get_variant_info`](#get_variant_info) | `variant` | Retrieves gnomAD variant information for a variant ID (e.g. 1-55051215-G-GA) from v2/v3/v4 datasets. |
| [`get_variant_liftover`](#get_variant_liftover) | `variant` | Retrieves gnomAD v2 liftover mapping for a variant between GRCh37 and GRCh38. |
| [`search_for_genes`](#search_for_genes) | `genome_annotation` | Searches gnomAD genes by query string in v2/v3 datasets. |
| [`search_for_variants`](#search_for_variants) | `variant` | Searches gnomAD variants by variant ID query in v2/v3/v4 datasets. |

## `get_clinvar_variant_info` {#get_clinvar_variant_info}

Retrieves ClinVar information for a variant ID from a gnomAD dataset (v2/v3/v4) on GRCh37 or GRCh38.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`dataset`, `reference_genome`, `variant_id`
- **关键词：**clinvar, clinical significance, gnomad, variant id, pathogenicity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset` | string | 是 |  |
| `reference_genome` | string | 是 |  |
| `variant_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve ClinVar variant info (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r4/gnomad_r3/gnomad_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variant_id (str): Variant ID
    Returns:
    dict: ClinVar info

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_clinvar_variant_info',
    {'dataset': '...', 'reference_genome': '...', 'variant_id': '...'},
)
```

## `get_copy_number_variant_info` {#get_copy_number_variant_info}

Retrieves copy number variant (CNV) information by CNV ID from gnomAD v4 (GRCh38).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`reference_genome`, `variantId`
- **关键词：**cnv, copy number, duplication, deletion, gnomad

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `reference_genome` | string | 是 |  |
| `variantId` | string | 是 |  |
| `dataset` | string | 否 | 默认 `gnomad_cnv_r4` |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve copy number variant info (v4 only)
    Args:
    reference_genome (str): Reference genome (GRCh38)
    variantId (str): CNV ID (e.g. 18714__DUP)
    Returns:
    dict: CNV info
    Note:
    Not supported in v2/v3.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_copy_number_variant_info',
    {'reference_genome': '...', 'variantId': '...'},
)
```

## `get_gene_info` {#get_gene_info}

Retrieves gnomAD v4 gene information by Ensembl gene ID or gene symbol.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**gnomad, gene constraint, ensembl id, gene symbol, loeuf

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string \| null | 否 | 默认 `None` |
| `gene_symbol` | string \| null | 否 | 默认 `None` |
| `dataset` | string | 否 | 默认 `gnomad_r4` |
| `reference_genome` | string | 否 | 默认 `GRCh38` |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve gene information (v4 only)
    Args:
    gene_id (str, optional): Ensembl gene ID (e.g. ENSG00000139618)
    gene_symbol (str, optional): Gene symbol (e.g. PCSK9)
    Returns:
    dict: gene info
    Note:
    Not supported in v2/v3. Use search_for_genes instead.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_gene_info',
    {},
)
```

## `get_mitochondrial_variant_info` {#get_mitochondrial_variant_info}

Retrieves mitochondrial variant information (e.g. M-8602-T-C) from gnomAD v4 on GRCh38.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`reference_genome`, `variant_id`
- **关键词：**mitochondrial, mtdna, heteroplasmy, gnomad, allele frequency

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `reference_genome` | string | 是 |  |
| `variant_id` | string | 是 |  |
| `dataset` | string | 否 | 默认 `gnomad_r4` |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve mitochondrial variant info (v4 only)
    Args:
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variant_id (str): Mitochondrial variant ID (e.g. M-8602-T-C). Build must be GRCh38.
    Returns:
    dict: mito variant info
    Note:
    Not supported in v2/v3.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_mitochondrial_variant_info',
    {'reference_genome': '...', 'variant_id': '...'},
)
```

## `get_region_info` {#get_region_info}

Retrieves gnomAD v4 information for a GRCh38 chromosomal region given chromosome, start and stop.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`reference_genome`, `chrom`, `start`, `stop`
- **关键词：**genomic region, gnomad, coordinates, population variants, allele frequency

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `reference_genome` | string | 是 |  |
| `chrom` | string | 是 |  |
| `start` | integer | 是 |  |
| `stop` | integer | 是 |  |
| `dataset` | string | 否 | 默认 `gnomad_r4` |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve region information (v4 only)
    Args:
    reference_genome (str): Reference genome (GRCh38)
    chrom (str): Chromosome
    start (int): Start position. Build must be GRCh38.
    stop (int): End position. Build must be GRCh38.
    Returns:
    dict: region info
    Note:
    Not supported in v2/v3.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_region_info',
    {'reference_genome': '...', 'chrom': '...', 'start': 0, 'stop': 0},
)
```

## `get_str_info` {#get_str_info}

Retrieves short tandem repeat (STR) information, e.g. ATXN1, from gnomAD v4.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`reference_genome`, `id`
- **关键词：**str, short tandem repeat, repeat expansion, gnomad

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `reference_genome` | string | 是 |  |
| `id` | string | 是 |  |
| `dataset` | string | 否 | 默认 `gnomad_r4` |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve STR info (v4 only)
    Args:
    reference_genome (str): Reference genome (GRCh38)
    id (str): STR ID (e.g. ATXN1)
    Returns:
    dict: STR info
    Note:
    Not supported in v2/v3.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_str_info',
    {'reference_genome': '...', 'id': '...'},
)
```

## `get_structural_variant_info` {#get_structural_variant_info}

Retrieves structural variant information by SV ID from gnomAD SV datasets (v2/v4).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`dataset`, `reference_genome`, `variantId`
- **关键词：**structural variant, sv, gnomad sv, population frequency

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset` | string | 是 |  |
| `reference_genome` | string | 是 |  |
| `variantId` | string | 是 |  |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve structural variant info (v2/v4)
    Args:
    dataset (str): SV dataset ID (gnomad_sv_r4/gnomad_sv_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variantId (str): Structural variant ID
    Returns:
    dict: SV info

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_structural_variant_info',
    {'dataset': '...', 'reference_genome': '...', 'variantId': '...'},
)
```

## `get_variant_info` {#get_variant_info}

Retrieves gnomAD variant information for a variant ID (e.g. 1-55051215-G-GA) from v2/v3/v4 datasets.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`dataset`, `reference_genome`, `variantId`
- **关键词：**allele frequency, population frequency, gnomad, snv, indel

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset` | string | 是 |  |
| `reference_genome` | string | 是 |  |
| `variantId` | string | 是 |  |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve variant information (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r4/gnomad_r3/gnomad_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variantId (str): Variant ID (e.g. 1-55051215-G-GA)
    Returns:
    dict: variant info

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_variant_info',
    {'dataset': '...', 'reference_genome': '...', 'variantId': '...'},
)
```

## `get_variant_liftover` {#get_variant_liftover}

Retrieves gnomAD v2 liftover mapping for a variant between GRCh37 and GRCh38.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`reference_genome`
- **关键词：**liftover, grch37, grch38, coordinate conversion, gnomad

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `reference_genome` | string | 是 |  |
| `source_variant_id` | string \| null | 否 | 默认 `None` |
| `liftover_variant_id` | string \| null | 否 | 默认 `None` |
| `dataset` | string | 否 | 默认 `gnomad_r2_1` |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Retrieve liftover info (v2 only)
    Args:
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    source_variant_id (str, optional): Source variant ID (e.g. 12-112241766-G-A on GRCh37)
    liftover_variant_id (str, optional): Lifted over variant ID (e.g. 12-111803962-G-A on GRCh38)
    Returns:
    dict: liftover info
    Note:
    Not supported in v3/v4.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'get_variant_liftover',
    {'reference_genome': '...'},
)
```

## `search_for_genes` {#search_for_genes}

Searches gnomAD genes by query string in v2/v3 datasets.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`dataset`, `reference_genome`, `query`
- **关键词：**gnomad, gene search, gene symbol, ensembl id

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset` | string | 是 |  |
| `reference_genome` | string | 是 |  |
| `query` | string | 是 |  |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Search for genes (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r3/gnomad_r2_1)
    reference_genome (str, optional): Reference genome (GRCh37 or GRCh38)
    query (str): Search string
    Returns:
    dict: search results

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'search_for_genes',
    {'dataset': '...', 'reference_genome': '...', 'query': '...'},
)
```

## `search_for_variants` {#search_for_variants}

Searches gnomAD variants by variant ID query in v2/v3/v4 datasets.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`dataset`, `reference_genome`, `query`
- **关键词：**gnomad, variant search, variant id, rsid

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset` | string | 是 |  |
| `reference_genome` | string | 是 |  |
| `query` | string | 是 |  |

??? note "上游说明（英文原文）"

    \[gnomAD API\] Search for variants (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r4/gnomad_r3/gnomad_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    query (str): Search string (variant_id)
    Returns:
    dict: search results (variant_id)

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'search_for_variants',
    {'dataset': '...', 'reference_genome': '...', 'query': '...'},
)
```
