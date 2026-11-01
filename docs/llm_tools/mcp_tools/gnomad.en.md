# gnomAD MCP (`gnomad`)

DNAKit's MCP tool index keeps **11** DNA/gene tools of this server (12 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [koido/gnomad-mcp @ `63f09d2`](https://github.com/koido/gnomad-mcp/tree/63f09d21385751430c968ed770bb214f3067d786) |
| License | Apache-2.0 |
| Command | `python server.py` |
| Transports | stdio |
| Installation | Clone the pinned koido/gnomad-mcp checkout (not on PyPI, no console script) and install its dependencies (fastmcp, gql, aiohttp, httpx) in Python &gt;=3.13; set workdir to the checkout and executable_path to that Python. No environment variables or credentials are needed. |
| Effects | Sends GraphQL queries over HTTPS to the public gnomAD browser API at https://gnomad.broadinstitute.org/api (gene, variant, region and metadata lookups); the server tools do not write local files or modify remote data. |

| Tool | Category | Summary |
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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `dataset`, `reference_genome`, `variant_id`
- **Keywords:** clinvar, clinical significance, gnomad, variant id, pathogenicity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset` | string | yes |  |
| `reference_genome` | string | yes |  |
| `variant_id` | string | yes |  |

??? note "Upstream description"

    \[gnomAD API\] Retrieve ClinVar variant info (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r4/gnomad_r3/gnomad_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variant_id (str): Variant ID
    Returns:
    dict: ClinVar info

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `reference_genome`, `variantId`
- **Keywords:** cnv, copy number, duplication, deletion, gnomad

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `reference_genome` | string | yes |  |
| `variantId` | string | yes |  |
| `dataset` | string | no | Default `gnomad_cnv_r4` |

??? note "Upstream description"

    \[gnomAD API\] Retrieve copy number variant info (v4 only)
    Args:
    reference_genome (str): Reference genome (GRCh38)
    variantId (str): CNV ID (e.g. 18714__DUP)
    Returns:
    dict: CNV info
    Note:
    Not supported in v2/v3.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** gnomad, gene constraint, ensembl id, gene symbol, loeuf

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string \| null | no | Default `None` |
| `gene_symbol` | string \| null | no | Default `None` |
| `dataset` | string | no | Default `gnomad_r4` |
| `reference_genome` | string | no | Default `GRCh38` |

??? note "Upstream description"

    \[gnomAD API\] Retrieve gene information (v4 only)
    Args:
    gene_id (str, optional): Ensembl gene ID (e.g. ENSG00000139618)
    gene_symbol (str, optional): Gene symbol (e.g. PCSK9)
    Returns:
    dict: gene info
    Note:
    Not supported in v2/v3. Use search_for_genes instead.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `reference_genome`, `variant_id`
- **Keywords:** mitochondrial, mtdna, heteroplasmy, gnomad, allele frequency

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `reference_genome` | string | yes |  |
| `variant_id` | string | yes |  |
| `dataset` | string | no | Default `gnomad_r4` |

??? note "Upstream description"

    \[gnomAD API\] Retrieve mitochondrial variant info (v4 only)
    Args:
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variant_id (str): Mitochondrial variant ID (e.g. M-8602-T-C). Build must be GRCh38.
    Returns:
    dict: mito variant info
    Note:
    Not supported in v2/v3.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `reference_genome`, `chrom`, `start`, `stop`
- **Keywords:** genomic region, gnomad, coordinates, population variants, allele frequency

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `reference_genome` | string | yes |  |
| `chrom` | string | yes |  |
| `start` | integer | yes |  |
| `stop` | integer | yes |  |
| `dataset` | string | no | Default `gnomad_r4` |

??? note "Upstream description"

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

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `reference_genome`, `id`
- **Keywords:** str, short tandem repeat, repeat expansion, gnomad

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `reference_genome` | string | yes |  |
| `id` | string | yes |  |
| `dataset` | string | no | Default `gnomad_r4` |

??? note "Upstream description"

    \[gnomAD API\] Retrieve STR info (v4 only)
    Args:
    reference_genome (str): Reference genome (GRCh38)
    id (str): STR ID (e.g. ATXN1)
    Returns:
    dict: STR info
    Note:
    Not supported in v2/v3.

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `dataset`, `reference_genome`, `variantId`
- **Keywords:** structural variant, sv, gnomad sv, population frequency

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset` | string | yes |  |
| `reference_genome` | string | yes |  |
| `variantId` | string | yes |  |

??? note "Upstream description"

    \[gnomAD API\] Retrieve structural variant info (v2/v4)
    Args:
    dataset (str): SV dataset ID (gnomad_sv_r4/gnomad_sv_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variantId (str): Structural variant ID
    Returns:
    dict: SV info

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `dataset`, `reference_genome`, `variantId`
- **Keywords:** allele frequency, population frequency, gnomad, snv, indel

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset` | string | yes |  |
| `reference_genome` | string | yes |  |
| `variantId` | string | yes |  |

??? note "Upstream description"

    \[gnomAD API\] Retrieve variant information (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r4/gnomad_r3/gnomad_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    variantId (str): Variant ID (e.g. 1-55051215-G-GA)
    Returns:
    dict: variant info

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `reference_genome`
- **Keywords:** liftover, grch37, grch38, coordinate conversion, gnomad

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `reference_genome` | string | yes |  |
| `source_variant_id` | string \| null | no | Default `None` |
| `liftover_variant_id` | string \| null | no | Default `None` |
| `dataset` | string | no | Default `gnomad_r2_1` |

??? note "Upstream description"

    \[gnomAD API\] Retrieve liftover info (v2 only)
    Args:
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    source_variant_id (str, optional): Source variant ID (e.g. 12-112241766-G-A on GRCh37)
    liftover_variant_id (str, optional): Lifted over variant ID (e.g. 12-111803962-G-A on GRCh38)
    Returns:
    dict: liftover info
    Note:
    Not supported in v3/v4.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `dataset`, `reference_genome`, `query`
- **Keywords:** gnomad, gene search, gene symbol, ensembl id

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset` | string | yes |  |
| `reference_genome` | string | yes |  |
| `query` | string | yes |  |

??? note "Upstream description"

    \[gnomAD API\] Search for genes (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r3/gnomad_r2_1)
    reference_genome (str, optional): Reference genome (GRCh37 or GRCh38)
    query (str): Search string
    Returns:
    dict: search results

Call it through DNAKit once the server is installed:

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

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `dataset`, `reference_genome`, `query`
- **Keywords:** gnomad, variant search, variant id, rsid

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset` | string | yes |  |
| `reference_genome` | string | yes |  |
| `query` | string | yes |  |

??? note "Upstream description"

    \[gnomAD API\] Search for variants (v2/v3/v4)
    Args:
    dataset (str): gnomAD dataset ID (gnomad_r4/gnomad_r3/gnomad_r2_1)
    reference_genome (str): Reference genome (GRCh37 or GRCh38)
    query (str): Search string (variant_id)
    Returns:
    dict: search results (variant_id)

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gnomad',
    'search_for_variants',
    {'dataset': '...', 'reference_genome': '...', 'query': '...'},
)
```
