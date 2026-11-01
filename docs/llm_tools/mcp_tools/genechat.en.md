# GeneChat MCP (`genechat`)

DNAKit's MCP tool index keeps **10** DNA/gene tools of this server (10 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [natecostello/genechat-mcp @ `6f40f93`](https://github.com/natecostello/genechat-mcp/tree/6f40f93c1ce6a0f914a4db34c759f733264e6cdb) |
| License | MIT |
| Command | `genechat serve` |
| Transports | stdio |
| Installation | Install and initialize GeneChat and its local genome/annotation databases separately. GENECHAT_CONFIG optionally selects an existing TOML configuration. Pin 'mcp&lt;2' (upstream declares mcp&gt;=1.0.0; under mcp 2.x the server fails to import), and configure at least one genome: serve exits without one. |
| Effects | Reads configured local genomes; returned tool results can include individual genotypes. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`calculate_prs`](#calculate_prs) | `variant` | Calculates a polygenic risk score for a trait or PGS Catalog ID on a registered personal genome, optionally comparing two genomes. |
| [`genome_summary`](#genome_summary) | `variant` | Returns a high-level summary of a registered personal genome: variant counts, ClinVar annotation summary and pharmacogenomics overview. |
| [`list_genomes`](#list_genomes) | `variant` | Lists registered personal genomes with labels, VCF paths and annotation state. |
| [`query_clinvar`](#query_clinvar) | `variant` | Finds variants in a registered genome with a given ClinVar significance, optionally filtered by gene or condition. |
| [`query_gene`](#query_gene) | `variant` | Queries variants in one gene of a registered genome, filtered by SnpEff impact level (HIGH, MODERATE, LOW, MODIFIER). |
| [`query_genes`](#query_genes) | `variant` | Queries variants across a comma-separated list of genes in a registered genome in one VCF pass, with impact filtering. |
| [`query_gwas`](#query_gwas) | `association` | Searches GWAS Catalog findings by trait, gene or rsID, ordered by p-value, with optional check against a registered genome. |
| [`query_pgx`](#query_pgx) | `variant` | Looks up pharmacogenomic information for a drug or gene and reports relevant genotypes from a registered genome. |
| [`query_variant`](#query_variant) | `variant` | Looks up a single variant by rsID or genomic position, returning genotype, functional annotation, ClinVar significance and population frequency. |
| [`query_variants`](#query_variants) | `variant` | Looks up multiple variants by comma-separated rsIDs in one VCF scan of a registered genome. |

## `calculate_prs` {#calculate_prs}

Calculates a polygenic risk score for a trait or PGS Catalog ID on a registered personal genome, optionally comparing two genomes.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** polygenic risk score, pgs catalog, personal genome, vcf, trait risk

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `trait` | string \| null | no | Default `None` |
| `prs_id` | string \| null | no | Default `None` |
| `genome` | string \| null | no | Default `None` |
| `genome2` | string \| null | no | Default `None` |

??? note "Upstream description"

    Calculate a polygenic risk score (PRS) for a trait.

    PRS aggregates the effects of many common genetic variants to estimate
    genetic predisposition to a trait. Provide either a trait name or a PGS Catalog ID.

    Optional: 'genome' selects which registered genome to query (default: primary genome).
    'genome2' calculates the score for a second genome for comparison.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'calculate_prs',
    {},
)
```

## `genome_summary` {#genome_summary}

Returns a high-level summary of a registered personal genome: variant counts, ClinVar annotation summary and pharmacogenomics overview.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** personal genome, vcf, variant counts, clinvar, pharmacogenomics

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string \| null | no | Default `None` |
| `genome2` | string \| null | no | Default `None` |

??? note "Upstream description"

    Get a high-level summary of your genome data.

    Returns variant counts, ClinVar annotation summary, and pharmacogenomics overview.
    Use this as a starting point when a user wants a general overview of their genome.

    Optional: 'genome' selects which registered genome to query (default: primary genome).
    'genome2' produces a side-by-side summary of a second genome.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'genome_summary',
    {},
)
```

## `list_genomes` {#list_genomes}

Lists registered personal genomes with labels, VCF paths and annotation state.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** personal genome, vcf, registered genomes, annotation state

??? note "Upstream description"

    List all registered genomes with their labels and basic info.

    Use this tool FIRST to discover which genomes are available before
    querying variants, genes, or other genomic data. Returns genome labels,
    VCF paths, and annotation state for each registered genome.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'list_genomes',
    {},
)
```

## `query_clinvar` {#query_clinvar}

Finds variants in a registered genome with a given ClinVar significance, optionally filtered by gene or condition.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `significance`
- **Keywords:** clinvar, pathogenic, personal genome, clinical significance, disease risk

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `significance` | string | yes |  |
| `gene` | string \| null | no | Default `None` |
| `condition` | string \| null | no | Default `None` |
| `max_results` | integer | no | Default `50` |
| `genome` | string \| null | no | Default `None` |
| `genome2` | string \| null | no | Default `None` |

??? note "Upstream description"

    Search your genome for variants with a specific ClinVar clinical significance.

    Use this when a user asks about pathogenic variants, disease risk, or clinical findings.
    Significance values: Pathogenic, Likely_pathogenic, Benign, Likely_benign,
    Uncertain_significance, risk_factor, drug_response, protective.

    Optional: 'genome' selects which registered genome to query (default: primary genome).
    'genome2' queries a second genome for side-by-side comparison.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'query_clinvar',
    {'significance': '...'},
)
```

## `query_gene` {#query_gene}

Queries variants in one gene of a registered genome, filtered by SnpEff impact level (HIGH, MODERATE, LOW, MODIFIER).

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene`
- **Keywords:** snpeff, impact, personal genome, vcf, gene variants

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes |  |
| `impact_filter` | string | no | Default `HIGH,MODERATE` |
| `max_results` | integer | no | Default `50` |
| `smart_filter` | boolean | no | Default `True` |
| `genome` | string \| null | no | Default `None` |
| `genome2` | string \| null | no | Default `None` |

??? note "Upstream description"

    Query variants in a specific gene from your genome.

    Use this when a user asks about variants in a gene (e.g. "What variants do I have in BRCA1?").
    Filters by SnpEff impact level (HIGH, MODERATE, LOW, MODIFIER) by default.
    If functional annotation is not available, all variants are included regardless
    of impact filter.

    Optional: 'genome' selects which registered genome to query (default: primary genome).
    'genome2' queries a second genome for side-by-side comparison.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'query_gene',
    {'gene': '...'},
)
```

## `query_genes` {#query_genes}

Queries variants across a comma-separated list of genes in a registered genome in one VCF pass, with impact filtering.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `genes`
- **Keywords:** gene panel, snpeff, personal genome, vcf, batch

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | string | yes |  |
| `impact_filter` | string | no | Default `HIGH,MODERATE` |
| `max_results_per_gene` | integer | no | Default `20` |
| `genome` | string \| null | no | Default `None` |

??? note "Upstream description"

    Query variants across multiple genes in a single call.

    Accepts a comma-separated list of gene symbols (e.g. "BRCA1,BRCA2,TP53").
    Much more efficient than calling query_gene repeatedly — opens the VCF once.
    Use this when investigating a pathway, gene panel, or related group of genes.
    For example: "Check APOB,LDLR,PCSK9 for cardiovascular risk" or
    "Look at CYP2D6,CYP2C19,CYP2C9 for drug metabolism".

    Optional: 'genome' selects which registered genome to query (default: primary genome).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'query_genes',
    {'genes': '...'},
)
```

## `query_gwas` {#query_gwas}

Searches GWAS Catalog findings by trait, gene or rsID, ordered by p-value, with optional check against a registered genome.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** none
- **Keywords:** gwas catalog, trait, rsid, p-value, genetic association

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `trait` | string \| null | no | Default `None` |
| `gene` | string \| null | no | Default `None` |
| `rsid` | string \| null | no | Default `None` |
| `max_results` | integer | no | Default `30` |
| `deduplicate` | boolean | no | Default `True` |
| `check_vcf` | boolean | no | Default `False` |
| `genome` | string \| null | no | Default `None` |

??? note "Upstream description"

    Search the GWAS Catalog for genome-wide association study findings.

    Use this when a user asks about genetic associations with a disease or trait,
    or to discover what conditions are linked to a gene or variant.
    Provide at least one of: trait (e.g. "type 2 diabetes"), gene (e.g. "FTO"),
    or rsid (e.g. "rs9939609").
    Results are ordered by statistical significance (lowest p-value first).

    Optional: 'genome' selects which registered genome to use for VCF cross-reference
    when check_vcf=true (default: primary genome).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'query_gwas',
    {},
)
```

## `query_pgx` {#query_pgx}

Looks up pharmacogenomic information for a drug or gene and reports relevant genotypes from a registered genome.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** pharmacogenomics, cyp2d6, drug metabolism, genotype, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `drug` | string \| null | no | Default `None` |
| `gene` | string \| null | no | Default `None` |
| `include_all_variants` | boolean | no | Default `False` |
| `genome` | string \| null | no | Default `None` |
| `genome2` | string \| null | no | Default `None` |

??? note "Upstream description"

    Look up pharmacogenomic information for a drug or gene.

    Use this when a user asks about drug interactions, medication safety, or
    how their genetics affect drug metabolism. Provide either a drug name
    (e.g. "simvastatin") or a gene symbol (e.g. "CYP2D6").

    Optional: 'genome' selects which registered genome to query (default: primary genome).
    'genome2' queries a second genome for side-by-side comparison.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'query_pgx',
    {},
)
```

## `query_variant` {#query_variant}

Looks up a single variant by rsID or genomic position, returning genotype, functional annotation, ClinVar significance and population frequency.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** rsid, snp, genotype, clinvar, allele frequency, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string \| null | no | Default `None` |
| `position` | string \| null | no | Default `None` |
| `genome` | string \| null | no | Default `None` |
| `genome2` | string \| null | no | Default `None` |

??? note "Upstream description"

    Look up a specific genetic variant by rsID (e.g. rs4149056) or genomic position (e.g. chr22:42127941).

    Use this when a user asks about a specific SNP, rsID, or genomic coordinate.
    Returns genotype, functional annotation, ClinVar significance, and population frequency.

    Optional: 'genome' selects which registered genome to query (default: primary genome).
    'genome2' queries a second genome for side-by-side comparison (e.g. carrier screening for couples).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'query_variant',
    {},
)
```

## `query_variants` {#query_variants}

Looks up multiple variants by comma-separated rsIDs in one VCF scan of a registered genome.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `rsids`
- **Keywords:** rsid, snp, genotype, batch, personal genome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsids` | string | yes |  |
| `genome` | string \| null | no | Default `None` |
| `genome2` | string \| null | no | Default `None` |

??? note "Upstream description"

    Look up multiple genetic variants by rsID in a single query.

    Accepts a comma-separated list of rsIDs (e.g. "rs4149056,rs1801133,rs762551").
    Much more efficient than calling query_variant repeatedly — scans the VCF once
    for all requested variants.
    Use this when a user asks about multiple specific SNPs at once, or when you
    need to check several variants for a comprehensive answer.

    Optional: 'genome' selects which registered genome to query (default: primary genome).
    'genome2' queries a second genome for side-by-side comparison.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'genechat',
    'query_variants',
    {'rsids': '...'},
)
```
