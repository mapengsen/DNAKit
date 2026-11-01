# AWS Kiro Life Sciences: genomics (`kiro-genomics`)

DNAKit's MCP tool index keeps **17** DNA/gene tools of this server (18 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| License | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| Command | `life-sciences-genomics` |
| Transports | stdio |
| Installation | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-genomics component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-genomics console script. No uvx or automatic install is used. |
| Effects | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`clinvar_get_variation`](#clinvar_get_variation) | `variant` | Fetches the full ClinVar variation record for a given ClinVar variation ID. |
| [`clinvar_search`](#clinvar_search) | `variant` | Searches ClinVar by rsID, HGVS notation or gene name and returns matching variation records. |
| [`cosmic_search`](#cosmic_search) | `variant` | Searches the COSMIC database for somatic cancer mutations in a given gene symbol (requires COSMIC_API_KEY). |
| [`dbsnp_lookup`](#dbsnp_lookup) | `variant` | Retrieves the dbSNP record for a variant given its rsID. |
| [`ddbj_fetch_sequence`](#ddbj_fetch_sequence) | `sequence` | Fetches a sequence record from DDBJ by accession number in FASTA format. |
| [`ddbj_search`](#ddbj_search) | `sequence` | Searches DDBJ (DNA Data Bank of Japan) sequence entries by keyword. |
| [`encode_search`](#encode_search) | `sequencing_data` | Searches ENCODE experiments filtered by biosample, assay and/or target label. |
| [`ensembl_gene_lookup`](#ensembl_gene_lookup) | `genome_annotation` | Looks up an Ensembl gene record (ID, location, biotype) by gene symbol and species via the Ensembl REST API. |
| [`ensembl_sequence`](#ensembl_sequence) | `sequence` | Retrieves the nucleotide sequence for a genomic region of a species via the Ensembl REST API. |
| [`ensembl_variants`](#ensembl_variants) | `variant` | Lists known variants overlapping a genomic region (e.g. 7:140424943-140624564) for a species via Ensembl REST. |
| [`geo_get_dataset`](#geo_get_dataset) | `sequencing_data` | Retrieves a GEO dataset/series record by accession (e.g. GSE12345). |
| [`geo_search`](#geo_search) | `sequencing_data` | Searches NCBI GEO datasets by keyword and returns matching dataset summaries. |
| [`gnomad_variant`](#gnomad_variant) | `variant` | Retrieves gnomAD allele frequency data for a variant given in chrom-pos-ref-alt form. |
| [`ncbi_fetch_sequence`](#ncbi_fetch_sequence) | `sequence` | Fetches a nucleotide sequence record by accession via NCBI Entrez efetch in FASTA or GenBank format. |
| [`ncbi_search`](#ncbi_search) | `genome_annotation` | Runs an NCBI Entrez esearch over a chosen database (e.g. gene, nucleotide, protein) and returns matching IDs. |
| [`sra_search`](#sra_search) | `sequencing_data` | Searches NCBI SRA runs by keyword or BioProject accession. |
| [`thousand_genomes_frequency`](#thousand_genomes_frequency) | `variant` | Retrieves 1000 Genomes population allele frequencies for a variant (e.g. rsID) via Ensembl. |

## `clinvar_get_variation` {#clinvar_get_variation}

Fetches the full ClinVar variation record for a given ClinVar variation ID.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variation_id`
- **Keywords:** clinical significance, pathogenicity, variation id, germline variant, hgvs

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variation_id` | string | yes | ClinVar variation ID. |

??? note "Upstream description"

    Get a full ClinVar variation record by variation ID.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'clinvar_get_variation',
    {'variation_id': '...'},
)
```

## `clinvar_search` {#clinvar_search}

Searches ClinVar by rsID, HGVS notation or gene name and returns matching variation records.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query`
- **Keywords:** rsid, hgvs, pathogenic, clinical significance, gene variants

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search term (rsID, HGVS, or gene name). |
| `max_results` | integer | no | Maximum results. Default `10` |

??? note "Upstream description"

    Search ClinVar by rsID, HGVS notation, or gene name.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'clinvar_search',
    {'query': '...'},
)
```

## `cosmic_search` {#cosmic_search}

Searches the COSMIC database for somatic cancer mutations in a given gene symbol (requires COSMIC_API_KEY).

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene`
- **Keywords:** somatic mutation, cancer, tumor, mutation catalogue, oncogene

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol (e.g. 'BRAF'). |

??? note "Upstream description"

    Search COSMIC somatic mutations by gene name. Requires COSMIC_API_KEY.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'cosmic_search',
    {'gene': '...'},
)
```

## `dbsnp_lookup` {#dbsnp_lookup}

Retrieves the dbSNP record for a variant given its rsID.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `rsid`
- **Keywords:** rsid, snp, single nucleotide polymorphism, refsnp, allele

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | yes | dbSNP rsID (e.g. 'rs328'). |

??? note "Upstream description"

    Lookup a variant by rsID from dbSNP.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'dbsnp_lookup',
    {'rsid': '...'},
)
```

## `ddbj_fetch_sequence` {#ddbj_fetch_sequence}

Fetches a sequence record from DDBJ by accession number in FASTA format.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** fasta, nucleotide, accession, insdc, dna sequence

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | DDBJ accession number. |

??? note "Upstream description"

    Fetch a sequence from DDBJ in FASTA format.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ddbj_fetch_sequence',
    {'accession': '...'},
)
```

## `ddbj_search` {#ddbj_search}

Searches DDBJ (DNA Data Bank of Japan) sequence entries by keyword.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `query`
- **Keywords:** nucleotide, insdc, accession, sequence records, dna databank

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search keyword. |
| `max_results` | integer | no | Maximum results. Default `10` |

??? note "Upstream description"

    Search DDBJ by keyword.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ddbj_search',
    {'query': '...'},
)
```

## `encode_search` {#encode_search}

Searches ENCODE experiments filtered by biosample, assay and/or target label.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** chip-seq, atac-seq, biosample, assay, functional genomics

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `biosample` | string | no | Biosample term name. Default `` |
| `assay` | string | no | Assay title. Default `` |
| `target` | string | no | Target label. Default `` |

??? note "Upstream description"

    Search ENCODE experiments by biosample, assay, and/or target.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'encode_search',
    {},
)
```

## `ensembl_gene_lookup` {#ensembl_gene_lookup}

Looks up an Ensembl gene record (ID, location, biotype) by gene symbol and species via the Ensembl REST API.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `symbol`
- **Keywords:** gene symbol, ensembl id, genomic location, biotype, species

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `symbol` | string | yes | Gene symbol (e.g. 'BRCA2'). |
| `species` | string | no | Species name. Default `homo_sapiens` |

??? note "Upstream description"

    Lookup a gene by symbol and species via Ensembl REST.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ensembl_gene_lookup',
    {'symbol': '...'},
)
```

## `ensembl_sequence` {#ensembl_sequence}

Retrieves the nucleotide sequence for a genomic region of a species via the Ensembl REST API.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `region`
- **Keywords:** genomic region, dna sequence, reference genome, coordinates, fasta

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes | Genomic region. |
| `species` | string | no | Species name. Default `homo_sapiens` |

??? note "Upstream description"

    Get nucleotide sequence for a genomic region via Ensembl REST.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ensembl_sequence',
    {'region': '...'},
)
```

## `ensembl_variants` {#ensembl_variants}

Lists known variants overlapping a genomic region (e.g. 7:140424943-140624564) for a species via Ensembl REST.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `region`
- **Keywords:** genomic region, snp, overlap, variation, coordinates

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes | Genomic region (e.g. '7:140424943-140624564'). |
| `species` | string | no | Species name. Default `homo_sapiens` |

??? note "Upstream description"

    Get variants in a genomic region via Ensembl REST.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ensembl_variants',
    {'region': '...'},
)
```

## `geo_get_dataset` {#geo_get_dataset}

Retrieves a GEO dataset/series record by accession (e.g. GSE12345).

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** gse, expression dataset, series, microarray, rna-seq

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | GEO accession number. |

??? note "Upstream description"

    Get a GEO dataset by accession (e.g. GSE12345).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'geo_get_dataset',
    {'accession': '...'},
)
```

## `geo_search` {#geo_search}

Searches NCBI GEO datasets by keyword and returns matching dataset summaries.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `term`
- **Keywords:** gse, expression data, rna-seq, microarray, datasets

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `term` | string | yes | Search keyword. |
| `max_results` | integer | no | Maximum results. Default `10` |

??? note "Upstream description"

    Search GEO datasets by keyword.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'geo_search',
    {'term': '...'},
)
```

## `gnomad_variant` {#gnomad_variant}

Retrieves gnomAD allele frequency data for a variant given in chrom-pos-ref-alt form.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant_id`
- **Keywords:** allele frequency, population frequency, gnomad, chrom-pos-ref-alt, maf

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes | Variant identifier (e.g. '1-55516888-G-A'). |

??? note "Upstream description"

    Get allele frequency data for a variant from gnomAD.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'gnomad_variant',
    {'variant_id': '...'},
)
```

## `ncbi_fetch_sequence` {#ncbi_fetch_sequence}

Fetches a nucleotide sequence record by accession via NCBI Entrez efetch in FASTA or GenBank format.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** efetch, genbank, fasta, nucleotide, accession

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | NCBI accession number. |
| `format` | string | no | Output format ('fasta' or 'gb'). Default `fasta` |

??? note "Upstream description"

    Fetch a nucleotide sequence by accession via Entrez efetch.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ncbi_fetch_sequence',
    {'accession': '...'},
)
```

## `ncbi_search` {#ncbi_search}

Runs an NCBI Entrez esearch over a chosen database (e.g. gene, nucleotide, protein) and returns matching IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `database`, `term`
- **Keywords:** entrez, esearch, gene database, nucleotide, uid

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `database` | string | yes | NCBI database name (e.g. 'gene', 'nucleotide', 'protein'). |
| `term` | string | yes | Search query string. |
| `max_results` | integer | no | Maximum results to return. Default `20` |

??? note "Upstream description"

    Search any NCBI database via Entrez esearch.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ncbi_search',
    {'database': '...', 'term': '...'},
)
```

## `sra_search` {#sra_search}

Searches NCBI SRA runs by keyword or BioProject accession.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `term`
- **Keywords:** sequencing runs, bioproject, reads, srr, raw data

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `term` | string | yes | Search keyword or BioProject accession. |
| `max_results` | integer | no | Maximum results. Default `10` |

??? note "Upstream description"

    Search SRA runs by keyword or BioProject.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'sra_search',
    {'term': '...'},
)
```

## `thousand_genomes_frequency` {#thousand_genomes_frequency}

Retrieves 1000 Genomes population allele frequencies for a variant (e.g. rsID) via Ensembl.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant`
- **Keywords:** 1000 genomes, allele frequency, population, rsid, superpopulation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | string | yes | Variant identifier (e.g. 'rs56116432'). |
| `species` | string | no | Species name. Default `homo_sapiens` |

??? note "Upstream description"

    Get population allele frequencies from 1000 Genomes.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'thousand_genomes_frequency',
    {'variant': '...'},
)
```
