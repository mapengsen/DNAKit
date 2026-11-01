# BioThings MCP (`biothings`)

DNAKit's MCP tool index keeps **10** DNA/gene tools of this server (19 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [longevity-genie/biothings-mcp @ `bbfb870`](https://github.com/longevity-genie/biothings-mcp/tree/bbfb870f55ee18482200db8492915e59d70f5fe6) |
| License | MIT |
| Command | `biothings-mcp` |
| Transports | stdio |
| Installation | Install the pinned longevity-genie/biothings-mcp project separately in Python &gt;=3.10; use its installed biothings-mcp entry point. |
| Effects | Queries BioThings services and may write cached/exported result data under the upstream output directory. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`biothings_download_entrez_data`](#biothings_download_entrez_data) | `sequence` | Downloads records from an NCBI Entrez database (e.g. nucleotide, protein, gene) by ID in FASTA or GenBank format. |
| [`biothings_download_entrez_data_local`](#biothings_download_entrez_data_local) | `sequence` | Downloads NCBI Entrez records by ID in FASTA or GenBank format and saves them to a local file. |
| [`biothings_get_gene`](#biothings_get_gene) | `genome_annotation` | Retrieves the MyGene.info annotation record for one gene by Entrez or Ensembl gene ID. |
| [`biothings_get_genes`](#biothings_get_genes) | `genome_annotation` | Retrieves MyGene.info annotation records for multiple genes by comma-separated Entrez or Ensembl IDs. |
| [`biothings_get_variant`](#biothings_get_variant) | `variant` | Retrieves MyVariant.info annotation (dbSNP, ClinVar, dbNSFP, gnomAD, etc.) for one variant by HGVS ID or rsID. |
| [`biothings_get_variants`](#biothings_get_variants) | `variant` | Retrieves MyVariant.info annotations for multiple variants by comma-separated HGVS IDs or rsIDs. |
| [`biothings_query_genes`](#biothings_query_genes) | `genome_annotation` | Searches MyGene.info genes with fielded Lucene queries (e.g. symbol, name, taxid) returning gene annotations. |
| [`biothings_query_many_genes`](#biothings_query_many_genes) | `genome_annotation` | Batch-queries MyGene.info by multiple gene symbols or IDs against chosen scopes, useful for ID mapping. |
| [`biothings_query_many_variants`](#biothings_query_many_variants) | `variant` | Batch-queries MyVariant.info by multiple rsIDs or HGVS IDs against chosen scopes. |
| [`biothings_query_variants`](#biothings_query_variants) | `variant` | Searches MyVariant.info variants with Lucene queries over fields such as rsID, gene, variant type and PolyPhen/SIFT scores. |

## `biothings_download_entrez_data` {#biothings_download_entrez_data}

Downloads records from an NCBI Entrez database (e.g. nucleotide, protein, gene) by ID in FASTA or GenBank format.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `ids`, `db`, `reftype`
- **Keywords:** ncbi entrez, fasta, genbank, nucleotide, efetch, accession

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | array of string | yes |  |
| `db` | string | yes | One of `pubmed`, `protein`, `nuccore`, `ipg`, `nucleotide`, `structure`, `genome`, `annotinfo`, `assembly`, `bioproject`, `biosample`, `blastdbinfo` |
| `reftype` | string | yes | One of `fasta`, `gb` |

??? note "Upstream description"

    Download data from NCBI Entrez databases using Bio.Entrez.

    Downloads data records from specified NCBI Entrez databases. This tool is designed to be called
    by automated agents (like LLMs) or other services.

    **Critical Configuration:**
    The server hosting this API *must* have the `ENTREZ_EMAIL` environment variable set
    to a valid email address. NCBI requires this for Entrez queries to monitor usage
    and prevent abuse. Without it, NCBI may block requests.

    **Parameters:**
    - `ids` (List\[str\], required): A list of unique identifiers for the records to fetch
    from the specified Entrez database. Example: `\["NM_000546.6", "AY123456.1"\]`
    - `db` (DB_LITERAL, required): The target NCBI Entrez database.
    Common choices for sequences: 'nucleotide', 'protein'.
    Other examples: 'gene', 'pubmed', 'taxonomy'.
    Ensure the `ids` provided are appropriate for the selected `db`.
    - `reftype` (Literal\["fasta", "gb"\], required): The desired format for the
    downloaded data.
    - "fasta": Returns data in FASTA format.
    - "gb": Returns data in GenBank format.
    Ensure the chosen `reftype` is compatible with the selected `db`.

    **Returns:**
    On success: Returns the downloaded data as a single raw string with the
    data fetched from Entrez in the specified `reftype`.

    **Example Usage:**
    To fetch the FASTA sequence for human TP53 mRNA (NM_000546.6):
    ```
    download_entrez_data(
    ids=\["NM_000546.6"\],
    db="nucleotide",
    reftype="fasta"
    )
    ```

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_download_entrez_data',
    {'ids': [], 'db': 'pubmed', 'reftype': 'fasta'},
)
```

## `biothings_download_entrez_data_local` {#biothings_download_entrez_data_local}

Downloads NCBI Entrez records by ID in FASTA or GenBank format and saves them to a local file.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `ids`, `db`, `reftype`
- **Keywords:** ncbi entrez, fasta, genbank, nucleotide, download file

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | array of string | yes |  |
| `db` | string | yes | One of `pubmed`, `protein`, `nuccore`, `ipg`, `nucleotide`, `structure`, `genome`, `annotinfo`, `assembly`, `bioproject`, `biosample`, `blastdbinfo` |
| `reftype` | string | yes | One of `fasta`, `gb` |
| `output_path` | string \| null | no | Default `None` |

??? note "Upstream description"

    Download data from NCBI Entrez databases and save to local file.

    Same as download_entrez_data but saves the result to a local file instead of returning the content.
    This is useful for large downloads or when you want to persist the data.

    **Parameters:**
    - `ids` (List\[str\], required): A list of unique identifiers for the records to fetch
    - `db` (DB_LITERAL, required): The target NCBI Entrez database
    - `reftype` (Literal\["fasta", "gb"\], required): The desired format for the downloaded data
    - `output_path` (Optional\[str\]): Custom output path. If None, generates unique filename

    **Returns:**
    LocalFileResult containing:
    - `path`: Path to the saved file
    - `format`: File format used
    - `success`: Whether the operation succeeded
    - `error`: Error message if failed

    **Example Usage:**
    ```
    download_entrez_data_local(
    ids=\["NM_000546.6"\],
    db="nucleotide",
    reftype="fasta",
    output_path="tp53_sequence.fasta"
    )
    ```

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_download_entrez_data_local',
    {'ids': [], 'db': 'pubmed', 'reftype': 'fasta'},
)
```

## `biothings_get_gene` {#biothings_get_gene}

Retrieves the MyGene.info annotation record for one gene by Entrez or Ensembl gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** mygene, entrez gene id, ensembl id, gene annotation, symbol

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes |  |
| `fields` | string \| null | no | Default `all` |
| `species` | string \| null | no | Default `None` |
| `email` | string \| null | no | Default `None` |

??? note "Upstream description"

    Fetch a specific gene by Entrez or Ensembl ID.

    Retrieves detailed information for a **single, specific gene** using its exact known identifier.
    **This is the preferred tool over `query_genes` for fetching a specific gene when you already know
    its standard ID (Entrez or Ensembl) and don't need complex search filters.**

    **Supported Identifiers:**
    - Entrez Gene ID: e.g., `1017`
    - Ensembl Gene ID: e.g., `ENSG00000123374`

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_get_gene',
    {'gene_id': '...'},
)
```

## `biothings_get_genes` {#biothings_get_genes}

Retrieves MyGene.info annotation records for multiple genes by comma-separated Entrez or Ensembl IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_ids`
- **Keywords:** mygene, entrez gene id, ensembl id, batch, gene annotation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_ids` | string | yes |  |
| `fields` | string \| null | no | Default `all` |
| `species` | string \| null | no | Default `None` |
| `email` | string \| null | no | Default `None` |

??? note "Upstream description"

    Fetch multiple genes by a comma-separated list of Entrez or Ensembl IDs.

    Retrieves detailed information for **multiple specific genes** in a single request using their exact known identifiers.
    **This is the preferred tool over `query_many_genes` for fetching multiple specific genes when you already know
    their standard IDs (Entrez, Ensembl) and don't need complex search filters.**

    **Input Format:** Accepts a comma-separated list of gene IDs (Entrez or Ensembl).
    **Example:** `gene_ids=1017,1018` or `gene_ids=ENSG00000123374,ENSG00000134057`

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_get_genes',
    {'gene_ids': '...'},
)
```

## `biothings_get_variant` {#biothings_get_variant}

Retrieves MyVariant.info annotation (dbSNP, ClinVar, dbNSFP, gnomAD, etc.) for one variant by HGVS ID or rsID.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant_id`
- **Keywords:** myvariant, hgvs, rsid, dbsnp, clinvar, dbnsfp

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_id` | string | yes |  |
| `fields` | string \| null | no | Default `all` |
| `email` | string \| null | no | Default `None` |

??? note "Upstream description"

    Fetch a specific variant by HGVS or rsID.

    Retrieves detailed annotation data for a **single, specific variant** using its identifier.
    **This is the preferred tool over `query_variants` for fetching a specific variant when you already know
    its standard ID (HGVS or rsID) and don't need complex search filters.**

    **Supported Identifiers:**
    - HGVS ID (e.g., `chr7:g.140453134T&gt;C`). *Note: MyVariant.info primarily uses hg19-based HGVS IDs.*
    - dbSNP rsID (e.g., `rs58991260`)

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_get_variant',
    {'variant_id': '...'},
)
```

## `biothings_get_variants` {#biothings_get_variants}

Retrieves MyVariant.info annotations for multiple variants by comma-separated HGVS IDs or rsIDs.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `variant_ids`
- **Keywords:** myvariant, hgvs, rsid, batch, variant annotation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant_ids` | string | yes |  |
| `fields` | string \| null | no | Default `all` |
| `email` | string \| null | no | Default `None` |

??? note "Upstream description"

    Fetch multiple variants by a comma-separated list of HGVS or rsIDs.

    Retrieves annotation data for **multiple specific variants** in a single request using their identifiers.
    **This is the preferred tool over `query_many_variants` for fetching multiple specific variants when you already know
    their standard IDs (HGVS or rsID).**

    **Input Format:** Accepts a comma-separated list of variant IDs (HGVS or dbSNP rsIDs).
    **Examples:** `variant_ids=chr7:g.140453134T&gt;C,chr1:g.69511A&gt;G` or `variant_ids=rs58991260,rs2500`

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_get_variants',
    {'variant_ids': '...'},
)
```

## `biothings_query_genes` {#biothings_query_genes}

Searches MyGene.info genes with fielded Lucene queries (e.g. symbol, name, taxid) returning gene annotations.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `q`
- **Keywords:** mygene, lucene, gene search, symbol, taxid

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes |  |
| `fields` | string \| null | no | Default `all` |
| `size` | integer | no | Default `10` |
| `skip` | integer | no | Default `0` |
| `sort` | string \| null | no | Default `None` |
| `species` | string \| null | no | Default `None` |
| `email` | string \| null | no | Default `None` |

??? note "Upstream description"

    Search genes via Lucene query, returning gene details and query metadata.

    **IMPORTANT:** This endpoint requires structured queries using specific field names.
    Simple natural language queries like "CDK2 gene" or "human kinase" will **NOT** work.
    You **MUST** specify the field you are querying, e.g., `symbol:CDK2`, `name:"cyclin-dependent kinase 2"`, `taxid:9606`.
    Use this tool when you need to *search* for genes based on criteria, not when you already know the specific gene ID.
    If you know the exact Entrez or Ensembl ID, use the `get_gene` tool instead for faster retrieval.

    **Supported Query Features (based on Lucene syntax):**
    1. Simple Term Queries: `q=cdk2` (Searches across default fields)
    2. Fielded Queries: `q=symbol:CDK2`, `q=name:"cyclin-dependent kinase 2"`
    3. Range Queries: `q=taxid:\[9606 TO 10090\]`
    4. Boolean Queries: `q=symbol:CDK2 AND taxid:9606`
    5. Wildcard Queries: `q=symbol:CDK*`

    Returns gene details including symbol, name, taxid, and entrezgene.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_query_genes',
    {'q': '...'},
)
```

## `biothings_query_many_genes` {#biothings_query_many_genes}

Batch-queries MyGene.info by multiple gene symbols or IDs against chosen scopes, useful for ID mapping.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query_list`
- **Keywords:** mygene, id mapping, symbol, entrez, ensembl, batch

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query_list` | string | yes |  |
| `scopes` | string \| null | no | Default `entrezgene,ensemblgene,retired` |
| `fields` | string \| null | no | Default `all` |
| `species` | string \| null | no | Default `None` |
| `email` | string \| null | no | Default `None` |
| `size` | integer | no | Default `10` |

??? note "Upstream description"

    Batch query genes by multiple terms, returning a list of gene details.

    Perform multiple gene searches in a single request using a comma-separated list of query terms.
    Unlike `query_genes`, the `query_list` parameter takes multiple **terms** (like gene IDs, symbols, names)
    rather than full query strings. The `scopes` parameter defines which fields these terms should be searched against.

    **Endpoint Usage:**
    - Query multiple symbols: `query_list=CDK2,BRCA1` with `scopes=symbol`
    - Query multiple Entrez IDs: `query_list=1017,672` with `scopes=entrezgene`
    - Query mixed IDs/symbols: `query_list=CDK2,672` with `scopes=symbol,entrezgene`

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_query_many_genes',
    {'query_list': '...'},
)
```

## `biothings_query_many_variants` {#biothings_query_many_variants}

Batch-queries MyVariant.info by multiple rsIDs or HGVS IDs against chosen scopes.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `query_list`
- **Keywords:** myvariant, rsid, hgvs, batch, variant annotation

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query_list` | string | yes |  |
| `scopes` | string \| null | no | Default `None` |
| `fields` | string \| null | no | Default `all` |
| `email` | string \| null | no | Default `None` |

??? note "Upstream description"

    Batch query variants by multiple identifiers (e.g., rsIDs, HGVS IDs).

    Perform multiple variant queries in a single request using a comma-separated list of variant identifiers.
    This tool takes multiple **terms** (like rsIDs, HGVS IDs) in `query_list` and searches for them within the specified `scopes`.

    **Endpoint Usage:**
    - Query multiple rsIDs: `query_list=rs58991260,rs2500` with `scopes=dbsnp.rsid`
    - Query multiple HGVS IDs: `query_list=chr7:g.140453134T&gt;C,chr1:g.69511A&gt;G`
    - Query mixed IDs: `query_list=rs58991260,chr1:g.69511A&gt;G` with `scopes=dbsnp.rsid,_id`

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_query_many_variants',
    {'query_list': '...'},
)
```

## `biothings_query_variants` {#biothings_query_variants}

Searches MyVariant.info variants with Lucene queries over fields such as rsID, gene, variant type and PolyPhen/SIFT scores.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `q`
- **Keywords:** myvariant, lucene, polyphen, dbnsfp, variant search

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `q` | string | yes |  |
| `fields` | string \| null | no | Default `all` |
| `size` | integer | no | Default `10` |
| `skip` | integer | no | Default `0` |
| `sort` | string \| null | no | Default `None` |
| `email` | string \| null | no | Default `None` |

??? note "Upstream description"

    Search variants via Lucene query (e.g., rsID, gene name), returning variant details and query metadata.

    Search for variants using a query string with various filtering options, leveraging the MyVariant.info API.
    **Use this tool for *searching* variants based on criteria.**
    If you already know the exact variant ID (HGVS, rsID), use the `get_variant` tool for faster direct retrieval.

    **Supported Query Features (Lucene syntax):**
    1. Simple Queries: `q=rs58991260` (Find by rsID)
    2. Fielded Queries: `q=dbsnp.vartype:snp`, `q=dbnsfp.polyphen2.hdiv.pred:(D P)`
    3. Range Queries: `q=dbnsfp.polyphen2.hdiv.score:&gt;0.99`
    4. Wildcard Queries: `q=dbnsfp.genename:CDK*`
    5. Boolean Queries: `q=_exists_:dbsnp AND dbsnp.vartype:snp`
    6. Genomic Interval Queries: `q=chr1:69000-70000`

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'biothings',
    'biothings_query_variants',
    {'q': '...'},
)
```
