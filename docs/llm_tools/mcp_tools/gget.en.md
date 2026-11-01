# gget MCP (`gget`)

DNAKit's MCP tool index keeps **15** DNA/gene tools of this server (20 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [longevity-genie/gget-mcp @ `f3c6312`](https://github.com/longevity-genie/gget-mcp/tree/f3c63125e29954c20b3a81d1d5df70405ddf5896) |
| License | MIT |
| Command | `gget-mcp stdio` |
| Transports | stdio |
| Installation | Install gget-mcp and its separately specified gget dependencies in Python &gt;=3.10; the verified command is gget-mcp stdio. |
| Effects | Queries gget-supported services; selected tools may run external searches, load models, download assets or write outputs. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`gget_archs4`](#gget_archs4) | `expression` | Finds the most correlated genes or the tissue expression atlas for a gene using ARCHS4 RNA-seq data. |
| [`gget_bgee`](#gget_bgee) | `expression` | Retrieves orthologs or anatomical expression data for an Ensembl gene ID from the Bgee database. |
| [`gget_blast`](#gget_blast) | `sequence` | Runs NCBI BLAST (blastn, blastp, blastx, tblastn, tblastx) for a sequence against nt, nr, RefSeq or SwissProt. |
| [`gget_blat`](#gget_blat) | `sequence` | Aligns a nucleotide or amino acid sequence to a UCSC genome assembly (hg38, mm39, taeGut2) with BLAT to get genomic location. |
| [`gget_cellxgene`](#gget_cellxgene) | `expression` | Queries CZ CELLxGENE Census single-cell RNA-seq data filtered by genes, tissues and cell types. |
| [`gget_cosmic`](#gget_cosmic) | `variant` | Searches the COSMIC database for cancer somatic mutations by gene symbol or name. |
| [`gget_enrichr`](#gget_enrichr) | `genome_annotation` | Performs Enrichr enrichment analysis on a gene list against KEGG, ChEA, GO, GWAS, PanglaoDB or KEA libraries. |
| [`gget_info`](#gget_info) | `genome_annotation` | Retrieves gene and transcript metadata (from Ensembl, UniProt, NCBI) for one or more Ensembl IDs. |
| [`gget_muscle`](#gget_muscle) | `sequence` | Aligns multiple nucleotide or amino acid sequences with MUSCLE and saves the result as FASTA/aligned FASTA. |
| [`gget_mutate`](#gget_mutate) | `sequence` | Mutates nucleotide sequences according to mutations in standard annotation (e.g. c.2C&gt;T). |
| [`gget_opentargets`](#gget_opentargets) | `association` | Queries Open Targets for a gene's associated diseases, drugs, tractability, pharmacogenetics, expression, DepMap or interactions. |
| [`gget_ref`](#gget_ref) | `genome_annotation` | Fetches FTP download links for Ensembl reference genome DNA, cDNA, CDS, GTF or peptide files for a species. |
| [`gget_search`](#gget_search) | `genome_annotation` | Searches Ensembl for genes or transcripts whose names or descriptions match free-text terms in a species. |
| [`gget_search_genes`](#gget_search_genes) | `genome_annotation` | Searches Ensembl by gene symbol(s) with descriptive terms to return Ensembl gene or transcript IDs. |
| [`gget_seq`](#gget_seq) | `sequence` | Fetches nucleotide sequences (or translated amino acid sequences) for Ensembl IDs and saves them as FASTA. |

## `gget_archs4` {#gget_archs4}

Finds the most correlated genes or the tissue expression atlas for a gene using ARCHS4 RNA-seq data.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `gene`
- **Keywords:** rna-seq, coexpression, tissue expression, correlated genes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes | Gene symbol (e.g., 'STAT4') or Ensembl ID if ensembl=True (e.g., 'ENSG00000138378') |
| `which` | string | no | Analysis type - 'correlation' (most correlated genes) or 'tissue' (tissue expression atlas) Default `correlation` |
| `species` | string | no | Target species - 'human' or 'mouse' (only for tissue expression atlas) Default `human` |

??? note "Upstream description"

    Find correlated genes or tissue expression atlas using ARCHS4 RNA-seq database.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_archs4',
    {'gene': '...'},
)
```

## `gget_bgee` {#gget_bgee}

Retrieves orthologs or anatomical expression data for an Ensembl gene ID from the Bgee database.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `gene_id`
- **Keywords:** bgee, orthologs, tissue expression, ensembl id, anatomy

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes | Ensembl gene ID (e.g., 'ENSG00000012048' for BRCA1) |
| `type` | string | no | Type of data to retrieve - 'orthologs' or 'expression' Default `orthologs` |

??? note "Upstream description"

    Get orthologs or expression data for a gene from Bgee database.

    PREREQUISITE: Use search_genes to get Ensembl ID first.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_bgee',
    {'gene_id': '...'},
)
```

## `gget_blast` {#gget_blast}

Runs NCBI BLAST (blastn, blastp, blastx, tblastn, tblastx) for a sequence against nt, nr, RefSeq or SwissProt.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** blastn, similarity search, nt, refseq, alignment

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Nucleotide or amino acid sequence (string) or path to FASTA file |
| `program` | string | no | BLAST program - 'blastn', 'blastp', 'blastx', 'tblastn', 'tblastx', or 'default' (auto-detect) Default `default` |
| `database` | string | no | BLAST database - 'nt', 'nr', 'refseq_rna', 'refseq_protein', 'swissprot', or 'default' (auto-detect) Default `default` |

??? note "Upstream description"

    BLAST a nucleotide or amino acid sequence against any BLAST database.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_blast',
    {'sequence': '...'},
)
```

## `gget_blat` {#gget_blat}

Aligns a nucleotide or amino acid sequence to a UCSC genome assembly (hg38, mm39, taeGut2) with BLAT to get genomic location.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequence`
- **Keywords:** ucsc, genome mapping, hg38, mm39, genomic coordinates

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Nucleotide or amino acid sequence (string) or path to FASTA file containing one sequence |
| `assembly` | string | no | Genome assembly - 'human' (hg38), 'mouse' (mm39), 'zebrafinch' (taeGut2) Default `human` |

??? note "Upstream description"

    BLAT a nucleotide or amino acid sequence against any BLAT UCSC assembly.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_blat',
    {'sequence': '...'},
)
```

## `gget_cellxgene` {#gget_cellxgene}

Queries CZ CELLxGENE Census single-cell RNA-seq data filtered by genes, tissues and cell types.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** none
- **Keywords:** single-cell, scrna-seq, census, cell type, anndata

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | array \| string \| null | no | Gene name(s) or Ensembl ID(s) (e.g., \['ACE2', 'SLC5A1'\]) Default `None` |
| `tissue` | array \| string \| null | no | Tissue(s) to query (e.g., \['lung', 'blood'\]) Default `None` |
| `cell_type` | array \| string \| null | no | Cell type(s) to query (e.g., \['mucus secreting cell'\]) Default `None` |
| `species` | string | no | Target species - 'homo_sapiens' or 'mus_musculus' Default `homo_sapiens` |

??? note "Upstream description"

    Query single-cell RNA-seq data from CZ CELLxGENE Discover using Census.

    NOTE: Querying large datasets requires &gt;16 GB RAM and &gt;5 Mbps internet connection.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_cellxgene',
    {},
)
```

## `gget_cosmic` {#gget_cosmic}

Searches the COSMIC database for cancer somatic mutations by gene symbol or name.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `searchterm`
- **Keywords:** cosmic, somatic mutation, cancer, mutations

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `searchterm` | string | yes | Gene symbol or name to search for (e.g., 'PIK3CA', 'BRCA1') |
| `cosmic_tsv_path` | string \| null | no | Path to COSMIC TSV file (optional, uses default if None) Default `None` |

??? note "Upstream description"

    Search COSMIC database for cancer mutations or download COSMIC databases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_cosmic',
    {'searchterm': '...'},
)
```

## `gget_enrichr` {#gget_enrichr}

Performs Enrichr enrichment analysis on a gene list against KEGG, ChEA, GO, GWAS, PanglaoDB or KEA libraries.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genes`
- **Keywords:** enrichment, gene set, kegg, go, pathway

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genes` | array of string | yes | List of gene symbols (e.g., \['PHF14', 'RBM3'\]) or Ensembl IDs if ensembl=True |
| `database` | string | no | Reference database shortcuts: 'pathway' (KEGG), 'transcription' (ChEA), 'ontology' (GO), 'diseases_drugs' (GWAS), 'celltypes' (PanglaoDB), 'kinase_interactions' (KEA) Default `pathway` |
| `species` | string | no | Species database - 'human', 'mouse', 'fly', 'yeast', 'worm', 'fish' Default `human` |

??? note "Upstream description"

    Perform functional enrichment analysis on gene list using Enrichr.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_enrichr',
    {'genes': []},
)
```

## `gget_info` {#gget_info}

Retrieves gene and transcript metadata (from Ensembl, UniProt, NCBI) for one or more Ensembl IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ensembl_ids`
- **Keywords:** ensembl id, transcript, gene metadata, uniprot, ncbi

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_ids` | string \| array | yes | One or more Ensembl gene IDs (e.g., 'ENSG00000141510' or \['ENSG00000141510'\]) Also supports WormBase and FlyBase IDs |

??? note "Upstream description"

    Get detailed gene and transcript metadata using Ensembl IDs.

    PREREQUISITE: Use search_genes first to get Ensembl IDs from gene names/symbols.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_info',
    {'ensembl_ids': '...'},
)
```

## `gget_muscle` {#gget_muscle}

Aligns multiple nucleotide or amino acid sequences with MUSCLE and saves the result as FASTA/aligned FASTA.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequences`
- **Keywords:** multiple sequence alignment, msa, fasta, afa

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequences` | array \| string | yes | List of sequences or ABSOLUTE path to FASTA file containing sequences to be aligned |
| `output_path` | string \| null | no | ABSOLUTE path to output file (e.g., '/home/user/alignment.fasta'). AVOID relative paths as they cause file location issues. Auto-generated if not provided. Default `None` |
| `format` | string | no | Output format ('fasta' for FASTA format, 'afa' for aligned FASTA format) One of `fasta`, `afa` Default `fasta` |

??? note "Upstream description"

    Align sequences using MUSCLE and save to local file.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_muscle',
    {'sequences': '...'},
)
```

## `gget_mutate` {#gget_mutate}

Mutates nucleotide sequences according to mutations in standard annotation (e.g. c.2C&gt;T).

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `sequences`, `mutations`
- **Keywords:** mutation, hgvs, nucleotide sequence, substitution, in silico mutagenesis

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequences` | string \| array | yes | Path to FASTA file or sequences as string/list (e.g., 'AGCTAGCT' or \['ACTG', 'ATCG'\]) |
| `mutations` | string \| array | yes | Path to CSV/TSV file or mutations as string/list (e.g., 'c.2C&gt;T' or \['c.2C&gt;T', 'c.1A&gt;C'\]) Standard mutation annotation format required |

??? note "Upstream description"

    Mutate nucleotide sequences according to provided mutations in standard annotation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_mutate',
    {'sequences': '...', 'mutations': '...'},
)
```

## `gget_opentargets` {#gget_opentargets}

Queries Open Targets for a gene's associated diseases, drugs, tractability, pharmacogenetics, expression, DepMap or interactions.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `ensembl_id`
- **Keywords:** open targets, gene-disease, target, depmap, tractability

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_id` | string | yes | Ensembl gene ID (e.g., 'ENSG00000169194') |
| `resource` | string | no | Type of information - 'diseases', 'drugs', 'tractability', 'pharmacogenetics', 'expression', 'depmap', 'interactions' Default `diseases` |

??? note "Upstream description"

    Query OpenTargets for diseases, drugs, and other data associated with a gene.

    PREREQUISITE: Use search_genes to get Ensembl ID first.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_opentargets',
    {'ensembl_id': '...'},
)
```

## `gget_ref` {#gget_ref}

Fetches FTP download links for Ensembl reference genome DNA, cDNA, CDS, GTF or peptide files for a species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** reference genome, gtf, ftp, cdna, ensembl release

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species in format "genus_species" (e.g., "homo_sapiens"). Shortcuts supported: "human", "mouse", "human_grch37" Default `homo_sapiens` |
| `which` | string \| array | no | Which results to return. Options: 'gtf', 'cdna', 'dna', 'cds', 'cdrna', 'pep', 'all' Default `all` |

??? note "Upstream description"

    Fetch FTPs for reference genomes and annotations by species from Ensembl.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_ref',
    {},
)
```

## `gget_search` {#gget_search}

Searches Ensembl for genes or transcripts whose names or descriptions match free-text terms in a species.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `search_terms`
- **Keywords:** ensembl, gene search, free text, synonyms, species

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search_terms` | string \| array | yes | Search terms, names, or synonyms (e.g., 'cancer' or \['apoptosis', 'death'\]) |
| `species` | string | no | Target species (e.g., 'homo_sapiens', 'mus_musculus') Default `homo_sapiens` |

??? note "Upstream description"

    General search for any biological terms using gene symbols, names, or synonyms.

    This is a general search that looks broadly across gene names and descriptions.
    For specific gene symbol searches, use search_genes_simple instead.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_search',
    {'search_terms': '...'},
)
```

## `gget_search_genes` {#gget_search_genes}

Searches Ensembl by gene symbol(s) with descriptive terms to return Ensembl gene or transcript IDs.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `search_terms`
- **Keywords:** ensembl id, gene symbol, id resolution, batch

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search_terms` | string \| array | yes | SINGLE gene symbol OR LIST of gene symbols with descriptive terms Single: 'APP amyloid precursor' Batch: \['BACE1 beta secretase', 'MAPT tau', 'APOE apolipoprotein'\] |
| `species` | string | no | Target species (e.g., 'homo_sapiens', 'mus_musculus') Default `homo_sapiens` |
| `id_type` | string | no | "gene" (default) or "transcript" - whether to return genes or transcripts Default `gene` |

??? note "Upstream description"

    Search for specific genes using gene symbols with enhanced search strategy.

    🚀 **BATCH PROCESSING SUPPORTED**: This function can process multiple genes in a single call!
    Use this tool FIRST when you have gene names/symbols and need to find their Ensembl IDs.
    Returns Ensembl IDs which are required for get_gene_info and get_sequences tools.

    IMPORTANT: Due to limitations in Ensembl search, short gene names often fail to find results.
    For best results, provide descriptive terms along with gene symbols:

    RECOMMENDED FORMAT: "GENE_SYMBOL descriptive_terms"
    Examples:
    - Instead of: "APP"
    - Use: "APP amyloid precursor" or "APP amyloid beta precursor protein"
    - Instead of: \["BACE1", "MAPT"\]
    - Use: \["BACE1 beta secretase", "MAPT microtubule tau"\]

    This function uses AND search for multi-word terms and OR search for single words.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_search_genes',
    {'search_terms': '...'},
)
```

## `gget_seq` {#gget_seq}

Fetches nucleotide sequences (or translated amino acid sequences) for Ensembl IDs and saves them as FASTA.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `ensembl_ids`
- **Keywords:** fasta, ensembl id, nucleotide sequence, transcript sequence

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ensembl_ids` | string \| array | yes | One or more Ensembl gene IDs (e.g., 'ENSG00000141510' or \['ENSG00000141510'\]) |
| `translate` | boolean | no | If True, returns amino acid sequences; if False, returns nucleotide sequences Default `False` |
| `output_path` | string \| null | no | ABSOLUTE path to output file (e.g., '/home/user/sequences.fasta'). AVOID relative paths as they cause file location issues. Auto-generated if not provided. Default `None` |
| `format` | string | no | Output format (currently supports 'fasta') Default `fasta` |

??? note "Upstream description"

    Fetch sequences and save to local file in stdio mode.

    PREREQUISITE: Use search_genes first to get Ensembl IDs from gene names/symbols.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_seq',
    {'ensembl_ids': '...'},
)
```
