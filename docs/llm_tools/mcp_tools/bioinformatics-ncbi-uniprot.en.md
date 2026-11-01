# Bioinformatics MCP Server (NCBI and UniProt) (`bioinformatics-ncbi-uniprot`)

DNAKit's MCP tool index keeps **4** DNA/gene tools of this server (10 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [AiAgentKarl/bioinformatics-mcp-server @ `b073742`](https://github.com/AiAgentKarl/bioinformatics-mcp-server/tree/b073742a7059f5cad816d0cf41f9fcdf49586933) |
| License | MIT |
| Command | `bioinformatics-server` |
| Transports | stdio |
| Installation | pip install bioinformatics-mcp-server (PyPI 0.1.0) into a Python &gt;=3.11 environment; no environment variables or API keys are needed. The wheel installs a top-level package named src, so use a dedicated environment. Pin the dependency 'mcp&lt;2': the upstream package leaves mcp unpinned and its server fails to import under mcp 2.x (FastMCP was renamed). |
| Effects | Read-only HTTPS queries to NCBI E-utilities (Gene, ClinVar, Taxonomy, Nucleotide, PubMed) and the UniProt REST API; every request sends tool and a fixed contact email parameter to NCBI, nothing is written locally. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`get_gene_details`](#get_gene_details) | `genome_annotation` | Retrieves detailed NCBI Gene information (name, summary, location) for an NCBI Gene ID. |
| [`get_ncbi_sequence`](#get_ncbi_sequence) | `sequence` | Fetches a nucleotide or protein sequence from NCBI by accession in FASTA or GenBank format. |
| [`search_clinvar_variants`](#search_clinvar_variants) | `variant` | Searches NCBI ClinVar for variants in a gene filtered by clinical significance (pathogenic, benign, uncertain, all). |
| [`search_genes`](#search_genes) | `genome_annotation` | Searches the NCBI Gene database by symbol or keyword for a given organism. |

## `get_gene_details` {#get_gene_details}

Retrieves detailed NCBI Gene information (name, summary, location) for an NCBI Gene ID.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_id`
- **Keywords:** ncbi gene, entrez, gene id, gene summary, genomic location

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_id` | string | yes |  |

??? note "Upstream description"

    Ruft detaillierte Informationen zu einem Gen via NCBI Gene ID ab.

    Args:
    gene_id: NCBI Gene ID (z.B. "672" für BRCA1, "7157" für TP53)

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'get_gene_details',
    {'gene_id': '...'},
)
```

## `get_ncbi_sequence` {#get_ncbi_sequence}

Fetches a nucleotide or protein sequence from NCBI by accession in FASTA or GenBank format.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `accession`
- **Keywords:** ncbi, accession, fasta, genbank, nucleotide, mrna

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes |  |
| `database` | string | no | Default `nucleotide` |
| `sequence_type` | string | no | Default `fasta` |

??? note "Upstream description"

    Ruft eine DNA- oder Protein-Sequenz aus NCBI ab.

    Args:
    accession: NCBI Accession-Nummer (z.B. "NM_007294" für BRCA1 mRNA, "NP_009225" für BRCA1 Protein)
    database: Datenbank — "nucleotide" (DNA/RNA) oder "protein". Standard: nucleotide
    sequence_type: Format — "fasta" (Sequenz) oder "gb" (GenBank-Eintrag). Standard: fasta

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'get_ncbi_sequence',
    {'accession': '...'},
)
```

## `search_clinvar_variants` {#search_clinvar_variants}

Searches NCBI ClinVar for variants in a gene filtered by clinical significance (pathogenic, benign, uncertain, all).

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `gene_symbol`
- **Keywords:** clinvar, pathogenic, clinical significance, mutation, gene variants

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes |  |
| `significance` | string | no | Default `pathogenic` |
| `max_results` | integer | no | Default `10` |

??? note "Upstream description"

    Sucht klinische Varianten (Mutationen) für ein Gen via NCBI ClinVar.

    Args:
    gene_symbol: Gen-Symbol (z.B. "BRCA1", "TP53", "CFTR")
    significance: Klinische Signifikanz — "pathogenic", "benign", "uncertain" oder "all". Standard: pathogenic
    max_results: Maximale Ergebnisanzahl (1-20). Standard: 10

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'search_clinvar_variants',
    {'gene_symbol': '...'},
)
```

## `search_genes` {#search_genes}

Searches the NCBI Gene database by symbol or keyword for a given organism.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** ncbi gene, gene search, gene symbol, entrez, organism

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes |  |
| `organism` | string | no | Default `human` |
| `max_results` | integer | no | Default `10` |

??? note "Upstream description"

    Sucht nach Genen in der NCBI Gene Datenbank.

    Args:
    query: Suchbegriff (z.B. "BRCA1", "TP53", "insulin receptor")
    organism: Organismus (z.B. "human", "mouse", "rat"). Standard: human
    max_results: Maximale Ergebnisanzahl (1-20). Standard: 10

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'search_genes',
    {'query': '...'},
)
```
