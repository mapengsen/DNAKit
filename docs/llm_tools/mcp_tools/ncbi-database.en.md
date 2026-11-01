# NCBI Database MCP (hpend2373) (`ncbi-database`)

DNAKit's MCP tool index keeps **4** DNA/gene tools of this server (4 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [hpend2373/NCBI-Database-MCP @ `03a6b80`](https://github.com/hpend2373/NCBI-Database-MCP/tree/03a6b80490753ba47166903daa7d7bce35848af3) |
| License | MIT declared in pyproject.toml; standalone license text not found |
| Command | `python src/gene_to_genomic_server.py` |
| Transports | stdio |
| Installation | Prepare the pinned NCBI-Database-MCP checkout and its requirements, set workdir to the checkout and executable_path to its Python. Uses the README-supported standard MCP script rather than executing its shell installer/launcher. |
| Effects | Queries NCBI gene/sequence/GEO resources; selected results can contain supplied gene or region queries. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`gene_to_genomic_sequence`](#gene_to_genomic_sequence) | `sequence` | Retrieves the genomic, CDS, mRNA or protein sequence of a gene by symbol and organism via NCBI E-utilities. |
| [`get_genomic_sequence`](#get_genomic_sequence) | `sequence` | Retrieves genomic DNA sequence for a chromosome accession (e.g. NC_000017.11) between start and end positions. |
| [`search_gene_info`](#search_gene_info) | `genome_annotation` | Searches NCBI Gene for a gene by name and organism and returns its information and genomic location. |
| [`search_geo_datasets`](#search_geo_datasets) | `sequencing_data` | Searches NCBI GEO datasets by disease and organism, filtered by study type such as RNA-seq, array or ChIP-seq. |

## `gene_to_genomic_sequence` {#gene_to_genomic_sequence}

Retrieves the genomic, CDS, mRNA or protein sequence of a gene by symbol and organism via NCBI E-utilities.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `gene_name`
- **Keywords:** gene symbol, cds, mrna, fasta, genbank

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string | yes | Gene symbol or name (e.g., TP53, BRCA1) |
| `organism` | string | no | Organism name (default: human) Default `human` |
| `sequence_type` | string | no | Type of sequence to retrieve (default: genomic) One of `genomic`, `cds`, `mrna`, `protein` Default `genomic` |
| `output_format` | string | no | Output format (default: fasta) One of `fasta`, `genbank`, `json` Default `fasta` |

??? note "Upstream description"

    Convert gene name to genomic DNA sequence using NCBI E-utilities

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'gene_to_genomic_sequence',
    {'gene_name': '...'},
)
```

## `get_genomic_sequence` {#get_genomic_sequence}

Retrieves genomic DNA sequence for a chromosome accession (e.g. NC_000017.11) between start and end positions.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `chromosome`, `start`, `end`
- **Keywords:** chromosome, coordinates, refseq, fasta, genbank

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome accession (e.g., NC_000017.11) |
| `start` | integer | yes | Start position |
| `end` | integer | yes | End position |
| `output_format` | string | no | Output format (default: fasta) One of `fasta`, `genbank` Default `fasta` |

??? note "Upstream description"

    Get genomic sequence from chromosome coordinates

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'get_genomic_sequence',
    {'chromosome': '...', 'start': 0, 'end': 0},
)
```

## `search_gene_info` {#search_gene_info}

Searches NCBI Gene for a gene by name and organism and returns its information and genomic location.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_name`
- **Keywords:** ncbi gene, genomic location, gene id, chromosome, organism

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string | yes | Gene symbol or name |
| `organism` | string | no | Organism name (default: human) Default `human` |

??? note "Upstream description"

    Search for gene information including genomic location

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'search_gene_info',
    {'gene_name': '...'},
)
```

## `search_geo_datasets` {#search_geo_datasets}

Searches NCBI GEO datasets by disease and organism, filtered by study type such as RNA-seq, array or ChIP-seq.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `disease`
- **Keywords:** gse, expression profiling, rna-seq, microarray, disease

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `disease` | string | yes | Disease or condition name (e.g., 'cancer', 'diabetes', 'Alzheimer') |
| `organism` | string | no | Organism (default: 'Homo sapiens') One of `Homo sapiens`, `Mus musculus`, `Rattus norvegicus` Default `Homo sapiens` |
| `study_type` | string | no | Type of expression study (default: Expression profiling by high throughput sequencing) One of `Expression profiling by array`, `Expression profiling by high throughput sequencing`, `Genome binding/occupancy profiling by high throughput sequencing` Default `Expression profiling by high throughput sequencing` |
| `max_results` | integer | no | Maximum number of results to return (default: 10) Default `10` |

??? note "Upstream description"

    Search GEO datasets by disease/condition and organism

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'search_geo_datasets',
    {'disease': '...'},
)
```
