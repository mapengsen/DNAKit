# UCSC Genome MCP (`ucsc-genome`)

DNAKit's MCP tool index keeps **11** DNA/gene tools of this server (12 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [hlydecker/ucsc-genome-mcp @ `7f0f83a`](https://github.com/hlydecker/ucsc-genome-mcp/tree/7f0f83a2f1935916f8cd35ddfb8dd3c36dac2e84) |
| License | MIT |
| Command | `python ucsc-genome-mcp.py` |
| Transports | stdio |
| Installation | Prepare the pinned source checkout and its requirements in Python &gt;=3.10, set workdir to the checkout and executable_path to that Python. Use the existing hyphenated ucsc-genome-mcp.py script; upstream README/packaging also mention nonexistent underscore filenames. Pin 'mcp&lt;2': under mcp 2.x the script fails to import. |
| Effects | Sends genome, track, region and sequence queries to UCSC and caller-selected public assembly hubs. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`find_genome`](#find_genome) | `genome_annotation` | Searches UCSC Genome Browser genome assemblies by name or accession with filters for assembly level, category and year. |
| [`get_sequence`](#get_sequence) | `sequence` | Retrieves DNA sequence for a chromosome or coordinate range from a UCSC or hub genome assembly, optionally reverse complemented. |
| [`get_track_data`](#get_track_data) | `genome_annotation` | Retrieves records from a UCSC or hub data track, optionally restricted to a chromosome region. |
| [`list_chromosomes`](#list_chromosomes) | `genome_annotation` | Lists chromosomes (and sizes) of a UCSC database genome, assembly hub or track. |
| [`list_files`](#list_files) | `genome_annotation` | Lists download files available for a UCSC genome assembly such as hg38 or mm10. |
| [`list_genark_genomes`](#list_genark_genomes) | `genome_annotation` | Lists UCSC GenArk assembly hub genomes or tests whether a specific genome exists. |
| [`list_hub_genomes`](#list_hub_genomes) | `genome_annotation` | Lists all genomes available in a given track hub or assembly hub URL. |
| [`list_public_hubs`](#list_public_hubs) | `genome_annotation` | Lists all public track hubs registered with the UCSC Genome Browser. |
| [`list_tracks`](#list_tracks) | `genome_annotation` | Lists all data tracks available for a UCSC database genome or hub genome. |
| [`list_ucsc_genomes`](#list_ucsc_genomes) | `genome_annotation` | Lists all UCSC Genome Browser database genome assemblies. |
| [`search_genome`](#search_genome) | `genome_annotation` | Searches within a UCSC genome assembly across track data, help docs and public hubs for a term. |

## `find_genome` {#find_genome}

Searches UCSC Genome Browser genome assemblies by name or accession with filters for assembly level, category and year.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `query`
- **Keywords:** genome assembly, gcf accession, species, assembly level, grch38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search string to find genomes (e.g., 'dog', 'GRCh38', 'GCF_028858775.2') |
| `browser` | string | no | Filter by browser availability (default: mustExist) One of `mustExist`, `mayExist`, `notExist` |
| `stats_only` | boolean | no | Only show statistics about search results |
| `year` | integer | no | Filter results by year |
| `category` | string | no | Filter by NCBI category One of `reference`, `representative` |
| `status` | string | no | Filter by NCBI status One of `reference`, `representative` |
| `level` | string | no | Filter by NCBI assembly level One of `complete`, `chromosome`, `scaffold`, `contig` |
| `max_items` | integer | no | Maximum number of items to return (default: 1000000, use -1 for max) |

??? note "Upstream description"

    Search for a genome in the UCSC browser using a search string. Supports advanced search with +word (force inclusion), -word (exclusion), and word* (wildcard).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'find_genome',
    {'query': '...'},
)
```

## `get_sequence` {#get_sequence}

Retrieves DNA sequence for a chromosome or coordinate range from a UCSC or hub genome assembly, optionally reverse complemented.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `genome`, `chrom`
- **Keywords:** dna sequence, chromosome, coordinates, reverse complement, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | Genome assembly name (e.g., 'hg38') |
| `chrom` | string | yes | Chromosome name (e.g., 'chr1', 'chrM') |
| `start` | integer | no | Start coordinate (0-based, optional, requires end) |
| `end` | integer | no | End coordinate (1-based, optional, requires start) |
| `hub_url` | string | no | URL of assembly hub (optional) |
| `reverse_complement` | boolean | no | Return reverse complement of sequence |

??? note "Upstream description"

    Retrieve DNA sequence from a specified genome assembly. Can retrieve entire chromosome or specific coordinates.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'get_sequence',
    {'genome': '...', 'chrom': '...'},
)
```

## `get_track_data` {#get_track_data}

Retrieves records from a UCSC or hub data track, optionally restricted to a chromosome region.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome`, `track`
- **Keywords:** annotation track, genomic region, bed, coordinates, ucsc table

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | Genome assembly name |
| `track` | string | yes | Track name |
| `chrom` | string | no | Chromosome name (optional) |
| `start` | integer | no | Start coordinate (0-based, optional, requires end) |
| `end` | integer | no | End coordinate (1-based, optional, requires start) |
| `hub_url` | string | no | URL of track/assembly hub (optional) |
| `max_items` | integer | no | Maximum number of items to return (default: 1000000) |
| `json_output_arrays` | boolean | no | Return data as JSON arrays instead of objects |

??? note "Upstream description"

    Retrieve data from a specified track in a hub or UCSC database genome. Can be filtered by chromosome and coordinates.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'get_track_data',
    {'genome': '...', 'track': '...'},
)
```

## `list_chromosomes` {#list_chromosomes}

Lists chromosomes (and sizes) of a UCSC database genome, assembly hub or track.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome`
- **Keywords:** chromosome sizes, assembly, contigs, hub, chrom

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | Genome assembly name |
| `hub_url` | string | no | URL of track/assembly hub (optional) |
| `track` | string | no | Specific track name to list chromosomes from (optional) |

??? note "Upstream description"

    List chromosomes in an assembly hub, track hub, or UCSC database genome. Optionally filter by specific track.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_chromosomes',
    {'genome': '...'},
)
```

## `list_files` {#list_files}

Lists download files available for a UCSC genome assembly such as hg38 or mm10.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome`
- **Keywords:** downloads, assembly files, 2bit, annotation files, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | Genome assembly name (e.g., 'hg38', 'mm10') |
| `format` | string | no | Output format (default: json) One of `json`, `text` |
| `max_items` | integer | no | Maximum number of items to return |

??? note "Upstream description"

    List download files available for a specified UCSC genome assembly.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_files',
    {'genome': '...'},
)
```

## `list_genark_genomes` {#list_genark_genomes}

Lists UCSC GenArk assembly hub genomes or tests whether a specific genome exists.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** genark, assembly hub, genome assembly, gca, species

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | no | Specific genome to test for existence (optional) |
| `max_items` | integer | no | Maximum number of items to return (default: 1000000) |

??? note "Upstream description"

    List UCSC Genome Browser database genomes from assembly hub host (GenArk). Can also test for existence of a specific genome.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_genark_genomes',
    {},
)
```

## `list_hub_genomes` {#list_hub_genomes}

Lists all genomes available in a given track hub or assembly hub URL.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `hub_url`
- **Keywords:** track hub, assembly hub, genome assembly, hub url, genomes

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `hub_url` | string | yes | URL of the track hub or assembly hub |

??? note "Upstream description"

    List all genomes available in a specified track or assembly hub.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_hub_genomes',
    {'hub_url': '...'},
)
```

## `list_public_hubs` {#list_public_hubs}

Lists all public track hubs registered with the UCSC Genome Browser.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** track hubs, public hubs, annotation collections, genome browser, datasets

??? note "Upstream description"

    List all available public track hubs in the UCSC Genome Browser.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_public_hubs',
    {},
)
```

## `list_tracks` {#list_tracks}

Lists all data tracks available for a UCSC database genome or hub genome.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `genome`
- **Keywords:** annotation tracks, track list, genome browser, hub, composite

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | yes | Genome assembly name |
| `hub_url` | string | no | URL of track/assembly hub (optional, required with genome for hub tracks) |
| `track_leaves_only` | boolean | no | Only show tracks without composite container information |

??? note "Upstream description"

    List all data tracks available in a specified hub or UCSC database genome.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_tracks',
    {'genome': '...'},
)
```

## `list_ucsc_genomes` {#list_ucsc_genomes}

Lists all UCSC Genome Browser database genome assemblies.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** none
- **Keywords:** genome assembly, hg38, mm10, database genomes, species

??? note "Upstream description"

    List all UCSC Genome Browser database genomes available on the database host.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'list_ucsc_genomes',
    {},
)
```

## `search_genome` {#search_genome}

Searches within a UCSC genome assembly across track data, help docs and public hubs for a term.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `search`, `genome`
- **Keywords:** gene search, position, track search, hg38, genome browser

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search` | string | yes | Search term |
| `genome` | string | yes | Genome assembly to search in |
| `categories` | string | no | Restrict search to specific categories One of `helpDocs`, `publicHubs`, `trackDb` |

??? note "Upstream description"

    Search for matches within a UCSC Genome Browser genome assembly across tracks, help docs, and public hubs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ucsc-genome',
    'search_genome',
    {'search': '...', 'genome': '...'},
)
```
