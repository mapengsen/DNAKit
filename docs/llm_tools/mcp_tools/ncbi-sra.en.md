# NCBI SRA MCP Server (`ncbi-sra`)

DNAKit's MCP tool index keeps **4** DNA/gene tools of this server (6 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [CSI-Genomics-and-Data-Analytics-Core/ncbi-sra-mcp-server @ `5569215`](https://github.com/CSI-Genomics-and-Data-Analytics-Core/ncbi-sra-mcp-server/tree/5569215722ba8554c80aa66f8d971b5bfe09a984) |
| License | MIT |
| Command | `python -c from sra_mcp_server import mcp; mcp.run()` |
| Transports | stdio |
| Installation | Not on PyPI: install separately from the pinned source checkout in a Python &gt;=3.11 environment (pip install -r requirements.txt, which provides fastmcp), set workdir to the checkout and executable_path to that Python; optionally set NCBI_EMAIL and NCBI_API_KEY, and install sra-tools (prefetch, fasterq-dump) for toolkit downloads. |
| Effects | Queries NCBI Entrez E-utilities and the ENA portal API; sra_download starts background jobs that create the destination directory (default ./downloads) and write .sra or FASTQ files via urllib or the sra-tools prefetch/fasterq-dump subprocesses. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`sra_download`](#sra_download) | `sequencing_data` | Launches a download/conversion job for an SRA run accession to SRA or FASTQ using SRA Toolkit or direct HTTP/S3. |
| [`sra_files`](#sra_files) | `sequencing_data` | Retrieves direct FASTQ/SRA file URLs, MD5 checksums and sizes for run accessions via the ENA mirror. |
| [`sra_metadata`](#sra_metadata) | `sequencing_data` | Retrieves SRA run/study metadata (platform, library strategy, spots, bases, organism) via Entrez esummary. |
| [`sra_search`](#sra_search) | `sequencing_data` | Searches SRA runs and studies with Entrez query syntax (organism, platform, strategy, dates). |

## `sra_download` {#sra_download}

Launches a download/conversion job for an SRA run accession to SRA or FASTQ using SRA Toolkit or direct HTTP/S3.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `id`
- **Keywords:** fastq, sra toolkit, fasterq-dump, srr, paired-end

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `id` | string | yes | SRA run accession (e.g., "SRR1234567") |
| `format` | string | no | Output format - "sra" or "fastq" (default: "fastq") Default `fastq` |
| `split` | string | no | Split mode - "none", "spots", or "files" (default: "files") Default `files` |
| `threads` | integer | no | Number of threads (default: 4) Default `4` |
| `dest` | string | no | Destination directory (default: "./downloads") Default `./downloads` |
| `prefer_direct_http` | boolean | no | Use direct HTTP download when possible (default: True) Default `True` |

??? note "Upstream description"

    Launch download/conversion job using SRA Toolkit or direct HTTP/S3.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_download',
    {'id': '...'},
)
```

## `sra_files` {#sra_files}

Retrieves direct FASTQ/SRA file URLs, MD5 checksums and sizes for run accessions via the ENA mirror.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `ids`
- **Keywords:** fastq, ena, file urls, md5, srr

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | array of string | yes | List of run accessions (e.g., \["SRR1234567", "SRR7890123"\]) Note: Must be run accessions (SRR/ERR/DRR), not study accessions |

??? note "Upstream description"

    Retrieve direct file endpoints via ENA mirror.

    RECOMMENDED LIMITS:
    - For download planning: ≤ 50 IDs per call
    - For bulk file listing: ≤ 100 IDs per call

    FILE FORMATS PROVIDED:
    - FASTQ: Compressed FASTQ files (.fastq.gz)
    - SRA: Original SRA format files (.sra)
    - MD5: Checksums for verification
    - File sizes: In bytes for download planning

    DOWNLOAD OPTIONS:
    - Direct HTTP: Use provided URLs with wget/curl
    - ENA mirror: Fast, reliable European server
    - Multiple FASTQ files: For paired-end reads (fastq_1, fastq_2)

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_files',
    {'ids': []},
)
```

## `sra_metadata` {#sra_metadata}

Retrieves SRA run/study metadata (platform, library strategy, spots, bases, organism) via Entrez esummary.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `ids`
- **Keywords:** esummary, library strategy, platform, spots, bases

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | array of string | yes | List of SRA IDs (UIDs or accessions) Examples: \["123456789", "987654321"\] or \["SRR1234567", "SRR7890123"\] |

??? note "Upstream description"

    Retrieve metadata for SRA runs/studies via Entrez esummary.

    RECOMMENDED LIMITS:
    - For detailed analysis: ≤ 50 IDs per call
    - For bulk metadata: ≤ 200 IDs per call
    - Use pagination if you have more IDs

    METADATA FIELDS RETURNED:
    - accession: SRA run/study accession (e.g., SRR1234567)
    - title: Experiment title
    - organism: Scientific name of organism
    - platform: Sequencing platform (e.g., ILLUMINA, PACBIO)
    - instrument: Specific instrument model
    - library_strategy: Library strategy (RNA-Seq, WGS, etc.)
    - library_source: Library source (GENOMIC, TRANSCRIPTOMIC, etc.)
    - spots: Total number of spots/reads
    - bases: Total number of bases
    - study_accession: Associated study accession

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_metadata',
    {'ids': []},
)
```

## `sra_search` {#sra_search}

Searches SRA runs and studies with Entrez query syntax (organism, platform, strategy, dates).

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `term`
- **Keywords:** esearch, entrez, rna-seq, wgs, illumina

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `term` | string | yes | Search term using Entrez query syntax Examples: - "illumina\[platform\] AND human\[orgn\]" - "SARS-CoV-2\[orgn\] AND RNA-Seq\[strategy\]" - "breast cancer\[all\] AND illumina\[platform\]" |
| `retmax` | integer | no | Maximum results to return (default: 200, recommended ≤ 200) Default `200` |
| `retstart` | integer | no | Starting position for pagination (default: 0) Default `0` |
| `sort` | string | no | Sort order - "relevance" (default) or "date" Default `relevance` |
| `datetype` | string \| null | no | Date type filter - "pdat" (publication), "mdat" (modification), "edat" (entry) Default `None` |
| `mindate` | string \| null | no | Minimum date in YYYY/MM/DD format (e.g., "2020/01/01") Default `None` |
| `maxdate` | string \| null | no | Maximum date in YYYY/MM/DD format (e.g., "2023/12/31") Default `None` |

??? note "Upstream description"

    Search SRA runs/studies using Entrez esearch.

    RECOMMENDED LIMITS:
    - For exploration: retmax ≤ 50
    - For bulk retrieval: Use pagination with retmax ≤ 200
    - Maximum supported: 100,000 records per query

    SEARCH SYNTAX:
    - Organism: "Homo sapiens\[orgn\]" or "human\[orgn\]"
    - Platform: "illumina\[platform\]" or "pacbio\[platform\]"
    - Strategy: "RNA-Seq\[strategy\]" or "WGS\[strategy\]"
    - Combine: Use AND, OR, NOT operators
    - Date ranges: Use datetype, mindate, maxdate parameters

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-sra',
    'sra_search',
    {'term': '...'},
)
```
