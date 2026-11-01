# Rewire Genomics MCP (`rewire-genomics`)

DNAKit's MCP tool index keeps **20** DNA/gene tools of this server (24 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [rewire-bio/genomics-mcp @ `c419069`](https://github.com/rewire-bio/genomics-mcp/tree/c4190699201b8a2f447dc52075c5977d7b32e287) |
| License | MIT |
| Command | `genomics-mcp` |
| Transports | stdio |
| Installation | Install the pinned repository checkout into a Python &gt;=3.12,&lt;3.13 environment (pip install of the checkout, with --no-binary pybigwig; needs a C compiler plus libcurl and zlib development files for remote bigWig/bigBed); the upstream docs say there is no PyPI package although PyPI lists rewire-genomics-mcp 0.1.0 by the same author, so prefer the checkout. Optionally set GENOMICS_MCP_ALLOWED_ROOTS to the folders the server may read, and EGA credentials only for controlled EGA data. |
| Effects | Calls public services (EGA, ENA, ENCODE, GEO, NCBI, Ensembl, HGNC, ClinVar, gnomAD, UniProt, Open Targets) over HTTPS and reads local files under the allowed roots; fetch_file, cancel_transfer and convert_sra_run write downloaded files into the server work directory (default under the XDG cache), values from private files go to external APIs only when a call sets allow_external_annotation. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`compare_samples`](#compare_samples) | `variant` | Compares genotypes and coverage at a genomic interval across multiple files or samples. |
| [`convert_sra_run`](#convert_sra_run) | `sequencing_data` | Downloads a public SRR/ERR/DRR run with SRA Toolkit prefetch and converts it to FASTQ via fasterq-dump as a job. |
| [`describe_dataset`](#describe_dataset) | `sequencing_data` | Describes a study or dataset by native accession with linked entities from EGA, ENA, ENCODE, GEO or NCBI Datasets. |
| [`fetch_file`](#fetch_file) | `sequencing_data` | Starts a bounded download of a dataset file and its index to the local work dir, with checksum verification. |
| [`get_coverage`](#get_coverage) | `sequencing_data` | Computes read depth over an interval from indexed BAM/CRAM with mapping/base quality and flag filters. |
| [`get_features`](#get_features) | `genome_annotation` | Returns features overlapping an interval from indexed BED, GFF3, GTF or bigBed files. |
| [`get_pileup`](#get_pileup) | `sequencing_data` | Returns per-position base counts over an interval from indexed BAM/CRAM alignments. |
| [`get_reads`](#get_reads) | `sequencing_data` | Returns alignments overlapping an interval from indexed BAM/CRAM with samtools-style flag filters. |
| [`get_sample_metadata`](#get_sample_metadata) | `sequencing_data` | Returns one sample's metadata and phenotype values as supplied by EGA, ENA, ENCODE, GEO or NCBI Datasets. |
| [`get_sequence`](#get_sequence) | `sequence` | Returns reference sequence for a genomic interval from an indexed FASTA file. |
| [`get_signal`](#get_signal) | `regulation_epigenome` | Returns signal values or binned summaries over an interval from a bigWig file. |
| [`get_variants`](#get_variants) | `variant` | Returns VCF/BCF variant records and optional sample genotypes overlapping an interval. |
| [`inspect_locus`](#inspect_locus) | `genome_annotation` | Combines data from several files at one locus with optional public reference evidence. |
| [`list_files`](#list_files) | `sequencing_data` | Lists files of a study/dataset (EGA, ENA, ENCODE, GEO, NCBI Datasets, local, S3) with format, size and index status. |
| [`list_samples`](#list_samples) | `sequencing_data` | Lists samples of a study/dataset as supplied by EGA, ENA, ENCODE, GEO or NCBI Datasets. |
| [`lookup_gene`](#lookup_gene) | `genome_annotation` | Returns gene identifiers, transcripts and source-attributed annotations for a gene on a given assembly. |
| [`lookup_variant`](#lookup_variant) | `variant` | Returns source-attributed evidence (ClinVar, gnomAD, Ensembl) for a variant given as VCF-style, HGVS or rsID. |
| [`normalize_variant`](#normalize_variant) | `variant` | Normalizes a VCF-style, HGVS or rsID variant on an explicit assembly with a trace. |
| [`resolve_identifier`](#resolve_identifier) | `genome_annotation` | Resolves a gene, transcript, protein or variant identifier and reports ambiguity. |
| [`search_datasets`](#search_datasets) | `sequencing_data` | Searches one archive (EGA, ENA, ENCODE, GEO, NCBI Datasets) for studies or datasets by text, assembly and organism. |

## `compare_samples` {#compare_samples}

Compares genotypes and coverage at a genomic interval across multiple files or samples.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `interval`, `files`
- **Keywords:** genotype comparison, coverage, interval, multi-sample, vcf, bam
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `interval` | any | yes |  |
| `files` | array | yes |  |
| `samples` | array \| null | no | Default `None` |
| `reference` | ? \| null | no | Default `None` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Compare genotypes/coverage at an interval across files or samples.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'compare_samples',
    {'interval': '...', 'files': []},
)
```

## `convert_sra_run` {#convert_sra_run}

Downloads a public SRR/ERR/DRR run with SRA Toolkit prefetch and converts it to FASTQ via fasterq-dump as a job.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** sra, fastq, fasterq-dump, prefetch, sequencing reads, srr

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | Public SRR, ERR or DRR run accession. |
| `budget_bytes` | integer \| null | no | Default `None` |
| `timeout_s` | number \| null | no | Default `None` |

??? note "Upstream description"

    Optional: download a public SRA run with SRA Toolkit prefetch and convert the whole run to local FASTQ files with fasterq-dump (paired and unpaired reads kept). Starts a job; poll get_transfer_status. Needs SRA Toolkit; ENA FASTQ files via list_files/fetch_file do not. FASTQ is not locus-ready.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'convert_sra_run',
    {'accession': '...'},
)
```

## `describe_dataset` {#describe_dataset}

Describes a study or dataset by native accession with linked entities from EGA, ENA, ENCODE, GEO or NCBI Datasets.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `source`, `accession`
- **Keywords:** ena, geo, encode, ega, study accession, dataset
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes |  |
| `accession` | string | yes |  |

??? note "Upstream description"

    Describe a study/dataset by native accession, with linked entities.
    Available for: ega, ena, encode, geo, ncbi_datasets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'describe_dataset',
    {'source': '...', 'accession': '...'},
)
```

## `fetch_file` {#fetch_file}

Starts a bounded download of a dataset file and its index to the local work dir, with checksum verification.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `file`
- **Keywords:** download, bam, vcf, index, checksum, transfer

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `include_index` | boolean | no | Default `True` |
| `budget_bytes` | integer \| null | no | Default `None` |
| `verify_checksum` | boolean | no | Default `True` |
| `prepare` | boolean | no | Default `False` |

??? note "Upstream description"

    Start a bounded download of a file (and index) to the local work dir. Returns a transfer job; files above the default limit need budget_bytes.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'fetch_file',
    {'file': '...'},
)
```

## `get_coverage` {#get_coverage}

Computes read depth over an interval from indexed BAM/CRAM with mapping/base quality and flag filters.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `file`, `interval`
- **Keywords:** read depth, bam, cram, samtools depth, coverage, mapping quality
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `interval` | any | yes |  |
| `reference` | ? \| null | no | Default `None` |
| `min_mapping_quality` | integer | no | Default `0` |
| `min_base_quality` | integer | no | Default `0` |
| `exclude_flags` | integer | no | Default `1796` |
| `bin_size` | integer \| null | no | Default `None` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Read depth over an interval from indexed BAM/CRAM (default excludes UNMAP, SECONDARY, QCFAIL, DUP, like samtools depth).
    Available for: bam, cram.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_coverage',
    {'file': '...', 'interval': '...'},
)
```

## `get_features` {#get_features}

Returns features overlapping an interval from indexed BED, GFF3, GTF or bigBed files.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `file`, `interval`
- **Keywords:** bed, gff3, gtf, bigbed, interval, gene features
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `interval` | any | yes |  |
| `feature_types` | array \| null | no | Default `None` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Features overlapping an interval from indexed BED/GFF3/GTF or bigBed.
    Available for: bed, bigbed, gff3, gtf.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_features',
    {'file': '...', 'interval': '...'},
)
```

## `get_pileup` {#get_pileup}

Returns per-position base counts over an interval from indexed BAM/CRAM alignments.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `file`, `interval`
- **Keywords:** pileup, base counts, bam, cram, mpileup
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `interval` | any | yes |  |
| `reference` | ? \| null | no | Default `None` |
| `min_mapping_quality` | integer | no | Default `0` |
| `min_base_quality` | integer | no | Default `13` |
| `exclude_flags` | integer | no | Default `1796` |
| `max_depth` | integer | no | Default `8000` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Per-position base counts over an interval from indexed BAM/CRAM.
    Available for: bam, cram.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_pileup',
    {'file': '...', 'interval': '...'},
)
```

## `get_reads` {#get_reads}

Returns alignments overlapping an interval from indexed BAM/CRAM with samtools-style flag filters.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `file`, `interval`
- **Keywords:** alignments, bam, cram, reads, samtools, flags
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `interval` | any | yes |  |
| `reference` | ? \| null | no | Default `None` |
| `min_mapping_quality` | integer | no | Default `0` |
| `require_flags` | integer | no | Default `0` |
| `exclude_flags` | integer | no | Default `0` |
| `include_sequence` | boolean | no | Default `False` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Alignments overlapping an interval from indexed BAM/CRAM. Flag filters follow samtools -f/-F. CRAM needs a matching reference unless self-contained.
    Available for: bam, cram.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_reads',
    {'file': '...', 'interval': '...'},
)
```

## `get_sample_metadata` {#get_sample_metadata}

Returns one sample's metadata and phenotype values as supplied by EGA, ENA, ENCODE, GEO or NCBI Datasets.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `source`, `accession`
- **Keywords:** sample metadata, biosample, phenotype, ena, geo
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes |  |
| `accession` | string | yes |  |

??? note "Upstream description"

    Get one sample's metadata and phenotype values exactly as supplied by the source.
    Available for: ega, ena, encode, geo, ncbi_datasets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_sample_metadata',
    {'source': '...', 'accession': '...'},
)
```

## `get_sequence` {#get_sequence}

Returns reference sequence for a genomic interval from an indexed FASTA file.

- **Category:** `sequence` (DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation)
- **Required arguments:** `file`, `interval`
- **Keywords:** fasta, reference sequence, interval, faidx
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `interval` | any | yes |  |

??? note "Upstream description"

    Reference sequence for an interval from an indexed FASTA.
    Available for: fasta.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_sequence',
    {'file': '...', 'interval': '...'},
)
```

## `get_signal` {#get_signal}

Returns signal values or binned summaries over an interval from a bigWig file.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `file`, `interval`
- **Keywords:** bigwig, signal track, bins, chip-seq, atac-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `interval` | any | yes |  |
| `bins` | integer \| null | no | Default `None` |
| `summary` | string | no | One of `mean`, `min`, `max`, `coverage`, `std`, `sum` Default `mean` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Signal values or binned summaries over an interval from bigWig.
    Available for: bigwig.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_signal',
    {'file': '...', 'interval': '...'},
)
```

## `get_variants` {#get_variants}

Returns VCF/BCF variant records and optional sample genotypes overlapping an interval.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `file`, `interval`
- **Keywords:** vcf, bcf, genotypes, interval, pass filter
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file` | any | yes |  |
| `interval` | any | yes |  |
| `samples` | array \| null | no | Default `None` |
| `include_genotypes` | boolean | no | Default `True` |
| `pass_only` | boolean | no | Default `False` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Variant records and optional genotypes overlapping an interval (VCF/BCF).
    Available for: bcf, vcf.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'get_variants',
    {'file': '...', 'interval': '...'},
)
```

## `inspect_locus` {#inspect_locus}

Combines data from several files at one locus with optional public reference evidence.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `interval`
- **Keywords:** locus, multi-file, interval, annotation, reference sources
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `interval` | any | yes |  |
| `files` | array \| null | no | Default `None` |
| `reference` | ? \| null | no | Default `None` |
| `reference_sources` | array \| null | no | Source names; default all available. Default `None` |
| `allow_external_annotation` | boolean | no | Default `False` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |

??? note "Upstream description"

    Combine data from several files at one locus, plus public evidence. Values derived from private files are only sent to external sources if allow_external_annotation is true.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'inspect_locus',
    {'interval': '...'},
)
```

## `list_files` {#list_files}

Lists files of a study/dataset (EGA, ENA, ENCODE, GEO, NCBI Datasets, local, S3) with format, size and index status.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `source`, `accession`
- **Keywords:** file listing, ena, encode, geo, checksum, fastq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes |  |
| `accession` | string | yes |  |
| `formats` | array \| null | no | Default `None` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |
| `cursor` | string \| null | no | Default `None` |
| `storage_profile` | string \| null | no | Default `None` |

??? note "Upstream description"

    List files of a study/dataset with format, size, checksums, index and readiness.
    Available for: ega, ena, encode, geo, local, ncbi_datasets, s3.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'list_files',
    {'source': '...', 'accession': '...'},
)
```

## `list_samples` {#list_samples}

Lists samples of a study/dataset as supplied by EGA, ENA, ENCODE, GEO or NCBI Datasets.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `source`, `accession`
- **Keywords:** samples, biosample, ena, geo, encode
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes |  |
| `accession` | string | yes |  |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |
| `cursor` | string \| null | no | Default `None` |

??? note "Upstream description"

    List samples of a study/dataset as supplied by the source.
    Available for: ega, ena, encode, geo, ncbi_datasets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'list_samples',
    {'source': '...', 'accession': '...'},
)
```

## `lookup_gene` {#lookup_gene}

Returns gene identifiers, transcripts and source-attributed annotations for a gene on a given assembly.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene`
- **Keywords:** gene identifiers, transcripts, ensembl, assembly, annotation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string | yes |  |
| `assembly` | string \| null | no | Default `None` |
| `sources` | array \| null | no | Source names; default all available. Default `None` |
| `include` | array \| null | no | Default `None` |

??? note "Upstream description"

    Gene identifiers, transcripts and source-attributed annotations.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'lookup_gene',
    {'gene': '...'},
)
```

## `lookup_variant` {#lookup_variant}

Returns source-attributed evidence (ClinVar, gnomAD, Ensembl) for a variant given as VCF-style, HGVS or rsID.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** clinvar, gnomad, ensembl vep, hgvs, rsid
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | ? \| null | no | Default `None` |
| `hgvs` | string \| null | no | Default `None` |
| `rsid` | string \| null | no | Default `None` |
| `assembly` | string \| null | no | Default `None` |
| `sources` | array \| null | no | Source names; default all available. Default `None` |
| `include` | array \| null | no | Default `None` |

??? note "Upstream description"

    Source-attributed evidence for a variant (ClinVar, gnomAD, Ensembl, optional Atlas). No consensus or clinical verdict.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'lookup_variant',
    {},
)
```

## `normalize_variant` {#normalize_variant}

Normalizes a VCF-style, HGVS or rsID variant on an explicit assembly with a trace.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** none
- **Keywords:** variant normalization, hgvs, rsid, left alignment, assembly
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | ? \| null | no | Default `None` |
| `hgvs` | string \| null | no | Default `None` |
| `rsid` | string \| null | no | Default `None` |
| `assembly` | string \| null | no | Default `None` |
| `reference` | ? \| null | no | Default `None` |
| `sources` | array \| null | no | Source names; default all available. Default `None` |

??? note "Upstream description"

    Normalize a variant (VCF-style, HGVS or rsID) on an explicit assembly with a trace.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'normalize_variant',
    {},
)
```

## `resolve_identifier` {#resolve_identifier}

Resolves a gene, transcript, protein or variant identifier and reports ambiguity.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `identifier`
- **Keywords:** id mapping, gene symbol, transcript, ensembl, ambiguity
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `identifier` | string | yes |  |
| `target_types` | array \| null | no | Default `None` |
| `assembly` | string \| null | no | Default `None` |
| `sources` | array \| null | no | Source names; default all available. Default `None` |

??? note "Upstream description"

    Resolve a gene/transcript/protein/variant identifier; reports ambiguity.
    Available for: default.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'resolve_identifier',
    {'identifier': '...'},
)
```

## `search_datasets` {#search_datasets}

Searches one archive (EGA, ENA, ENCODE, GEO, NCBI Datasets) for studies or datasets by text, assembly and organism.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `source`, `query`
- **Keywords:** ena, geo, encode, ega, dataset search, study
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `source` | string | yes |  |
| `query` | string | yes |  |
| `assembly` | string \| null | no | Default `None` |
| `organism` | string \| null | no | Default `None` |
| `max_records` | integer \| null | no | Lower the record limit. Default `None` |
| `cursor` | string \| null | no | Default `None` |

??? note "Upstream description"

    Search one archive/catalog for studies or datasets by text query.
    Available for: ega, ena, encode, geo, ncbi_datasets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'search_datasets',
    {'source': '...', 'query': '...'},
)
```
