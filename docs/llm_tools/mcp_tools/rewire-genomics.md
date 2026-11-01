# Rewire Genomics MCP (`rewire-genomics`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **20** 个 DNA/基因相关工具（服务器共列出 24 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [rewire-bio/genomics-mcp @ `c419069`](https://github.com/rewire-bio/genomics-mcp/tree/c4190699201b8a2f447dc52075c5977d7b32e287) |
| 许可证 | MIT |
| 启动命令 | `genomics-mcp` |
| 传输方式 | stdio |
| 安装 | Install the pinned repository checkout into a Python &gt;=3.12,&lt;3.13 environment (pip install of the checkout, with --no-binary pybigwig; needs a C compiler plus libcurl and zlib development files for remote bigWig/bigBed); the upstream docs say there is no PyPI package although PyPI lists rewire-genomics-mcp 0.1.0 by the same author, so prefer the checkout. Optionally set GENOMICS_MCP_ALLOWED_ROOTS to the folders the server may read, and EGA credentials only for controlled EGA data. |
| 副作用 | Calls public services (EGA, ENA, ENCODE, GEO, NCBI, Ensembl, HGNC, ClinVar, gnomAD, UniProt, Open Targets) over HTTPS and reads local files under the allowed roots; fetch_file, cancel_transfer and convert_sra_run write downloaded files into the server work directory (default under the XDG cache), values from private files go to external APIs only when a call sets allow_external_annotation. |

| 工具 | 类别 | 简介 |
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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`interval`, `files`
- **关键词：**genotype comparison, coverage, interval, multi-sample, vcf, bam
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `interval` | any | 是 |  |
| `files` | array | 是 |  |
| `samples` | array \| null | 否 | 默认 `None` |
| `reference` | ? \| null | 否 | 默认 `None` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Compare genotypes/coverage at an interval across files or samples.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`accession`
- **关键词：**sra, fastq, fasterq-dump, prefetch, sequencing reads, srr

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | Public SRR, ERR or DRR run accession. |
| `budget_bytes` | integer \| null | 否 | 默认 `None` |
| `timeout_s` | number \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Optional: download a public SRA run with SRA Toolkit prefetch and convert the whole run to local FASTQ files with fasterq-dump (paired and unpaired reads kept). Starts a job; poll get_transfer_status. Needs SRA Toolkit; ENA FASTQ files via list_files/fetch_file do not. FASTQ is not locus-ready.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`source`, `accession`
- **关键词：**ena, geo, encode, ega, study accession, dataset
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 是 |  |
| `accession` | string | 是 |  |

??? note "上游说明（英文原文）"

    Describe a study/dataset by native accession, with linked entities.
    Available for: ega, ena, encode, geo, ncbi_datasets.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`file`
- **关键词：**download, bam, vcf, index, checksum, transfer

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `include_index` | boolean | 否 | 默认 `True` |
| `budget_bytes` | integer \| null | 否 | 默认 `None` |
| `verify_checksum` | boolean | 否 | 默认 `True` |
| `prepare` | boolean | 否 | 默认 `False` |

??? note "上游说明（英文原文）"

    Start a bounded download of a file (and index) to the local work dir. Returns a transfer job; files above the default limit need budget_bytes.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`file`, `interval`
- **关键词：**read depth, bam, cram, samtools depth, coverage, mapping quality
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `interval` | any | 是 |  |
| `reference` | ? \| null | 否 | 默认 `None` |
| `min_mapping_quality` | integer | 否 | 默认 `0` |
| `min_base_quality` | integer | 否 | 默认 `0` |
| `exclude_flags` | integer | 否 | 默认 `1796` |
| `bin_size` | integer \| null | 否 | 默认 `None` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Read depth over an interval from indexed BAM/CRAM (default excludes UNMAP, SECONDARY, QCFAIL, DUP, like samtools depth).
    Available for: bam, cram.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`file`, `interval`
- **关键词：**bed, gff3, gtf, bigbed, interval, gene features
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `interval` | any | 是 |  |
| `feature_types` | array \| null | 否 | 默认 `None` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Features overlapping an interval from indexed BED/GFF3/GTF or bigBed.
    Available for: bed, bigbed, gff3, gtf.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`file`, `interval`
- **关键词：**pileup, base counts, bam, cram, mpileup
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `interval` | any | 是 |  |
| `reference` | ? \| null | 否 | 默认 `None` |
| `min_mapping_quality` | integer | 否 | 默认 `0` |
| `min_base_quality` | integer | 否 | 默认 `13` |
| `exclude_flags` | integer | 否 | 默认 `1796` |
| `max_depth` | integer | 否 | 默认 `8000` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Per-position base counts over an interval from indexed BAM/CRAM.
    Available for: bam, cram.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`file`, `interval`
- **关键词：**alignments, bam, cram, reads, samtools, flags
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `interval` | any | 是 |  |
| `reference` | ? \| null | 否 | 默认 `None` |
| `min_mapping_quality` | integer | 否 | 默认 `0` |
| `require_flags` | integer | 否 | 默认 `0` |
| `exclude_flags` | integer | 否 | 默认 `0` |
| `include_sequence` | boolean | 否 | 默认 `False` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Alignments overlapping an interval from indexed BAM/CRAM. Flag filters follow samtools -f/-F. CRAM needs a matching reference unless self-contained.
    Available for: bam, cram.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`source`, `accession`
- **关键词：**sample metadata, biosample, phenotype, ena, geo
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 是 |  |
| `accession` | string | 是 |  |

??? note "上游说明（英文原文）"

    Get one sample's metadata and phenotype values exactly as supplied by the source.
    Available for: ega, ena, encode, geo, ncbi_datasets.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`file`, `interval`
- **关键词：**fasta, reference sequence, interval, faidx
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `interval` | any | 是 |  |

??? note "上游说明（英文原文）"

    Reference sequence for an interval from an indexed FASTA.
    Available for: fasta.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`file`, `interval`
- **关键词：**bigwig, signal track, bins, chip-seq, atac-seq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `interval` | any | 是 |  |
| `bins` | integer \| null | 否 | 默认 `None` |
| `summary` | string | 否 | 可选值：`mean`, `min`, `max`, `coverage`, `std`, `sum` 默认 `mean` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Signal values or binned summaries over an interval from bigWig.
    Available for: bigwig.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`file`, `interval`
- **关键词：**vcf, bcf, genotypes, interval, pass filter
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file` | any | 是 |  |
| `interval` | any | 是 |  |
| `samples` | array \| null | 否 | 默认 `None` |
| `include_genotypes` | boolean | 否 | 默认 `True` |
| `pass_only` | boolean | 否 | 默认 `False` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Variant records and optional genotypes overlapping an interval (VCF/BCF).
    Available for: bcf, vcf.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`interval`
- **关键词：**locus, multi-file, interval, annotation, reference sources
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `interval` | any | 是 |  |
| `files` | array \| null | 否 | 默认 `None` |
| `reference` | ? \| null | 否 | 默认 `None` |
| `reference_sources` | array \| null | 否 | Source names; default all available. 默认 `None` |
| `allow_external_annotation` | boolean | 否 | 默认 `False` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |

??? note "上游说明（英文原文）"

    Combine data from several files at one locus, plus public evidence. Values derived from private files are only sent to external sources if allow_external_annotation is true.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`source`, `accession`
- **关键词：**file listing, ena, encode, geo, checksum, fastq
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 是 |  |
| `accession` | string | 是 |  |
| `formats` | array \| null | 否 | 默认 `None` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |
| `cursor` | string \| null | 否 | 默认 `None` |
| `storage_profile` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    List files of a study/dataset with format, size, checksums, index and readiness.
    Available for: ega, ena, encode, geo, local, ncbi_datasets, s3.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`source`, `accession`
- **关键词：**samples, biosample, ena, geo, encode
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 是 |  |
| `accession` | string | 是 |  |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |
| `cursor` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    List samples of a study/dataset as supplied by the source.
    Available for: ega, ena, encode, geo, ncbi_datasets.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene`
- **关键词：**gene identifiers, transcripts, ensembl, assembly, annotation
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 |  |
| `assembly` | string \| null | 否 | 默认 `None` |
| `sources` | array \| null | 否 | Source names; default all available. 默认 `None` |
| `include` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Gene identifiers, transcripts and source-attributed annotations.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**clinvar, gnomad, ensembl vep, hgvs, rsid
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | ? \| null | 否 | 默认 `None` |
| `hgvs` | string \| null | 否 | 默认 `None` |
| `rsid` | string \| null | 否 | 默认 `None` |
| `assembly` | string \| null | 否 | 默认 `None` |
| `sources` | array \| null | 否 | Source names; default all available. 默认 `None` |
| `include` | array \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Source-attributed evidence for a variant (ClinVar, gnomAD, Ensembl, optional Atlas). No consensus or clinical verdict.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**无
- **关键词：**variant normalization, hgvs, rsid, left alignment, assembly
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | ? \| null | 否 | 默认 `None` |
| `hgvs` | string \| null | 否 | 默认 `None` |
| `rsid` | string \| null | 否 | 默认 `None` |
| `assembly` | string \| null | 否 | 默认 `None` |
| `reference` | ? \| null | 否 | 默认 `None` |
| `sources` | array \| null | 否 | Source names; default all available. 默认 `None` |

??? note "上游说明（英文原文）"

    Normalize a variant (VCF-style, HGVS or rsID) on an explicit assembly with a trace.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`identifier`
- **关键词：**id mapping, gene symbol, transcript, ensembl, ambiguity
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `identifier` | string | 是 |  |
| `target_types` | array \| null | 否 | 默认 `None` |
| `assembly` | string \| null | 否 | 默认 `None` |
| `sources` | array \| null | 否 | Source names; default all available. 默认 `None` |

??? note "上游说明（英文原文）"

    Resolve a gene/transcript/protein/variant identifier; reports ambiguity.
    Available for: default.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`source`, `query`
- **关键词：**ena, geo, encode, ega, dataset search, study
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `source` | string | 是 |  |
| `query` | string | 是 |  |
| `assembly` | string \| null | 否 | 默认 `None` |
| `organism` | string \| null | 否 | 默认 `None` |
| `max_records` | integer \| null | 否 | Lower the record limit. 默认 `None` |
| `cursor` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Search one archive/catalog for studies or datasets by text query.
    Available for: ega, ena, encode, geo, ncbi_datasets.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'rewire-genomics',
    'search_datasets',
    {'source': '...', 'query': '...'},
)
```
