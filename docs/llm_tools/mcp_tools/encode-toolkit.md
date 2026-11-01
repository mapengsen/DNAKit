# ENCODE Toolkit (`encode-toolkit`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **9** 个 DNA/基因相关工具（服务器共列出 20 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [ammawla/encode-toolkit @ `1c70476`](https://github.com/ammawla/encode-toolkit/tree/1c7047617d112b45db1d33927d9cc11902dc38d2) |
| 许可证 | AGPL-3.0 |
| 启动命令 | `python -m encode_connector.server.main` |
| 传输方式 | stdio |
| 安装 | Install the Python encode-toolkit distribution separately and set executable_path to its interpreter. DNAKit invokes the verified Python module directly; it does not use the npm shim which can bootstrap uv/Python. Pin 'mcp&lt;2': under mcp 2.x the server fails to import. |
| 副作用 | Queries ENCODE experiments and files; selected upstream tools may download data, create workflows or manage local credential state. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`encode_batch_download`](#encode_batch_download) | `sequencing_data` | Searches ENCODE files by assay, organ, target, assembly and format, then downloads them in batch (dry-run preview by default). |
| [`encode_compare_experiments`](#encode_compare_experiments) | `sequencing_data` | Compares two tracked ENCODE experiments on organism, assembly, assay, biosample, target and replication to report compatibility issues. |
| [`encode_download_files`](#encode_download_files) | `sequencing_data` | Downloads specific ENCODE files by accession to a local directory with optional MD5 verification and folder organization. |
| [`encode_get_experiment`](#encode_get_experiment) | `sequencing_data` | Retrieves full ENCODE experiment details by accession, including files, quality metrics, controls, replicates and audit status. |
| [`encode_get_facets`](#encode_get_facets) | `sequencing_data` | Returns live faceted counts of ENCODE experiments or files per filter value to explore available data. |
| [`encode_get_file_info`](#encode_get_file_info) | `sequencing_data` | Retrieves detailed metadata for an ENCODE file accession, including format, size, download URL, MD5 and assembly. |
| [`encode_list_files`](#encode_list_files) | `sequencing_data` | Lists files of an ENCODE experiment with filters for format, output type, category, assembly and status. |
| [`encode_search_experiments`](#encode_search_experiments) | `sequencing_data` | Searches ENCODE experiments by assay, organism, organ, biosample, target, lab and assembly or free text. |
| [`encode_search_files`](#encode_search_files) | `sequencing_data` | Searches files across all ENCODE experiments using combined file and experiment filters (format, output type, assay, target). |

## `encode_batch_download` {#encode_batch_download}

Searches ENCODE files by assay, organ, target, assembly and format, then downloads them in batch (dry-run preview by default).

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`download_dir`
- **关键词：**fastq, bigwig, peaks, chip-seq, bulk download, assembly

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `download_dir` | string | 是 |  |
| `file_format` | string \| null | 否 | 默认 `None` |
| `output_type` | string \| null | 否 | 默认 `None` |
| `output_category` | string \| null | 否 | 默认 `None` |
| `assembly` | string \| null | 否 | 默认 `None` |
| `assay_title` | string \| null | 否 | 默认 `None` |
| `organism` | string | 否 | 默认 `Homo sapiens` |
| `organ` | string \| null | 否 | 默认 `None` |
| `biosample_type` | string \| null | 否 | 默认 `None` |
| `target` | string \| null | 否 | 默认 `None` |
| `preferred_default` | boolean \| null | 否 | 默认 `None` |
| `organize_by` | string | 否 | 可选值：`flat`, `experiment`, `format`, `experiment_format` 默认 `experiment` |
| `verify_md5` | boolean | 否 | 默认 `True` |
| `limit` | integer | 否 | 默认 `100` |
| `dry_run` | boolean | 否 | 默认 `True` |

??? note "上游说明（英文原文）"

    Search for files and download them all in batch.

    First searches for files matching the criteria, then downloads them.
    By default runs in dry_run mode to preview what would be downloaded.
    Set dry_run=False to actually download.

    WHEN TO USE: Use for searching and downloading files in one step. Always use
    dry_run=True first to preview. For specific file accessions, use encode_download_files.
    RELATED TOOLS: encode_download_files, encode_search_files

    Examples:
    - Download all BED files from human pancreas ChIP-seq:
    file_format="bed", assay_title="Histone ChIP-seq", organ="pancreas",
    download_dir="/data/encode", dry_run=False
    - Preview FASTQ downloads for mouse brain RNA-seq:
    file_format="fastq", assay_title="total RNA-seq", organ="brain",
    organism="Mus musculus", download_dir="/data/encode"
    - Download IDR peaks for H3K27me3 in GRCh38:
    output_type="IDR thresholded peaks", target="H3K27me3", assembly="GRCh38",
    download_dir="/data/encode", dry_run=False

    Args:
    download_dir: Local directory to save files
    file_format: File format filter ("fastq", "bam", "bed", "bigWig", etc.)
    output_type: Output type filter ("reads", "peaks", "signal", etc.)
    output_category: Output category ("raw data", "alignment", "annotation", etc.)
    assembly: Genome assembly ("GRCh38", "mm10", etc.)
    assay_title: Assay type ("Histone ChIP-seq", "ATAC-seq", "total RNA-seq", etc.)
    organism: Organism (default: "Homo sapiens")
    organ: Organ/tissue ("pancreas", "brain", "liver", etc.)
    biosample_type: Biosample type ("tissue", "cell line", "primary cell", etc.)
    target: ChIP/CUT&amp;RUN target ("H3K27me3", "CTCF", etc.)
    preferred_default: If True, only download default/recommended files
    organize_by: File organization ("flat", "experiment", "format", "experiment_format")
    verify_md5: Verify downloads with MD5 checksums (default True)
    limit: Max files to download (default 100, safety limit)
    dry_run: If True (default), only preview what would be downloaded. Set False to download.

    Returns:
    JSON with download preview (dry_run=True) or download results (dry_run=False).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_batch_download',
    {'download_dir': '...'},
)
```

## `encode_compare_experiments` {#encode_compare_experiments}

Compares two tracked ENCODE experiments on organism, assembly, assay, biosample, target and replication to report compatibility issues.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`accession1`, `accession2`
- **关键词：**compatibility, assembly, replicates, combined analysis, assay
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession1` | string | 是 |  |
| `accession2` | string | 是 |  |

??? note "上游说明（英文原文）"

    Analyze whether two ENCODE experiments are compatible for combined analysis.

    Compares organism, genome assembly, assay type, biosample, organ, target,
    replication strategy, and lab to identify potential issues.

    Both experiments must be tracked first (use encode_track_experiment).

    WHEN TO USE: Use to check if two experiments are compatible for combined analysis.
    Both must be tracked first.
    RELATED TOOLS: encode_track_experiment, encode_list_tracked

    Args:
    accession1: First experiment accession (e.g., "ENCSR133RZO")
    accession2: Second experiment accession (e.g., "ENCSR000AKS")

    Returns:
    JSON compatibility report with verdict, issues, warnings, and recommendations.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_compare_experiments',
    {'accession1': '...', 'accession2': '...'},
)
```

## `encode_download_files` {#encode_download_files}

Downloads specific ENCODE files by accession to a local directory with optional MD5 verification and folder organization.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`file_accessions`, `download_dir`
- **关键词：**encff, md5, download, fastq, bam

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_accessions` | array of string | 是 |  |
| `download_dir` | string | 是 |  |
| `organize_by` | string | 否 | 可选值：`flat`, `experiment`, `format`, `experiment_format` 默认 `flat` |
| `verify_md5` | boolean | 否 | 默认 `True` |

??? note "上游说明（英文原文）"

    Download specific ENCODE files by accession to a local directory.

    Downloads files from ENCODE to your local machine. Supports MD5 verification,
    concurrent downloads, and skip-if-already-downloaded.

    WHEN TO USE: Use for downloading specific files by accession. For bulk downloads,
    prefer encode_batch_download.
    RELATED TOOLS: encode_batch_download, encode_search_files, encode_log_derived_file

    Args:
    file_accessions: List of file accessions to download (e.g., \["ENCFF635JIA", "ENCFF388RZD"\])
    download_dir: Local directory path to save files (e.g., "./data/encode")
    organize_by: How to organize downloaded files:
    - "flat": All files in download_dir (default)
    - "experiment": download_dir/ENCSR.../filename
    - "format": download_dir/bed/filename
    - "experiment_format": download_dir/ENCSR.../bed/filename
    verify_md5: Verify file integrity with MD5 checksum (default True)

    Returns:
    JSON with download results for each file (path, size, success/error, MD5 status).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_download_files',
    {'file_accessions': [], 'download_dir': '...'},
)
```

## `encode_get_experiment` {#encode_get_experiment}

Retrieves full ENCODE experiment details by accession, including files, quality metrics, controls, replicates and audit status.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`accession`
- **关键词：**encsr, quality metrics, replicates, audit, experiment metadata
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 |  |

??? note "上游说明（英文原文）"

    Get full details for a specific ENCODE experiment by accession ID.

    Returns complete experiment metadata including all associated files,
    quality metrics, controls, replicate information, and audit status.

    WHEN TO USE: Use when you have a specific accession and need full details
    including files, quality metrics, and audit status.
    RELATED TOOLS: encode_list_files, encode_track_experiment, encode_compare_experiments

    Args:
    accession: ENCODE experiment accession (e.g., "ENCSR133RZO", "ENCSR000AKS")

    Returns:
    JSON with full experiment details and file listing.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_get_experiment',
    {'accession': '...'},
)
```

## `encode_get_facets` {#encode_get_facets}

Returns live faceted counts of ENCODE experiments or files per filter value to explore available data.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**无
- **关键词：**facets, counts, assay, organ, data discovery
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `search_type` | string | 否 | 默认 `Experiment` |
| `assay_title` | string \| null | 否 | 默认 `None` |
| `organism` | string \| null | 否 | 默认 `None` |
| `organ` | string \| null | 否 | 默认 `None` |
| `biosample_type` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Get live filter counts from ENCODE to discover what data is available.

    Returns faceted counts showing how many experiments/files exist for each
    filter value. Useful for exploring what's available before searching.

    WHEN TO USE: Use to explore what data exists before searching. Shows counts
    per filter value. Best first step for unknown datasets.
    RELATED TOOLS: encode_get_metadata, encode_search_experiments

    Examples:
    - What assays are available for pancreas?
    organ="pancreas"
    - What organs have Histone ChIP-seq data?
    assay_title="Histone ChIP-seq"
    - What targets are available for mouse brain ChIP-seq?
    assay_title="Histone ChIP-seq", organism="Mus musculus", organ="brain"

    Args:
    search_type: Object type ("Experiment" or "File")
    assay_title: Pre-filter by assay type
    organism: Pre-filter by organism
    organ: Pre-filter by organ
    biosample_type: Pre-filter by biosample type

    Returns:
    JSON with facet names and their term counts.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_get_facets',
    {},
)
```

## `encode_get_file_info` {#encode_get_file_info}

Retrieves detailed metadata for an ENCODE file accession, including format, size, download URL, MD5 and assembly.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`accession`
- **关键词：**encff, file metadata, assembly, md5, download url
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 |  |

??? note "上游说明（英文原文）"

    Get detailed information about a specific ENCODE file.

    WHEN TO USE: Use when you need detailed metadata for a specific file
    (size, md5, assembly, biological replicate info).
    RELATED TOOLS: encode_download_files, encode_list_files

    Args:
    accession: File accession ID (e.g., "ENCFF635JIA")

    Returns:
    JSON with file metadata including format, size, download URL, MD5, assembly, etc.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_get_file_info',
    {'accession': '...'},
)
```

## `encode_list_files` {#encode_list_files}

Lists files of an ENCODE experiment with filters for format, output type, category, assembly and status.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`experiment_accession`
- **关键词：**fastq, bam, bed, bigwig, peaks, gene quantifications
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `experiment_accession` | string | 是 |  |
| `file_format` | string \| null | 否 | 默认 `None` |
| `file_type` | string \| null | 否 | 默认 `None` |
| `output_type` | string \| null | 否 | 默认 `None` |
| `output_category` | string \| null | 否 | 默认 `None` |
| `assembly` | string \| null | 否 | 默认 `None` |
| `status` | string \| null | 否 | 默认 `None` |
| `preferred_default` | boolean \| null | 否 | 默认 `None` |
| `limit` | integer | 否 | 默认 `200` |

??? note "上游说明（英文原文）"

    List all files for a specific ENCODE experiment, with optional filters.

    Examples:
    - All BED files: experiment_accession="ENCSR133RZO", file_format="bed"
    - FASTQs only: experiment_accession="ENCSR133RZO", file_format="fastq"
    - Signal tracks: experiment_accession="ENCSR133RZO", output_category="signal"
    - Default/recommended files: preferred_default=True
    - Peaks from GRCh38: file_format="bed", output_type="IDR thresholded peaks", assembly="GRCh38"

    Common file_format values: "fastq", "bam", "bed", "bigWig", "bigBed", "tsv", "hic"

    Common output_type values: "reads", "alignments", "signal of unique reads",
    "signal of all reads", "fold change over control", "IDR thresholded peaks",
    "pseudoreplicated peaks", "replicated peaks", "gene quantifications",
    "transcript quantifications", "contact matrix"

    WHEN TO USE: Use to browse files within a known experiment. Use encode_search_files
    instead to find files across experiments.
    RELATED TOOLS: encode_search_files, encode_get_file_info, encode_download_files

    Args:
    experiment_accession: ENCODE experiment accession (e.g., "ENCSR133RZO")
    file_format: Filter by format ("fastq", "bam", "bed", "bigWig", "bigBed", etc.)
    file_type: Filter by specific type ("bed narrowPeak", "bed broadPeak", etc.)
    output_type: Filter by output type ("reads", "peaks", "signal", etc.)
    output_category: Filter by category ("raw data", "alignment", "signal", "annotation")
    assembly: Filter by genome assembly ("GRCh38", "hg19", "mm10")
    status: Filter by status ("released", "archived", "in progress")
    preferred_default: If True, return only default/recommended files
    limit: Max files to return (default 200)

    Returns:
    JSON list of files with accession, format, size, download URL, and metadata.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_list_files',
    {'experiment_accession': '...'},
)
```

## `encode_search_experiments` {#encode_search_experiments}

Searches ENCODE experiments by assay, organism, organ, biosample, target, lab and assembly or free text.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**无
- **关键词：**chip-seq, atac-seq, rna-seq, histone, biosample, target
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `assay_title` | string \| null | 否 | 默认 `None` |
| `organism` | string | 否 | 默认 `Homo sapiens` |
| `organ` | string \| null | 否 | 默认 `None` |
| `biosample_type` | string \| null | 否 | 默认 `None` |
| `biosample_term_name` | string \| null | 否 | 默认 `None` |
| `target` | string \| null | 否 | 默认 `None` |
| `status` | string | 否 | 默认 `released` |
| `lab` | string \| null | 否 | 默认 `None` |
| `award` | string \| null | 否 | 默认 `None` |
| `assembly` | string \| null | 否 | 默认 `None` |
| `replication_type` | string \| null | 否 | 默认 `None` |
| `life_stage` | string \| null | 否 | 默认 `None` |
| `sex` | string \| null | 否 | 默认 `None` |
| `treatment` | string \| null | 否 | 默认 `None` |
| `genetic_modification` | string \| null | 否 | 默认 `None` |
| `perturbed` | boolean \| null | 否 | 默认 `None` |
| `search_term` | string \| null | 否 | 默认 `None` |
| `date_released_from` | string \| null | 否 | 默认 `None` |
| `date_released_to` | string \| null | 否 | 默认 `None` |
| `limit` | integer | 否 | 默认 `25` |
| `offset` | integer | 否 | 默认 `0` |

??? note "上游说明（英文原文）"

    Search ENCODE experiments with comprehensive filters.

    Examples:
    - Find all Histone ChIP-seq on human pancreas tissue:
    assay_title="Histone ChIP-seq", organ="pancreas", biosample_type="tissue"
    - Find ATAC-seq on human brain:
    assay_title="ATAC-seq", organ="brain"
    - Find RNA-seq on GM12878 cell line:
    assay_title="total RNA-seq", biosample_term_name="GM12878"
    - Find ChIP-seq targeting H3K27me3:
    assay_title="Histone ChIP-seq", target="H3K27me3"
    - Find all mouse liver experiments:
    organism="Mus musculus", organ="liver"
    - Free text search:
    search_term="CRISPR screen pancreatic"

    Common assay_title values: "Histone ChIP-seq", "TF ChIP-seq", "ATAC-seq",
    "DNase-seq", "total RNA-seq", "polyA plus RNA-seq", "WGBS", "intact Hi-C",
    "CUT&amp;RUN", "CUT&amp;Tag", "STARR-seq", "MPRA", "eCLIP", "CRISPR screen"

    Common organ values: "pancreas", "liver", "brain", "heart", "kidney",
    "lung", "intestine", "skin of body", "blood", "spleen", "thymus"

    biosample_type values: "tissue", "cell line", "primary cell",
    "in vitro differentiated cells", "organoid"

    WHEN TO USE: Use as the primary entry point when users want to find experiments.
    Start with encode_get_facets if unsure what filters to use.
    RELATED TOOLS: encode_get_facets, encode_get_metadata, encode_search_files

    Args:
    assay_title: Assay type (e.g., "Histone ChIP-seq", "ATAC-seq", "total RNA-seq")
    organism: Species (default: "Homo sapiens"). Also: "Mus musculus"
    organ: Organ/tissue system (e.g., "pancreas", "brain", "liver")
    biosample_type: Sample classification ("tissue", "cell line", "primary cell", "organoid")
    biosample_term_name: Specific cell/tissue name (e.g., "GM12878", "HepG2", "pancreas")
    target: ChIP/CUT&amp;RUN target (e.g., "H3K27me3", "CTCF", "p300")
    status: Data status (default: "released"). Also: "archived", "revoked"
    lab: Submitting lab name
    award: Funding project
    assembly: Genome assembly (e.g., "GRCh38", "mm10")
    replication_type: "isogenic", "anisogenic", or "unreplicated"
    life_stage: "embryonic", "postnatal", "child", "adult"
    sex: "male", "female", "mixed"
    treatment: Treatment name if perturbation experiment
    genetic_modification: Modification type ("CRISPR", "RNAi")
    perturbed: True for perturbation experiments only
    search_term: Free text search across all fields
    date_released_from: Start date (YYYY-MM-DD) for date range filter
    date_released_to: End date (YYYY-MM-DD) for date range filter
    limit: Max results to return (default 25, use larger for comprehensive searches)
    offset: Skip first N results (for pagination)

    Returns:
    JSON with experiment results, total count, and pagination info.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_search_experiments',
    {},
)
```

## `encode_search_files` {#encode_search_files}

Searches files across all ENCODE experiments using combined file and experiment filters (format, output type, assay, target).

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**无
- **关键词：**fastq, bigwig, peaks, narrowpeak, assembly, cross-experiment
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_format` | string \| null | 否 | 默认 `None` |
| `file_type` | string \| null | 否 | 默认 `None` |
| `output_type` | string \| null | 否 | 默认 `None` |
| `output_category` | string \| null | 否 | 默认 `None` |
| `assembly` | string \| null | 否 | 默认 `None` |
| `assay_title` | string \| null | 否 | 默认 `None` |
| `organism` | string \| null | 否 | 默认 `None` |
| `organ` | string \| null | 否 | 默认 `None` |
| `biosample_type` | string \| null | 否 | 默认 `None` |
| `target` | string \| null | 否 | 默认 `None` |
| `status` | string | 否 | 默认 `released` |
| `preferred_default` | boolean \| null | 否 | 默认 `None` |
| `search_term` | string \| null | 否 | 默认 `None` |
| `limit` | integer | 否 | 默认 `25` |
| `offset` | integer | 否 | 默认 `0` |

??? note "上游说明（英文原文）"

    Search files across ALL experiments with combined experiment + file filters.

    This is powerful for finding specific file types across many experiments.

    Examples:
    - All BED files from human pancreas ChIP-seq:
    file_format="bed", assay_title="Histone ChIP-seq", organ="pancreas"
    - FASTQs from mouse liver RNA-seq:
    file_format="fastq", assay_title="total RNA-seq", organ="liver", organism="Mus musculus"
    - All IDR peak files for H3K27me3:
    output_type="IDR thresholded peaks", target="H3K27me3"
    - BigWig signal tracks from ATAC-seq on brain tissue:
    file_format="bigWig", assay_title="ATAC-seq", organ="brain", biosample_type="tissue"

    WHEN TO USE: Use to find specific file types across ALL experiments. More powerful
    than encode_list_files for cross-experiment file discovery.
    RELATED TOOLS: encode_list_files, encode_batch_download, encode_get_file_info

    Args:
    file_format: File format ("fastq", "bam", "bed", "bigWig", etc.)
    file_type: Specific file type ("bed narrowPeak", "bed broadPeak", etc.)
    output_type: Output type ("reads", "peaks", "signal", etc.)
    output_category: Output category ("raw data", "alignment", "signal", "annotation")
    assembly: Genome assembly ("GRCh38", "hg19", "mm10")
    assay_title: Filter by assay type of parent experiment
    organism: Filter by organism of parent experiment
    organ: Filter by organ of parent experiment
    biosample_type: Filter by biosample type ("tissue", "cell line", etc.)
    target: Filter by ChIP/CUT&amp;RUN target
    status: File status (default: "released")
    preferred_default: If True, only default/recommended files
    search_term: Free text search
    limit: Max results (default 25)
    offset: Skip first N results (pagination)

    Returns:
    JSON with file results, total count, and pagination info.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'encode-toolkit',
    'encode_search_files',
    {},
)
```
