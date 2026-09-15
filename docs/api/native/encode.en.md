# ENCODE queries, downloads and data management

This page documents **20 newly added public functions**. [Complete index](index.md) · [Conventions and examples](usage.md)

Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.

ENCODE returns text: mostly JSON, but get_citations can return BibTeX/RIS and export_data can return CSV/TSV. Use json.loads() only for JSON. Tracking and credentials use ~/.dnakit-encode/; persistent credential storage needs dnakit[encode-credentials]. Names starting with encode_ in source descriptions are upstream names; use this page's dnakit.encode.* functions in DNAKit.

## Function index

- [`dnakit.encode.batch_download`](#dnakit.encode.batch_download)
- [`dnakit.encode.compare_experiments`](#dnakit.encode.compare_experiments)
- [`dnakit.encode.download_files`](#dnakit.encode.download_files)
- [`dnakit.encode.export_data`](#dnakit.encode.export_data)
- [`dnakit.encode.get_citations`](#dnakit.encode.get_citations)
- [`dnakit.encode.get_experiment`](#dnakit.encode.get_experiment)
- [`dnakit.encode.get_facets`](#dnakit.encode.get_facets)
- [`dnakit.encode.get_file_info`](#dnakit.encode.get_file_info)
- [`dnakit.encode.get_metadata`](#dnakit.encode.get_metadata)
- [`dnakit.encode.get_provenance`](#dnakit.encode.get_provenance)
- [`dnakit.encode.get_references`](#dnakit.encode.get_references)
- [`dnakit.encode.link_reference`](#dnakit.encode.link_reference)
- [`dnakit.encode.list_files`](#dnakit.encode.list_files)
- [`dnakit.encode.list_tracked`](#dnakit.encode.list_tracked)
- [`dnakit.encode.log_derived_file`](#dnakit.encode.log_derived_file)
- [`dnakit.encode.manage_credentials`](#dnakit.encode.manage_credentials)
- [`dnakit.encode.search_experiments`](#dnakit.encode.search_experiments)
- [`dnakit.encode.search_files`](#dnakit.encode.search_files)
- [`dnakit.encode.summarize_collection`](#dnakit.encode.summarize_collection)
- [`dnakit.encode.track_experiment`](#dnakit.encode.track_experiment)

## `dnakit.encode.batch_download` {#dnakit.encode.batch_download}

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


```text
dnakit.encode.batch_download(download_dir, file_format=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, target=None, preferred_default=None, organize_by='experiment', verify_md5=True, limit=100, dry_run=True, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `download_dir` | str | Required | Local directory to save files |
| `file_format` | str &#124; None | None | File format filter (&quot;fastq&quot;, &quot;bam&quot;, &quot;bed&quot;, &quot;bigWig&quot;, etc.) |
| `output_type` | str &#124; None | None | Output type filter (&quot;reads&quot;, &quot;peaks&quot;, &quot;signal&quot;, etc.) |
| `output_category` | str &#124; None | None | Output category (&quot;raw data&quot;, &quot;alignment&quot;, &quot;annotation&quot;, etc.) |
| `assembly` | str &#124; None | None | Genome assembly (&quot;GRCh38&quot;, &quot;mm10&quot;, etc.) |
| `assay_title` | str &#124; None | None | Assay type (&quot;Histone ChIP-seq&quot;, &quot;ATAC-seq&quot;, &quot;total RNA-seq&quot;, etc.) |
| `organism` | str | &#x27;Homo sapiens&#x27; | Organism (default: &quot;Homo sapiens&quot;) |
| `organ` | str &#124; None | None | Organ/tissue (&quot;pancreas&quot;, &quot;brain&quot;, &quot;liver&quot;, etc.) |
| `biosample_type` | str &#124; None | None | Biosample type (&quot;tissue&quot;, &quot;cell line&quot;, &quot;primary cell&quot;, etc.) |
| `target` | str &#124; None | None | ChIP/CUT&amp;RUN target (&quot;H3K27me3&quot;, &quot;CTCF&quot;, etc.) |
| `preferred_default` | bool &#124; None | None | If True, only download default/recommended files |
| `organize_by` | Literal&#91;&#x27;flat&#x27;, &#x27;experiment&#x27;, &#x27;format&#x27;, &#x27;experiment_format&#x27;&#93; | &#x27;experiment&#x27; | File organization (&quot;flat&quot;, &quot;experiment&quot;, &quot;format&quot;, &quot;experiment_format&quot;) |
| `verify_md5` | bool | True | Verify downloads with MD5 checksums (default True) |
| `limit` | int | 100 | Max files to download (default 100, safety limit) |
| `dry_run` | bool | True | If True (default), only preview what would be downloaded. Set False to download. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with download preview (dry_run=True) or download results (dry_run=False).

## `dnakit.encode.compare_experiments` {#dnakit.encode.compare_experiments}

Analyze whether two ENCODE experiments are compatible for combined analysis.

Compares organism, genome assembly, assay type, biosample, organ, target,
replication strategy, and lab to identify potential issues.

Both experiments must be tracked first (use encode_track_experiment).

WHEN TO USE: Use to check if two experiments are compatible for combined analysis.
Both must be tracked first.
RELATED TOOLS: encode_track_experiment, encode_list_tracked


```text
dnakit.encode.compare_experiments(accession1, accession2, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession1` | str | Required | First experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `accession2` | str | Required | Second experiment accession (e.g., &quot;ENCSR000AKS&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON compatibility report with verdict, issues, warnings, and recommendations.

## `dnakit.encode.download_files` {#dnakit.encode.download_files}

Download specific ENCODE files by accession to a local directory.

Downloads files from ENCODE to your local machine. Supports MD5 verification,
concurrent downloads, and skip-if-already-downloaded.

WHEN TO USE: Use for downloading specific files by accession. For bulk downloads,
prefer encode_batch_download.
RELATED TOOLS: encode_batch_download, encode_search_files, encode_log_derived_file


```text
dnakit.encode.download_files(file_accessions, download_dir, organize_by='flat', verify_md5=True, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file_accessions` | list&#91;str&#93; | Required | List of file accessions to download (e.g., &#91;&quot;ENCFF635JIA&quot;, &quot;ENCFF388RZD&quot;&#93;) |
| `download_dir` | str | Required | Local directory path to save files (e.g., &quot;./data/encode&quot;) |
| `organize_by` | Literal&#91;&#x27;flat&#x27;, &#x27;experiment&#x27;, &#x27;format&#x27;, &#x27;experiment_format&#x27;&#93; | &#x27;flat&#x27; | How to organize downloaded files: - &quot;flat&quot;: All files in download_dir (default) - &quot;experiment&quot;: download_dir/ENCSR.../filename - &quot;format&quot;: download_dir/bed/filename - &quot;experiment_format&quot;: download_dir/ENCSR.../bed/filename |
| `verify_md5` | bool | True | Verify file integrity with MD5 checksum (default True) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with download results for each file (path, size, success/error, MD5 status).

## `dnakit.encode.export_data` {#dnakit.encode.export_data}

Export tracked experiments as a table (CSV, TSV, or JSON).

Creates a tabular export of all tracked experiments with metadata,
publication counts, PMIDs, and derived file counts. Useful for loading
into Excel, R, pandas, or sharing with collaborators.

PMIDs in the output can be directly used with PubMed MCP tools for
further literature analysis.

WHEN TO USE: Use to create shareable tables of tracked experiments (CSV, TSV, JSON).
Good for manuscripts and reports.
RELATED TOOLS: encode_list_tracked, encode_summarize_collection


```text
dnakit.encode.export_data(format='csv', assay_title=None, organism=None, organ=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `format` | Literal&#91;&#x27;csv&#x27;, &#x27;tsv&#x27;, &#x27;json&#x27;&#93; | &#x27;csv&#x27; | Output format: - &quot;csv&quot;: Comma-separated values (default, for Excel/spreadsheets) - &quot;tsv&quot;: Tab-separated values (for R, pandas) - &quot;json&quot;: JSON array (for programmatic use) |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  Formatted table data in the requested format.

## `dnakit.encode.get_citations` {#dnakit.encode.get_citations}

Get publications and citations for tracked experiments.

Returns publication data with authors, journal, DOI, PMID.
Can export as BibTeX or RIS (Endnote/Zotero/Mendeley compatible).

WHEN TO USE: Use to get publication data for tracked experiments. Supports
BibTeX and RIS export for reference managers.
RELATED TOOLS: encode_track_experiment, encode_link_reference


```text
dnakit.encode.get_citations(accession=None, export_format='json', *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str &#124; None | None | Specific experiment accession. If None, returns all publications. |
| `export_format` | Literal&#91;&#x27;json&#x27;, &#x27;bibtex&#x27;, &#x27;ris&#x27;&#93; | &#x27;json&#x27; | Output format: - &quot;json&quot;: Structured data (default) - &quot;bibtex&quot;: BibTeX format for LaTeX - &quot;ris&quot;: RIS format (Endnote, Zotero, Mendeley) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  Publications in the requested format.

## `dnakit.encode.get_experiment` {#dnakit.encode.get_experiment}

Get full details for a specific ENCODE experiment by accession ID.

Returns complete experiment metadata including all associated files,
quality metrics, controls, replicate information, and audit status.

WHEN TO USE: Use when you have a specific accession and need full details
including files, quality metrics, and audit status.
RELATED TOOLS: encode_list_files, encode_track_experiment, encode_compare_experiments


```text
dnakit.encode.get_experiment(accession, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;, &quot;ENCSR000AKS&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with full experiment details and file listing.

## `dnakit.encode.get_facets` {#dnakit.encode.get_facets}

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


```text
dnakit.encode.get_facets(search_type='Experiment', assay_title=None, organism=None, organ=None, biosample_type=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `search_type` | str | &#x27;Experiment&#x27; | Object type (&quot;Experiment&quot; or &quot;File&quot;) |
| `assay_title` | str &#124; None | None | Pre-filter by assay type |
| `organism` | str &#124; None | None | Pre-filter by organism |
| `organ` | str &#124; None | None | Pre-filter by organ |
| `biosample_type` | str &#124; None | None | Pre-filter by biosample type |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with facet names and their term counts.

## `dnakit.encode.get_file_info` {#dnakit.encode.get_file_info}

Get detailed information about a specific ENCODE file.

WHEN TO USE: Use when you need detailed metadata for a specific file
(size, md5, assembly, biological replicate info).
RELATED TOOLS: encode_download_files, encode_list_files


```text
dnakit.encode.get_file_info(accession, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | File accession ID (e.g., &quot;ENCFF635JIA&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with file metadata including format, size, download URL, MD5, assembly, etc.

## `dnakit.encode.get_metadata` {#dnakit.encode.get_metadata}

Get available filter values for ENCODE searches.

Use this to discover valid values for search parameters.

WHEN TO USE: Use to discover valid filter values before searching. Helps prevent
typos in assay_title, organ, biosample_type etc.
RELATED TOOLS: encode_get_facets, encode_search_experiments


```text
dnakit.encode.get_metadata(metadata_type, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `metadata_type` | Literal&#91;&#x27;assays&#x27;, &#x27;organisms&#x27;, &#x27;organs&#x27;, &#x27;biosample_types&#x27;, &#x27;file_formats&#x27;, &#x27;output_types&#x27;, &#x27;output_categories&#x27;, &#x27;assemblies&#x27;, &#x27;life_stages&#x27;, &#x27;replication_types&#x27;, &#x27;statuses&#x27;, &#x27;file_statuses&#x27;&#93; | Required | Type of metadata to retrieve. Options: - &quot;assays&quot;: Available assay types (Histone ChIP-seq, ATAC-seq, total RNA-seq, etc.) - &quot;organisms&quot;: Available organisms (Homo sapiens, Mus musculus, etc.) - &quot;organs&quot;: Available organ/tissue systems (pancreas, brain, liver, etc.) - &quot;biosample_types&quot;: Biosample classifications (tissue, cell line, primary cell, etc.) - &quot;file_formats&quot;: File format types (fastq, bam, bed, bigWig, etc.) - &quot;output_types&quot;: Output data types (reads, peaks, signal, etc.) - &quot;output_categories&quot;: Output categories (raw data, alignment, signal, etc.) - &quot;assemblies&quot;: Genome assemblies (GRCh38, hg19, mm10, etc.) - &quot;life_stages&quot;: Life stages (embryonic, adult, child, etc.) - &quot;replication_types&quot;: Replication types (isogenic, anisogenic, unreplicated) - &quot;statuses&quot;: Experiment statuses (released, archived, etc.) - &quot;file_statuses&quot;: File statuses (released, archived, in progress, etc.) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON list of valid values for the specified metadata type.

## `dnakit.encode.get_provenance` {#dnakit.encode.get_provenance}

Get provenance information for derived files.

Shows the chain from your derived files back to original ENCODE data,
including what tools and parameters were used.

WHEN TO USE: Use to trace a derived file back to original ENCODE data.
Shows tools and parameters used.
RELATED TOOLS: encode_log_derived_file


```text
dnakit.encode.get_provenance(file_path=None, source_accession=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file_path` | str &#124; None | None | Get provenance for a specific derived file |
| `source_accession` | str &#124; None | None | List all files derived from a specific ENCODE accession |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON provenance chain or list of derived files.

## `dnakit.encode.get_references` {#dnakit.encode.get_references}

Get external references linked to tracked experiments.

Returns PMIDs, DOIs, NCT IDs, GEO accessions and other identifiers
linked to experiments. These identifiers can be directly passed to
PubMed, bioRxiv, ClinicalTrials.gov, or other MCP tools.

WHEN TO USE: Use to retrieve external references linked to experiments.
PMIDs can be passed to PubMed MCP tools.
RELATED TOOLS: encode_link_reference, encode_get_citations


```text
dnakit.encode.get_references(experiment_accession=None, reference_type=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `experiment_accession` | str &#124; None | None | Filter by specific experiment (optional) |
| `reference_type` | Literal&#91;&#x27;pmid&#x27;, &#x27;doi&#x27;, &#x27;nct_id&#x27;, &#x27;preprint_doi&#x27;, &#x27;geo_accession&#x27;, &#x27;other&#x27;&#93; &#124; None | None | Filter by reference type (optional): &quot;pmid&quot;, &quot;doi&quot;, &quot;nct_id&quot;, &quot;preprint_doi&quot;, &quot;geo_accession&quot;, &quot;other&quot; |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with linked external references.

## `dnakit.encode.link_reference` {#dnakit.encode.link_reference}

Link an external reference to a tracked ENCODE experiment.

This is the cross-server bridge. Attach PubMed IDs, bioRxiv DOIs,
ClinicalTrials.gov NCT IDs, GEO accessions, or any external identifier
to your tracked experiments for provenance and cross-referencing.

After finding a relevant paper with PubMed MCP or a preprint on bioRxiv,
link it to the ENCODE experiment for future reference.

WHEN TO USE: Use to attach external IDs (PMID, DOI, GEO, NCT) to tracked
experiments for cross-referencing.
RELATED TOOLS: encode_get_references, encode_get_citations


```text
dnakit.encode.link_reference(experiment_accession, reference_type, reference_id, description='', *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `experiment_accession` | str | Required | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `reference_type` | Literal&#91;&#x27;pmid&#x27;, &#x27;doi&#x27;, &#x27;nct_id&#x27;, &#x27;preprint_doi&#x27;, &#x27;geo_accession&#x27;, &#x27;other&#x27;&#93; | Required | Type of external reference: - &quot;pmid&quot;: PubMed ID (e.g., &quot;32728249&quot;) - &quot;doi&quot;: DOI (e.g., &quot;10.1038/s41586-020-2493-4&quot;) - &quot;nct_id&quot;: ClinicalTrials.gov ID (e.g., &quot;NCT04567890&quot;) - &quot;preprint_doi&quot;: bioRxiv/medRxiv DOI - &quot;geo_accession&quot;: GEO accession (e.g., &quot;GSE123456&quot;) - &quot;other&quot;: Any other identifier |
| `reference_id` | str | Required | The actual identifier value |
| `description` | str | &#x27;&#x27; | Optional description of why this reference is linked |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with linking result.

## `dnakit.encode.list_files` {#dnakit.encode.list_files}

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


```text
dnakit.encode.list_files(experiment_accession, file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, status=None, preferred_default=None, limit=200, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `experiment_accession` | str | Required | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `file_format` | str &#124; None | None | Filter by format (&quot;fastq&quot;, &quot;bam&quot;, &quot;bed&quot;, &quot;bigWig&quot;, &quot;bigBed&quot;, etc.) |
| `file_type` | str &#124; None | None | Filter by specific type (&quot;bed narrowPeak&quot;, &quot;bed broadPeak&quot;, etc.) |
| `output_type` | str &#124; None | None | Filter by output type (&quot;reads&quot;, &quot;peaks&quot;, &quot;signal&quot;, etc.) |
| `output_category` | str &#124; None | None | Filter by category (&quot;raw data&quot;, &quot;alignment&quot;, &quot;signal&quot;, &quot;annotation&quot;) |
| `assembly` | str &#124; None | None | Filter by genome assembly (&quot;GRCh38&quot;, &quot;hg19&quot;, &quot;mm10&quot;) |
| `status` | str &#124; None | None | Filter by status (&quot;released&quot;, &quot;archived&quot;, &quot;in progress&quot;) |
| `preferred_default` | bool &#124; None | None | If True, return only default/recommended files |
| `limit` | int | 200 | Max files to return (default 200) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON list of files with accession, format, size, download URL, and metadata.

## `dnakit.encode.list_tracked` {#dnakit.encode.list_tracked}

List all experiments you've tracked locally, with optional filters.

Shows your local library of tracked ENCODE experiments, their metadata,
publication counts, and derived file counts.

WHEN TO USE: Use to see all experiments in your local library. Filter by assay,
organism, or organ.
RELATED TOOLS: encode_summarize_collection, encode_export_data


```text
dnakit.encode.list_tracked(assay_title=None, organism=None, organ=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with tracked experiments metadata table and tracker stats.

## `dnakit.encode.log_derived_file` {#dnakit.encode.log_derived_file}

Log a file you've derived from ENCODE data for provenance tracking.

Use this when you create new files from ENCODE data (e.g., running a pipeline,
filtering peaks, merging samples). This creates a provenance record linking
your derived file back to the original ENCODE source data.

WHEN TO USE: Use after creating files from ENCODE data (filtered peaks, merged
signals). Creates provenance chain back to source.
RELATED TOOLS: encode_get_provenance, encode_download_files


```text
dnakit.encode.log_derived_file(file_path, source_accessions, description='', file_type='', tool_used='', parameters='', *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file_path` | str | Required | Path to the derived file you created |
| `source_accessions` | list&#91;str&#93; | Required | List of ENCODE accessions this file was derived from (experiment or file accessions, e.g., &#91;&quot;ENCSR133RZO&quot;, &quot;ENCFF635JIA&quot;&#93;) |
| `description` | str | &#x27;&#x27; | What this derived file contains |
| `file_type` | str | &#x27;&#x27; | Type of file (e.g., &quot;filtered_peaks&quot;, &quot;merged_signal&quot;, &quot;differential&quot;) |
| `tool_used` | str | &#x27;&#x27; | Tool/software used to create it (e.g., &quot;bedtools intersect&quot;, &quot;DESeq2&quot;) |
| `parameters` | str | &#x27;&#x27; | Parameters or command used |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with the provenance record ID.

## `dnakit.encode.manage_credentials` {#dnakit.encode.manage_credentials}

Manage ENCODE API credentials for accessing restricted/unreleased data.

Most ENCODE data is public and requires no authentication.
Credentials are only needed for unreleased or restricted datasets.

Credentials are stored securely in your OS keyring (macOS Keychain,
Linux Secret Service, Windows Credential Locker) and never in plaintext.

WHEN TO USE: Use only for accessing unreleased/restricted ENCODE data.
Public data requires no authentication.
RELATED TOOLS: encode_search_experiments


```text
dnakit.encode.manage_credentials(action, access_key=None, secret_key=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `action` | Literal&#91;&#x27;store&#x27;, &#x27;check&#x27;, &#x27;clear&#x27;&#93; | Required | What to do: - &quot;store&quot;: Save new credentials (requires access_key and secret_key) - &quot;check&quot;: Check if credentials are configured - &quot;clear&quot;: Remove stored credentials |
| `access_key` | str &#124; None | None | Your ENCODE access key (only for action=&quot;store&quot;) |
| `secret_key` | str &#124; None | None | Your ENCODE secret key (only for action=&quot;store&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with action result.

## `dnakit.encode.search_experiments` {#dnakit.encode.search_experiments}

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
"CUT&RUN", "CUT&Tag", "STARR-seq", "MPRA", "eCLIP", "CRISPR screen"

Common organ values: "pancreas", "liver", "brain", "heart", "kidney",
"lung", "intestine", "skin of body", "blood", "spleen", "thymus"

biosample_type values: "tissue", "cell line", "primary cell",
"in vitro differentiated cells", "organoid"

WHEN TO USE: Use as the primary entry point when users want to find experiments.
Start with encode_get_facets if unsure what filters to use.
RELATED TOOLS: encode_get_facets, encode_get_metadata, encode_search_files


```text
dnakit.encode.search_experiments(assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, biosample_term_name=None, target=None, status='released', lab=None, award=None, assembly=None, replication_type=None, life_stage=None, sex=None, treatment=None, genetic_modification=None, perturbed=None, search_term=None, date_released_from=None, date_released_to=None, limit=25, offset=0, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `assay_title` | str &#124; None | None | Assay type (e.g., &quot;Histone ChIP-seq&quot;, &quot;ATAC-seq&quot;, &quot;total RNA-seq&quot;) |
| `organism` | str | &#x27;Homo sapiens&#x27; | Species (default: &quot;Homo sapiens&quot;). Also: &quot;Mus musculus&quot; |
| `organ` | str &#124; None | None | Organ/tissue system (e.g., &quot;pancreas&quot;, &quot;brain&quot;, &quot;liver&quot;) |
| `biosample_type` | str &#124; None | None | Sample classification (&quot;tissue&quot;, &quot;cell line&quot;, &quot;primary cell&quot;, &quot;organoid&quot;) |
| `biosample_term_name` | str &#124; None | None | Specific cell/tissue name (e.g., &quot;GM12878&quot;, &quot;HepG2&quot;, &quot;pancreas&quot;) |
| `target` | str &#124; None | None | ChIP/CUT&amp;RUN target (e.g., &quot;H3K27me3&quot;, &quot;CTCF&quot;, &quot;p300&quot;) |
| `status` | str | &#x27;released&#x27; | Data status (default: &quot;released&quot;). Also: &quot;archived&quot;, &quot;revoked&quot; |
| `lab` | str &#124; None | None | Submitting lab name |
| `award` | str &#124; None | None | Funding project |
| `assembly` | str &#124; None | None | Genome assembly (e.g., &quot;GRCh38&quot;, &quot;mm10&quot;) |
| `replication_type` | str &#124; None | None | &quot;isogenic&quot;, &quot;anisogenic&quot;, or &quot;unreplicated&quot; |
| `life_stage` | str &#124; None | None | &quot;embryonic&quot;, &quot;postnatal&quot;, &quot;child&quot;, &quot;adult&quot; |
| `sex` | str &#124; None | None | &quot;male&quot;, &quot;female&quot;, &quot;mixed&quot; |
| `treatment` | str &#124; None | None | Treatment name if perturbation experiment |
| `genetic_modification` | str &#124; None | None | Modification type (&quot;CRISPR&quot;, &quot;RNAi&quot;) |
| `perturbed` | bool &#124; None | None | True for perturbation experiments only |
| `search_term` | str &#124; None | None | Free text search across all fields |
| `date_released_from` | str &#124; None | None | Start date (YYYY-MM-DD) for date range filter |
| `date_released_to` | str &#124; None | None | End date (YYYY-MM-DD) for date range filter |
| `limit` | int | 25 | Max results to return (default 25, use larger for comprehensive searches) |
| `offset` | int | 0 | Skip first N results (for pagination) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with experiment results, total count, and pagination info.

## `dnakit.encode.search_files` {#dnakit.encode.search_files}

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


```text
dnakit.encode.search_files(file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism=None, organ=None, biosample_type=None, target=None, status='released', preferred_default=None, search_term=None, limit=25, offset=0, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file_format` | str &#124; None | None | File format (&quot;fastq&quot;, &quot;bam&quot;, &quot;bed&quot;, &quot;bigWig&quot;, etc.) |
| `file_type` | str &#124; None | None | Specific file type (&quot;bed narrowPeak&quot;, &quot;bed broadPeak&quot;, etc.) |
| `output_type` | str &#124; None | None | Output type (&quot;reads&quot;, &quot;peaks&quot;, &quot;signal&quot;, etc.) |
| `output_category` | str &#124; None | None | Output category (&quot;raw data&quot;, &quot;alignment&quot;, &quot;signal&quot;, &quot;annotation&quot;) |
| `assembly` | str &#124; None | None | Genome assembly (&quot;GRCh38&quot;, &quot;hg19&quot;, &quot;mm10&quot;) |
| `assay_title` | str &#124; None | None | Filter by assay type of parent experiment |
| `organism` | str &#124; None | None | Filter by organism of parent experiment |
| `organ` | str &#124; None | None | Filter by organ of parent experiment |
| `biosample_type` | str &#124; None | None | Filter by biosample type (&quot;tissue&quot;, &quot;cell line&quot;, etc.) |
| `target` | str &#124; None | None | Filter by ChIP/CUT&amp;RUN target |
| `status` | str | &#x27;released&#x27; | File status (default: &quot;released&quot;) |
| `preferred_default` | bool &#124; None | None | If True, only default/recommended files |
| `search_term` | str &#124; None | None | Free text search |
| `limit` | int | 25 | Max results (default 25) |
| `offset` | int | 0 | Skip first N results (pagination) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with file results, total count, and pagination info.

## `dnakit.encode.summarize_collection` {#dnakit.encode.summarize_collection}

Summarize your tracked experiment collection with grouped statistics.

Provides an overview of tracked experiments grouped by assay type, target,
organism, organ, biosample type, and lab. Shows total counts for publications,
derived files, and external references.

Useful when tracking 10+ experiments and needing a bird's-eye view of your
research data collection.

WHEN TO USE: Use for a bird's-eye view of tracked experiments grouped by assay,
target, organ. Best for 10+ tracked experiments.
RELATED TOOLS: encode_list_tracked, encode_export_data


```text
dnakit.encode.summarize_collection(assay_title=None, organism=None, organ=None, *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON summary with experiment counts grouped by multiple dimensions.

## `dnakit.encode.track_experiment` {#dnakit.encode.track_experiment}

Track an ENCODE experiment locally with its publications, methods, and pipeline info.

Fetches full experiment metadata from ENCODE and stores it in a local SQLite
database along with any associated publications (PMIDs, DOIs, authors, journal)
and pipeline/analysis information (software versions, methods).

This is like adding an experiment to your "library" - similar to Endnote for papers.

WHEN TO USE: Use to save an experiment to your local library with publications
and pipeline info. Required before compare or citations.
RELATED TOOLS: encode_compare_experiments, encode_get_citations, encode_export_data


```text
dnakit.encode.track_experiment(accession, fetch_publications=True, fetch_pipelines=True, notes='', *, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `fetch_publications` | bool | True | Also fetch and store publications/citations (default True) |
| `fetch_pipelines` | bool | True | Also fetch and store pipeline/analysis info (default True) |
| `notes` | str | &#x27;&#x27; | Optional notes to attach to this experiment |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |

**Returns**

Python return type: str.  JSON with tracking result including publications and pipeline info found.

