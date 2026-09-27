# Download: parameters and notes

See the [tutorial](../../api/features/15_download.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

ENCODE returns text: mostly JSON, but get_citations can return BibTeX/RIS and export_data can return CSV/TSV. Use json.loads() only for JSON. Tracking and credentials use ~/.dnakit-encode/; persistent credential storage needs dnakit[encode-credentials].

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

## Download a file over HTTP or HTTPS {#dnakit.databases.download_file-parameters}

`dnakit.databases.download_file(*, url, output_path=None, timeout=None, return_content=None, chunk_size=None, follow_redirects=None, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `url` | string | Required | HTTP or HTTPS URL to download from (e.g., https://example.com/file.txt) {&quot;format&quot;: &quot;uri&quot;} |
| `output_path` | string | None | Optional path to save the file. If not specified, file will be saved to system temp directory. {&quot;example&quot;: &quot;/tmp/downloaded_file.txt or C:\\Users\\Downloads\\file.txt&quot;} |
| `timeout` | integer | None | Request timeout in seconds {&quot;minimum&quot;: 1, &quot;maximum&quot;: 300, &quot;default&quot;: 30} |
| `return_content` | boolean | None | If true, return file content as text instead of saving to disk (default: false) {&quot;default&quot;: false} |
| `chunk_size` | integer | None | Download chunk size in bytes (default: 8192) {&quot;minimum&quot;: 1024, &quot;maximum&quot;: 10485760, &quot;default&quot;: 8192} |
| `follow_redirects` | boolean | None | Follow HTTP redirects (default: true) {&quot;default&quot;: true} |
| `verbose` | bool | True | Progress flag. |


## Download an ENA nucleotide sequence with a checksum manifest {#dnakit.download.ena.sequence-parameters}

`dnakit.download.ena.sequence(accession, output_path, *, format='fasta', config=None, search_config=None, api_base_url='https://www.ebi.ac.uk/ena/browser/api', verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | Versioned sequence accession; ENA uses INSDC, NCBI also accepts RefSeq. |
| `output_path` | str &#124; Path | Required | Destination sequence file path. |
| `format` | Literal&#91;&#x27;fasta&#x27;, &#x27;embl&#x27;, &#x27;genbank&#x27;&#93; | &#x27;fasta&#x27; | fasta, genbank or embl; see the download guide for provider and conversion support. |
| `config` | DownloadConfig &#124; None | None | SearchConfig for queries; DownloadConfig for downloads. |
| `search_config` | SearchConfig &#124; None | None | SearchConfig for query timeout, response size and NCBI credentials. |
| `api_base_url` | str | &#x27;https://www.ebi.ac.uk/ena/browser/api&#x27; | Service API base URL; see the signature for its default. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Download an NCBI nucleotide or protein sequence {#dnakit.download.ncbi.sequence-parameters}

`dnakit.download.ncbi.sequence(accession, output_path, *, database='nuccore', format='fasta', seq_start=None, seq_stop=None, strand=1, config=None, search_config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | Versioned sequence accession; ENA uses INSDC, NCBI also accepts RefSeq. |
| `output_path` | str &#124; Path | Required | Destination sequence file path. |
| `database` | Literal&#91;&#x27;nuccore&#x27;, &#x27;protein&#x27;&#93; | &#x27;nuccore&#x27; | NCBI sequence database: nuccore or protein. |
| `format` | Literal&#91;&#x27;fasta&#x27;, &#x27;genbank&#x27;, &#x27;embl&#x27;&#93; | &#x27;fasta&#x27; | fasta, genbank or embl; see the download guide for provider and conversion support. |
| `seq_start` | int &#124; None | None | NCBI nucleotide start, 1-based inclusive; requires seq_stop. |
| `seq_stop` | int &#124; None | None | NCBI nucleotide end, 1-based inclusive; requires seq_start. |
| `strand` | Literal&#91;1, 2&#93; | 1 | NCBI nucleotide strand: 1 forward, 2 reverse complement. |
| `config` | DownloadConfig &#124; None | None | SearchConfig for queries; DownloadConfig for downloads. |
| `search_config` | SearchConfig &#124; None | None | SearchConfig for query timeout, response size and NCBI credentials. |
| `api_base_url` | str | &#x27;https://eutils.ncbi.nlm.nih.gov/entrez/eutils&#x27; | Service API base URL; see the signature for its default. |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Find and batch-download ENCODE files {#dnakit.encode.batch_download-parameters}

`dnakit.encode.batch_download(download_dir, file_format=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, target=None, preferred_default=None, organize_by='experiment', verify_md5=True, limit=100, dry_run=True, *, verbose=True)`

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


## Download ENCODE files by accession {#dnakit.encode.download_files-parameters}

`dnakit.encode.download_files(file_accessions, download_dir, organize_by='flat', verify_md5=True, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file_accessions` | list&#91;str&#93; | Required | List of file accessions to download (e.g., &#91;&quot;ENCFF635JIA&quot;, &quot;ENCFF388RZD&quot;&#93;) |
| `download_dir` | str | Required | Local directory path to save files (e.g., &quot;./data/encode&quot;) |
| `organize_by` | Literal&#91;&#x27;flat&#x27;, &#x27;experiment&#x27;, &#x27;format&#x27;, &#x27;experiment_format&#x27;&#93; | &#x27;flat&#x27; | How to organize downloaded files: - &quot;flat&quot;: All files in download_dir (default) - &quot;experiment&quot;: download_dir/ENCSR.../filename - &quot;format&quot;: download_dir/bed/filename - &quot;experiment_format&quot;: download_dir/ENCSR.../bed/filename |
| `verify_md5` | bool | True | Verify file integrity with MD5 checksum (default True) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |
