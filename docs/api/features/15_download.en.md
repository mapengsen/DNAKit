# Download

Download sequences, reference genomes and annotations; [FAQ](../../faq/tutorial_integrations.md#downloads) covers formats, retries and species names.

<span id="accession-sequences"></span>By accession: [NCBI](#dnakit.download.ncbi.sequence) · [ENA](#dnakit.download.ena.sequence). Formats and coordinates: [FAQ](../../faq/tutorial_integrations.md#downloads).

## 1) complete genome {#1-dbd-001-complete-genome}

- **Purpose:** Download the complete reference genome FASTA by assembly accession or common species name, and save the assembly report, source URL and verification information for establishing traceable local reference sequences.
- **Method:** Resolve the assembly, download its NCBI package and validate the FASTA and report.
- **API:** `dnakit.references.resolve_genome_assembly(query[required], timeout[optional], api_key[optional], api_base_url[optional])`, `dnakit.references.download_genome(query[required], output_dir[required], overwrite[optional], keep_package[optional], timeout[optional], api_key[optional], api_base_url[optional], chunk_size[optional], max_download_bytes[optional], progress[optional])`.
- **Input:** NCBI versioned accession, fixed alias (such as `hg38`), or species name; requires local output directory.
- **Output:** `resolve_genome_assembly` returns assembly metadata; `download_genome` returns `GenomeAssembly` with the FASTA path and checksums.
- **Example code:**

```python
from pathlib import Path

from dnakit.references import download_genome

result = download_genome("hg38", Path("references/hg38"))
print(result.fasta_path)
print(result.fasta_md5, result.fasta_sha256)
```

- **Example output:**

```text
references/hg38/<genomic FASTA>
<NCBI MD5> <local SHA-256>
```

## 2) Partial genome data package {#2-dbd-002-partial-genome-data-package}

- **Purpose:** Download only the required chromosomes, sequence categories or file types from the specified assembly and generate a manifest to reduce the amount of transmission and storage of large genome data packages.
- **Method:** Select NCBI assembly package content by chromosome and file type.
- **API:** `dnakit.download.genome_package(accessions[required], output_dir[required], include[optional], chromosomes[optional], hydrated[optional], config[optional], progress[optional])`.
- **Input:** One or more `GCA_`/`GCF_` accessions with version numbers; optional chromosome name and `include` file type.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import genome_package

result = genome_package(
    "GCF_000001405.40",
    "downloads/human_chrM",
    chromosomes=("MT",),
    include=("GENOME_FASTA", "GENOME_GFF", "SEQUENCE_REPORT"),
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/human_chrM/ncbi_genome_package_manifest.json
```

## 3) region sequence {#3-dbd-003-region-sequence}

- **Purpose:** Download the target region sequence by assembly, chromosome and start and end coordinates, and return FASTA and coordinate source records for extracting loci or verifying local analysis results.
- **Method:** Submit region and strand parameters to Ensembl and save the returned sequence.
- **API:** `dnakit.download.sequence(species[required], region[required], output_path[required], strand[optional], upstream[optional], downstream[optional], mask[optional], config[optional], query_progress[optional])`.
- **Input:** Ensembl species name and 0-based half-open region; supports chain direction, upstream and downstream expansion, and multiple regions.
- **Output:** Returns `DatasetDownloadResult` with the region FASTA file and manifest path.
- **Example code:**

```python
from dnakit.download import sequence

result = sequence(
    "homo_sapiens",
    "7:140424943-140424950",
    "downloads/region.fa",
)
print(result.files[0].path)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/region.fa
downloads/region.fa.manifest.json
```

## 4) Genome annotation {#4-dbd-004-genome-annotation}

- **Purpose:** Download genome annotation files such as GFF, GTF, and GBFF that match the specified assembly, and record the version and verification values for gene, transcript, and functional region analysis.
- **Method:** Retrieve selected annotation formats from the NCBI assembly package and save a manifest.
- **API:** `dnakit.download.annotation(accession[required], output_dir[required], formats[optional], chromosomes[optional], config[optional], progress[optional])`.
- **Input:** Assembly accession with version number; optional annotation format and chromosome filtering.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import annotation

result = annotation(
    "GCF_000001405.40",
    "downloads/hg38_annotation",
    formats=("GENOME_GFF", "GENOME_GTF", "SEQUENCE_REPORT"),
)
print(result.files)
print(result.manifest_path)
```

- **Example output:**

```text
(DownloadedFile(...),)
downloads/hg38_annotation/ncbi_genome_package_manifest.json
```

## 5) Gene sequence {#5-dbd-005-gene-sequence}

- **Purpose:** Selectively download genome segments, RNA, CDS, UTR and protein FASTA according to GeneID, retain the source accession, and use it to build a local sequence collection of the specified gene.
- **Method:** Request an NCBI gene package using GeneIDs and include options.
- **API:** `dnakit.download.gene(gene_ids[required], output_dir[required], include[optional], accession_filter[optional], include_product_report[optional], config[optional], progress[optional])`.
- **Input:** One or more numeric NCBI GeneIDs; optional FASTA type and accession filter.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import gene

result = gene(
    672,
    "downloads/gene_672",
    include=("FASTA_GENE", "FASTA_RNA", "FASTA_CDS"),
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/gene_672/ncbi_gene_package_manifest.json
```

## 6) protein sequence {#6-dbd-006-protein-sequence}

- **Purpose:** Download public protein FASTA associated with gene or virus records, retaining protein accession and source information for translation result verification or protein reference library preparation.
- **Method:** Retrieve protein sequences by accession and save files with provenance.
- **API:** `dnakit.download.gene(..., include=("FASTA_PROTEIN",))`; viral proteins use `dnakit.download.virus_package(..., include=("PROTEIN",))`.
- **Input:** Gene downloads use a numeric GeneID; virus downloads use a virus accession or a single taxon.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import gene

result = gene(
    672,
    "downloads/gene_672_protein",
    include=("FASTA_PROTEIN",),
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/gene_672_protein/ncbi_gene_package_manifest.json
```

## 7) Assembly information {#7-dbd-007-assembly-information}

- **Purpose:** Download metadata files such as assembly reports, sequence reports, and annotation statistics, which are used to check chromosome names, assembly levels, file versions, and annotation completeness.
- **Method:** Retrieve structured assembly metadata and export it in the chosen format.
- **API:** `dnakit.download.genome_package(accessions[required], output_dir[required], include=("SEQUENCE_REPORT",), config[optional], progress[optional])`.
- **Input:** Assembly accession with version number; usually just select `SEQUENCE_REPORT` or other report file.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import genome_package

result = genome_package(
    "GCF_000001405.40",
    "downloads/hg38_assembly_report",
    include=("SEQUENCE_REPORT",),
)
print(result.metadata)
print(result.manifest_path)
```

- **Example output:**

```text
{'accessions': ('GCF_000001405.40',), 'include': ('SEQUENCE_REPORT',), ...}
downloads/hg38_assembly_report/ncbi_genome_package_manifest.json
```

## 8) Classification information {#8-dbd-008-classification-information}

- **Purpose:** Download the scientific name, alias and complete taxonomic pedigree corresponding to the specified TaxID, and save the structured summary, which is used to unify species identification and organize data by taxonomic level.
- **Method:** Retrieve taxonomy records and lineages and save structured results.
- **API:** `dnakit.download.taxonomy(tax_ids[required], output_dir[required], include_names[optional], include_summary[optional], config[optional], progress[optional])`.
- **Input:** One or more numeric NCBI TaxIDs; optionally named and summary reports.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import taxonomy

result = taxonomy(
    (9606, 10090),
    "downloads/taxonomy",
    include_names=True,
    include_summary=True,
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/taxonomy/ncbi_taxonomy_package_manifest.json
```

## 9) Gene metadata {#query-result-export}

- **Purpose:** Export genes, samples or run metadata obtained from public database queries to CSV, TSV or JSON to facilitate offline screening, auditing and association with downloaded files.
- **Method:** Query gene metadata and export the result table to a local file.
- **API:** `dnakit.download.metadata(query[required], output_path[required], format[optional], config[optional])`.
- **Input:** `dnakit.search` Returned `QueryResult`; format supports JSON, JSONL, CSV, TSV and XML.
- **Output:** Returns `DatasetDownloadResult` with the exported metadata file and manifest path.
- **Example code:**

```python
from dnakit.download import metadata
from dnakit.search import sample

samples = sample("Homo sapiens[Organism]", retmax=20)
result = metadata(samples, "downloads/samples.csv", format="csv")
print(result.manifest_path)
```

- **Example output:**

```text
downloads/samples.csv
downloads/samples.csv.manifest.json
```

## 10) Virus packet {#10-dbd-010-virus-packet}

- **Purpose:** Download the data package consisting of genome, CDS, protein and annotation by virus name, TaxID or accession, and save the source list for establishing a clear version of the virus reference set.
- **Method:** Retrieve the selected genome and optional annotations from an NCBI virus package.
- **API:** `dnakit.download.virus_package(query[required], output_dir[required], by[optional], include[optional], aux_reports[optional], refseq_only[optional], annotated_only[optional], complete_only[optional], host[optional], geo_location[optional], released_since[optional], updated_since[optional], config[optional], progress[optional])`.
- **Input:** Viral accession or individual taxon; optional RefSeq, annotated, complete genome, host, region and time filters.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import virus_package

result = virus_package(
    "NC_045512.2",
    "downloads/sars_cov_2",
    by="accession",
    include=("GENOME", "CDS", "PROTEIN"),
    aux_reports=("ANNOTATION", "BIOSAMPLE_REPORT"),
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/sars_cov_2/ncbi_virus_NC_045512.2_manifest.json
```

## 11) Raw sequencing data {#11-dbd-011-raw-sequencing-data}

- **Purpose:** Press ENA Run accession to download single-ended or double-ended FASTQ, verify the integrity based on the MD5 provided by ENA, and record the corresponding relationship between the sample and the file.
- **Method:** Resolve ENA Run FASTQ URLs, download files and verify provider MD5 values.
- **API:** `dnakit.download.reads(query[required], output_dir[required], file_kind="fastq", config[optional], progress[optional])`.
- **Input:** ENA Run accession, ENA query expression, or `QueryResult` returned by `dnakit.search.reads()`; `file_kind` is `fastq`.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import reads

result = reads(
    "SRR390728",
    "downloads/SRR390728_fastq",
    file_kind="fastq",
)
print(result.files)
print(result.manifest_path)
```

- **Example output:**

```text
(DownloadedFile(...),)
downloads/SRR390728_fastq/ena_fastq_manifest.json
```

## 12) Compare data {#12-dbd-012-compare-data}

- **Purpose:** Download ENA's public BAM, CRAM or SAM files and available verification information, which can be used to reuse existing reads-to-reference alignment results without re-alignment.
- **Method:** Read ENA submitted file lists and download available alignment files.
- **API:** `dnakit.download.reads(query[required], output_dir[required], file_kind="submitted", config[optional], progress[optional])`.
- **Input:** ENA Run accession, query expression, or ENA `QueryResult`; `file_kind` is `submitted`.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import reads

result = reads(
    "SRR390728",
    "downloads/SRR390728_submitted",
    file_kind="submitted",
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/SRR390728_submitted/ena_submitted_manifest.json
```

## 13) Sequencing metadata {#13-dbd-013-sequencing-metadata}

- **Purpose:** Query and save structured information such as Study, Sample, Experiment, Run and file links of ENA, which is used to filter sequencing data and track sample sources.
- **Method:** Query ENA metadata and serialize it to the selected tabular or structured format.
- **API:** `dnakit.search.reads(query[required], fields[optional], result[optional], limit[optional], offset[optional], config[optional])` plus `dnakit.download.metadata(query[required], output_path[required], format[optional], config[optional])`.
- **Input:** ENA query expression and optional fields; query results can be exported to JSON, JSONL, CSV, TSV, or XML.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.download import metadata
from dnakit.search import reads as search_reads

runs = search_reads('run_accession="SRR390728"', limit=20)
result = metadata(runs, "downloads/ena_runs.json", format="json")
print(result.manifest_path)
```

- **Example output:**

```text
downloads/ena_runs.json
downloads/ena_runs.json.manifest.json
```

## 14) Mutation data {#14-dbd-014-mutation-data}

- **Purpose:** Download the dbSNP VCF, index and verification files of the specified reference version for local query of included variants and annotation of sequence variant results.
- **Method:** Select dbSNP release files for the assembly and download VCF, indexes and checksums.
- **API:** `dnakit.download.variants(output_dir[required], source="dbsnp", format="vcf", assembly[optional], include_index[optional], include_checksum_files[optional], config[optional], progress[optional])`.
- **Input:** `source="dbsnp"`, `format="vcf"` and `GRCh37`/`GRCh38` assembly versions.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import variants

result = variants(
    "downloads/dbsnp_grch38",
    source="dbsnp",
    format="vcf",
    assembly="GRCh38",
)
print(result.files)
print(result.manifest_path)
```

- **Example output:**

```text
(dbsnp_GRCh38.vcf.gz, dbsnp_GRCh38.vcf.gz.tbi, ...)
downloads/dbsnp_grch38/dbsnp_vcf_manifest.json
```

## 15) ClinVar data {#15-dbd-015-clinvar-data}

- **Purpose:** Download ClinVar's VCF, XML or TSV release files and version information for locally associated variants and public clinical significance annotations.
- **Method:** Select the ClinVar release format and assembly and download the corresponding files.
- **API:** `dnakit.download.variants(output_dir[required], source="clinvar", format[optional], assembly[optional], include_index[optional], config[optional], progress[optional])`.
- **Input:** `source="clinvar"`; format `vcf`, `xml` or `tsv`, VCF optional `GRCh37`/`GRCh38` and index.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import variants

result = variants(
    "downloads/clinvar",
    source="clinvar",
    format="vcf",
    assembly="GRCh38",
    include_index=True,
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/clinvar/clinvar_vcf_manifest.json
```

## 16) Expression data {#16-dbd-016-expression-data}

- **Purpose:** Press GEO accession to download public expression matrices or supplementary files, and save sample and source information for offline expression data analysis.
- **Method:** Locate and download expression matrices or supplementary files by GEO accession.
- **API:** `dnakit.download.expression(accession[required], output_dir[required], format[optional], platform[optional], config[optional], progress[optional])`.
- **Input:** `GSE`, `GDS`, `GPL` or `GSM` accession; Series Matrix and common RAW typically require `GSE`, multi-platform matrices may specify `GPL`.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
from dnakit.download import expression

result = expression(
    "GSE100",
    "downloads/GSE100",
    format="matrix",
    platform="GPL96",
)
print(result.manifest_path)
```

- **Example output:**

```text
downloads/GSE100/geo_matrix_manifest.json
```

## 17) Control data {#17-dbd-017-control-data}

- **Purpose:** Download ENCODE/UCSC regulatory data by species, assembly, experiment or track conditions, retaining file metadata for annotating open chromatin, binding sites and other regions.
- **Method:** Filter ENCODE files or explicit track resources, then download the selected results.
- **API:** `dnakit.search.encode_search(...)` plus `dnakit.download.encode_files(query[required], output_dir[required], config[optional], progress[optional])`; explicit resources can also use `dnakit.download.tracks(resources[required], output_dir[required], source[optional], config[optional], progress[optional])`.
- **Input:** Query results for ENCODE `object_type="File"`, or a list of HTTPS `RemoteFile`s that the caller has confirmed.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.download import encode_files
from dnakit.search import encode_search

found = encode_search(
    object_type="File",
    filters={"file_format": "bigWig", "status": "released"},
    limit=10,
)
result = encode_files(found, "downloads/encode_regulation")
print(result.manifest_path)
```

- **Example output:**

```text
downloads/encode_regulation/encode_files_manifest.json
```

## 18) Repeating sequence data {#18-dbd-018-repeating-sequence-data}

- **Purpose:** Download repeat sequence tracks such as UCSC RepeatMasker that match the specified assembly and are used to mark repeat categories and coordinate intervals in the genome.
- **Method:** Filter UCSC RepeatMasker resources by filename pattern and download matches.
- **API:** `dnakit.search.ucsc_files(genome[required], pattern[optional], limit[optional], config[optional])` plus `dnakit.download.ucsc_files(query[required], output_dir[required], config[optional], progress[optional])`.
- **Input:** UCSC assembly and file glob; the download interface receives the corresponding `QueryResult`.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*repeat*", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_repeats")
print(result.manifest_path)
```

- **Example output:**

```text
downloads/hg38_repeats/ucsc_files_manifest.json
```

## 19) Conservative data {#19-dbd-019-conservative-data}

- **Purpose:** Download the conserved tracks and description information of the specified assembly such as phyloP and phastCons, which are used to associate cross-species conservation scores for genome positions.
- **Method:** Select UCSC phyloP or phastCons files and download existing conservation tracks.
- **API:** `dnakit.search.ucsc_files(...)` plus `dnakit.download.ucsc_files(...)`; individual tracks with specified URLs are also available `dnakit.download.tracks(...)`.
- **Input:** UCSC assembly and conservation file globs, or explicit HTTPS files confirmed by the caller.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*phastCons*", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_conservation")
print(result.manifest_path)
```

- **Example output:**

```text
downloads/hg38_conservation/ucsc_files_manifest.json
```

## 20) Multi-species comparison {#20-dbd-020-multi-species-comparison}

- **Purpose:** Download the UCSC multi-species MAF alignment file of the specified assembly and species collection for offline viewing of homologous regions and comparative genome analysis.
- **Method:** Query UCSC multiple-species alignment resources and download selected files.
- **API:** `dnakit.search.ucsc_files(...)` plus `dnakit.download.ucsc_files(...)`; query result metadata can be exported with `dnakit.download.metadata(...)`.
- **Input:** UCSC assembly and MAF/multi-sequence file globs, or explicit HTTPS files.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*.maf.gz", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_maf")
print(result.manifest_path)
```

- **Example output:**

```text
downloads/hg38_maf/ucsc_files_manifest.json
```

## 21) coordinate conversion file {#21-dbd-021-coordinate-conversion-file}

- **Purpose:** Download a UCSC chain file between two genome assembly versions and record the directions, used to map region, variant and annotation coordinates from one assembly to the other.
- **Method:** Select a UCSC chain file for the source and target assemblies and download it.
- **API:** `dnakit.search.ucsc_files(...)` plus `dnakit.download.ucsc_files(...)`.
- **Input:** UCSC assembly and `liftOver`/`chain` file globs.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*liftOver*chain.gz", limit=1000)
result = download_ucsc_files(catalog, "downloads/ucsc_chain")
print(result.manifest_path)
```

- **Example output:**

```text
downloads/ucsc_chain/ucsc_files_manifest.json
```

## 22) Download from a URL {#22-dbd-022-download-from-a-url}

- **Purpose:** Download a file from a direct URL.
- **Method:** Save response bytes unchanged and optionally verify MD5 or SHA-256.
- **API:** `dnakit.download.download_url(url, output_dir, *, filename=None, config=None, ...)`
- **Input:** A direct HTTPS file URL and output directory; replace the example URL. Existing files are protected, and archives are not extracted.
- **Output:** Returns `DownloadedFile` with the local path, size, and checksums.
- **Example code:**

```python
from dnakit.download import download_url

result = download_url(
    "https://example.org/files/data.fa.gz",  # Replace with your direct file URL
    "downloads/my_data",
)
print(result.path)       # Absolute path to the downloaded file
print(result.byte_size)  # File size in bytes
```

- **Example output:** Absolute file path and downloaded byte count; values depend on the URL.

## 23) Search and download a sequence {#sequence-download-example}

- **Purpose:** Find sequences by organism and download the first match.
- **Method:** Search NCBI Entrez, resolve record IDs to accessions, then download the selected sequence.
- **API:** `dnakit.search.ncbi.entrez(...)`; `dnakit.search.ncbi.accessions(...)`; `dnakit.download.ncbi.sequence(...)`
- **Input:** Organism and search filters, plus an output path; requires network access and writes local files.
- **Output:** Search returns `QueryResult`; download returns `DatasetDownloadResult` with the sequence file, length, and manifest.
- **Example code:**

```python
import dnakit

found = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(found.metadata["uids"])
if accessions:
    downloaded = dnakit.download.ncbi.sequence(accessions[0], "ecoli_k12.fasta", format="fasta")
    print(downloaded.files[0].path)
    print(downloaded.metadata["sequence_length"])
    print(downloaded.manifest_path)
else:
    print("No sequence matches; review the query.")
```

- **Example output:** The downloaded path, sequence length and manifest; prints a no-match message if the query returns no records.

<span id="common-genome-names"></span>Species names, ENA examples and options: [FAQ](../../faq/tutorial_integrations.md#downloads).

<span id="encode-download-example"></span>[ENCODE batch downloads and dry-run preview](#dnakit.encode.batch_download).

<!-- BEGIN GENERATED API: downloads -->

<span id="downloads"></span>

## 24) Download a file over HTTP or HTTPS {#dnakit.databases.download_file}

- **Purpose:** Download a file over HTTP or HTTPS.
- **Method:** Stream an HTTP(S) response to the destination, or return text when return_content=True.
- **API：** `dnakit.databases.download_file(*, url, output_path=None, timeout=None, return_content=None, chunk_size=None, follow_redirects=None, verbose=True)`
- **Input:** Provide `url`. Prepare the local files and external programs used by the example, and replace paths with your input or output locations.
- **Output:** Returns the file path when saving; with `return_content=True`, returns the file contents as text.
- **Example code:**

```python
import dnakit

result = dnakit.databases.download_file(
    url='https://www.wikipedia.org',
    output_path='./downloads/output.txt',
    timeout=30,
    return_content=False,
    chunk_size=8192,
    follow_redirects=True,
    verbose=False,
)
print(result)
```

- **Example output:** Returns the file path when saving; with `return_content=True`, returns the file contents as text.

[Full parameters and usage notes](../../faq/api/15_download.md#dnakit.databases.download_file-parameters)

## 25) Download an ENA nucleotide sequence with a checksum manifest {#dnakit.download.ena.sequence}

- **Purpose:** Download an ENA nucleotide sequence with a checksum manifest.
- **Method:** Resolve the accession, retrieve sequence data, and write the files and provenance manifest.
- **API：** `dnakit.download.ena.sequence(accession, output_path, *, format='fasta', config=None, search_config=None, api_base_url='https://www.ebi.ac.uk/ena/browser/api', verbose=True)`
- **Input:** Provide `accession`, `output_path`.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
import dnakit

result = dnakit.download.ena.sequence(
    accession='U00096.3',
    output_path='ecoli_ena.fasta',
    format='fasta',
    verbose=False,
)
print(result)
```

- **Example output:** DatasetDownloadResult with files, metadata (source, resolved accession and length), and manifest_path. See the [sequence download guide](15_download.md#accession-sequences) for formats and coordinates.

[Full parameters and usage notes](../../faq/api/15_download.md#dnakit.download.ena.sequence-parameters)

## 26) Download an NCBI nucleotide or protein sequence {#dnakit.download.ncbi.sequence}

- **Purpose:** Download an NCBI nucleotide or protein sequence.
- **Method:** Resolve the accession, retrieve sequence data, and write the files and provenance manifest.
- **API：** `dnakit.download.ncbi.sequence(accession, output_path, *, database='nuccore', format='fasta', seq_start=None, seq_stop=None, strand=1, config=None, search_config=None, api_base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils', verbose=True)`
- **Input:** Provide `accession`, `output_path`.
- **Output:** Returns `DatasetDownloadResult` with downloaded `files`, source `metadata`, and `manifest_path`.
- **Example code:**

```python
import dnakit

result = dnakit.download.ncbi.sequence(
    accession='NC_000913.3',
    output_path='ecoli_ncbi.fasta',
    format='fasta',
    verbose=False,
)
print(result)
```

- **Example output:** DatasetDownloadResult with files, metadata (source, resolved accession and length), and manifest_path. See the [sequence download guide](15_download.md#accession-sequences) for formats and coordinates.

[Full parameters and usage notes](../../faq/api/15_download.md#dnakit.download.ncbi.sequence-parameters)

<!-- END GENERATED API: downloads -->

<!-- BEGIN GENERATED API: encode-downloads -->

<span id="encode-downloads"></span>

## 27) Find and batch-download ENCODE files {#dnakit.encode.batch_download}

- **Purpose:** Find and batch-download ENCODE files.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.batch_download(download_dir, file_format=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, target=None, preferred_default=None, organize_by='experiment', verify_md5=True, limit=100, dry_run=True, *, verbose=True)`
- **Input:** Provide `download_dir`.
- **Output:** Python return type: str.  JSON with download preview (dry_run=True) or download results (dry_run=False).
- **Example code:**

```python
import dnakit

result = dnakit.encode.batch_download(
    download_dir='downloads/encode',
    assay_title='ATAC-seq',
    assembly='GRCh38',
    file_format='bed',
    limit=3,
    dry_run=True,
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with download preview (dry_run=True) or download results (dry_run=False).

[Full parameters and usage notes](../../faq/api/15_download.md#dnakit.encode.batch_download-parameters)

## 28) Download ENCODE files by accession {#dnakit.encode.download_files}

- **Purpose:** Download ENCODE files by accession.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.download_files(file_accessions, download_dir, organize_by='flat', verify_md5=True, *, verbose=True)`
- **Input:** Provide `file_accessions`, `download_dir`. Retrieve a real file accession with the file query API and replace `ENCFF_REPLACE_ME`.
- **Output:** Python return type: str.  JSON with download results for each file (path, size, success/error, MD5 status).
- **Example code:**

```python
import dnakit

result = dnakit.encode.download_files(
    file_accessions=['ENCFF_REPLACE_ME'],
    download_dir='downloads/encode',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with download results for each file (path, size, success/error, MD5 status).

[Full parameters and usage notes](../../faq/api/15_download.md#dnakit.encode.download_files-parameters)

<!-- END GENERATED API: encode-downloads -->
