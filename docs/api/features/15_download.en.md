# Download

Download reference genomes, sequences, annotations, and other public database files based on species name, database number, or query results.

<span id="common-genome-names"></span>**Common species names**

This section lists items that can be copied directly to `download_genome(query, ...)` or
Common NCBI species names in `dnakit download-genome` commands.

- **API**: `dnakit.references.download_genome(query[required], output_dir[required], overwrite[optional], keep_package[optional], timeout[optional], api_key[optional], api_base_url[optional], chunk_size[optional], max_download_bytes[optional], progress[optional])`, `dnakit download-genome query[required] output_dir[required] --overwrite[optional] --keep-package[optional] --api-key[optional] --progress/--no-progress[optional]`.

<span id="1"></span>**People and commonly used model animals**

| Chinese name | `query` that can be copied directly | Current reference assembly accession |
| -------- | -------------------------- | -------------------------- |
| Human | `Homo sapiens` | `GCF_000001405.40` |
| Mouse | `Mus musculus` | `GCF_000001635.27` |
| Rat | `Rattus norvegicus` | `GCF_036323735.1` |
| Pig | `Sus scrofa` | `GCF_054392235.1` |
| Cow | `Bos taurus` | `GCF_002263795.3` |
| Sheep | `Ovis aries` | `GCF_016772045.2` |
| Goat | `Capra hircus` | `GCF_001704415.2` |
| Horse | `Equus caballus` | `GCF_041296265.1` |
| Dog | `Canis lupus familiaris` | `GCF_011100685.1` |
| Macaque | `Macaca mulatta` | `GCF_049350105.2` |
| Chicken | `Gallus gallus` | `GCF_016699485.2` |
| Zebrafish | `Danio rerio` | `GCF_049306965.2` |
| Xenopus | `Xenopus tropicalis` | `GCF_000004195.4` |

<span id="2"></span>**Common crops and plants**

| Chinese name | `query` that can be copied directly | Current reference assembly accession |
| -------------------------- | ------------------------------- | ----------------------- |
| Arabidopsis | `Arabidopsis thaliana` | `GCF_000001735.4` |
| Rice | `Oryza sativa Japonica Group` | `GCF_034140825.1` |
| Corn | `Zea mays` | `GCF_902167145.1` |
| Wheat | `Triticum aestivum` | `GCF_018294505.1` |
| Soybeans | `Glycine max` | `GCF_000004515.6` |
| Tomato | `Solanum lycopersicum` | `GCF_036512215.1` |
| Potatoes | `Solanum tuberosum` | `GCF_000226075.1` |
| Upland cotton | `Gossypium hirsutum` | `GCF_007990345.1` |
| Tobacco | `Nicotiana tabacum` | `GCF_000715075.1` |
| Populus trichocarpa | `Populus trichocarpa` | `GCF_000002775.5` |
| Arabica coffee | `Coffea arabica` | `GCF_036785885.1` |

<span id="3"></span>**Other commonly used model creatures**

| Chinese name | `query` that can be copied directly | Current reference assembly accession |
| --------------------- | ---------------------------------- | ----------------------- |
| Drosophila | `Drosophila melanogaster` | `GCF_000001215.4` |
| Caenorhabditis elegans | `Caenorhabditis elegans` | `GCF_000002985.6` |
| Bombyx mori | `Bombyx mori` | `GCF_030269925.1` |
| Bee | `Apis mellifera` | `GCF_003254395.2` |
| Saccharomyces cerevisiae | `Saccharomyces cerevisiae S288C` | `GCF_000146045.2` |

## 1) `DBD-001` complete genome

- **Function:** Download the complete reference genome FASTA by assembly accession or common species name, and save the assembly report, source URL and verification information for establishing traceable local reference sequences.
- **API:** `dnakit.references.resolve_genome_assembly(query[required], timeout[optional], api_key[optional], api_base_url[optional])`, `dnakit.references.download_genome(query[required], output_dir[required], overwrite[optional], keep_package[optional], timeout[optional], api_key[optional], api_base_url[optional], chunk_size[optional], max_download_bytes[optional], progress[optional])`.
- **Input:** NCBI versioned accession, fixed alias (such as `hg38`), or species name; requires local output directory.
- **Sample code:**

```python
from pathlib import Path

from dnakit.references import download_genome

result = download_genome("hg38", Path("references/hg38"))
print(result.fasta_path)
print(result.fasta_md5, result.fasta_sha256)
```

- **Example results:**

```text
references/hg38/<genomic FASTA>
<NCBI MD5> <local SHA-256>
```

## 2) `DBD-002` Partial genome data package

- **Function:** Download only the required chromosomes, sequence categories or file types from the specified assembly and generate a manifest to reduce the amount of transmission and storage of large genome data packages.
- **API:** `dnakit.download.genome_package(accessions[required], output_dir[required], include[optional], chromosomes[optional], hydrated[optional], config[optional], progress[optional])`.
- **Input:** One or more `GCA_`/`GCF_` accessions with version numbers; optional chromosome name and `include` file type.
- **Sample code:**

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

- **Example results:**

```text
downloads/human_chrM/ncbi_genome_package_manifest.json
```

## 3) `DBD-003` region sequence

- **Function:** Download the target region sequence by assembly, chromosome and start and end coordinates, and return FASTA and coordinate source records for extracting loci or verifying local analysis results.
- **API:** `dnakit.download.sequence(species[required], region[required], output_path[required], strand[optional], upstream[optional], downstream[optional], mask[optional], config[optional], query_progress[optional])`.
- **Input:** Ensembl species name and 0-based half-open region; supports chain direction, upstream and downstream expansion, and multiple regions.
- **Sample code:**

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

- **Example results:**

```text
downloads/region.fa
downloads/region.fa.manifest.json
```

## 4) `DBD-004` Genome annotation

- **Function:** Download genome annotation files such as GFF, GTF, and GBFF that match the specified assembly, and record the version and verification values for gene, transcript, and functional region analysis.
- **API:** `dnakit.download.annotation(accession[required], output_dir[required], formats[optional], chromosomes[optional], config[optional], progress[optional])`.
- **Input:** Assembly accession with version number; optional annotation format and chromosome filtering.
- **Sample code:**

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

- **Example results:**

```text
(DownloadedFile(...),)
downloads/hg38_annotation/ncbi_genome_package_manifest.json
```

## 5) `DBD-005` Gene sequence

- **Function:** Selectively download genome segments, RNA, CDS, UTR and protein FASTA according to GeneID, retain the source accession, and use it to build a local sequence collection of the specified gene.
- **API:** `dnakit.download.gene(gene_ids[required], output_dir[required], include[optional], accession_filter[optional], include_product_report[optional], config[optional], progress[optional])`.
- **Input:** One or more numeric NCBI GeneIDs; optional FASTA type and accession filter.
- **Sample code:**

```python
from dnakit.download import gene

result = gene(
    672,
    "downloads/gene_672",
    include=("FASTA_GENE", "FASTA_RNA", "FASTA_CDS"),
)
print(result.manifest_path)
```

- **Example results:**

```text
downloads/gene_672/ncbi_gene_package_manifest.json
```

## 6) `DBD-006` protein sequence

- **Function:** Download public protein FASTA associated with gene or virus records, retaining protein accession and source information for translation result verification or protein reference library preparation.
- **API:** `dnakit.download.gene(..., include=("FASTA_PROTEIN",))`; viral proteins use `dnakit.download.virus_package(..., include=("PROTEIN",))`.
- **Input:** Gene downloads use a numeric GeneID; virus downloads use a virus accession or a single taxon.
- **Sample code:**

```python
from dnakit.download import gene

result = gene(
    672,
    "downloads/gene_672_protein",
    include=("FASTA_PROTEIN",),
)
print(result.manifest_path)
```

- **Example results:**

```text
downloads/gene_672_protein/ncbi_gene_package_manifest.json
```

## 7) `DBD-007` Assembly information

- **Function:** Download metadata files such as assembly reports, sequence reports, and annotation statistics, which are used to check chromosome names, assembly levels, file versions, and annotation completeness.
- **API:** `dnakit.download.genome_package(accessions[required], output_dir[required], include=("SEQUENCE_REPORT",), config[optional], progress[optional])`.
- **Input:** Assembly accession with version number; usually just select `SEQUENCE_REPORT` or other report file.
- **Sample code:**

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

- **Example results:**

```text
{'accessions': ('GCF_000001405.40',), 'include': ('SEQUENCE_REPORT',), ...}
downloads/hg38_assembly_report/ncbi_genome_package_manifest.json
```

## 8) `DBD-008` Classification information

- **Function:** Download the scientific name, alias and complete taxonomic pedigree corresponding to the specified TaxID, and save the structured summary, which is used to unify species identification and organize data by taxonomic level.
- **API:** `dnakit.download.taxonomy(tax_ids[required], output_dir[required], include_names[optional], include_summary[optional], config[optional], progress[optional])`.
- **Input:** One or more numeric NCBI TaxIDs; optionally named and summary reports.
- **Sample code:**

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

- **Example results:**

```text
downloads/taxonomy/ncbi_taxonomy_package_manifest.json
```

## 9) `DBD-009` Gene metadata {#query-result-export}

- **Function:** Export genes, samples or run metadata obtained from public database queries to CSV, TSV or JSON to facilitate offline screening, auditing and association with downloaded files.
- **API:** `dnakit.download.metadata(query[required], output_path[required], format[optional], config[optional])`.
- **Input:** `dnakit.search` Returned `QueryResult`; format supports JSON, JSONL, CSV, TSV and XML.
- **Sample code:**

```python
from dnakit.download import metadata
from dnakit.search import sample

samples = sample("Homo sapiens[Organism]", retmax=20)
result = metadata(samples, "downloads/samples.csv", format="csv")
print(result.manifest_path)
```

- **Example results:**

```text
downloads/samples.csv
downloads/samples.csv.manifest.json
```

## 10) `DBD-010` Virus packet

- **Function:** Download the data package consisting of genome, CDS, protein and annotation by virus name, TaxID or accession, and save the source list for establishing a clear version of the virus reference set.
- **API:** `dnakit.download.virus_package(query[required], output_dir[required], by[optional], include[optional], aux_reports[optional], refseq_only[optional], annotated_only[optional], complete_only[optional], host[optional], geo_location[optional], released_since[optional], updated_since[optional], config[optional], progress[optional])`.
- **Input:** Viral accession or individual taxon; optional RefSeq, annotated, complete genome, host, region and time filters.
- **Sample code:**

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

- **Example results:**

```text
downloads/sars_cov_2/ncbi_virus_NC_045512.2_manifest.json
```

## 11) `DBD-011` Raw sequencing data

- **Function:** Press ENA Run accession to download single-ended or double-ended FASTQ, verify the integrity based on the MD5 provided by ENA, and record the corresponding relationship between the sample and the file.
- **API:** `dnakit.download.reads(query[required], output_dir[required], file_kind="fastq", config[optional], progress[optional])`.
- **Input:** ENA Run accession, ENA query expression, or `QueryResult` returned by `dnakit.search.reads()`; `file_kind` is `fastq`.
- **Sample code:**

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

- **Example results:**

```text
(DownloadedFile(...),)
downloads/SRR390728_fastq/ena_fastq_manifest.json
```

## 12) `DBD-012` Compare data

- **Function:** Download ENA's public BAM, CRAM or SAM files and available verification information, which can be used to reuse existing reads-to-reference alignment results without re-alignment.
- **API:** `dnakit.download.reads(query[required], output_dir[required], file_kind="submitted", config[optional], progress[optional])`.
- **Input:** ENA Run accession, query expression, or ENA `QueryResult`; `file_kind` is `submitted`.
- **Sample code:**

```python
from dnakit.download import reads

result = reads(
    "SRR390728",
    "downloads/SRR390728_submitted",
    file_kind="submitted",
)
print(result.manifest_path)
```

- **Example results:**

```text
downloads/SRR390728_submitted/ena_submitted_manifest.json
```

## 13) `DBD-013` Sequencing metadata

- **Function:** Query and save structured information such as Study, Sample, Experiment, Run and file links of ENA, which is used to filter sequencing data and track sample sources.
- **API:** `dnakit.search.reads(query[required], fields[optional], result[optional], limit[optional], offset[optional], config[optional])` plus `dnakit.download.metadata(query[required], output_path[required], format[optional], config[optional])`.
- **Input:** ENA query expression and optional fields; query results can be exported to JSON, JSONL, CSV, TSV, or XML.
- **Sample code:**

```python
from dnakit.download import metadata
from dnakit.search import reads as search_reads

runs = search_reads('run_accession="SRR390728"', limit=20)
result = metadata(runs, "downloads/ena_runs.json", format="json")
print(result.manifest_path)
```

- **Example results:**

```text
downloads/ena_runs.json
downloads/ena_runs.json.manifest.json
```

## 14) `DBD-014` Mutation data

- **Function:** Download the dbSNP VCF, index and verification files of the specified reference version for local query of included variants and annotation of sequence variant results.
- **API:** `dnakit.download.variants(output_dir[required], source="dbsnp", format="vcf", assembly[optional], include_index[optional], include_checksum_files[optional], config[optional], progress[optional])`.
- **Input:** `source="dbsnp"`, `format="vcf"` and `GRCh37`/`GRCh38` assembly versions.
- **Sample code:**

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

- **Example results:**

```text
(dbsnp_GRCh38.vcf.gz, dbsnp_GRCh38.vcf.gz.tbi, ...)
downloads/dbsnp_grch38/dbsnp_vcf_manifest.json
```

## 15) `DBD-015` ClinVar data

- **What it does:** Download ClinVar's VCF, XML or TSV release files and version information for locally associated variants and public clinical significance annotations.
- **API:** `dnakit.download.variants(output_dir[required], source="clinvar", format[optional], assembly[optional], include_index[optional], config[optional], progress[optional])`.
- **Input:** `source="clinvar"`; format `vcf`, `xml` or `tsv`, VCF optional `GRCh37`/`GRCh38` and index.
- **Sample code:**

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

- **Example results:**

```text
downloads/clinvar/clinvar_vcf_manifest.json
```

## 16) `DBD-016` Expression data

- **Function:** Press GEO accession to download public expression matrices or supplementary files, and save sample and source information for offline expression data analysis.
- **API:** `dnakit.download.expression(accession[required], output_dir[required], format[optional], platform[optional], config[optional], progress[optional])`.
- **Input:** `GSE`, `GDS`, `GPL` or `GSM` accession; Series Matrix and common RAW typically require `GSE`, multi-platform matrices may specify `GPL`.
- **Sample code:**

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

- **Example results:**

```text
downloads/GSE100/geo_matrix_manifest.json
```

## 17) `DBD-017` Control data

- **Function:** Download ENCODE/UCSC regulatory data by species, assembly, experiment or track conditions, retaining file metadata for annotating open chromatin, binding sites and other regions.
- **API:** `dnakit.search.encode_search(...)` plus `dnakit.download.encode_files(query[required], output_dir[required], config[optional], progress[optional])`; explicit resources can also use `dnakit.download.tracks(resources[required], output_dir[required], source[optional], config[optional], progress[optional])`.
- **Input:** Query results for ENCODE `object_type="File"`, or a list of HTTPS `RemoteFile`s that the caller has confirmed.
- **Sample code:**

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

- **Example results:**

```text
downloads/encode_regulation/encode_files_manifest.json
```

## 18) `DBD-018` Repeating sequence data

- **Function:** Download repeat sequence tracks such as UCSC RepeatMasker that match the specified assembly and are used to mark repeat categories and coordinate intervals in the genome.
- **API:** `dnakit.search.ucsc_files(genome[required], pattern[optional], limit[optional], config[optional])` plus `dnakit.download.ucsc_files(query[required], output_dir[required], config[optional], progress[optional])`.
- **Input:** UCSC assembly and file glob; the download interface receives the corresponding `QueryResult`.
- **Sample code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*repeat*", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_repeats")
print(result.manifest_path)
```

- **Example results:**

```text
downloads/hg38_repeats/ucsc_files_manifest.json
```

## 19) `DBD-019` Conservative data

- **Function:** Download the conserved tracks and description information of the specified assembly such as phyloP and phastCons, which are used to associate cross-species conservation scores for genome positions.
- **API:** `dnakit.search.ucsc_files(...)` plus `dnakit.download.ucsc_files(...)`; individual tracks with specified URLs are also available `dnakit.download.tracks(...)`.
- **Input:** UCSC assembly and conservation file globs, or explicit HTTPS files confirmed by the caller.
- **Sample code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*phastCons*", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_conservation")
print(result.manifest_path)
```

- **Example results:**

```text
downloads/hg38_conservation/ucsc_files_manifest.json
```

## 20) `DBD-020` Multi-species comparison

- **Function:** Download the UCSC multi-species MAF alignment file of the specified assembly and species collection for offline viewing of homologous regions and comparative genome analysis.
- **API:** `dnakit.search.ucsc_files(...)` plus `dnakit.download.ucsc_files(...)`; query result metadata can be exported with `dnakit.download.metadata(...)`.
- **Input:** UCSC assembly and MAF/multi-sequence file globs, or explicit HTTPS files.
- **Sample code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*.maf.gz", limit=100)
result = download_ucsc_files(catalog, "downloads/hg38_maf")
print(result.manifest_path)
```

- **Example results:**

```text
downloads/hg38_maf/ucsc_files_manifest.json
```

## 21) `DBD-021` coordinate conversion file

- **What it does:** Download a UCSC chain file between two genome assembly versions and record the directions, used to map region, variant and annotation coordinates from one assembly to the other.
- **API:** `dnakit.search.ucsc_files(...)` plus `dnakit.download.ucsc_files(...)`.
- **Input:** UCSC assembly and `liftOver`/`chain` file globs.
- **Sample code:**

```python
from dnakit.download import ucsc_files as download_ucsc_files
from dnakit.search import ucsc_files as list_ucsc_files

catalog = list_ucsc_files("hg38", pattern="*liftOver*chain.gz", limit=1000)
result = download_ucsc_files(catalog, "downloads/ucsc_chain")
print(result.manifest_path)
```

- **Example results:**

```text
downloads/ucsc_chain/ucsc_files_manifest.json
```
