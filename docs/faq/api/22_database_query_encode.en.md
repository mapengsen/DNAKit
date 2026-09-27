# ENCODE experiment and file queries: parameters and notes

See the [tutorial](../../api/features/22_database_query_encode.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

ENCODE returns text: mostly JSON, but get_citations can return BibTeX/RIS and export_data can return CSV/TSV. Use json.loads() only for JSON. Tracking and credentials use ~/.dnakit-encode/; persistent credential storage needs dnakit[encode-credentials].

## Retrieve ENCODE experiment details by accession {#dnakit.encode.get_experiment-parameters}

`dnakit.encode.get_experiment(accession, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;, &quot;ENCSR000AKS&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Retrieve live ENCODE filter counts {#dnakit.encode.get_facets-parameters}

`dnakit.encode.get_facets(search_type='Experiment', assay_title=None, organism=None, organ=None, biosample_type=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `search_type` | str | &#x27;Experiment&#x27; | Object type (&quot;Experiment&quot; or &quot;File&quot;) |
| `assay_title` | str &#124; None | None | Pre-filter by assay type |
| `organism` | str &#124; None | None | Pre-filter by organism |
| `organ` | str &#124; None | None | Pre-filter by organ |
| `biosample_type` | str &#124; None | None | Pre-filter by biosample type |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Retrieve ENCODE file details {#dnakit.encode.get_file_info-parameters}

`dnakit.encode.get_file_info(accession, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | File accession ID (e.g., &quot;ENCFF635JIA&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Retrieve available ENCODE search filters {#dnakit.encode.get_metadata-parameters}

`dnakit.encode.get_metadata(metadata_type, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `metadata_type` | Literal&#91;&#x27;assays&#x27;, &#x27;organisms&#x27;, &#x27;organs&#x27;, &#x27;biosample_types&#x27;, &#x27;file_formats&#x27;, &#x27;output_types&#x27;, &#x27;output_categories&#x27;, &#x27;assemblies&#x27;, &#x27;life_stages&#x27;, &#x27;replication_types&#x27;, &#x27;statuses&#x27;, &#x27;file_statuses&#x27;&#93; | Required | Type of metadata to retrieve. Options: - &quot;assays&quot;: Available assay types (Histone ChIP-seq, ATAC-seq, total RNA-seq, etc.) - &quot;organisms&quot;: Available organisms (Homo sapiens, Mus musculus, etc.) - &quot;organs&quot;: Available organ/tissue systems (pancreas, brain, liver, etc.) - &quot;biosample_types&quot;: Biosample classifications (tissue, cell line, primary cell, etc.) - &quot;file_formats&quot;: File format types (fastq, bam, bed, bigWig, etc.) - &quot;output_types&quot;: Output data types (reads, peaks, signal, etc.) - &quot;output_categories&quot;: Output categories (raw data, alignment, signal, etc.) - &quot;assemblies&quot;: Genome assemblies (GRCh38, hg19, mm10, etc.) - &quot;life_stages&quot;: Life stages (embryonic, adult, child, etc.) - &quot;replication_types&quot;: Replication types (isogenic, anisogenic, unreplicated) - &quot;statuses&quot;: Experiment statuses (released, archived, etc.) - &quot;file_statuses&quot;: File statuses (released, archived, in progress, etc.) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## List files belonging to an ENCODE experiment {#dnakit.encode.list_files-parameters}

`dnakit.encode.list_files(experiment_accession, file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, status=None, preferred_default=None, limit=200, *, verbose=True)`

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


## Search ENCODE experiments by experimental conditions {#dnakit.encode.search_experiments-parameters}

`dnakit.encode.search_experiments(assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, biosample_term_name=None, target=None, status='released', lab=None, award=None, assembly=None, replication_type=None, life_stage=None, sex=None, treatment=None, genetic_modification=None, perturbed=None, search_term=None, date_released_from=None, date_released_to=None, limit=25, offset=0, *, verbose=True)`

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


## Search ENCODE files using experiment and file filters {#dnakit.encode.search_files-parameters}

`dnakit.encode.search_files(file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism=None, organ=None, biosample_type=None, target=None, status='released', preferred_default=None, search_term=None, limit=25, offset=0, *, verbose=True)`

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
