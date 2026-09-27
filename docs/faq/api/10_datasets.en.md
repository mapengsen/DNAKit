# Data partitioning: parameters and notes

See the [tutorial](../../api/features/10_datasets.md) for usage and examples.

Each function includes its purpose, API, inputs, outputs, example code and result description. Prepare any required credentials, local files and optional dependencies; remote results vary with database releases and query time.

ENCODE returns text: mostly JSON, but get_citations can return BibTeX/RIS and export_data can return CSV/TSV. Use json.loads() only for JSON. Tracking and credentials use ~/.dnakit-encode/; persistent credential storage needs dnakit[encode-credentials].

## Compare ENCODE experiments for joint analysis compatibility {#dnakit.encode.compare_experiments-parameters}

`dnakit.encode.compare_experiments(accession1, accession2, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession1` | str | Required | First experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `accession2` | str | Required | Second experiment accession (e.g., &quot;ENCSR000AKS&quot;) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Export tracked ENCODE experiments as a table {#dnakit.encode.export_data-parameters}

`dnakit.encode.export_data(format='csv', assay_title=None, organism=None, organ=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `format` | Literal&#91;&#x27;csv&#x27;, &#x27;tsv&#x27;, &#x27;json&#x27;&#93; | &#x27;csv&#x27; | Output format: - &quot;csv&quot;: Comma-separated values (default, for Excel/spreadsheets) - &quot;tsv&quot;: Tab-separated values (for R, pandas) - &quot;json&quot;: JSON array (for programmatic use) |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Retrieve publications and citations for tracked experiments {#dnakit.encode.get_citations-parameters}

`dnakit.encode.get_citations(accession=None, export_format='json', *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str &#124; None | None | Specific experiment accession. If None, returns all publications. |
| `export_format` | Literal&#91;&#x27;json&#x27;, &#x27;bibtex&#x27;, &#x27;ris&#x27;&#93; | &#x27;json&#x27; | Output format: - &quot;json&quot;: Structured data (default) - &quot;bibtex&quot;: BibTeX format for LaTeX - &quot;ris&quot;: RIS format (Endnote, Zotero, Mendeley) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Retrieve derived-file provenance {#dnakit.encode.get_provenance-parameters}

`dnakit.encode.get_provenance(file_path=None, source_accession=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file_path` | str &#124; None | None | Get provenance for a specific derived file |
| `source_accession` | str &#124; None | None | List all files derived from a specific ENCODE accession |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Retrieve external references linked to experiments {#dnakit.encode.get_references-parameters}

`dnakit.encode.get_references(experiment_accession=None, reference_type=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `experiment_accession` | str &#124; None | None | Filter by specific experiment (optional) |
| `reference_type` | Literal&#91;&#x27;pmid&#x27;, &#x27;doi&#x27;, &#x27;nct_id&#x27;, &#x27;preprint_doi&#x27;, &#x27;geo_accession&#x27;, &#x27;other&#x27;&#93; &#124; None | None | Filter by reference type (optional): &quot;pmid&quot;, &quot;doi&quot;, &quot;nct_id&quot;, &quot;preprint_doi&quot;, &quot;geo_accession&quot;, &quot;other&quot; |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Link an external reference to a tracked experiment {#dnakit.encode.link_reference-parameters}

`dnakit.encode.link_reference(experiment_accession, reference_type, reference_id, description='', *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `experiment_accession` | str | Required | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `reference_type` | Literal&#91;&#x27;pmid&#x27;, &#x27;doi&#x27;, &#x27;nct_id&#x27;, &#x27;preprint_doi&#x27;, &#x27;geo_accession&#x27;, &#x27;other&#x27;&#93; | Required | Type of external reference: - &quot;pmid&quot;: PubMed ID (e.g., &quot;32728249&quot;) - &quot;doi&quot;: DOI (e.g., &quot;10.1038/s41586-020-2493-4&quot;) - &quot;nct_id&quot;: ClinicalTrials.gov ID (e.g., &quot;NCT04567890&quot;) - &quot;preprint_doi&quot;: bioRxiv/medRxiv DOI - &quot;geo_accession&quot;: GEO accession (e.g., &quot;GSE123456&quot;) - &quot;other&quot;: Any other identifier |
| `reference_id` | str | Required | The actual identifier value |
| `description` | str | &#x27;&#x27; | Optional description of why this reference is linked |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## List locally tracked ENCODE experiments {#dnakit.encode.list_tracked-parameters}

`dnakit.encode.list_tracked(assay_title=None, organism=None, organ=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Record a derived file and its provenance {#dnakit.encode.log_derived_file-parameters}

`dnakit.encode.log_derived_file(file_path, source_accessions, description='', file_type='', tool_used='', parameters='', *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `file_path` | str | Required | Path to the derived file you created |
| `source_accessions` | list&#91;str&#93; | Required | List of ENCODE accessions this file was derived from (experiment or file accessions, e.g., &#91;&quot;ENCSR133RZO&quot;, &quot;ENCFF635JIA&quot;&#93;) |
| `description` | str | &#x27;&#x27; | What this derived file contains |
| `file_type` | str | &#x27;&#x27; | Type of file (e.g., &quot;filtered_peaks&quot;, &quot;merged_signal&quot;, &quot;differential&quot;) |
| `tool_used` | str | &#x27;&#x27; | Tool/software used to create it (e.g., &quot;bedtools intersect&quot;, &quot;DESeq2&quot;) |
| `parameters` | str | &#x27;&#x27; | Parameters or command used |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Summarize tracked experiment collections by field {#dnakit.encode.summarize_collection-parameters}

`dnakit.encode.summarize_collection(assay_title=None, organism=None, organ=None, *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `assay_title` | str &#124; None | None | Filter by assay type (partial match) |
| `organism` | str &#124; None | None | Filter by organism (partial match) |
| `organ` | str &#124; None | None | Filter by organ (partial match) |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |


## Track an ENCODE experiment and its provenance locally {#dnakit.encode.track_experiment-parameters}

`dnakit.encode.track_experiment(accession, fetch_publications=True, fetch_pipelines=True, notes='', *, verbose=True)`

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `accession` | str | Required | ENCODE experiment accession (e.g., &quot;ENCSR133RZO&quot;) |
| `fetch_publications` | bool | True | Also fetch and store publications/citations (default True) |
| `fetch_pipelines` | bool | True | Also fetch and store pipeline/analysis info (default True) |
| `notes` | str | &#x27;&#x27; | Optional notes to attach to this experiment |
| `verbose` | bool | True | Show DNAKit progress messages without changing results. |
