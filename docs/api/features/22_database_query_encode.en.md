# ENCODE experiment and file queries

Find ENCODE experiments, samples and downloadable files.

## 1) Filter ENCODE objects and files {#search-encode-search}

- **Purpose:** Find ENCODE records by object type, keyword and field filters, such as released bigWig files.
- **Method:** Submit ENCODE search filters and normalize matching objects and file metadata into QueryResult.
- **API:** `dnakit.search.encode_search(object_type[optional], search_term[optional], filters[optional], limit[optional], config[optional], verbose[optional])`
- **Input:** Optional object type, search term, filter dictionary, record limit and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import encode_search

result = encode_search(
    object_type="File",
    filters={"file_format": "bigWig", "status": "released"},
    limit=20,
    verbose=False,
)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
regulation ENCODE <record count returned>
```

<!-- BEGIN GENERATED API: encode -->

<span id="encode"></span>

## 2) Retrieve ENCODE experiment details by accession {#dnakit.encode.get_experiment}

- **Purpose:** Retrieve ENCODE experiment details by accession.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.get_experiment(accession, *, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Python return type: str.  JSON with full experiment details and file listing.
- **Example code:**

```python
import dnakit

result = dnakit.encode.get_experiment(
    accession='ENCSR133RZO',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with full experiment details and file listing.

[Full parameters and usage notes](../../faq/api/22_database_query_encode.md#dnakit.encode.get_experiment-parameters)

## 3) Retrieve live ENCODE filter counts {#dnakit.encode.get_facets}

- **Purpose:** Retrieve live ENCODE filter counts.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.get_facets(search_type='Experiment', assay_title=None, organism=None, organ=None, biosample_type=None, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Python return type: str.  JSON with facet names and their term counts.
- **Example code:**

```python
import dnakit

result = dnakit.encode.get_facets(
    assay_title='ATAC-seq',
    organism='Homo sapiens',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with facet names and their term counts.

[Full parameters and usage notes](../../faq/api/22_database_query_encode.md#dnakit.encode.get_facets-parameters)

## 4) Retrieve ENCODE file details {#dnakit.encode.get_file_info}

- **Purpose:** Retrieve ENCODE file details.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.get_file_info(accession, *, verbose=True)`
- **Input:** Provide `accession`. Retrieve a real file accession with the file query API and replace `ENCFF_REPLACE_ME`.
- **Output:** Python return type: str.  JSON with file metadata including format, size, download URL, MD5, assembly, etc.
- **Example code:**

```python
import dnakit

result = dnakit.encode.get_file_info(
    accession='ENCFF_REPLACE_ME',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with file metadata including format, size, download URL, MD5, assembly, etc.

[Full parameters and usage notes](../../faq/api/22_database_query_encode.md#dnakit.encode.get_file_info-parameters)

## 5) Retrieve available ENCODE search filters {#dnakit.encode.get_metadata}

- **Purpose:** Retrieve available ENCODE search filters.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.get_metadata(metadata_type, *, verbose=True)`
- **Input:** Provide `metadata_type`.
- **Output:** Python return type: str.  JSON list of valid values for the specified metadata type.
- **Example code:**

```python
import dnakit

result = dnakit.encode.get_metadata(
    metadata_type='assays',
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON list of valid values for the specified metadata type.

[Full parameters and usage notes](../../faq/api/22_database_query_encode.md#dnakit.encode.get_metadata-parameters)

## 6) List files belonging to an ENCODE experiment {#dnakit.encode.list_files}

- **Purpose:** List files belonging to an ENCODE experiment.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.list_files(experiment_accession, file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, status=None, preferred_default=None, limit=200, *, verbose=True)`
- **Input:** Provide `experiment_accession`.
- **Output:** Python return type: str.  JSON list of files with accession, format, size, download URL, and metadata.
- **Example code:**

```python
import dnakit

result = dnakit.encode.list_files(
    experiment_accession='ENCSR133RZO',
    file_format='bed',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON list of files with accession, format, size, download URL, and metadata.

[Full parameters and usage notes](../../faq/api/22_database_query_encode.md#dnakit.encode.list_files-parameters)

## 7) Search ENCODE experiments by experimental conditions {#dnakit.encode.search_experiments}

- **Purpose:** Search ENCODE experiments by experimental conditions.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.search_experiments(assay_title=None, organism='Homo sapiens', organ=None, biosample_type=None, biosample_term_name=None, target=None, status='released', lab=None, award=None, assembly=None, replication_type=None, life_stage=None, sex=None, treatment=None, genetic_modification=None, perturbed=None, search_term=None, date_released_from=None, date_released_to=None, limit=25, offset=0, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Python return type: str.  JSON with experiment results, total count, and pagination info.
- **Example code:**

```python
import dnakit

result = dnakit.encode.search_experiments(
    assay_title='ATAC-seq',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with experiment results, total count, and pagination info.

[Full parameters and usage notes](../../faq/api/22_database_query_encode.md#dnakit.encode.search_experiments-parameters)

## 8) Search ENCODE files using experiment and file filters {#dnakit.encode.search_files}

- **Purpose:** Search ENCODE files using experiment and file filters.
- **Method:** Query ENCODE with the supplied filters and parse experiment, file or download results.
- **API：** `dnakit.encode.search_files(file_format=None, file_type=None, output_type=None, output_category=None, assembly=None, assay_title=None, organism=None, organ=None, biosample_type=None, target=None, status='released', preferred_default=None, search_term=None, limit=25, offset=0, *, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Python return type: str.  JSON with file results, total count, and pagination info.
- **Example code:**

```python
import dnakit

result = dnakit.encode.search_files(
    file_format='bed',
    assembly='GRCh38',
    assay_title='ATAC-seq',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** Python return type: str.  JSON with file results, total count, and pagination info.

[Full parameters and usage notes](../../faq/api/22_database_query_encode.md#dnakit.encode.search_files-parameters)

<!-- END GENERATED API: encode -->
