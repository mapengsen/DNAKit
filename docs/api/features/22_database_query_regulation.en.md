# Expression, regulation and cancer queries

Query gene expression, regulation and cancer datasets.

## 1) Query GEO expression data by keyword {#search-expression}

- **Purpose:** Find GEO expression records matching a keyword query and retrieve summaries and source metadata.
- **Method:** Query GEO through NCBI Entrez and return a normalized QueryResult.
- **API:** `dnakit.search.expression(term[required], retmax[optional], config[optional], verbose[optional])`
- **Input:** An Entrez search term; optional retmax record limit and SearchConfig. Requires network access.
- **Output:** Returns `QueryResult` with query `records`, source `provider`, and pagination information.
- **Example code:**

```python
from dnakit.search import expression

result = expression("single cell RNA sequencing[All Fields]", retmax=20, verbose=False)
print(result.query_type, result.provider, len(result.records))
```

- **Example output:** Return format only; the network query was not run. The live database determines the record count.

```text
expression NCBI Entrez <record count returned>
```

<!-- BEGIN GENERATED API: regulation -->

<span id="regulation"></span>

## 2) Browse the cBioPortal cancer study catalogue {#dnakit.databases.c_bio_portal_get_cancer_studies}

- **Purpose:** Browse the cBioPortal cancer study catalogue.
- **Method:** Use `CBioPortal` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.c_bio_portal_get_cancer_studies(*, limit=None, offset=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `studyId`, `name`, `description`, `cancerTypeId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_cancer_studies(
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "studyId": "<string>",
    "name": "<string>",
    "description": "<string>",
    "cancerTypeId": "<string>",
    "publicStudy": "<boolean>",
    "pmid": "<string>",
    "citation": "<string>",
    "groups": [
      "<value>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_cancer_studies-parameters)

## 3) Retrieve cBioPortal gene information by Entrez ID {#dnakit.databases.c_bio_portal_get_gene_info}

- **Purpose:** Retrieve cBioPortal gene information by Entrez ID.
- **Method:** Use `CBioPortal` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.c_bio_portal_get_gene_info(*, entrez_gene_id, verbose=True)`
- **Input:** Provide `entrez_gene_id`.
- **Output:** Returns a mapping with `entrezGeneId`, `hugoGeneSymbol`, `type`, `cytoband` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_gene_info(
    entrez_gene_id=672,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "entrezGeneId": "<integer>",
  "hugoGeneSymbol": "<string>",
  "type": "<string>",
  "cytoband": "<string>",
  "length": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_gene_info-parameters)

## 4) Retrieve molecular profiles for a cancer study {#dnakit.databases.c_bio_portal_get_molecular_profiles}

- **Purpose:** Retrieve molecular profiles for a cancer study.
- **Method:** Use `CBioPortal` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.c_bio_portal_get_molecular_profiles(*, study_id, verbose=True)`
- **Input:** Provide `study_id`.
- **Output:** Returns a list of records; each record includes `molecularProfileId`, `molecularAlterationType`, `datatype`, `name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_molecular_profiles(
    study_id='brca_tcga',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "molecularProfileId": "<string>",
    "molecularAlterationType": "<string>",
    "datatype": "<string>",
    "name": "<string>",
    "description": "<string>",
    "showProfileInAnalysisTab": "<boolean>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_molecular_profiles-parameters)

## 5) Retrieve gene mutations from a cancer study {#dnakit.databases.c_bio_portal_get_mutations}

- **Purpose:** Retrieve gene mutations from a cancer study.
- **Method:** Use `CBioPortal` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.c_bio_portal_get_mutations(*, study_id, gene_list, sample_list_id=None, verbose=True)`
- **Input:** Provide `study_id`, `gene_list`.
- **Output:** Returns a list of records; each record includes `molecularProfileId`, `sampleId`, `patientId`, `studyId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_mutations(
    study_id='brca_tcga',
    gene_list='BRCA2',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "molecularProfileId": "<string>",
    "sampleId": "<string>",
    "patientId": "<string>",
    "studyId": "<string>",
    "entrezGeneId": "<integer>",
    "gene": {
      "entrezGeneId": "<integer>",
      "hugoGeneSymbol": "<string>",
      "type": "<string>"
    },
    "mutationType": "<string>",
    "proteinChange": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_mutations-parameters)

## 6) Retrieve samples in a cancer study {#dnakit.databases.c_bio_portal_get_samples}

- **Purpose:** Retrieve samples in a cancer study.
- **Method:** Use `CBioPortal` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.c_bio_portal_get_samples(*, study_id, page_size=None, page_number=None, verbose=True)`
- **Input:** Provide `study_id`.
- **Output:** Returns a list of records; each record includes `sampleId`, `sampleType`, `patientId`, `studyId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.c_bio_portal_get_samples(
    study_id='brca_tcga',
    page_size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "sampleId": "<string>",
    "sampleType": "<string>",
    "patientId": "<string>",
    "studyId": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.c_bio_portal_get_samples-parameters)

## 7) Retrieve original CELLxGENE H5AD files {#dnakit.databases.cellx_gene_download_h5ad}

- **Purpose:** Retrieve original CELLxGENE H5AD files.
- **Method:** Locate the dataset through CELLxGENE Census and download its H5AD file.
- **API：** `dnakit.databases.cellx_gene_download_h5ad(*, operation=None, dataset_id, output_path=None, census_version=None, verbose=True)`
- **Input:** Provide `dataset_id`. Requires `cellxgene-census`. Prepare the local files and external programs used by the example, and replace paths with your input or output locations.
- **Output:** Returns a mapping with `dataset_id`, `output_path`, `uri`, `message` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cellx_gene_download_h5ad(
    operation='download_h5ad',
    dataset_id='4eb29386-de81-452f-b3c0-e00844e8c7fd',
    census_version='stable',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "dataset_id": "<string>",
  "output_path": "<string>",
  "uri": "<string>",
  "message": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_download_h5ad-parameters)

## 8) Query filtered single-cell metadata {#dnakit.databases.cellx_gene_get_cell_metadata}

- **Purpose:** Query filtered single-cell metadata.
- **Method:** Use `CELLxGENECensus` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cellx_gene_get_cell_metadata(*, operation=None, organism=None, obs_value_filter=None, column_names=None, census_version=None, gene=None, tissue=None, max_cells=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Requires `cellxgene-census`.
- **Output:** Returns a list of records; each record includes `cell_type`, `tissue_general`, `disease`, `assay` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cellx_gene_get_cell_metadata(
    operation='get_obs_metadata',
    organism='Homo sapiens',
    obs_value_filter="tissue_general == 'brain' and disease == 'normal'",
    census_version='stable',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "cell_type": "<string>",
    "tissue_general": "<string>",
    "disease": "<string>",
    "assay": "<string>",
    "sex": "<string>",
    "donor_id": "<string>",
    "dataset_id": "<string>",
    "soma_joinid": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_get_cell_metadata-parameters)

## 9) Retrieve precomputed cell embeddings {#dnakit.databases.cellx_gene_get_embeddings}

- **Purpose:** Retrieve precomputed cell embeddings.
- **Method:** Use `CELLxGENECensus` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cellx_gene_get_embeddings(*, operation=None, organism=None, embedding_name=None, census_version=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Requires `cellxgene-census`.
- **Output:** Returns a mapping with `message`, `organism`, `embedding_name`, `available` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cellx_gene_get_embeddings(
    operation='get_embeddings',
    organism='Homo sapiens',
    census_version='stable',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "message": "<string>",
  "organism": "<string>",
  "embedding_name": "<string>",
  "available": "<boolean>",
  "shape": [
    "<value>"
  ],
  "available_embeddings": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_get_embeddings-parameters)

## 10) Query expression data with cell and gene filters {#dnakit.databases.cellx_gene_get_expression_data}

- **Purpose:** Query expression data with cell and gene filters.
- **Method:** Use `CELLxGENECensus` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cellx_gene_get_expression_data(*, operation=None, organism=None, obs_value_filter=None, var_value_filter=None, obs_column_names=None, var_column_names=None, census_version=None, gene=None, tissue=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Requires `cellxgene-census`.
- **Output:** Returns a mapping with `n_obs`, `n_vars`, `obs_columns`, `var_columns` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cellx_gene_get_expression_data(
    operation='get_anndata',
    organism='Homo sapiens',
    obs_value_filter="tissue_general == 'brain' and disease == 'normal'",
    var_value_filter="feature_name in ['SOX2', 'GFAP', 'OLIG2']",
    census_version='stable',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "n_obs": "<integer>",
  "n_vars": "<integer>",
  "obs_columns": [
    "<value>"
  ],
  "var_columns": [
    "<value>"
  ],
  "message": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.cellx_gene_get_expression_data-parameters)

## 11) Prepare ChIP-Atlas enrichment submission instructions {#dnakit.databases.ch_ipatlas_enrichment_analysis}

- **Purpose:** Prepare ChIP-Atlas enrichment submission instructions from BED, motif or gene-list input; this function does not calculate enrichment statistics.
- **Method:** Use `ChIPAtlas` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_ipatlas_enrichment_analysis(*, operation=None, bed_data=None, motif=None, gene_list=None, genome=None, antigen_class=None, cell_type_class=None, threshold=None, distance=None, gene=None, cell_type=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Mapping with `status=requires_submission`, `analysis_completed=False` and task submission instructions; the enrichment test has not been run.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_enrichment_analysis(
    operation='enrichment_analysis',
    gene_list=['TP53', 'MDM2', 'CDKN1A'],
    genome='hg38',
    distance='5000',
    verbose=False,
)
print(result)
```

- **Example output:** Mapping with `status=requires_submission`, `analysis_completed=False` and task submission instructions; the enrichment test has not been run.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "message": "<string>",
  "genes": [
    "<string>"
  ],
  "distance_from_tss": "<string>",
  "url": "<string>",
  "note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_enrichment_analysis-parameters)

## 12) Search ChIP-Atlas experiment metadata {#dnakit.databases.ch_ipatlas_get_experiments}

- **Purpose:** Search ChIP-Atlas experiment metadata.
- **Method:** Use `ChIPAtlas` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_ipatlas_get_experiments(*, operation=None, genome=None, antigen=None, cell_type=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `message`, `metadata_file_url`, `web_search_url`, `file_size` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_get_experiments(
    operation='get_experiment_list',
    genome='hg38',
    antigen='CTCF',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "message": "<string>",
  "metadata_file_url": "<string>",
  "web_search_url": "<string>",
  "file_size": "<string>",
  "total_experiments": "<string>",
  "filters_requested": "<object>",
  "recommendation": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_get_experiments-parameters)

## 13) Retrieve ChIP-Atlas peak resources and local overlaps {#dnakit.databases.ch_ipatlas_get_peak_data}

- **Purpose:** Query ChIP-Atlas experiment peak resources. The gene mode resolves a GRCh38 locus and intersects a supplied peak_file BED, or returns a data retrieval plan when no BED is supplied.
- **Method:** Use `ChIPAtlas` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_ipatlas_get_peak_data(*, operation=None, experiment_id=None, genome=None, format=None, threshold=None, gene=None, experiment_type=None, limit=None, peak_file=None, flank_bp=None, max_file_bytes=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `experiment_id`, `genome`, `format`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_get_peak_data(
    operation='get_peak_data',
    experiment_id='SRX097088',
    genome='hg19',
    format='bigwig',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "experiment_id": "<string>",
  "genome": "<string>",
  "format": "<string>",
  "url": "<string>",
  "message": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_get_peak_data-parameters)

## 14) Search ChIP-Atlas datasets by antigen or cell type {#dnakit.databases.ch_ipatlas_search_datasets}

- **Purpose:** Search ChIP-Atlas datasets by antigen or cell type.
- **Method:** Use `ChIPAtlas` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_ipatlas_search_datasets(*, operation=None, antigen=None, cell_type=None, genome=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `search_key`, `search_value`, `genome`, `num_results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_ipatlas_search_datasets(
    operation='search_datasets',
    antigen='CTCF',
    genome='hg38',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "search_key": "<string>",
  "search_value": "<string>",
  "genome": "<string>",
  "num_results": "<integer>",
  "results": [
    {
      "genome": "<string>",
      "class": "<string>",
      "name": "<string>",
      "num_experiments": "<string>",
      "experiment_ids": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.ch_ipatlas_search_datasets-parameters)

## 15) Retrieve gene dependency information from CRISPR screens {#dnakit.databases.dep_map_get_gene_dependencies}

- **Purpose:** Retrieve gene dependency information from CRISPR screens.
- **Method:** Use `DepMap` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dep_map_get_gene_dependencies(*, gene_symbol, model_id=None, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** Returns a mapping with `status`, `data`, `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dep_map_get_gene_dependencies(
    gene_symbol='EGFR',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "gene_symbol": "<string>",
    "genes": [
      "<value>"
    ]
  },
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.dep_map_get_gene_dependencies-parameters)

## 16) Search cancer cell lines by name {#dnakit.databases.dep_map_search_cell_lines}

- **Purpose:** Search cancer cell lines by name.
- **Method:** Use `DepMap` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dep_map_search_cell_lines(*, query, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `status`, `data`, `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dep_map_search_cell_lines(
    query='A549',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": {
    "query": "<string>",
    "cell_lines": [
      "<value>"
    ],
    "count": "<integer>"
  },
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.dep_map_search_cell_lines-parameters)

## 17) Retrieve variant–expression associations from an eQTL dataset {#dnakit.databases.e_qtl_get_associations}

- **Purpose:** Retrieve variant–expression associations from an eQTL dataset.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.e_qtl_get_associations(*, dataset_id, size=None, gene_id=None, variant=None, molecular_trait_id=None, verbose=True)`
- **Input:** Provide `dataset_id`.
- **Output:** Returns a list of records; each record includes `molecular_trait_id`, `gene_id`, `variant`, `rsid` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.e_qtl_get_associations(
    dataset_id='QTD000001',
    size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "molecular_trait_id": "<string>",
    "gene_id": "<string>",
    "variant": "<string>",
    "rsid": "<string>",
    "chromosome": "<string>",
    "position": "<integer>",
    "ref": "<string>",
    "alt": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.e_qtl_get_associations-parameters)

## 18) List available eQTL study and tissue datasets {#dnakit.databases.e_qtl_list_datasets}

- **Purpose:** List available eQTL study and tissue datasets.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.e_qtl_list_datasets(*, size=None, quant_method=None, tissue_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `dataset_id`, `study_id`, `study_label`, `quant_method` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.e_qtl_list_datasets(
    size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "dataset_id": "<string>",
    "study_id": "<string>",
    "study_label": "<string>",
    "quant_method": "<string>",
    "tissue_id": "<string>",
    "tissue_label": "<string>",
    "condition_label": "<string>",
    "sample_size": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.e_qtl_list_datasets-parameters)

## 19) Retrieve ENCODE biosample metadata {#dnakit.databases.encode_get_biosample}

- **Purpose:** Retrieve ENCODE biosample metadata.
- **Method:** Use `ENCODE` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_get_biosample(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a mapping with `accession`, `url`, `aliases`, `schema_version` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_get_biosample(
    accession='ENCBS000AAA',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "accession": "<string>",
  "url": "<string>",
  "aliases": [
    "<string>"
  ],
  "schema_version": "<string>",
  "status": "<string>",
  "lab": "<string>",
  "award": "<string>",
  "date_created": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_get_biosample-parameters)

## 20) Search chromatin state segmentation annotations {#dnakit.databases.encode_get_chromatin_state}

- **Purpose:** Search chromatin state segmentation annotations.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_get_chromatin_state(*, biosample_term_name=None, organism=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `annotations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_get_chromatin_state(
    biosample_term_name='K562',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "annotations": [
    {
      "accession": "<string>",
      "annotation_type": "<string>",
      "description": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_get_chromatin_state-parameters)

## 21) Retrieve ENCODE experiment metadata and analysis details {#dnakit.databases.encode_get_experiment}

- **Purpose:** Retrieve ENCODE experiment metadata and analysis details.
- **Method:** Use `ENCODE` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_get_experiment(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a mapping with `assay_term_name`, `biosample_ontology`, `documents`, `references` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_get_experiment(
    accession='ENCSR000AKS',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "assay_term_name": "<string>",
  "biosample_ontology": "<string>",
  "documents": [
    "<value>"
  ],
  "references": [
    "<string>"
  ],
  "schema_version": "<string>",
  "accession": "<string>",
  "alternate_accessions": [
    "<value>"
  ],
  "analyses": [
    "<string>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_get_experiment-parameters)

## 22) List ENCODE files by format and assay {#dnakit.databases.encode_list_files}

- **Purpose:** List ENCODE files by format and assay.
- **Method:** Use `ENCODEFiles` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_list_files(*, file_type=None, assay_title=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `@graph`, `@type`, `all`, `notification` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_list_files(
    file_type='fastq',
    limit=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "@graph": [
    {
      "@id": "<string>",
      "@type": "<array>",
      "accession": "<string>",
      "assay_term_name": "<string>",
      "award": "<object>",
      "biological_replicates": "<array>",
      "biosample_ontology": "<object>",
      "dataset": "<string>"
    }
  ],
  "@type": [
    "<string>"
  ],
  "all": "<string>",
  "notification": "<string>",
  "title": "<string>",
  "total": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_list_files-parameters)

## 23) Search cis-regulatory and epigenomic annotations {#dnakit.databases.encode_search_annotations}

- **Purpose:** Search cis-regulatory and epigenomic annotations.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_search_annotations(*, annotation_type=None, biosample_term_name=None, organism=None, assembly=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `annotations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_search_annotations(
    annotation_type='candidate Cis-Regulatory Elements',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "annotations": [
    {
      "accession": "<string>",
      "annotation_type": "<string>",
      "description": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_annotations-parameters)

## 24) Search ENCODE biosamples by organism and treatment {#dnakit.databases.encode_search_biosamples}

- **Purpose:** Search ENCODE biosamples by organism and treatment.
- **Method:** Use `ENCODESearchBiosamples` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_search_biosamples(*, organism=None, biosample_type=None, treatment=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `@graph`, `@type`, `notification`, `title` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_search_biosamples(
    organism='Homo sapiens',
    limit=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "@graph": [
    {
      "@id": "<string>",
      "@type": "<array>",
      "accession": "<string>",
      "age": "<string>",
      "applied_modifications": "<array>",
      "biosample_ontology": "<object>",
      "lab": "<object>",
      "life_stage": "<string>"
    }
  ],
  "@type": [
    "<string>"
  ],
  "notification": "<string>",
  "title": "<string>",
  "total": "<integer>",
  "all": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_biosamples-parameters)

## 25) Search ATAC-seq and DNase-seq experiments {#dnakit.databases.encode_search_chromatin_accessibility}

- **Purpose:** Search ATAC-seq and DNase-seq experiments.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_search_chromatin_accessibility(*, assay_type=None, biosample_term_name=None, organism=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `experiments` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_search_chromatin_accessibility(
    assay_type='ATAC-seq',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "experiments": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>",
      "lab": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_chromatin_accessibility-parameters)

## 26) Search ENCODE experiments through the database module {#dnakit.databases.encode_search_experiments}

- **Purpose:** Search ENCODE experiments through the database module.
- **Method:** Use `ENCODESearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_search_experiments(*, assay_title=None, target=None, organism=None, status=None, limit=None, biosample=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns ENCODE experiment search results; ordinary queries preserve the response in data, while biosample queries return the parsed ENCODE search result.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_search_experiments(
    assay_title='ATAC-seq',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** Returns ENCODE experiment search results; ordinary queries preserve the response in data, while biosample queries return the parsed ENCODE search result.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_experiments-parameters)

## 27) Search histone ChIP-seq experiments by modification {#dnakit.databases.encode_search_histone_experiments}

- **Purpose:** Search histone ChIP-seq experiments by modification.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_search_histone_experiments(*, histone_mark=None, target=None, biosample_term_name=None, biosample=None, organism=None, limit=None, biosample_term=None, cell_type=None, tissue=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `experiments` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_search_histone_experiments(
    histone_mark='H3K27ac',
    organism='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "experiments": [
    {
      "accession": "<string>",
      "histone_mark": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>",
      "lab": "<string>",
      "date_released": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_histone_experiments-parameters)

## 28) Search ENCODE RNA-seq experiments by sample and assay {#dnakit.databases.encode_search_rnaseq_experiments}

- **Purpose:** Search ENCODE RNA-seq experiments by sample and assay.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.encode_search_rnaseq_experiments(*, biosample_term_name=None, biosample=None, cell_type=None, tissue=None, organism=None, assay_type=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `experiments` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.encode_search_rnaseq_experiments(
    biosample_term_name='K562',
    assay_type='total RNA-seq',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "experiments": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "biosample_summary": "<string>",
      "status": "<string>",
      "lab": "<string>",
      "date_released": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.encode_search_rnaseq_experiments-parameters)

## 29) Search differential expression experiment records {#dnakit.databases.expression_atlas_search_differential}

- **Purpose:** Search differential expression experiment records.
- **Method:** Use `ExpressionAtlas` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.expression_atlas_search_differential(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `gene`, `condition`, `species`, `experiments` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.expression_atlas_search_differential(
    gene='TP53',
    species='homo sapiens',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene": "<string>",
  "condition": "<string>",
  "species": "<string>",
  "experiments": [
    {
      "experiment_accession": "<string>",
      "experiment_type": "<string>",
      "experiment_description": "<string>",
      "species": "<string>",
      "num_assays": "<integer>",
      "gene_mentioned": "<boolean>",
      "factors": "<array>"
    }
  ],
  "experiment_count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.expression_atlas_search_differential-parameters)

## 30) Search baseline and differential expression experiments {#dnakit.databases.expression_atlas_search_experiments}

- **Purpose:** Search baseline and differential expression experiments.
- **Method:** Use `ExpressionAtlas` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.expression_atlas_search_experiments(*, gene=None, condition=None, species=None, limit=None, offset=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `gene`, `condition`, `species`, `experiments` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.expression_atlas_search_experiments(
    gene='TP53',
    species='homo sapiens',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene": "<string>",
  "condition": "<string>",
  "species": "<string>",
  "experiments": [
    {
      "experiment_accession": "<string>",
      "experiment_type": "<string>",
      "experiment_description": "<string>",
      "species": "<string>",
      "num_assays": "<integer>",
      "gene_mentioned": "<boolean>"
    }
  ],
  "total_count": "<integer>",
  "gene_specific_count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.expression_atlas_search_experiments-parameters)

## 31) Search Hi-C and chromatin conformation data {#dnakit.databases.four_dn_search_data}

- **Purpose:** Search Hi-C and chromatin conformation data.
- **Method:** Use `FourDN` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.four_dn_search_data(*, operation=None, query=None, item_type=None, file_type=None, assay_title=None, biosource_name=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `num_results`, `total`, `results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.four_dn_search_data(
    operation='search',
    query='Hi-C',
    item_type='File',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "num_results": "<integer>",
  "total": "<integer>",
  "results": [
    "<value>"
  ],
  "search_url": "<string>",
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.four_dn_search_data-parameters)

## 32) Retrieve GDC cancer case clinical data {#dnakit.databases.gdc_get_clinical_data}

- **Purpose:** Retrieve GDC cancer case clinical data.
- **Method:** Use `GDCClinicalData` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gdc_get_clinical_data(*, project_id=None, primary_site=None, disease_type=None, vital_status=None, gender=None, size=None, offset=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** data.cases contains case, diagnosis and treatment records; data.pagination contains paging information.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gdc_get_clinical_data(
    project_id='TCGA-BRCA',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data.cases contains case, diagnosis and treatment records; data.pagination contains paging information.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_clinical_data-parameters)

## 33) Retrieve pan-cancer gene mutation frequencies {#dnakit.databases.gdc_get_mutation_frequency}

- **Purpose:** Retrieve pan-cancer gene mutation frequencies.
- **Method:** Use `GDCMutationFrequency` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gdc_get_mutation_frequency(*, gene_symbol=None, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `source`, `gene`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gdc_get_mutation_frequency(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "success",
  "source": "GDC",
  "gene": "<str>",
  "data": {
    "gene_info": "<dict>",
    "total_ssm_occurrences": "<int>",
    "is_cancer_gene_census": "<bool | None>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_mutation_frequency-parameters)

## 34) Retrieve GDC somatic mutations by gene {#dnakit.databases.gdc_get_ssm_by_gene}

- **Purpose:** Retrieve GDC somatic mutations by gene.
- **Method:** Use `GDCSSM` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gdc_get_ssm_by_gene(*, gene_symbol, project_id=None, size=None, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** data contains GDC somatic mutation records; gene identifies the queried gene.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gdc_get_ssm_by_gene(
    gene_symbol='TP53',
    project_id='TCGA-BRCA',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data contains GDC somatic mutation records; gene identifies the queried gene.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_ssm_by_gene-parameters)

## 35) Retrieve cancer cohort survival and censoring data {#dnakit.databases.gdc_get_survival}

- **Purpose:** Retrieve cancer cohort survival and censoring data.
- **Method:** Use `GDCSurvival` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gdc_get_survival(*, project_id, gene_symbol=None, verbose=True)`
- **Input:** Provide `project_id`.
- **Output:** data contains follow-up times, censoring, survival estimates and cohort summaries; donors includes at most 50 records.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gdc_get_survival(
    project_id='TCGA-BRCA',
    verbose=False,
)
print(result)
```

- **Example output:** data contains follow-up times, censoring, survival estimates and cohort summaries; donors includes at most 50 records.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_get_survival-parameters)

## 36) List GDC files by data type {#dnakit.databases.gdc_list_files}

- **Purpose:** List GDC files by data type.
- **Method:** Use `GDCFiles` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gdc_list_files(*, data_type=None, size=None, offset=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** data preserves the GDC file listing and pagination response, including file identifiers and metadata.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gdc_list_files(
    data_type='Gene Expression Quantification',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data preserves the GDC file listing and pagination response, including file identifiers and metadata.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_list_files-parameters)

## 37) List GDC cancer projects and summaries {#dnakit.databases.gdc_list_projects}

- **Purpose:** List GDC cancer projects and summaries.
- **Method:** Use `GDCProjects` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gdc_list_projects(*, program=None, size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** data preserves the GDC project listing and project statistics.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gdc_list_projects(
    program='TCGA',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data preserves the GDC project listing and project statistics.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_list_projects-parameters)

## 38) Search GDC cancer cases with filters {#dnakit.databases.gdc_search_cases}

- **Purpose:** Search GDC cancer cases with filters.
- **Method:** Use `GDCCases` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gdc_search_cases(*, project_id=None, size=None, offset=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** data preserves the GDC case search response and pagination information.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gdc_search_cases(
    project_id='TCGA-BRCA',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data preserves the GDC case search response and pagination information.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gdc_search_cases-parameters)

## 39) Search GEO ATAC-seq datasets {#dnakit.databases.geo_search_atacseq_datasets}

- **Purpose:** Search GEO ATAC-seq datasets.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.geo_search_atacseq_datasets(*, query=None, organism=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `datasets` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.geo_search_atacseq_datasets(
    query='T cell differentiation',
    organism='Homo sapiens',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "datasets": [
    {
      "accession": "<string>",
      "title": "<string>",
      "summary": "<string>",
      "organism": "<string>",
      "n_samples": "<integer>",
      "date_published": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.geo_search_atacseq_datasets-parameters)

## 40) Search GEO ChIP-seq datasets {#dnakit.databases.geo_search_chipseq_datasets}

- **Purpose:** Search GEO ChIP-seq datasets.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.geo_search_chipseq_datasets(*, query, organism=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `total`, `datasets` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.geo_search_chipseq_datasets(
    query='H3K27ac liver',
    organism='Homo sapiens',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "datasets": [
    {
      "accession": "<string>",
      "title": "<string>",
      "summary": "<string>",
      "organism": "<string>",
      "n_samples": "<integer>",
      "date_published": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.geo_search_chipseq_datasets-parameters)

## 41) Search GEO RNA-seq datasets {#dnakit.databases.geo_search_rnaseq_datasets}

- **Purpose:** Search GEO RNA-seq datasets.
- **Method:** Use `Epigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.geo_search_rnaseq_datasets(*, query=None, organism=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total`, `datasets` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.geo_search_rnaseq_datasets(
    query='breast cancer BRCA1',
    organism='Homo sapiens',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total": "<integer>",
  "datasets": [
    {
      "accession": "<string>",
      "title": "<string>",
      "summary": "<string>",
      "organism": "<string>",
      "n_samples": "<integer>",
      "date_published": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.geo_search_rnaseq_datasets-parameters)

## 42) Calculate a gene–variant eQTL association in a tissue {#dnakit.databases.gtex_calculate_eqtl}

- **Purpose:** Calculate a gene–variant eQTL association in a tissue.
- **Method:** Use `GTExV2` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtex_calculate_eqtl(*, operation=None, gencode_id, variant_id, tissue_site_detail_id, dataset_id=None, verbose=True)`
- **Input:** Provide `gencode_id`, `variant_id`, `tissue_site_detail_id`.
- **Output:** Returns a mapping with `data`, `error`, `gencodeId`, `geneSymbol` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtex_calculate_eqtl(
    operation='calculate_eqtl',
    gencode_id='ENSG00000227232.5',
    variant_id='chr1_64764_C_T_b38',
    tissue_site_detail_id='Whole_Blood',
    dataset_id='gtex_v8',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": [
    "<number>"
  ],
  "error": "<number>",
  "gencodeId": "<string>",
  "geneSymbol": "<string>",
  "genotypes": [
    "<integer>"
  ],
  "hetCount": "<integer>",
  "homoAltCount": "<integer>",
  "homoRefCount": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_calculate_eqtl-parameters)

## 43) Retrieve GTEx gene expression summaries {#dnakit.databases.gtex_get_expression_summary}

- **Purpose:** Retrieve GTEx gene expression summaries.
- **Method:** Use `GTExExpression` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtex_get_expression_summary(*, gene_symbol=None, ensembl_gene_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns the resolved gene identifier and tissue expression data; note explains missing expression records.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtex_get_expression_summary(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** Returns the resolved gene identifier and tissue expression data; note explains missing expression records.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_expression_summary-parameters)

## 44) Retrieve median gene expression across tissues {#dnakit.databases.gtex_get_median_gene_expression}

- **Purpose:** Retrieve median gene expression across tissues.
- **Method:** Use `GTExV2` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtex_get_median_gene_expression(*, operation=None, gencode_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `median`, `tissueSiteDetailId`, `ontologyId`, `datasetId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtex_get_median_gene_expression(
    operation='get_median_gene_expression',
    gencode_id='ENSG00000141510.18',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "median": "<number>",
    "tissueSiteDetailId": "<string>",
    "ontologyId": "<string>",
    "datasetId": "<string>",
    "gencodeId": "<string>",
    "geneSymbol": "<string>",
    "unit": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_median_gene_expression-parameters)

## 45) Retrieve multi-tissue eQTL meta-analysis results {#dnakit.databases.gtex_get_multi_tissue_eqtls}

- **Purpose:** Retrieve multi-tissue eQTL meta-analysis results.
- **Method:** Use `GTExV2` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtex_get_multi_tissue_eqtls(*, operation=None, gencode_id, variant_id=None, dataset_id=None, page=None, items_per_page=None, verbose=True)`
- **Input:** Provide `gencode_id`.
- **Output:** Returns a list of records; each record includes `gencodeId`, `datasetId`, `metaP`, `variantId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtex_get_multi_tissue_eqtls(
    operation='get_multi_tissue_eqtls',
    gencode_id='ENSG00000141510.16',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "gencodeId": "<string>",
    "datasetId": "<string>",
    "metaP": "<number>",
    "variantId": "<string>",
    "tissues": {
      "Thyroid": "<object>",
      "Testis": "<object>",
      "Small_Intestine_Terminal_Ileum": "<object>",
      "Brain_Frontal_Cortex_BA9": "<object>",
      "Skin_Not_Sun_Exposed_Suprapubic": "<object>",
      "Vagina": "<object>",
      "Whole_Blood": "<object>",
      "Breast_Mammary_Tissue": "<object>"
    }
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_multi_tissue_eqtls-parameters)

## 46) Retrieve significant single-tissue eQTL associations {#dnakit.databases.gtex_get_single_tissue_eqtls}

- **Purpose:** Retrieve significant single-tissue eQTL associations.
- **Method:** Use `GTExV2` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtex_get_single_tissue_eqtls(*, operation=None, gencode_id=None, variant_id=None, tissue_site_detail_id=None, dataset_id=None, page=None, items_per_page=None, gene_symbol=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `snpId`, `pos`, `snpIdUpper`, `variantId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtex_get_single_tissue_eqtls(
    operation='get_single_tissue_eqtls',
    gencode_id=['ENSG00000227232.5'],
    tissue_site_detail_id=['Whole_Blood'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "snpId": "<string>",
    "pos": "<integer>",
    "snpIdUpper": "<string>",
    "variantId": "<string>",
    "geneSymbol": "<string>",
    "pValue": "<number>",
    "geneSymbolUpper": "<string>",
    "datasetId": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_single_tissue_eqtls-parameters)

## 47) List GTEx tissue sites and sample statistics {#dnakit.databases.gtex_get_tissue_sites}

- **Purpose:** List GTEx tissue sites and sample statistics.
- **Method:** Use `GTExV2` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtex_get_tissue_sites(*, operation=None, dataset_id=None, page=None, items_per_page=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `tissueSiteDetailId`, `colorHex`, `colorRgb`, `datasetId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtex_get_tissue_sites(
    operation='get_tissue_sites',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "tissueSiteDetailId": "<string>",
    "colorHex": "<string>",
    "colorRgb": "<string>",
    "datasetId": "<string>",
    "eGeneCount": "<integer>",
    "expressedGeneCount": "<integer>",
    "hasEGenes": "<boolean>",
    "hasSGenes": "<boolean>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_get_tissue_sites-parameters)

## 48) Query GTEx single-tissue eQTL data {#dnakit.databases.gtex_query_eqtl}

- **Purpose:** Query GTEx single-tissue eQTL data.
- **Method:** Use `GTExEQTL` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gtex_query_eqtl(*, gene_symbol=None, ensembl_gene_id=None, page=None, size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** data.singleTissueEqtl contains single-tissue eQTL records; metadata records the source and query.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gtex_query_eqtl(
    gene_symbol='TP53',
    size=3,
    verbose=False,
)
print(result)
```

- **Example output:** data.singleTissueEqtl contains single-tissue eQTL records; metadata records the source and query.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.gtex_query_eqtl-parameters)

## 49) Aggregate gene information across genomics datasets {#dnakit.databases.harmonizome_get_gene}

- **Purpose:** Aggregate gene information across genomics datasets.
- **Method:** Use `Harmonizome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.harmonizome_get_gene(*, gene_symbol, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** Returns a mapping with `symbol`, `name`, `ncbi_entrez_gene_id`, `ncbi_entrez_gene_url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.harmonizome_get_gene(
    gene_symbol='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "symbol": "<string>",
  "name": "<string>",
  "ncbi_entrez_gene_id": "<integer>",
  "ncbi_entrez_gene_url": "<string>",
  "description": "<string>",
  "synonyms": [
    "<string>"
  ],
  "proteins": [
    {
      "symbol": "<string>",
      "href": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.harmonizome_get_gene-parameters)

## 50) Retrieve associations between gene expression and cancer survival {#dnakit.databases.hpa_get_cancer_prognostics_by_gene}

- **Purpose:** Retrieve associations between gene expression and cancer survival.
- **Method:** Use `HPAGetCancerPrognostics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpa_get_cancer_prognostics_by_gene(*, ensembl_id=None, gene_name=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `ensembl_id`, `gene`, `gene_synonym`, `cancers_assessed_count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpa_get_cancer_prognostics_by_gene(
    ensembl_id='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "ensembl_id": "<string>",
  "gene": "<string>",
  "gene_synonym": [
    "<string>"
  ],
  "cancers_assessed_count": "<integer>",
  "prognostic_cancers_count": "<integer>",
  "non_prognostic_cancers_count": "<integer>",
  "prognostic_summary": [
    {
      "cancer_type": "<string>",
      "prognostic_type": "<string>",
      "p_value": "<string>",
      "is_prognostic": "<boolean>"
    }
  ],
  "non_prognostic_summary": [
    {
      "cancer_type": "<string>",
      "prognostic_type": "<string>",
      "p_value": "<string>",
      "is_prognostic": "<boolean>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_cancer_prognostics_by_gene-parameters)

## 51) Inspect the retired HPA protein interaction query {#dnakit.databases.hpa_get_protein_interactions_by_gene}

- **Purpose:** Inspect the retired HPA protein interaction query.
- **Method:** Use `HPAGetProteinInteractions` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpa_get_protein_interactions_by_gene(*, gene_name, verbose=True)`
- **Input:** Provide `gene_name`. This data-source operation is currently unavailable and returns an error; the code only illustrates the call syntax.
- **Output:** This operation currently reports an unavailable service rather than valid interaction results.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpa_get_protein_interactions_by_gene(
    gene_name='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** This operation currently reports an unavailable service rather than valid interaction results.

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_protein_interactions_by_gene-parameters)

## 52) Retrieve RNA expression by biological source {#dnakit.databases.hpa_get_rna_expression_by_source}

- **Purpose:** Retrieve RNA expression by biological source.
- **Method:** Use `HPAGetRnaExpressionBySource` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpa_get_rna_expression_by_source(*, gene_name, source_type, source_name, verbose=True)`
- **Input:** Provide `gene_name`, `source_type`, `source_name`.
- **Output:** Returns a mapping with `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpa_get_rna_expression_by_source(
    gene_name='TP53',
    source_type='tissue',
    source_name='liver',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_rna_expression_by_source-parameters)

## 53) Retrieve gene RNA expression in selected tissues {#dnakit.databases.hpa_get_rna_expression_in_specific_tissues}

- **Purpose:** Retrieve gene RNA expression in selected tissues.
- **Method:** Use `HPAGetRnaExpressionByTissue` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpa_get_rna_expression_in_specific_tissues(*, ensembl_id=None, tissue_names=None, tissues=None, gene_name=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpa_get_rna_expression_in_specific_tissues(
    ensembl_id='ENSG00000254647',
    tissue_names=['pancreas'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_rna_expression_in_specific_tissues-parameters)

## 54) Retrieve protein subcellular localization annotations {#dnakit.databases.hpa_get_subcellular_location}

- **Purpose:** Retrieve protein subcellular localization annotations.
- **Method:** Use `HPAGetSubcellularLocation` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpa_get_subcellular_location(*, gene_name, verbose=True)`
- **Input:** Provide `gene_name`.
- **Output:** Returns a mapping with `gene_name`, `gene_synonym`, `main_locations`, `additional_locations` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpa_get_subcellular_location(
    gene_name='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_name": "<string>",
  "gene_synonym": [
    "<string>"
  ],
  "main_locations": [
    "<string>"
  ],
  "additional_locations": [
    "<string>"
  ],
  "total_locations": "<integer>",
  "location_summary": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_get_subcellular_location-parameters)

## 55) Search HPA genes and Ensembl identifiers {#dnakit.databases.hpa_search_genes_by_query}

- **Purpose:** Search HPA genes and Ensembl identifiers.
- **Method:** Use `HPASearchGenes` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpa_search_genes_by_query(*, search_query, max_results=None, verbose=True)`
- **Input:** Provide `search_query`.
- **Output:** Returns a mapping with `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpa_search_genes_by_query(
    search_query='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.hpa_search_genes_by_query-parameters)

## 56) Retrieve transcription factor binding matrices and metadata {#dnakit.databases.jaspar_get_matrix}

- **Purpose:** Retrieve transcription factor binding matrices and metadata.
- **Method:** Use `JASPAR` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.jaspar_get_matrix(*, matrix_id, verbose=True)`
- **Input:** Provide `matrix_id`.
- **Output:** Returns a mapping with `matrix_id`, `base_id`, `name`, `collection` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.jaspar_get_matrix(
    matrix_id='MA0002.2',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "matrix_id": "<string>",
  "base_id": "<string>",
  "name": "<string>",
  "collection": "<string>",
  "version": "<string>",
  "pfm": {
    "A": [
      "<number>"
    ],
    "C": [
      "<number>"
    ],
    "G": [
      "<number>"
    ],
    "T": [
      "<number>"
    ]
  },
  "tax_group": "<string>",
  "species": [
    {
      "tax_id": "<integer>",
      "name": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_get_matrix-parameters)

## 57) List available versions of a JASPAR matrix {#dnakit.databases.jaspar_get_matrix_versions}

- **Purpose:** List available versions of a JASPAR matrix.
- **Method:** Use `JASPAR` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.jaspar_get_matrix_versions(*, base_id, page=None, page_size=None, verbose=True)`
- **Input:** Provide `base_id`.
- **Output:** Returns a mapping with `count`, `next`, `previous`, `results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.jaspar_get_matrix_versions(
    base_id='MA0002',
    page=1,
    page_size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "count": "<integer>",
  "next": "<string>",
  "previous": "<string>",
  "results": [
    {
      "matrix_id": "<string>",
      "base_id": "<string>",
      "name": "<string>",
      "collection": "<string>",
      "version": "<string>",
      "url": "<string>",
      "sequence_logo": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_get_matrix_versions-parameters)

## 58) Browse JASPAR transcription factor binding profiles {#dnakit.databases.jaspar_get_transcription_factors}

- **Purpose:** Browse JASPAR transcription factor binding profiles.
- **Method:** Use `JASPAR` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.jaspar_get_transcription_factors(*, collection=None, page=None, page_size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `count`, `next`, `previous`, `results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.jaspar_get_transcription_factors(
    collection='CORE',
    page=1,
    page_size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "count": "<integer>",
  "next": "<string>",
  "previous": "<string>",
  "results": [
    {
      "matrix_id": "<string>",
      "base_id": "<string>",
      "name": "<string>",
      "collection": "<string>",
      "version": "<string>",
      "url": "<string>",
      "sequence_logo": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_get_transcription_factors-parameters)

## 59) Search binding matrices by name, species or taxonomic group {#dnakit.databases.jaspar_search_matrices}

- **Purpose:** Search binding matrices by name, species or taxonomic group.
- **Method:** Use `JASPAR` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.jaspar_search_matrices(*, search=None, name=None, collection=None, tax_group=None, species=None, page=None, page_size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `count`, `next`, `previous`, `results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.jaspar_search_matrices(
    search='CTCF',
    page=1,
    page_size=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "count": "<integer>",
  "next": "<string>",
  "previous": "<string>",
  "results": [
    {
      "matrix_id": "<string>",
      "base_id": "<string>",
      "name": "<string>",
      "collection": "<string>",
      "version": "<string>",
      "url": "<string>",
      "sequence_logo": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.jaspar_search_matrices-parameters)

## 60) Retrieve transcription factor binding sites by gene and cell type {#dnakit.databases.re_map_get_transcription_factor_binding}

- **Purpose:** Retrieve transcription factor binding sites by gene and cell type.
- **Method:** Use `ReMap` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.re_map_get_transcription_factor_binding(*, gene_name, cell_type=None, limit=None, verbose=True)`
- **Input:** Provide `gene_name`.
- **Output:** Returns a mapping with `experiments`, `gene_name`, `cell_type`, `count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.re_map_get_transcription_factor_binding(
    gene_name='BRCA1',
    cell_type='HepG2',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "experiments": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "target": "<object>",
      "biosample_ontology": "<object>",
      "description": "<string>",
      "status": "<string>"
    }
  ],
  "gene_name": "<string>",
  "cell_type": "<string>",
  "count": "<integer>",
  "url": "<string>",
  "note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.re_map_get_transcription_factor_binding-parameters)

## 61) Retrieve regulatory variant annotations by rsID {#dnakit.databases.regulome_db_query_variant}

- **Purpose:** Retrieve regulatory variant annotations by rsID.
- **Method:** Use `RegulomeDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.regulome_db_query_variant(*, rsid=None, variant=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `rsid`, `assembly`, `query_coordinates`, `regulome_score` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.regulome_db_query_variant(
    rsid='rs4994',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "rsid": "<string>",
  "assembly": "<string>",
  "query_coordinates": [
    "<string>"
  ],
  "regulome_score": {
    "probability": "<string>",
    "ranking": "<string>",
    "tissue_specific_scores": "<object>"
  },
  "variants": [
    {
      "chrom": "<string>",
      "start": "<integer>",
      "end": "<integer>",
      "rsids": "<array>"
    }
  ],
  "features": "<object>",
  "nearby_snps": [
    "<value>"
  ],
  "notifications": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.regulome_db_query_variant-parameters)

## 62) Retrieve candidate cis-regulatory elements in an interval {#dnakit.databases.screen_get_regulatory_elements}

- **Purpose:** Retrieve candidate cis-regulatory elements in an interval.
- **Method:** Use `SCREEN` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.screen_get_regulatory_elements(*, gene_name=None, query=None, element_type=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `gene_name`, `element_type`, `count`, `regulatory_elements` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.screen_get_regulatory_elements(
    gene_name='BRCA1',
    element_type='enhancer',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_name": "<string>",
  "element_type": "<string>",
  "count": "<integer>",
  "regulatory_elements": [
    {
      "accession": "<string>",
      "assay_title": "<string>",
      "target": "<string>",
      "biosample": "<string>",
      "status": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.screen_get_regulatory_elements-parameters)

## 63) Retrieve ENCODE candidate cis-regulatory elements from UCSC {#dnakit.databases.ucsc_get_encode_c_cres}

- **Purpose:** Retrieve ENCODE candidate cis-regulatory elements from UCSC.
- **Method:** Use `UCSCEpigenomics` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ucsc_get_encode_c_cres(*, genome=None, chrom, start, end, verbose=True)`
- **Input:** Provide `chrom`, `start`, `end`.
- **Output:** Returns a mapping with `genome`, `region`, `ccre_count`, `ccres` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ucsc_get_encode_c_cres(
    genome='hg38',
    chrom='chr17',
    start=7668421,
    end=7687490,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "genome": "<string>",
  "region": "<string>",
  "ccre_count": "<integer>",
  "ccres": [
    {
      "name": "<string>",
      "chrom": "<string>",
      "chromStart": "<integer>",
      "chromEnd": "<integer>",
      "cCRE_class": "<string>",
      "DNase_maxZ": "<number>",
      "H3K4me3_maxZ": "<number>",
      "H3K27ac_maxZ": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.ucsc_get_encode_c_cres-parameters)

## 64) Retrieve a direct transcription factor binding dataset {#dnakit.databases.uni_bind_get_dataset}

- **Purpose:** Retrieve a direct transcription factor binding dataset.
- **Method:** Use `UniBind` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_bind_get_dataset(*, operation=None, dataset_id, verbose=True)`
- **Input:** Provide `dataset_id`.
- **Output:** Returns a mapping with `dataset_id`, `tf_name`, `cell_line`, `biological_condition` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_bind_get_dataset(
    operation='get_dataset',
    dataset_id='EXP030726.neural_stem_cells.SMAD3',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "dataset_id": "<string>",
  "tf_name": "<string>",
  "cell_line": [
    "<string>"
  ],
  "biological_condition": [
    "<string>"
  ],
  "identifier": [
    "<string>"
  ],
  "jaspar_id": [
    "<string>"
  ],
  "prediction_models": [
    "<string>"
  ],
  "tfbs_models": [
    {
      "prediction_model": "<string>",
      "jaspar_id": "<string>",
      "jaspar_version": "<string>",
      "total_tfbs": "<string>",
      "score_threshold": "<string>",
      "distance_threshold": "<string>",
      "adj_centrimo_pvalue": "<number>",
      "bed_url": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.uni_bind_get_dataset-parameters)

## 65) List transcription factors in UniBind {#dnakit.databases.uni_bind_list_tfs}

- **Purpose:** List transcription factors in UniBind.
- **Method:** Use `UniBind` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_bind_list_tfs(*, operation=None, search=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list for “List transcription factors in UniBind”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_bind_list_tfs(
    operation='list_tfs',
    search='SMAD',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<string>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.uni_bind_list_tfs-parameters)

## 66) Search binding datasets by transcription factor and cell context {#dnakit.databases.uni_bind_search_datasets}

- **Purpose:** Search binding datasets by transcription factor and cell context.
- **Method:** Use `UniBind` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_bind_search_datasets(*, operation=None, tf_name=None, species=None, cell_line=None, collection=None, order=None, page=None, page_size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `tf_name`, `total_peaks`, `dataset_url`, `dataset_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_bind_search_datasets(
    operation='search_datasets',
    tf_name='SMAD3',
    page_size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "tf_name": "<string>",
    "total_peaks": "<string>",
    "dataset_url": "<string>",
    "dataset_id": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_regulation.md#dnakit.databases.uni_bind_search_datasets-parameters)

<!-- END GENERATED API: regulation -->
