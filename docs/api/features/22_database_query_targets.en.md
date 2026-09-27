# Target, drug and association queries

Query targets, drugs and association evidence.

<!-- BEGIN GENERATED API: targets -->

<span id="targets"></span>

## 1) Retrieve drug mechanisms by name or identifier {#dnakit.databases.ch_embl_get_drug_mechanisms}

- **Purpose:** Retrieve drug mechanisms by name or identifier.
- **Method:** Use `ChEMBL` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_embl_get_drug_mechanisms(*, limit=None, offset=None, drug_chembl_id=None, drug_name=None, molecule_chembl_id=None, chembl_id=None, drug_chembl_id__exact=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_embl_get_drug_mechanisms(
    drug_chembl_id='CHEMBL1201581',
    limit=5,
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
    "mechanisms": [
      "<object>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_get_drug_mechanisms-parameters)

## 2) Retrieve ChEMBL molecular structures and properties {#dnakit.databases.ch_embl_get_molecule}

- **Purpose:** Retrieve ChEMBL molecular structures and properties.
- **Method:** Use `ChEMBL` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_embl_get_molecule(*, chembl_id=None, format=None, molecule_chembl_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_embl_get_molecule(
    chembl_id='CHEMBL25',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "data": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_get_molecule-parameters)

## 3) Retrieve bioactivity measurements for a target {#dnakit.databases.ch_embl_get_target_activities}

- **Purpose:** Retrieve bioactivity measurements for a target.
- **Method:** Use `ChEMBL` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_embl_get_target_activities(*, target_chembl_id__exact=None, limit=None, offset=None, target_chembl_id=None, chembl_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select exactly one argument group: `target_chembl_id__exact` / `target_chembl_id`.
- **Output:** Returns a mapping with `status`, `data`, `count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_embl_get_target_activities(
    target_chembl_id__exact='CHEMBL2074',
    limit=20,
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
    "activities": [
      "<object>"
    ]
  },
  "count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_get_target_activities-parameters)

## 4) Search drugs by name and development phase {#dnakit.databases.ch_embl_search_drugs}

- **Purpose:** Search drugs by name and development phase.
- **Method:** Use `ChEMBL` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_embl_search_drugs(*, query=None, max_phase=None, limit=None, offset=None, molecule_chembl_id=None, q=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data`, `count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_embl_search_drugs(
    max_phase=4,
    limit=5,
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
    "drugs": [
      "<object>"
    ]
  },
  "count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_search_drugs-parameters)

## 5) Search drug targets by name, organism or type {#dnakit.databases.ch_embl_search_targets}

- **Purpose:** Search drug targets by name, organism or type.
- **Method:** Use `ChEMBL` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ch_embl_search_targets(*, target_chembl_id=None, pref_name__contains=None, organism=None, target_type=None, fields=None, limit=None, offset=None, q=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data`, `count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ch_embl_search_targets(
    pref_name__contains='kinase',
    organism='Homo sapiens',
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
  "data": {
    "targets": [
      "<object>"
    ]
  },
  "count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.ch_embl_search_targets-parameters)

## 6) Inspect the unavailable CTD chemical–gene interaction query {#dnakit.databases.ctd_get_chemical_gene_interactions}

- **Purpose:** Inspect the unavailable CTD chemical–gene interaction query.
- **Method:** Use `CTD` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ctd_get_chemical_gene_interactions(*, input_terms, verbose=True)`
- **Input:** Provide `input_terms`. This data-source operation is currently unavailable and returns an error; the code only illustrates the call syntax.
- **Output:** This operation currently reports an unavailable service rather than valid interaction results.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ctd_get_chemical_gene_interactions(
    input_terms='bisphenol A',
    verbose=False,
)
print(result)
```

- **Example output:** This operation currently reports an unavailable service rather than valid interaction results.

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.ctd_get_chemical_gene_interactions-parameters)

## 7) Retrieve drug–gene interaction evidence in batches {#dnakit.databases.dgidb_get_drug_gene_interactions}

- **Purpose:** Retrieve drug–gene interaction evidence in batches.
- **Method:** Use `DGIdb` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dgidb_get_drug_gene_interactions(*, genes=None, gene_name=None, gene=None, interaction_sources=None, interaction_types=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dgidb_get_drug_gene_interactions(
    genes=['EGFR'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "genes": {
      "nodes": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.dgidb_get_drug_gene_interactions-parameters)

## 8) Retrieve gene druggability categories {#dnakit.databases.dgidb_get_gene_druggability}

- **Purpose:** Retrieve gene druggability categories.
- **Method:** Use `DGIdb` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dgidb_get_gene_druggability(*, genes=None, gene_name=None, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dgidb_get_gene_druggability(
    genes=['EGFR', 'TP53', 'BRCA1'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "genes": {
      "nodes": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.dgidb_get_gene_druggability-parameters)

## 9) Retrieve basic drug information and approval status {#dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id}

- **Purpose:** Retrieve basic drug information and approval status.
- **Method:** Use `XML` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `error`, `query`, `total_matches`, `total_returned_results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id(
    query='Aspirin',
    case_sensitive=False,
    exact_match=False,
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "error": "<string>",
  "query": "<string>",
  "total_matches": "<integer>",
  "total_returned_results": "<integer>",
  "results": [
    {
      "_tag": "<string>",
      "_text": "<string>",
      "_attributes": "<object>",
      "drug_name": "<string>",
      "drugbank_id": "<string>",
      "secondary_ids": "<string>",
      "description": "<string>",
      "cas_number": "<string>"
    }
  ],
  "search_parameters": {
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id-parameters)

## 10) Retrieve pharmacodynamics, mechanisms and pharmacokinetics {#dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id}

- **Purpose:** Retrieve pharmacodynamics, mechanisms and pharmacokinetics.
- **Method:** Use `XML` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `error`, `query`, `total_matches`, `total_returned_results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id(
    query='Aspirin',
    case_sensitive=False,
    exact_match=False,
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "error": "<string>",
  "query": "<string>",
  "total_matches": "<integer>",
  "total_returned_results": "<integer>",
  "results": [
    {
      "_tag": "<string>",
      "_text": "<string>",
      "_attributes": "<object>",
      "drug_name": "<string>",
      "drugbank_id": "<string>",
      "synonyms": "<string>",
      "pharmacodynamics": "<string>",
      "mechanism_of_action": "<string>"
    }
  ],
  "search_parameters": {
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id-parameters)

## 11) Retrieve drug targets, enzymes, carriers and transporters {#dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id}

- **Purpose:** Retrieve drug targets, enzymes, carriers and transporters.
- **Method:** Use `XML` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id(*, query, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `error`, `query`, `total_matches`, `total_returned_results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id(
    query='Aspirin',
    case_sensitive=False,
    exact_match=False,
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "error": "<string>",
  "query": "<string>",
  "total_matches": "<integer>",
  "total_returned_results": "<integer>",
  "results": [
    {
      "_tag": "<string>",
      "_text": "<string>",
      "_attributes": "<object>",
      "drug_name": "<string>",
      "drugbank_id": "<string>",
      "synonyms": "<string>",
      "targets": "<array>",
      "enzymes": "<array>"
    }
  ],
  "search_parameters": {
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id-parameters)

## 12) Search the DrugBank vocabulary by name or synonym {#dnakit.databases.drugbank_vocab_search}

- **Purpose:** Search the DrugBank vocabulary by name or synonym.
- **Method:** Use `Dataset` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.drugbank_vocab_search(*, query, search_fields=None, case_sensitive=None, exact_match=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `query`, `total_results`, `results`, `search_parameters` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.drugbank_vocab_search(
    query='aspirin',
    search_fields=['Common name', 'Synonyms'],
    case_sensitive=False,
    exact_match=False,
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "query": "<string>",
  "total_results": "<integer>",
  "results": [
    {
      "DrugBank ID": "<string>",
      "Accession Numbers": "<string>",
      "Common name": "<string>",
      "CAS": "<string>",
      "UNII": "<string>",
      "Synonyms": "<string>",
      "Standard InChI Key": "<string>",
      "matched_fields": "<array>"
    }
  ],
  "search_parameters": {
    "search_fields": [
      "<string>"
    ],
    "case_sensitive": "<boolean>",
    "exact_match": "<boolean>",
    "limit": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.drugbank_vocab_search-parameters)

## 13) Retrieve literature-based disease–target association scores {#dnakit.databases.europepmc_disease_target_score}

- **Purpose:** Retrieve literature-based disease–target association scores.
- **Method:** Use `DiseaseTargetScore` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.europepmc_disease_target_score(*, efoId, pageSize=None, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `disease_info`, `datasource`, `total_targets_with_scores`, `target_scores` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.europepmc_disease_target_score(
    efoId='MONDO_0011996',
    pageSize=100,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "disease_info": {
    "disease_id": "<string>",
    "disease_name": "<string>"
  },
  "datasource": "<string>",
  "total_targets_with_scores": "<integer>",
  "target_scores": [
    {
      "target_symbol": "<string>",
      "target_id": "<string>",
      "datasource": "<string>",
      "score": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.europepmc_disease_target_score-parameters)

## 14) Retrieve expression-based disease–target association scores {#dnakit.databases.expression_atlas_disease_target_score}

- **Purpose:** Retrieve expression-based disease–target association scores.
- **Method:** Use `DiseaseTargetScore` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.expression_atlas_disease_target_score(*, efoId, pageSize=None, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `disease_info`, `datasource`, `total_targets_with_scores`, `target_scores` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.expression_atlas_disease_target_score(
    efoId='MONDO_0011996',
    pageSize=100,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "disease_info": {
    "disease_id": "<string>",
    "disease_name": "<string>"
  },
  "datasource": "<string>",
  "total_targets_with_scores": "<integer>",
  "target_scores": [
    {
      "target_symbol": "<string>",
      "target_id": "<string>",
      "datasource": "<string>",
      "score": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.expression_atlas_disease_target_score-parameters)

## 15) Retrieve active ingredients in a drug product {#dnakit.databases.fda_get_active_ingredient_info_by_drug_name}

- **Purpose:** Retrieve active ingredients in a drug product.
- **Method:** Use `FDADrugLabel` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fda_get_active_ingredient_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **Input:** Provide `drug_name`.
- **Output:** Returns a mapping with `meta`, `results`, `result_count`, `duplicates_removed` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fda_get_active_ingredient_info_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "active_ingredient": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_active_ingredient_info_by_drug_name-parameters)

## 16) Retrieve adverse reactions from drug labels {#dnakit.databases.fda_get_adverse_reactions_by_drug_name}

- **Purpose:** Retrieve adverse reactions from drug labels.
- **Method:** Use `FDADrugLabel` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fda_get_adverse_reactions_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **Input:** Provide `drug_name`.
- **Output:** Returns a mapping with `meta`, `results`, `result_count`, `duplicates_removed` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fda_get_adverse_reactions_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "adverse_reactions": "<array>",
      "warnings_and_cautions": "<array>",
      "spl_product_data_elements": "<array>",
      "warnings": "<array>",
      "boxed_warning": "<array>",
      "precautions": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>",
  "section_note": "<string>",
  "note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_adverse_reactions_by_drug_name-parameters)

## 17) Retrieve drug boxed warnings and safety information {#dnakit.databases.fda_get_boxed_warning_info_by_drug_name}

- **Purpose:** Retrieve drug boxed warnings and safety information.
- **Method:** Use `FDADrugLabel` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fda_get_boxed_warning_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **Input:** Provide `drug_name`.
- **Output:** Returns a mapping with `meta`, `results`, `result_count`, `duplicates_removed` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fda_get_boxed_warning_info_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "boxed_warning": "<array>",
      "spl_product_data_elements": "<array>",
      "warnings_and_cautions": "<array>",
      "warnings": "<array>",
      "related_sections_present": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>",
  "section_note": "<string>",
  "note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_boxed_warning_info_by_drug_name-parameters)

## 18) Retrieve clinical study information from drug labels {#dnakit.databases.fda_get_clinical_studies_info_by_drug_name}

- **Purpose:** Retrieve clinical study information from drug labels.
- **Method:** Use `FDADrugLabel` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fda_get_clinical_studies_info_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **Input:** Provide `drug_name`.
- **Output:** Returns a mapping with `meta`, `results`, `result_count`, `duplicates_removed` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fda_get_clinical_studies_info_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "clinical_studies": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_clinical_studies_info_by_drug_name-parameters)

## 19) Retrieve drug indications and usage information {#dnakit.databases.fda_get_indications_by_drug_name}

- **Purpose:** Retrieve drug indications and usage information.
- **Method:** Use `FDADrugLabel` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fda_get_indications_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **Input:** Provide `drug_name`.
- **Output:** Returns a mapping with `meta`, `results`, `result_count`, `duplicates_removed` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fda_get_indications_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "indications_and_usage": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_indications_by_drug_name-parameters)

## 20) Retrieve mechanisms of action from drug labels {#dnakit.databases.fda_get_mechanism_of_action_by_drug_name}

- **Purpose:** Retrieve mechanisms of action from drug labels.
- **Method:** Use `FDADrugLabel` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.fda_get_mechanism_of_action_by_drug_name(*, drug_name, limit=None, skip=None, verbose=True)`
- **Input:** Provide `drug_name`.
- **Output:** Returns a mapping with `meta`, `results`, `result_count`, `duplicates_removed` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.fda_get_mechanism_of_action_by_drug_name(
    drug_name='Aspirin',
    limit=1,
    skip=0,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "meta": {
    "skip": "<integer>",
    "limit": "<integer>",
    "total": "<integer>"
  },
  "results": [
    {
      "openfda.brand_name": "<array>",
      "openfda.generic_name": "<array>",
      "mechanism_of_action": "<array>",
      "spl_product_data_elements": "<array>"
    }
  ],
  "result_count": "<integer>",
  "duplicates_removed": "<integer>",
  "dedup_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.fda_get_mechanism_of_action_by_drug_name-parameters)

## 21) Retrieve metabolite–disease associations {#dnakit.databases.hmdb_get_diseases}

- **Purpose:** Retrieve metabolite–disease associations.
- **Method:** Use `Metabolite` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hmdb_get_diseases(*, operation=None, hmdb_id=None, compound_name=None, pubchem_cid=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `identifier`, `compound_name`, `pubchem_cid`, `ctd_query_term` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hmdb_get_diseases(
    hmdb_id='HMDB0000122',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "identifier": "<string>",
  "compound_name": "<string>",
  "pubchem_cid": "<integer>",
  "ctd_query_term": "<string>",
  "disease_count": "<integer>",
  "diseases": [
    {
      "disease_name": "<string>",
      "disease_id": "<string>",
      "disease_categories": "<string>",
      "direct_evidence": "<string>",
      "pubmed_ids": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.hmdb_get_diseases-parameters)

## 22) Retrieve metabolite information by identifier or name {#dnakit.databases.hmdb_get_metabolite}

- **Purpose:** Retrieve metabolite information by identifier or name.
- **Method:** Use `Metabolite` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hmdb_get_metabolite(*, operation=None, hmdb_id=None, compound_name=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `pubchem_cid`, `name`, `iupac_name`, `formula` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hmdb_get_metabolite(
    hmdb_id='HMDB0000122',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "pubchem_cid": "<integer>",
  "name": "<string>",
  "iupac_name": "<string>",
  "formula": "<string>",
  "molecular_weight": "<string>",
  "smiles": "<string>",
  "inchikey": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.hmdb_get_metabolite-parameters)

## 23) Search metabolites by name or molecular formula {#dnakit.databases.hmdb_search}

- **Purpose:** Search metabolites by name or molecular formula.
- **Method:** Use `Metabolite` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hmdb_search(*, operation=None, query, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `query`, `results`, `count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hmdb_search(
    query='glucose',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "query": "<string>",
  "results": [
    {
      "pubchem_cid": "<integer>",
      "name": "<string>",
      "formula": "<string>",
      "molecular_weight": "<string>",
      "smiles": "<string>"
    }
  ],
  "count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.hmdb_search-parameters)

## 24) Retrieve known drugs associated with a disease {#dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id}

- **Purpose:** Retrieve known drugs associated with a disease.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id(*, efoId, size=None, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id(
    efoId='MONDO_0004975',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "drugAndClinicalCandidates": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id-parameters)

## 25) Retrieve target-associated drugs and development phases {#dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id}

- **Purpose:** Retrieve target-associated drugs and development phases.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id(*, ensemblId, size=None, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id(
    ensemblId='ENSG00000121879',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "drugAndClinicalCandidates": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id-parameters)

## 26) Retrieve HPO phenotypes associated with a disease {#dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id}

- **Purpose:** Retrieve HPO phenotypes associated with a disease.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id(*, efoId, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "phenotypes": "<object>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id-parameters)

## 27) Retrieve disease-associated targets with pagination {#dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id}

- **Purpose:** Retrieve disease-associated targets with pagination.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id(*, efoId, size=None, index=None, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "associatedTargets": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id-parameters)

## 28) Retrieve GWAS credible set and locus-to-gene details {#dnakit.databases.open_targets_get_credible_set_detail}

- **Purpose:** Retrieve GWAS credible set and locus-to-gene details.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_credible_set_detail(*, studyLocusId, verbose=True)`
- **Input:** Provide `studyLocusId`.
- **Output:** Returns a mapping with `credibleSet` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_credible_set_detail(
    studyLocusId='b758d8fb10924f5338cbad8d27c7dee8',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "credibleSet": {
    "studyLocusId": "<string>",
    "studyId": "<string>",
    "region": "<string>",
    "chromosome": "<string>",
    "position": "<integer>",
    "finemappingMethod": "<string>",
    "pValueMantissa": "<number>",
    "pValueExponent": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_credible_set_detail-parameters)

## 29) Retrieve disease ontology parents and ancestors {#dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id}

- **Purpose:** Retrieve disease ontology parents and ancestors.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id(*, efoId, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "ancestors": "<array>",
      "parents": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id-parameters)

## 30) Retrieve disease ontology children and descendants {#dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id}

- **Purpose:** Retrieve disease ontology children and descendants.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id(*, efoId, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "descendants": "<array>",
      "children": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id-parameters)

## 31) Retrieve disease descriptions and cross-references {#dnakit.databases.open_targets_get_disease_description_by_efo_id}

- **Purpose:** Retrieve disease descriptions and cross-references.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_disease_description_by_efo_id(*, efoId, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_description_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "description": "<string>",
      "dbXRefs": "<array>",
      "isTherapeuticArea": "<boolean>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_description_by_efo_id-parameters)

## 32) Resolve disease names to ontology identifiers {#dnakit.databases.open_targets_get_disease_id_description_by_name}

- **Purpose:** Resolve disease names to ontology identifiers.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_disease_id_description_by_name(*, diseaseName=None, query=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_id_description_by_name(
    diseaseName='Diabetes Mellitus',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "search": {
      "hits": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_id_description_by_name-parameters)

## 33) Retrieve disease synonyms {#dnakit.databases.open_targets_get_disease_synonyms_by_efo_id}

- **Purpose:** Retrieve disease synonyms.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_disease_synonyms_by_efo_id(*, efoId, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_synonyms_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "synonyms": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_synonyms_by_efo_id-parameters)

## 34) Retrieve therapeutic areas associated with a disease {#dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id}

- **Purpose:** Retrieve therapeutic areas associated with a disease.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id(*, efoId, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id(
    efoId='MONDO_0005011',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "therapeuticAreas": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id-parameters)

## 35) Retrieve target-associated diseases and phenotypes with pagination {#dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl}

- **Purpose:** Retrieve target-associated diseases and phenotypes with pagination.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl(*, ensemblId, page=None, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "associatedDiseases": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl-parameters)

## 36) Retrieve significant drug adverse event records {#dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id}

- **Purpose:** Retrieve significant drug adverse event records.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id(*, chemblId, page=None, verbose=True)`
- **Input:** Provide `chemblId`.
- **Output:** Returns a mapping with `drug`, `metadata`, `error`, `error_details` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id(
    chemblId='CHEMBL25',
    page={'index': 0, 'size': 10},
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "drug": {
    "id": "<string>",
    "name": "<string>",
    "adverseEvents": {
      "count": "<integer>",
      "returned": "<integer>",
      "truncated": "<boolean>",
      "criticalValue": "<number>",
      "rows": "<array>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  },
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id-parameters)

## 37) Resolve a generic drug name to a ChEMBL identifier {#dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name}

- **Purpose:** Resolve a generic drug name to a ChEMBL identifier.
- **Method:** Use `OpentargetToolDrugNameMatch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name(*, drugName=None, query=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name(
    drugName='Aspirin',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "search": {
      "hits": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name-parameters)

## 38) Retrieve drug indications by ChEMBL identifier {#dnakit.databases.open_targets_get_drug_indications_by_chembl_id}

- **Purpose:** Retrieve drug indications by ChEMBL identifier.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_drug_indications_by_chembl_id(*, chemblId, verbose=True)`
- **Input:** Provide `chemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_indications_by_chembl_id(
    chemblId='CHEMBL25',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "drug": {
      "id": "<string>",
      "name": "<string>",
      "indications": "<object>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_indications_by_chembl_id-parameters)

## 39) Retrieve drug mechanisms by ChEMBL identifier {#dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id}

- **Purpose:** Retrieve drug mechanisms by ChEMBL identifier.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id(*, chemblId, verbose=True)`
- **Input:** Provide `chemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id(
    chemblId='CHEMBL25',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "drug": {
      "id": "<string>",
      "name": "<string>",
      "mechanismsOfAction": "<object>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id-parameters)

## 40) Retrieve drug warnings by ChEMBL identifier {#dnakit.databases.open_targets_get_drug_warnings_by_chembl_id}

- **Purpose:** Retrieve drug warnings by ChEMBL identifier.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_drug_warnings_by_chembl_id(*, chemblId, verbose=True)`
- **Input:** Provide `chemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_drug_warnings_by_chembl_id(
    chemblId='CHEMBL25',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "drug": {
      "id": "<string>",
      "name": "<string>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_drug_warnings_by_chembl_id-parameters)

## 41) Filter target–disease evidence by data source {#dnakit.databases.open_targets_get_evidence_by_datasource}

- **Purpose:** Filter target–disease evidence by data source.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_evidence_by_datasource(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, datasourceIds=None, size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_evidence_by_datasource(
    efoId='MONDO_0005011',
    ensemblId='ENSG00000141510',
    datasourceIds=['europepmc'],
    size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "evidences": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>",
    "datasource_rename_note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_evidence_by_datasource-parameters)

## 42) Retrieve GWAS cohort and study metadata {#dnakit.databases.open_targets_get_gwas_study}

- **Purpose:** Retrieve GWAS cohort and study metadata.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_gwas_study(*, studyId, verbose=True)`
- **Input:** Provide `studyId`.
- **Output:** Returns a mapping with `study` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_gwas_study(
    studyId='GCST000392',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "study": {
    "id": "<string>",
    "studyType": "<string>",
    "traitFromSource": "<string>",
    "publicationFirstAuthor": "<string>",
    "publicationTitle": "<string>",
    "publicationDate": "<string>",
    "publicationJournal": "<string>",
    "pubmedId": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_gwas_study-parameters)

## 43) Retrieve publications associated with a target {#dnakit.databases.open_targets_get_publications_by_target_ensembl_id}

- **Purpose:** Retrieve publications associated with a target.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_publications_by_target_ensembl_id(*, entityId, additionalIds=None, startYear=None, startMonth=None, endYear=None, endMonth=None, verbose=True)`
- **Input:** Provide `entityId`.
- **Output:** Returns a mapping with `data`, `metadata`, `error`, `error_details` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_publications_by_target_ensembl_id(
    entityId='ENSG00000141510',
    startYear=2023,
    startMonth=0,
    endYear=2023,
    endMonth=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "literatureOcurrences": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  },
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_publications_by_target_ensembl_id-parameters)

## 44) Retrieve biomedical entities similar to a disease {#dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id}

- **Purpose:** Retrieve biomedical entities similar to a disease.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id(*, efoId, threshold, size, verbose=True)`
- **Input:** Provide `efoId`, `threshold`, `size`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id(
    efoId='MONDO_0005011',
    threshold=1.0,
    size=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "similarEntities": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id-parameters)

## 45) Retrieve credible sets from a GWAS study {#dnakit.databases.open_targets_get_study_credible_sets}

- **Purpose:** Retrieve credible sets from a GWAS study.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_study_credible_sets(*, studyIds, size=None, index=None, verbose=True)`
- **Input:** Provide `studyIds`.
- **Output:** Returns a mapping with `credibleSets` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_study_credible_sets(
    studyIds=['GCST000392'],
    size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "credibleSets": {
    "count": "<integer>",
    "rows": [
      "<object>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_study_credible_sets-parameters)

## 46) Retrieve target functional classes {#dnakit.databases.open_targets_get_target_classes_by_ensembl_id}

- **Purpose:** Retrieve target functional classes.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_classes_by_ensembl_id(*, ensemblId, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_classes_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "targetClass": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_classes_by_ensembl_id-parameters)

## 47) Retrieve GO annotations for a target {#dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id}

- **Purpose:** Retrieve GO annotations for a target.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id(*, ensemblId, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "geneOntology": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id-parameters)

## 48) Retrieve target genomic coordinates and strand {#dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id}

- **Purpose:** Retrieve target genomic coordinates and strand.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id(*, ensemblId, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "genomicLocation": "<object>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id-parameters)

## 49) Retrieve homologues of a target {#dnakit.databases.open_targets_get_target_homologues_by_ensembl_id}

- **Purpose:** Retrieve homologues of a target.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_homologues_by_ensembl_id(*, ensemblId, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_homologues_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "homologues": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_homologues_by_ensembl_id-parameters)

## 50) Resolve target names to identifiers and descriptions {#dnakit.databases.open_targets_get_target_id_description_by_name}

- **Purpose:** Resolve target names to identifiers and descriptions.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_id_description_by_name(*, targetName, verbose=True)`
- **Input:** Provide `targetName`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_id_description_by_name(
    targetName='BRAF',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "search": {
      "hits": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_id_description_by_name-parameters)

## 51) Retrieve target interaction partners and evidence with pagination {#dnakit.databases.open_targets_get_target_interactions_by_ensembl_id}

- **Purpose:** Retrieve target interaction partners and evidence with pagination.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_interactions_by_ensembl_id(*, ensemblId, page=None, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data`, `metadata`, `error`, `error_details` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_interactions_by_ensembl_id(
    ensemblId='ENSG00000141510',
    page={'index': 0, 'size': 10},
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "interactions": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  },
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_interactions_by_ensembl_id-parameters)

## 52) Retrieve known target safety liabilities {#dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id}

- **Purpose:** Retrieve known target safety liabilities.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id(*, ensemblId, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "safetyLiabilities": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id-parameters)

## 53) Retrieve target tractability across treatment modalities {#dnakit.databases.open_targets_get_target_tractability_by_ensembl_id}

- **Purpose:** Retrieve target tractability across treatment modalities.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_target_tractability_by_ensembl_id(*, ensemblId, verbose=True)`
- **Input:** Provide `ensemblId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_target_tractability_by_ensembl_id(
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "target": {
      "id": "<string>",
      "approvedSymbol": "<string>",
      "tractability": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_target_tractability_by_ensembl_id-parameters)

## 54) Retrieve GWAS credible sets containing a variant {#dnakit.databases.open_targets_get_variant_credible_sets}

- **Purpose:** Retrieve GWAS credible sets containing a variant.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_variant_credible_sets(*, variantId=None, size=None, variant_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `variant` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_variant_credible_sets(
    variantId='10_112998590_C_T',
    size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": {
    "id": "<string>",
    "rsIds": [
      "<string>"
    ],
    "credibleSets": {
      "count": "<integer>",
      "rows": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_variant_credible_sets-parameters)

## 55) Retrieve variant coordinates, frequencies and consequences {#dnakit.databases.open_targets_get_variant_info}

- **Purpose:** Retrieve variant coordinates, frequencies and consequences.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_get_variant_info(*, variantId, verbose=True)`
- **Input:** Provide `variantId`.
- **Output:** Returns a mapping with `variant` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_get_variant_info(
    variantId='10_112998590_C_T',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": {
    "id": "<string>",
    "rsIds": [
      "<string>"
    ],
    "chromosome": "<string>",
    "position": "<integer>",
    "referenceAllele": "<string>",
    "alternateAllele": "<string>",
    "hgvsId": "<string>",
    "variantDescription": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_get_variant_info-parameters)

## 56) Map disease and phenotype identifiers across databases {#dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids}

- **Purpose:** Map disease and phenotype identifiers across databases.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids(*, inputId, verbose=True)`
- **Input:** Provide `inputId`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids(
    inputId='1',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "mapIds": {
      "mappings": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids-parameters)

## 57) Search disease, target and drug entities together {#dnakit.databases.open_targets_multi_entity_search_by_query_string}

- **Purpose:** Search disease, target and drug entities together.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_multi_entity_search_by_query_string(*, queryString, entityNames=None, page=None, verbose=True)`
- **Input:** Provide `queryString`.
- **Output:** Returns a mapping with `error`, `error_details` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_multi_entity_search_by_query_string(
    queryString='gene:TP53',
    entityNames=[],
    page={'index': 0, 'size': 10},
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "error": "<string>",
  "error_details": {
    "type": "<string>",
    "message": "<string>",
    "retriable": "<boolean>",
    "next_steps": [
      "<string>"
    ],
    "details": {
      "validation_error": "<string>",
      "path": "<array>",
      "schema": "<object>",
      "parameter": "<string>",
      "expected": "<string>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_multi_entity_search_by_query_string-parameters)

## 58) Find GWAS studies by disease or trait {#dnakit.databases.open_targets_search_gwas_studies_by_disease}

- **Purpose:** Find GWAS studies by disease or trait.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_search_gwas_studies_by_disease(*, diseaseIds=None, enableIndirect=None, disease_name=None, disease=None, trait=None, size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `studies` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_search_gwas_studies_by_disease(
    diseaseIds=['MONDO_0005148'],
    size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "studies": {
    "count": "<integer>",
    "rows": [
      "<object>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_search_gwas_studies_by_disease-parameters)

## 59) Retrieve IntOGen cancer driver association evidence {#dnakit.databases.open_targets_target_disease_evidence}

- **Purpose:** Retrieve IntOGen cancer driver association evidence.
- **Method:** Use `OpenTarget` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.open_targets_target_disease_evidence(*, efoId=None, ensemblId=None, gene_symbol=None, disease_name=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.open_targets_target_disease_evidence(
    efoId='MONDO_0005011',
    ensemblId='ENSG00000141510',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "disease": {
      "id": "<string>",
      "name": "<string>",
      "evidences": "<object>"
    }
  },
  "metadata": {
    "truncation_note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.open_targets_target_disease_evidence-parameters)

## 60) Resolve compound names to PubChem CIDs {#dnakit.databases.pub_chem_get_cid_by_compound_name}

- **Purpose:** Resolve compound names to PubChem CIDs.
- **Method:** Use `PubChem` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pub_chem_get_cid_by_compound_name(*, name=None, compound_name=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `IdentifierList` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pub_chem_get_cid_by_compound_name(
    name='aspirin',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "IdentifierList": {
    "CID": [
      "<integer>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.pub_chem_get_cid_by_compound_name-parameters)

## 61) Retrieve molecular properties by PubChem CID {#dnakit.databases.pub_chem_get_compound_properties_by_cid}

- **Purpose:** Retrieve molecular properties by PubChem CID.
- **Method:** Use `PubChem` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pub_chem_get_compound_properties_by_cid(*, cid, properties=None, verbose=True)`
- **Input:** Provide `cid`.
- **Output:** Returns a mapping with `PropertyTable` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pub_chem_get_compound_properties_by_cid(
    cid=2244,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "PropertyTable": {
    "Properties": [
      "<object>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.pub_chem_get_compound_properties_by_cid-parameters)

## 62) Retrieve pathway-based disease–target association scores {#dnakit.databases.reactome_disease_target_score}

- **Purpose:** Retrieve pathway-based disease–target association scores.
- **Method:** Use `DiseaseTargetScore` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.reactome_disease_target_score(*, efoId, pageSize=None, verbose=True)`
- **Input:** Provide `efoId`.
- **Output:** Returns a mapping with `disease_info`, `datasource`, `total_targets_with_scores`, `target_scores` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.reactome_disease_target_score(
    efoId='MONDO_0011996',
    pageSize=100,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "disease_info": {
    "disease_id": "<string>",
    "disease_name": "<string>"
  },
  "datasource": "<string>",
  "total_targets_with_scores": "<integer>",
  "target_scores": [
    {
      "target_symbol": "<string>",
      "target_id": "<string>",
      "datasource": "<string>",
      "score": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_targets.md#dnakit.databases.reactome_disease_target_score-parameters)

<!-- END GENERATED API: targets -->
