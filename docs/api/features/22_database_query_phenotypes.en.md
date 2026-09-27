# Disease, phenotype and ontology queries

Query diseases, phenotypes and associated ontology terms.

<!-- BEGIN GENERATED API: phenotypes -->

<span id="phenotypes"></span>

## 1) Retrieve disease-associated genes and association scores {#dnakit.databases.dis_ge_net_get_disease_genes}

- **Purpose:** Retrieve disease-associated genes and association scores.
- **Method:** Use `DisGeNET` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dis_ge_net_get_disease_genes(*, operation=None, disease, min_score=None, limit=None, verbose=True)`
- **Input:** Provide `disease`. Configure: `DISGENET_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dis_ge_net_get_disease_genes(
    disease='C0006142',
    min_score=0.3,
    limit=30,
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
    "disease": "<string>",
    "genes": [
      "<object>"
    ],
    "gene_count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_get_disease_genes-parameters)

## 2) Filter gene–disease associations by source and score {#dnakit.databases.dis_ge_net_get_gda}

- **Purpose:** Filter gene–disease associations by source and score.
- **Method:** Use `DisGeNET` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dis_ge_net_get_gda(*, operation=None, gene=None, disease=None, source=None, min_score=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Configure: `DISGENET_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dis_ge_net_get_gda(
    gene='TP53',
    min_score=0.3,
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
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_get_gda-parameters)

## 3) Retrieve variant–disease association evidence {#dnakit.databases.dis_ge_net_get_vda}

- **Purpose:** Retrieve variant–disease association evidence.
- **Method:** Use `DisGeNET` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dis_ge_net_get_vda(*, operation=None, variant=None, gene=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Configure: `DISGENET_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dis_ge_net_get_vda(
    gene='BRCA1',
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
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_get_vda-parameters)

## 4) Find disease-associated genes by UMLS concept ID {#dnakit.databases.dis_ge_net_search_disease}

- **Purpose:** Find disease-associated genes by UMLS concept ID.
- **Method:** Use `DisGeNET` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dis_ge_net_search_disease(*, operation=None, disease, limit=None, verbose=True)`
- **Input:** Provide `disease`. Configure: `DISGENET_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dis_ge_net_search_disease(
    disease='C0006142',
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
    "disease": "<string>",
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_search_disease-parameters)

## 5) Find diseases associated with a gene symbol {#dnakit.databases.dis_ge_net_search_gene}

- **Purpose:** Find diseases associated with a gene symbol.
- **Method:** Use `DisGeNET` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dis_ge_net_search_gene(*, operation=None, gene, limit=None, verbose=True)`
- **Input:** Provide `gene`. Configure: `DISGENET_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dis_ge_net_search_gene(
    gene='BRCA1',
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
    "gene": "<string>",
    "associations": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.dis_ge_net_search_gene-parameters)

## 6) Summarize GenCC gene–disease validity classifications {#dnakit.databases.gen_cc_get_classifications}

- **Purpose:** Summarize GenCC gene–disease validity classifications.
- **Method:** Use `GenCC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gen_cc_get_classifications(*, operation=None, gene_symbol=None, disease=None, submitter=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gen_cc_get_classifications(
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
    "classifications": [
      "<object>"
    ],
    "total_submissions": "<integer>",
    "unique_genes": "<integer>",
    "unique_diseases": "<integer>",
    "top_submitters": [
      "<object>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.gen_cc_get_classifications-parameters)

## 7) Search GenCC gene evidence for a disease {#dnakit.databases.gen_cc_search_disease}

- **Purpose:** Search GenCC gene evidence for a disease.
- **Method:** Use `GenCC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gen_cc_search_disease(*, operation=None, disease=None, disease_name=None, classification=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gen_cc_search_disease(
    disease='Marfan syndrome',
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
    "disease": "<string>",
    "submissions": [
      "<object>"
    ],
    "submission_count": "<integer>",
    "unique_genes": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.gen_cc_search_disease-parameters)

## 8) Retrieve GenCC gene–disease classifications {#dnakit.databases.gen_cc_search_gene}

- **Purpose:** Retrieve GenCC gene–disease classifications.
- **Method:** Use `GenCC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gen_cc_search_gene(*, operation=None, gene_symbol, classification=None, verbose=True)`
- **Input:** Provide `gene_symbol`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gen_cc_search_gene(
    gene_symbol='BRCA2',
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
    "submissions": [
      "<object>"
    ],
    "submission_count": "<integer>",
    "unique_diseases": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.gen_cc_search_gene-parameters)

## 9) Retrieve diseases associated with multiple HPO phenotypes {#dnakit.databases.get_joint_associated_diseases_by_hpo_id_list}

- **Purpose:** Retrieve diseases associated with multiple HPO phenotypes.
- **Method:** Use `MonarchDiseasesForMultiplePheno` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.get_joint_associated_diseases_by_hpo_id_list(*, HPO_ID_list, limit=None, offset=None, verbose=True)`
- **Input:** Provide `HPO_ID_list`.
- **Output:** Returns a list for “Retrieve diseases associated with multiple HPO phenotypes”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.get_joint_associated_diseases_by_hpo_id_list(
    HPO_ID_list=['HP:0001250'],
    limit=5,
    offset=1,
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

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.get_joint_associated_diseases_by_hpo_id_list-parameters)

## 10) Retrieve diseases associated with an HPO phenotype {#dnakit.databases.hpo_get_diseases_by_phenotype}

- **Purpose:** Retrieve diseases associated with an HPO phenotype.
- **Method:** Use `HPO` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpo_get_diseases_by_phenotype(*, term_id=None, limit=None, offset=None, hp_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `diseases` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpo_get_diseases_by_phenotype(
    term_id='HP:0001250',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "diseases": [
    {
      "id": "<string>",
      "name": "<string>",
      "mondo_id": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_get_diseases_by_phenotype-parameters)

## 11) Retrieve genes associated with an HPO phenotype {#dnakit.databases.hpo_get_genes_by_phenotype}

- **Purpose:** Retrieve genes associated with an HPO phenotype.
- **Method:** Use `HPO` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpo_get_genes_by_phenotype(*, term_id, limit=None, offset=None, verbose=True)`
- **Input:** Provide `term_id`.
- **Output:** Returns a mapping with `genes` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpo_get_genes_by_phenotype(
    term_id='HP:0001250',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "genes": [
    {
      "id": "<string>",
      "name": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_get_genes_by_phenotype-parameters)

## 12) Retrieve HPO term definitions and synonyms {#dnakit.databases.hpo_get_term}

- **Purpose:** Retrieve HPO term definitions and synonyms.
- **Method:** Use `HPO` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpo_get_term(*, term_id, verbose=True)`
- **Input:** Provide `term_id`.
- **Output:** Returns a mapping with `id`, `name`, `definition`, `comment` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpo_get_term(
    term_id='HP:0001250',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "id": "<string>",
  "name": "<string>",
  "definition": "<string>",
  "comment": "<string>",
  "synonyms": [
    "<string>"
  ],
  "descendant_count": "<integer>",
  "xrefs": [
    "<string>"
  ],
  "translations": [
    "<object>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_get_term-parameters)

## 13) Search HPO phenotype terms by keyword {#dnakit.databases.hpo_search_terms}

- **Purpose:** Search HPO phenotype terms by keyword.
- **Method:** Use `HPO` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.hpo_search_terms(*, query, max_results=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `id`, `name`, `definition`, `descendant_count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.hpo_search_terms(
    query='seizure',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "id": "<string>",
    "name": "<string>",
    "definition": "<string>",
    "descendant_count": "<integer>",
    "synonyms": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.hpo_search_terms-parameters)

## 14) Aggregate human gene identifiers and annotations {#dnakit.databases.marrvel_get_gene}

- **Purpose:** Aggregate human gene identifiers and annotations.
- **Method:** Use `MARRVELGene` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.marrvel_get_gene(*, symbol=None, gene_symbol=None, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select at least one argument group: `symbol` / `gene_symbol` / `gene`.
- **Output:** Returns gene identifier cross-references, names, location, aliases and summary; no OMIM identifier is included.
- **Example code:**

```python
import dnakit

result = dnakit.databases.marrvel_get_gene(
    symbol='CFTR',
    verbose=False,
)
print(result)
```

- **Example output:** Returns gene identifier cross-references, names, location, aliases and summary; no OMIM identifier is included.

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.marrvel_get_gene-parameters)

## 15) Retrieve OMIM phenotypes and inheritance for a gene {#dnakit.databases.marrvel_get_omim_phenotypes}

- **Purpose:** Retrieve OMIM phenotypes and inheritance for a gene.
- **Method:** Use `MARRVELOmim` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.marrvel_get_omim_phenotypes(*, symbol=None, gene_symbol=None, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select at least one argument group: `symbol` / `gene_symbol` / `gene`.
- **Output:** Returns a list for “Retrieve OMIM phenotypes and inheritance for a gene”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.marrvel_get_omim_phenotypes(
    symbol='CFTR',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<object>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.marrvel_get_omim_phenotypes-parameters)

## 16) Retrieve Monarch gene–disease associations {#dnakit.databases.monarch_get_gene_diseases}

- **Purpose:** Retrieve Monarch gene–disease associations.
- **Method:** Use `Monarch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_get_gene_diseases(*, subject, limit=None, offset=None, verbose=True)`
- **Input:** Provide `subject`.
- **Output:** Returns a mapping with `items`, `total` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_get_gene_diseases(
    subject='HGNC:11998',
    limit=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "items": [
    "<value>"
  ],
  "total": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_get_gene_diseases-parameters)

## 17) Retrieve cross-species gene–phenotype associations {#dnakit.databases.monarch_get_gene_phenotypes}

- **Purpose:** Retrieve cross-species gene–phenotype associations.
- **Method:** Use `Monarch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_get_gene_phenotypes(*, subject, limit=None, offset=None, verbose=True)`
- **Input:** Provide `subject`.
- **Output:** Returns a mapping with `items`, `total` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_get_gene_phenotypes(
    subject='HGNC:11998',
    limit=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "items": [
    "<value>"
  ],
  "total": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_get_gene_phenotypes-parameters)

## 18) Search Monarch biomedical entities by name {#dnakit.databases.monarch_search}

- **Purpose:** Search Monarch biomedical entities by name.
- **Method:** Use `MonarchV3` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_search(*, query, category=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `id`, `name`, `category`, `symbol` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_search(
    query='Alzheimer',
    category='biolink:Disease',
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
    "id": "<string>",
    "name": "<string>",
    "category": "<string>",
    "symbol": "<string>",
    "description": "<string>",
    "taxon": "<string>",
    "taxon_label": "<string>",
    "deprecated": "<boolean>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_search-parameters)

## 19) Search Monarch genes and cross-species identifiers {#dnakit.databases.monarch_search_gene}

- **Purpose:** Search Monarch genes and cross-species identifiers.
- **Method:** Use `Monarch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_search_gene(*, query, limit=None, offset=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `items`, `total` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_search_gene(
    query='TP53',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "items": [
    "<value>"
  ],
  "total": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_search_gene-parameters)

## 20) Query cross-species biomedical associations {#dnakit.databases.monarch_v3_get_associations}

- **Purpose:** Query cross-species biomedical associations.
- **Method:** Use `MonarchV3` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_v3_get_associations(*, subject=None, object=None, category, limit=None, verbose=True)`
- **Input:** Provide `category`.
- **Output:** Returns a list for “Query cross-species biomedical associations”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_associations(
    subject='HGNC:11998',
    category='biolink:GeneToPhenotypicFeatureAssociation',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_associations-parameters)

## 21) Retrieve biomedical entity details by CURIE {#dnakit.databases.monarch_v3_get_entity}

- **Purpose:** Retrieve biomedical entity details by CURIE.
- **Method:** Use `MonarchV3` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_v3_get_entity(*, entity_id, verbose=True)`
- **Input:** Provide `entity_id`.
- **Output:** data contains the entity name, category, description, synonyms and cross-references.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_entity(
    entity_id='HGNC:11998',
    verbose=False,
)
print(result)
```

- **Example output:** data contains the entity name, category, description, synonyms and cross-references.

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_entity-parameters)

## 22) Summarize disease phenotypes by body system {#dnakit.databases.monarch_v3_get_histopheno}

- **Purpose:** Summarize disease phenotypes by body system.
- **Method:** Use `MonarchV3` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_v3_get_histopheno(*, entity_id, verbose=True)`
- **Input:** Provide `entity_id`.
- **Output:** Returns a list for “Summarize disease phenotypes by body system”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_histopheno(
    entity_id='MONDO:0007947',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_histopheno-parameters)

## 23) Retrieve cross-ontology identifier mappings {#dnakit.databases.monarch_v3_get_mappings}

- **Purpose:** Retrieve cross-ontology identifier mappings.
- **Method:** Use `MonarchV3` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_v3_get_mappings(*, entity_id, limit=None, verbose=True)`
- **Input:** Provide `entity_id`.
- **Output:** Returns a list for “Retrieve cross-ontology identifier mappings”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_v3_get_mappings(
    entity_id='MONDO:0005148',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_get_mappings-parameters)

## 24) Find diseases by phenotype semantic similarity {#dnakit.databases.monarch_v3_phenotype_similarity_search}

- **Purpose:** Find diseases by phenotype semantic similarity.
- **Method:** Use `MonarchV3` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_v3_phenotype_similarity_search(*, phenotypes, group=None, limit=None, verbose=True)`
- **Input:** Provide `phenotypes`.
- **Output:** Returns a list for “Find diseases by phenotype semantic similarity”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_v3_phenotype_similarity_search(
    phenotypes=['HP:0001250', 'HP:0001249', 'HP:0002353'],
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_phenotype_similarity_search-parameters)

## 25) Search Monarch knowledge graph entities and categories {#dnakit.databases.monarch_v3_search}

- **Purpose:** Search Monarch knowledge graph entities and categories.
- **Method:** Use `MonarchV3` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.monarch_v3_search(*, query=None, limit=None, category=None, q=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list for “Search Monarch knowledge graph entities and categories”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.monarch_v3_search(
    query='BRCA1',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<value>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.monarch_v3_search-parameters)

## 26) Retrieve ontology term ancestors {#dnakit.databases.ols_get_term_ancestors}

- **Purpose:** Retrieve ontology term ancestors.
- **Method:** Use `OLS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ols_get_term_ancestors(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `page`, `numElements`, `totalPages`, `totalElements` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ols_get_term_ancestors(
    term_iri='http://www.ebi.ac.uk/efo/EFO_0000384',
    ontology='efo',
    include_obsolete=False,
    size=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "page": "<integer>",
  "numElements": "<integer>",
  "totalPages": "<integer>",
  "totalElements": "<integer>",
  "elements": [
    "<value>"
  ],
  "term_iri": "<string>",
  "ontology": "<string>",
  "filters": {
    "include_obsolete": "<boolean>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_get_term_ancestors-parameters)

## 27) Retrieve ontology term children {#dnakit.databases.ols_get_term_children}

- **Purpose:** Retrieve ontology term children.
- **Method:** Use `OLS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ols_get_term_children(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `page`, `numElements`, `totalPages`, `totalElements` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ols_get_term_children(
    term_iri='http://www.ebi.ac.uk/efo/EFO_0000384',
    ontology='efo',
    include_obsolete=False,
    size=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "page": "<integer>",
  "numElements": "<integer>",
  "totalPages": "<integer>",
  "totalElements": "<integer>",
  "elements": [
    "<value>"
  ],
  "term_iri": "<string>",
  "ontology": "<string>",
  "filters": {
    "include_obsolete": "<boolean>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_get_term_children-parameters)

## 28) Retrieve ontology term details {#dnakit.databases.ols_get_term_info}

- **Purpose:** Retrieve ontology term details.
- **Method:** Use `OLS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ols_get_term_info(*, operation=None, id=None, term_id=None, term_iri=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `iri`, `ontologyName`, `shortForm`, `label` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ols_get_term_info(
    id='HP:0001903',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "iri": "<string>",
  "ontologyName": "<string>",
  "shortForm": "<string>",
  "label": "<string>",
  "oboId": "<string>",
  "isObsolete": "<boolean>",
  "description": [
    "<value>"
  ],
  "synonyms": [
    "<value>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_get_term_info-parameters)

## 29) Search ontology terms by keyword {#dnakit.databases.ols_search_terms}

- **Purpose:** Search ontology terms by keyword.
- **Method:** Use `OLS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ols_search_terms(*, operation=None, query, rows=None, size=None, ontology=None, exact_match=None, include_obsolete=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `response`, `responseHeader`, `facet_counts`, `query` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ols_search_terms(
    query='gene:TP53',
    rows=10,
    ontology='efo',
    exact_match=False,
    include_obsolete=False,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "response": {
    "docs": [
      "<value>"
    ],
    "numFound": "<integer>",
    "start": "<integer>"
  },
  "responseHeader": {
    "QTime": "<integer>",
    "status": "<integer>"
  },
  "facet_counts": {
    "facet_fields": {
      "ontologyPreferredPrefix": "<array>",
      "isDefiningOntology": "<array>",
      "ontologyId": "<array>",
      "ontologyIri": "<array>",
      "isObsolete": "<array>",
      "type": "<array>"
    }
  },
  "query": "<string>",
  "filters": {
    "ontology": "<string>",
    "exact_match": "<boolean>",
    "exact_match_fields": "<string>",
    "include_obsolete": "<boolean>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.ols_search_terms-parameters)

## 30) Retrieve an OMIM clinical phenotype synopsis {#dnakit.databases.omim_get_clinical_synopsis}

- **Purpose:** Retrieve an OMIM clinical phenotype synopsis.
- **Method:** Use `OMIM` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.omim_get_clinical_synopsis(*, operation=None, mim_number, verbose=True)`
- **Input:** Provide `mim_number`. Configure: `OMIM_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.omim_get_clinical_synopsis(
    mim_number='219700',
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
    "mimNumber": "<integer>",
    "inheritance": "<string>",
    "neurologicCentralNervousSystem": [
      "<value>"
    ],
    "cardiovascular": [
      "<value>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_get_clinical_synopsis-parameters)

## 31) Retrieve an OMIM entry by MIM number {#dnakit.databases.omim_get_entry}

- **Purpose:** Retrieve an OMIM entry by MIM number.
- **Method:** Use `OMIM` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.omim_get_entry(*, operation=None, mim_number, include=None, verbose=True)`
- **Input:** Provide `mim_number`. Configure: `OMIM_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.omim_get_entry(
    mim_number='164730',
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
    "mimNumber": "<integer>",
    "titles": "<object>",
    "textSectionList": [
      "<value>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_get_entry-parameters)

## 32) Retrieve OMIM gene–disease and chromosome mappings {#dnakit.databases.omim_get_gene_map}

- **Purpose:** Retrieve OMIM gene–disease and chromosome mappings.
- **Method:** Use `OMIM` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.omim_get_gene_map(*, operation=None, mim_number=None, chromosome=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Configure: `OMIM_API_KEY`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.omim_get_gene_map(
    mim_number='164730',
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
    "gene_maps": [
      "<value>"
    ],
    "count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_get_gene_map-parameters)

## 33) Search Mendelian disorders and associated genes {#dnakit.databases.omim_search}

- **Purpose:** Search Mendelian disorders and associated genes.
- **Method:** Use `OMIM` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.omim_search(*, operation=None, query, limit=None, start=None, verbose=True)`
- **Input:** Provide `query`. Configure: `OMIM_API_KEY`.
- **Output:** Returns a mapping with `status`, `data`, `error` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.omim_search(
    query='BRCA1',
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
    "total_results": "<integer>",
    "entries": [
      "<value>"
    ]
  },
  "error": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.omim_search-parameters)

## 34) Retrieve a rare disease classification hierarchy {#dnakit.databases.orphanet_get_classification}

- **Purpose:** Retrieve a rare disease classification hierarchy.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_classification(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_classification(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "parent": "<object>",
    "children": [
      "<value>"
    ],
    "member_of_classifications": "<object>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_classification-parameters)

## 35) Retrieve rare disease details by ORPHA code {#dnakit.databases.orphanet_get_disease}

- **Purpose:** Retrieve rare disease details by ORPHA code.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_disease(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_disease(
    orpha_code='558',
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
    "ORPHAcode": "<string>",
    "Preferred term": "<string>",
    "Definition": "<string>",
    "Synonyms": [
      "<value>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_disease-parameters)

## 36) Retrieve reported gene associations for an ORPHA disorder {#dnakit.databases.orphanet_get_disease_genes}

- **Purpose:** Retrieve reported gene associations for an ORPHA disorder.
- **Method:** Use `API` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_disease_genes(*, orpha_code, verbose=True)`
- **Input:** Provide `orpha_code`.
- **Output:** Returns the Orphanet_get_genes association result; see [orphanet_get_genes](22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes).
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_disease_genes(
    orpha_code='558',
    verbose=False,
)
print(result)
```

- **Example output:** Returns the Orphanet_get_genes association result; see [orphanet_get_genes](22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes).

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_disease_genes-parameters)

## 37) Retrieve rare disease prevalence and incidence {#dnakit.databases.orphanet_get_epidemiology}

- **Purpose:** Retrieve rare disease prevalence and incidence.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_epidemiology(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_epidemiology(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "preferred_term": "<string>",
    "prevalences": [
      "<object>"
    ],
    "prevalence_count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_epidemiology-parameters)

## 38) Find rare diseases associated with a gene {#dnakit.databases.orphanet_get_gene_diseases}

- **Purpose:** Find rare diseases associated with a gene.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_gene_diseases(*, operation=None, gene_name=None, gene_symbol=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_gene_diseases(
    gene_name='FBN1',
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
    "gene_name": "<string>",
    "diseases": [
      "<object>"
    ],
    "disease_count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_gene_diseases-parameters)

## 39) Retrieve rare disease genes and association types {#dnakit.databases.orphanet_get_genes}

- **Purpose:** Retrieve rare disease genes and association types.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_genes(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_genes(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "genes": [
      "<object>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_genes-parameters)

## 40) Map rare disease codes to other medical coding systems {#dnakit.databases.orphanet_get_icd_mapping}

- **Purpose:** Map rare disease codes to other medical coding systems.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_icd_mapping(*, operation=None, orpha_code=None, coding_system=None, lang=None, orphacode=None, orpha_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_icd_mapping(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "mappings": "<object>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_icd_mapping-parameters)

## 41) Retrieve rare disease onset and inheritance information {#dnakit.databases.orphanet_get_natural_history}

- **Purpose:** Retrieve rare disease onset and inheritance information.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_natural_history(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_natural_history(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "preferred_term": "<string>",
    "average_age_of_onset": [
      "<string>"
    ],
    "type_of_inheritance": [
      "<string>"
    ],
    "disorder_group": "<string>",
    "typology": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_natural_history-parameters)

## 42) Retrieve rare disease phenotypes and frequencies {#dnakit.databases.orphanet_get_phenotypes}

- **Purpose:** Retrieve rare disease phenotypes and frequencies.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_get_phenotypes(*, operation=None, orpha_code=None, orpha_id=None, orphacode=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_get_phenotypes(
    orpha_code='558',
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
    "orpha_code": "<string>",
    "preferred_term": "<string>",
    "phenotypes": [
      "<object>"
    ],
    "phenotype_count": "<integer>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_get_phenotypes-parameters)

## 43) Search rare diseases by complete or partial name {#dnakit.databases.orphanet_search_by_name}

- **Purpose:** Search rare diseases by complete or partial name.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_search_by_name(*, operation=None, name=None, exact=None, lang=None, limit=None, query=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_search_by_name(
    name='Marfan syndrome',
    exact=True,
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
    "results": [
      "<value>"
    ],
    "count": "<integer>",
    "exact_match": "<boolean>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_search_by_name-parameters)

## 44) Search rare diseases by keyword {#dnakit.databases.orphanet_search_diseases}

- **Purpose:** Search rare diseases by keyword.
- **Method:** Use `Orphanet` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.orphanet_search_diseases(*, operation=None, query, limit=None, lang=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.orphanet_search_diseases(
    query='Marfan syndrome',
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
    "results": [
      "<value>"
    ],
    "count": "<integer>",
    "total_count": "<integer>",
    "query": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.orphanet_search_diseases-parameters)

## 45) Resolve a disease name to an EFO term {#dnakit.databases.osl_get_efo_id_by_disease_name}

- **Purpose:** Resolve a disease name to an EFO term.
- **Method:** Use `EFO` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.osl_get_efo_id_by_disease_name(*, disease, verbose=True)`
- **Input:** Provide `disease`.
- **Output:** Returns a mapping with `efo_id`, `name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.osl_get_efo_id_by_disease_name(
    disease='Diabetes Mellitus',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "efo_id": "<string>",
  "name": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_phenotypes.md#dnakit.databases.osl_get_efo_id_by_disease_name-parameters)

<!-- END GENERATED API: phenotypes -->
