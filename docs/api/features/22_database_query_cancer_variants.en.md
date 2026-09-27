# Cancer Variants and Treatment Evidence

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: cancer-variants -->

<span id="cancer-variants"></span>

## 1) Retrieve a CIViC molecular profile {#dnakit.databases.civic_get_molecular_profile}

- **Purpose:** Retrieve a CIViC molecular profile.
- **Method:** Use `CIViC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.civic_get_molecular_profile(*, molecular_profile_id, evidence_limit=None, verbose=True)`
- **Input:** Provide `molecular_profile_id`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.civic_get_molecular_profile(
    molecular_profile_id=12,
    evidence_limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "molecularProfile": {
      "id": "<integer>",
      "name": "<string>",
      "molecularProfileScore": "<number>",
      "description": "<string>",
      "evidenceItems": "<object>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_get_molecular_profile-parameters)

## 2) Retrieve a CIViC cancer variant {#dnakit.databases.civic_get_variant}

- **Purpose:** Retrieve a CIViC cancer variant.
- **Method:** Use `CIViC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.civic_get_variant(*, variant_id, verbose=True)`
- **Input:** Provide `variant_id`.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.civic_get_variant(
    variant_id=2,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "variant": {
      "id": "<integer>",
      "name": "<string>",
      "molecularProfiles": "<object>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_get_variant-parameters)

## 3) Find CIViC variants associated with a gene {#dnakit.databases.civic_get_variants_by_gene}

- **Purpose:** Find CIViC variants associated with a gene.
- **Method:** Use `CIViC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.civic_get_variants_by_gene(*, gene_id=None, gene_name=None, gene_symbol=None, gene=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.civic_get_variants_by_gene(
    gene_name='EGFR',
    limit=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "gene": {
      "id": "<integer>",
      "name": "<string>",
      "variants": "<object>"
    }
  },
  "metadata": {
    "source": "<string>",
    "format": "<string>",
    "note": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_get_variants_by_gene-parameters)

## 4) Search CIViC clinical assertions {#dnakit.databases.civic_search_assertions}

- **Purpose:** Search CIViC clinical assertions.
- **Method:** Use `CIViC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.civic_search_assertions(*, limit=None, therapy=None, therapy_name=None, disease=None, disease_name=None, variant_name=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.civic_search_assertions(
    limit=5,
    therapy='imatinib',
    disease='leukemia',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "assertions": {
      "totalCount": "<integer>",
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_assertions-parameters)

## 5) Search CIViC clinical evidence items {#dnakit.databases.civic_search_evidence_items}

- **Purpose:** Search CIViC clinical evidence items.
- **Method:** Use `CIViC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.civic_search_evidence_items(*, limit=None, status=None, therapy=None, therapy_name=None, disease=None, disease_name=None, molecular_profile=None, evidence_type=None, significance=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.civic_search_evidence_items(
    limit=5,
    therapy='imatinib',
    disease='leukemia',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "evidenceItems": {
      "totalCount": "<integer>",
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_evidence_items-parameters)

## 6) Search genes in CIViC {#dnakit.databases.civic_search_genes}

- **Purpose:** Search genes in CIViC.
- **Method:** Use `CIViC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.civic_search_genes(*, limit=None, name=None, query=None, gene_name=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.civic_search_genes(
    limit=5,
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
      "totalCount": "<integer>",
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  },
  "metadata": {
    "source": "<string>",
    "format": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_genes-parameters)

## 7) Search therapies in CIViC {#dnakit.databases.civic_search_therapies}

- **Purpose:** Search therapies in CIViC.
- **Method:** Use `CIViC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.civic_search_therapies(*, limit=None, name=None, query=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.civic_search_therapies(
    name='imatinib',
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "browseTherapies": {
      "pageInfo": "<object>",
      "nodes": "<array>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.civic_search_therapies-parameters)

## 8) Retrieve COSMIC somatic mutations for a gene {#dnakit.databases.cosmic_get_mutations_by_gene}

- **Purpose:** Retrieve COSMIC somatic mutations for a gene.
- **Method:** Use `COSMIC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cosmic_get_mutations_by_gene(*, operation=None, gene=None, gene_name=None, max_results=None, genome_build=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `gene`, `total_count`, `results`, `genome_build` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cosmic_get_mutations_by_gene(
    operation='get_by_gene',
    gene='BRAF',
    max_results=50,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene": "<string>",
  "total_count": "<integer>",
  "results": [
    {
      "mutation_id": "<string>",
      "GeneName": "<string>",
      "MutationCDS": "<string>",
      "MutationAA": "<string>",
      "MutationGenomePosition": "<string>",
      "PrimarySite": "<string>",
      "PrimaryHistology": "<string>",
      "Fathmm": "<string>"
    }
  ],
  "genome_build": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.cosmic_get_mutations_by_gene-parameters)

## 9) Search COSMIC cancer mutations {#dnakit.databases.cosmic_search_mutations}

- **Purpose:** Search COSMIC cancer mutations.
- **Method:** Use `COSMIC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.cosmic_search_mutations(*, operation=None, terms=None, query=None, max_results=None, genome_build=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total_count`, `results`, `genome_build` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.cosmic_search_mutations(
    operation='search',
    terms='BRAF V600E',
    max_results=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total_count": "<integer>",
  "results": [
    {
      "mutation_id": "<string>",
      "display": "<array>",
      "GeneName": "<string>",
      "MutationCDS": "<string>",
      "MutationAA": "<string>",
      "PrimarySite": "<string>",
      "PrimaryHistology": "<string>"
    }
  ],
  "genome_build": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.cosmic_search_mutations-parameters)

## 10) Retrieve OncoKB oncogenicity and treatment evidence {#dnakit.databases.onco_kb_annotate_variant}

- **Purpose:** Retrieve OncoKB oncogenicity and treatment evidence.
- **Method:** Use `OncoKB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.onco_kb_annotate_variant(*, operation=None, gene=None, gene_symbol=None, variant=None, tumor_type=None, alteration=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `query`, `geneExist`, `variantExist`, `mutationEffect` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.onco_kb_annotate_variant(
    gene='BRAF',
    variant='V600E',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "query": {
    "id": "<null>",
    "referenceGenome": "<string>",
    "hugoSymbol": "<string>",
    "entrezGeneId": "<integer>",
    "alteration": "<string>",
    "alterationType": "<null>",
    "svType": "<null>",
    "tumorType": "<string>"
  },
  "geneExist": "<boolean>",
  "variantExist": "<boolean>",
  "mutationEffect": {
    "knownEffect": "<string>",
    "description": "<string>",
    "citations": {
      "pmids": "<array>",
      "abstracts": "<array>"
    }
  },
  "highestSensitiveLevel": "<string>",
  "highestResistanceLevel": "<null>",
  "highestDiagnosticImplicationLevel": "<string>",
  "highestPrognosticImplicationLevel": "<null>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.onco_kb_annotate_variant-parameters)

## 11) Search cancer samples with regional copy number changes {#dnakit.databases.progenetix_cnv_search}

- **Purpose:** Search cancer samples with regional copy number changes.
- **Method:** Use `Progenetix` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.progenetix_cnv_search(*, reference_name, start, end, variant_type=None, filters=None, limit=None, verbose=True)`
- **Input:** Provide `reference_name`, `start`, `end`.
- **Output:** Returns a mapping with `region`, `variant_type`, `filters`, `total_count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.progenetix_cnv_search(
    reference_name='refseq:NC_000007.14',
    start=55019017,
    end=55211628,
    variant_type='DUP',
    filters='NCIT:C4017',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "region": "<string>",
  "variant_type": "<string>",
  "filters": "<string>",
  "total_count": "<integer>",
  "returned_count": "<integer>",
  "biosamples": [
    {
      "id": "<string>",
      "biosample_status": "<string>",
      "histological_diagnosis": "<string>",
      "histological_diagnosis_id": "<string>",
      "external_references": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_cnv_search-parameters)

## 12) List Progenetix cancer CNV cohorts {#dnakit.databases.progenetix_list_cohorts}

- **Purpose:** List Progenetix cancer CNV cohorts.
- **Method:** Use `Progenetix` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.progenetix_list_cohorts(*, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `total_cohorts`, `cohorts` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.progenetix_list_cohorts(
    limit=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total_cohorts": "<integer>",
  "cohorts": [
    {
      "id": "<string>",
      "name": "<string>",
      "cohort_size": "<integer>",
      "cohort_type": "<string>",
      "data_types": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_list_cohorts-parameters)

## 13) List Progenetix filtering terms {#dnakit.databases.progenetix_list_filtering_terms}

- **Purpose:** List Progenetix filtering terms.
- **Method:** Use `Progenetix` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.progenetix_list_filtering_terms(*, prefixes=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `prefixes`, `total_terms`, `terms` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.progenetix_list_filtering_terms(
    prefixes='NCIT',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "prefixes": "<string>",
  "total_terms": "<integer>",
  "terms": [
    {
      "id": "<string>",
      "label": "<string>",
      "count": "<integer>",
      "type": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_list_filtering_terms-parameters)

## 14) Search cancer biosamples by disease ontology {#dnakit.databases.progenetix_search_biosamples}

- **Purpose:** Search cancer biosamples by disease ontology.
- **Method:** Use `Progenetix` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.progenetix_search_biosamples(*, filters, limit=None, skip=None, verbose=True)`
- **Input:** Provide `filters`.
- **Output:** Returns a mapping with `filters`, `total_count`, `returned_count`, `biosamples` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.progenetix_search_biosamples(
    filters='NCIT:C3058',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "filters": "<string>",
  "total_count": "<integer>",
  "returned_count": "<integer>",
  "biosamples": [
    {
      "id": "<string>",
      "biosample_status": "<string>",
      "histological_diagnosis": "<string>",
      "histological_diagnosis_id": "<string>",
      "pathological_stage": "<string>",
      "pathological_tnm": "<string>",
      "sample_origin_type": "<string>",
      "external_references": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_cancer_variants.md#dnakit.databases.progenetix_search_biosamples-parameters)

<!-- END GENERATED API: cancer-variants -->
