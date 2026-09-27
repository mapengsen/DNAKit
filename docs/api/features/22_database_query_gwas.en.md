# GWAS Studies and Associations

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: gwas -->

<span id="gwas"></span>

## 1) Retrieve GWAS associations for a SNP {#dnakit.databases.gwas_get_associations_for_snp}

- **Purpose:** Retrieve GWAS associations for a SNP.
- **Method:** Use `GWASAssociationsForSNP` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_associations_for_snp(*, rs_id, sort=None, direction=None, size=None, page=None, verbose=True)`
- **Input:** Provide `rs_id`.
- **Output:** Returns a list of records; each record includes `association_id`, `p_value`, `beta`, `efo_traits` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_associations_for_snp(
    rs_id='rs7903146',
    sort='p_value',
    direction='asc',
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
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_associations_for_snp-parameters)

## 2) Retrieve associations from a GWAS study {#dnakit.databases.gwas_get_associations_for_study}

- **Purpose:** Retrieve associations from a GWAS study.
- **Method:** Use `GWASAssociationsForStudy` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_associations_for_study(*, accession_id, size=None, page=None, verbose=True)`
- **Input:** Provide `accession_id`.
- **Output:** Returns a list of records; each record includes `association_id`, `p_value`, `beta`, `efo_traits` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_associations_for_study(
    accession_id='GCST000392',
    size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_associations_for_study-parameters)

## 3) Retrieve GWAS associations for a trait {#dnakit.databases.gwas_get_associations_for_trait}

- **Purpose:** Retrieve GWAS associations for a trait.
- **Method:** Use `GWASAssociationsForTrait` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_associations_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, page=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select exactly one argument group: `disease_trait` / `efo_uri` / `efo_id` / `efo_trait`.
- **Output:** Returns a list of records; each record includes `association_id`, `p_value`, `beta`, `efo_traits` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_associations_for_trait(
    disease_trait='type 2 diabetes',
    size=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_associations_for_trait-parameters)

## 4) Retrieve GWAS SNP details by rsID {#dnakit.databases.gwas_get_snp_by_id}

- **Purpose:** Retrieve GWAS SNP details by rsID.
- **Method:** Use `GWASSNPByID` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_snp_by_id(*, rs_id, verbose=True)`
- **Input:** Provide `rs_id`.
- **Output:** Returns a mapping with `rs_id`, `merged`, `functional_class`, `last_update_date` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_snp_by_id(
    rs_id='rs429358',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "rs_id": "<string>",
  "merged": "<integer>",
  "functional_class": "<string>",
  "last_update_date": "<string>",
  "locations": [
    {
      "chromosome_name": "<string>",
      "chromosome_position": "<integer>",
      "region": "<object>"
    }
  ],
  "alleles": "<string>",
  "most_severe_consequence": "<string>",
  "mapped_genes": [
    "<string>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_snp_by_id-parameters)

## 5) Retrieve GWAS SNPs mapped to a gene {#dnakit.databases.gwas_get_snps_for_gene}

- **Purpose:** Retrieve GWAS SNPs mapped to a gene.
- **Method:** Use `GWASSNPsForGene` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_snps_for_gene(*, size=None, page=None, gene_symbol=None, gene=None, mapped_gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `rs_id`, `merged`, `functional_class`, `last_update_date` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_snps_for_gene(
    mapped_gene='BRCA1',
    size=10,
    page=1,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "rs_id": "<string>",
    "merged": "<integer>",
    "functional_class": "<string>",
    "last_update_date": "<string>",
    "locations": [
      "<object>"
    ],
    "alleles": "<string>",
    "most_severe_consequence": "<string>",
    "mapped_genes": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_snps_for_gene-parameters)

## 6) Retrieve GWAS studies for a trait {#dnakit.databases.gwas_get_studies_for_trait}

- **Purpose:** Retrieve GWAS studies for a trait.
- **Method:** Use `GWASStudiesForTrait` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_studies_for_trait(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select exactly one argument group: `disease_trait` / `efo_uri` / `efo_id` / `efo_trait`.
- **Output:** Returns a list of records; each record includes `accession_id`, `disease_trait`, `efo_traits`, `initial_sample_size` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_studies_for_trait(
    disease_trait='hypertension',
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
    "accession_id": "<string>",
    "disease_trait": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "initial_sample_size": "<string>",
    "replication_sample_size": "<string>",
    "discovery_ancestry": [
      "<string>"
    ],
    "replication_ancestry": [
      "<string>"
    ],
    "cohort": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_studies_for_trait-parameters)

## 7) Retrieve a GWAS study by accession {#dnakit.databases.gwas_get_study_by_id}

- **Purpose:** Retrieve a GWAS study by accession.
- **Method:** Use `GWASStudyByID` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_study_by_id(*, study_id=None, accession_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `accession_id`, `disease_trait`, `efo_traits`, `initial_sample_size` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_study_by_id(
    study_id='GCST000392',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "accession_id": "<string>",
  "disease_trait": "<string>",
  "efo_traits": [
    {
      "efo_id": "<string>",
      "efo_trait": "<string>"
    }
  ],
  "initial_sample_size": "<string>",
  "replication_sample_size": "<string>",
  "discovery_ancestry": [
    "<string>"
  ],
  "replication_ancestry": [
    "<string>"
  ],
  "cohort": [
    "<string>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_study_by_id-parameters)

## 8) Find GWAS variants associated with a disease or trait {#dnakit.databases.gwas_get_variants_for_trait}

- **Purpose:** Find GWAS variants associated with a disease or trait.
- **Method:** Use `GWASVariantsForTrait` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_get_variants_for_trait(*, disease_trait=None, trait=None, efo_uri=None, efo_id=None, efo_trait=None, size=None, limit=None, page=None, sort=None, direction=None, p_value_threshold=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select exactly one argument group: `disease_trait` / `trait` / `efo_uri` / `efo_id` / `efo_trait`.
- **Output:** Returns a list of records; each record includes `association_id`, `risk_frequency`, `pvalue_description`, `range` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_get_variants_for_trait(
    disease_trait='breast cancer',
    size=20,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "association_id": "<integer>",
    "risk_frequency": "<string>",
    "pvalue_description": "<string>",
    "range": "<string>",
    "locations": [
      "<object>"
    ],
    "rs_id": "<string>",
    "merged": "<integer>",
    "functional_class": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_get_variants_for_trait-parameters)

## 9) Search GWAS associations by keyword {#dnakit.databases.gwas_search_associations}

- **Purpose:** Search GWAS associations by keyword.
- **Method:** Use `GWASAssociationSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_search_associations(*, disease_trait=None, query=None, efo_uri=None, efo_id=None, efo_trait=None, rs_id=None, accession_id=None, p_value=None, p_value_threshold=None, sort=None, direction=None, size=None, page=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `association_id`, `p_value`, `beta`, `efo_traits` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_search_associations(
    disease_trait='coronary artery disease',
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
    "association_id": "<integer>",
    "p_value": "<number>",
    "beta": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "reported_trait": [
      "<string>"
    ],
    "accession_id": "<string>",
    "pubmed_id": "<string>",
    "first_author": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_associations-parameters)

## 10) Search gene GWAS associations and risk alleles {#dnakit.databases.gwas_search_associations_by_gene}

- **Purpose:** Search gene GWAS associations and risk alleles.
- **Method:** Use `GWASGeneSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_search_associations_by_gene(*, gene_name, size=None, verbose=True)`
- **Input:** Provide `gene_name`.
- **Output:** Returns a mapping with `gene_name`, `association_count`, `associations`, `total_found` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_search_associations_by_gene(
    gene_name='TP53',
    size=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_name": "<string>",
  "association_count": "<integer>",
  "associations": [
    "<object>"
  ],
  "total_found": "<integer>",
  "queried_gene_name": "<string>",
  "gene_name_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_associations_by_gene-parameters)

## 11) Search GWAS SNPs by rsID or gene {#dnakit.databases.gwas_search_snps}

- **Purpose:** Search GWAS SNPs by rsID or gene.
- **Method:** Use `GWASSNPSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_search_snps(*, rs_id=None, mapped_gene=None, size=None, page=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `rs_id`, `merged`, `functional_class`, `last_update_date` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_search_snps(
    mapped_gene='BRCA1',
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
    "rs_id": "<string>",
    "merged": "<integer>",
    "functional_class": "<string>",
    "last_update_date": "<string>",
    "locations": [
      "<object>"
    ],
    "alleles": "<string>",
    "most_severe_consequence": "<string>",
    "mapped_genes": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_snps-parameters)

## 12) Search GWAS studies by trait and cohort {#dnakit.databases.gwas_search_studies}

- **Purpose:** Search GWAS studies by trait and cohort.
- **Method:** Use `GWASStudySearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.gwas_search_studies(*, disease_trait=None, efo_uri=None, efo_id=None, efo_trait=None, cohort=None, gxe=None, full_pvalue_set=None, size=None, page=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `accession_id`, `disease_trait`, `efo_traits`, `initial_sample_size` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.gwas_search_studies(
    disease_trait='diabetes',
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
    "accession_id": "<string>",
    "disease_trait": "<string>",
    "efo_traits": [
      "<object>"
    ],
    "initial_sample_size": "<string>",
    "replication_sample_size": "<string>",
    "discovery_ancestry": [
      "<string>"
    ],
    "replication_ancestry": [
      "<string>"
    ],
    "cohort": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_gwas.md#dnakit.databases.gwas_search_studies-parameters)

<!-- END GENERATED API: gwas -->
