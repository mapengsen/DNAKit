# Pathway, interaction, protein and structure queries

Query pathways, molecular interactions, proteins and structures.

<!-- BEGIN GENERATED API: pathways -->

<span id="pathways"></span>

## 1) Retrieve stored AlphaFold protein structure predictions {#dnakit.databases.alphafold_get_prediction}

- **Purpose:** Retrieve stored AlphaFold protein structure predictions.
- **Method:** Use `AlphaFold` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.alphafold_get_prediction(*, qualifier=None, uniprot_id=None, uniprot_accession=None, accession=None, sequence_checksum=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `toolUsed`, `providerId`, `entityType`, `isUniProt` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.alphafold_get_prediction(
    qualifier='P69905',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "toolUsed": "<string>",
    "providerId": "<string>",
    "entityType": "<string>",
    "isUniProt": "<boolean>",
    "modelEntityId": "<string>",
    "modelCreatedDate": "<string>",
    "sequenceVersionDate": "<string>",
    "globalMetricValue": "<number>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.alphafold_get_prediction-parameters)

## 2) Retrieve experimentally supported protein and genetic interactions {#dnakit.databases.bio_grid_get_interactions}

- **Purpose:** Retrieve experimentally supported protein and genetic interactions.
- **Method:** Use `BioGRID` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.bio_grid_get_interactions(*, gene_names=None, gene_name=None, organism=None, interaction_type=None, evidence_types=None, limit=None, throughput=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Configure: `BIOGRID_ACCESS_KEY`.
- **Output:** Returns a mapping keyed by record; each record includes `BIOGRID_INTERACTION_ID`, `ENTREZ_GENE_A`, `ENTREZ_GENE_B`, `OFFICIAL_SYMBOL_A` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.bio_grid_get_interactions(
    gene_names=['TP53'],
    organism='9606',
    interaction_type='physical',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "<record key>": {
    "BIOGRID_INTERACTION_ID": "<integer>",
    "ENTREZ_GENE_A": "<string>",
    "ENTREZ_GENE_B": "<string>",
    "OFFICIAL_SYMBOL_A": "<string>",
    "OFFICIAL_SYMBOL_B": "<string>",
    "SYNONYMS_A": "<string>",
    "SYNONYMS_B": "<string>",
    "EXPERIMENTAL_SYSTEM": "<string>"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.bio_grid_get_interactions-parameters)

## 3) Search pathogen B-cell and T-cell epitopes {#dnakit.databases.bvbrc_search_epitopes}

- **Purpose:** Search pathogen B-cell and T-cell epitopes.
- **Method:** Use `BVBRC` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.bvbrc_search_epitopes(*, taxon_id=None, protein_name=None, epitope_type=None, organism=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `epitope_id`, `epitope_type`, `epitope_sequence`, `organism` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.bvbrc_search_epitopes(
    taxon_id='2697049',
    protein_name='Spike glycoprotein',
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
    "epitope_id": "<string>",
    "epitope_type": "<string>",
    "epitope_sequence": "<string>",
    "organism": "<string>",
    "protein_name": "<string>",
    "start": "<integer>",
    "end": "<integer>",
    "bcell_assays": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.bvbrc_search_epitopes-parameters)

## 4) Retrieve GO annotations for a gene or protein {#dnakit.databases.go_get_annotations_for_gene}

- **Purpose:** Retrieve GO annotations for a gene or protein.
- **Method:** Use `GeneOntology` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.go_get_annotations_for_gene(*, gene_id, rows=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a list of records; each record includes `document_category`, `source`, `bioentity`, `bioentity_internal_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.go_get_annotations_for_gene(
    gene_id='TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "document_category": "<string>",
    "source": "<string>",
    "bioentity": "<string>",
    "bioentity_internal_id": "<string>",
    "bioentity_label": "<string>",
    "bioentity_label_searchable": "<string>",
    "annotation_class": "<string>",
    "annotation_class_label": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.go_get_annotations_for_gene-parameters)

## 5) Search GO terms by keyword {#dnakit.databases.go_search_terms}

- **Purpose:** Search GO terms by keyword.
- **Method:** Use `GeneOntology` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.go_search_terms(*, query, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `responseHeader`, `response` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.go_search_terms(
    query='gene:TP53',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "responseHeader": {
    "status": "<integer>",
    "QTime": "<integer>",
    "params": {
      "q": "<string>",
      "fq": "<array>",
      "rows": "<string>",
      "wt": "<string>"
    }
  },
  "response": {
    "numFound": "<integer>",
    "start": "<integer>",
    "docs": [
      "<value>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.go_search_terms-parameters)

## 6) Retrieve tissue-specific protein interaction networks {#dnakit.databases.humanbase_ppi_analysis}

- **Purpose:** Retrieve tissue-specific protein interaction networks.
- **Method:** Use `HumanBase` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.humanbase_ppi_analysis(*, gene_list=None, tissue=None, max_node=None, interaction=None, string_mode=None, genes=None, top_n=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a string.
- **Example code:**

```python
import dnakit

result = dnakit.databases.humanbase_ppi_analysis(
    gene_list=['TP53', 'BRCA1'],
    tissue='brain',
    max_node=3,
    interaction='co-expression',
    string_mode=True,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
"<string>"
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.humanbase_ppi_analysis-parameters)

## 7) Retrieve antigens linked to an immune epitope {#dnakit.databases.iedb_get_epitope_antigens}

- **Purpose:** Retrieve antigens linked to an immune epitope.
- **Method:** Use `IEDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.iedb_get_epitope_antigens(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **Input:** Provide `structure_id`.
- **Output:** Returns a list of records; each record includes `structure_id`, `parent_source_antigen_iri` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.iedb_get_epitope_antigens(
    structure_id=2,
    limit=2,
    select=['structure_id', 'parent_source_antigen_iri'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "structure_id": "<integer>",
    "parent_source_antigen_iri": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_get_epitope_antigens-parameters)

## 8) Retrieve epitope MHC restriction associations {#dnakit.databases.iedb_get_epitope_mhc}

- **Purpose:** Retrieve epitope MHC restriction associations.
- **Method:** Use `IEDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.iedb_get_epitope_mhc(*, structure_id, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **Input:** Provide `structure_id`.
- **Output:** Returns a list of records; each record includes `structure_id`, `elution_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.iedb_get_epitope_mhc(
    structure_id=2,
    limit=2,
    select=['structure_id', 'elution_id'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "structure_id": "<integer>",
    "elution_id": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_get_epitope_mhc-parameters)

## 9) Search B-cell epitope assay records {#dnakit.databases.iedb_search_bcell}

- **Purpose:** Search B-cell epitope assay records.
- **Method:** Use `IEDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.iedb_search_bcell(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `bcell_id`, `structure_id`, `assay_description`, `reference_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.iedb_search_bcell(
    limit=2,
    select=['bcell_id', 'structure_id', 'assay_description', 'reference_id'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "bcell_id": "<integer>",
    "structure_id": "<integer>",
    "assay_description": "<string>",
    "reference_id": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_bcell-parameters)

## 10) Search immune epitopes and their structure identifiers {#dnakit.databases.iedb_search_epitopes}

- **Purpose:** Search immune epitopes and their structure identifiers.
- **Method:** Use `IEDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.iedb_search_epitopes(*, sequence_contains=None, structure_type=None, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `structure_id`, `structure_type`, `linear_sequence` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.iedb_search_epitopes(
    sequence_contains='KVF',
    structure_type='Linear peptide',
    limit=2,
    select=['structure_id', 'structure_type', 'linear_sequence'],
    order='structure_id.asc',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "structure_id": "<integer>",
    "structure_type": "<string>",
    "linear_sequence": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_epitopes-parameters)

## 11) Search MHC-related assay records {#dnakit.databases.iedb_search_mhc}

- **Purpose:** Search MHC-related assay records.
- **Method:** Use `IEDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.iedb_search_mhc(*, limit=None, offset=None, order=None, select=None, filters=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `structure_id`, `mhc_restriction`, `mhc_class`, `mhc_allele_resolution` and related fields.
- **Example code:**

```python
import dnakit

example_select = ['structure_id',
 'mhc_restriction',
 'mhc_class',
 'mhc_allele_resolution',
 'reference_id']

result = dnakit.databases.iedb_search_mhc(
    limit=2,
    select=example_select,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "structure_id": "<integer>",
    "mhc_restriction": "<string>",
    "mhc_class": "<string>",
    "mhc_allele_resolution": "<string>",
    "reference_id": "<integer>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_mhc-parameters)

## 12) Search T-cell response assays {#dnakit.databases.iedb_search_tcell_assays}

- **Purpose:** Search T-cell response assays.
- **Method:** Use `IEDB` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.iedb_search_tcell_assays(*, sequence=None, sequence_contains=None, mhc_class=None, qualitative_measure=None, limit=None, offset=None, select=None, filters=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `tcell_id`, `linear_sequence`, `structure_type`, `mhc_restriction` and related fields.
- **Example code:**

```python
import dnakit

example_select = ['tcell_id',
 'linear_sequence',
 'mhc_restriction',
 'qualitative_measure',
 'pubmed_id']

result = dnakit.databases.iedb_search_tcell_assays(
    sequence='GILGFVFTL',
    limit=3,
    select=example_select,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "tcell_id": "<integer>",
    "linear_sequence": "<string>",
    "structure_type": "<string>",
    "mhc_restriction": "<string>",
    "mhc_class": "<string>",
    "qualitative_measure": "<string>",
    "pubmed_id": "<string>",
    "journal_name": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.iedb_search_tcell_assays-parameters)

## 13) Retrieve immunogenetics resources and gene naming conventions {#dnakit.databases.imgt_get_gene_info}

- **Purpose:** Retrieve immunogenetics resources and gene naming conventions.
- **Method:** Use `IMGT` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.imgt_get_gene_info(*, operation=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.imgt_get_gene_info(
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.imgt_get_gene_info-parameters)

## 14) Locate immunoglobulin and TCR germline gene searches {#dnakit.databases.imgt_search_genes}

- **Purpose:** Locate immunoglobulin and TCR germline gene searches.
- **Method:** Use `IMGT` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.imgt_search_genes(*, operation=None, query=None, gene_type=None, species=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `data` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.imgt_search_genes(
    gene_type='IGHV',
    species='Homo sapiens',
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.imgt_search_genes-parameters)

## 15) Search immunology studies by disease and assay {#dnakit.databases.imm_port_search_studies}

- **Purpose:** Search immunology studies by disease and assay.
- **Method:** Use `ImmPort` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.imm_port_search_studies(*, query, condition_or_disease=None, assay_method=None, research_focus=None, species=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `study_accession`, `title`, `brief_description`, `condition_or_disease` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.imm_port_search_studies(
    query='influenza vaccine',
    species='Homo sapiens',
    limit=3,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "study_accession": "<string>",
    "title": "<string>",
    "brief_description": "<string>",
    "condition_or_disease": [
      "<value>"
    ],
    "research_focus": [
      "<value>"
    ],
    "species": [
      "<value>"
    ],
    "assay_methods": [
      "<value>"
    ],
    "assay_method_counts": [
      "<value>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.imm_port_search_studies-parameters)

## 16) Retrieve protein complex composition and function {#dnakit.databases.intact_get_complex_details}

- **Purpose:** Retrieve protein complex composition and function.
- **Method:** Use `IntAct` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.intact_get_complex_details(*, complex_ac, format=None, verbose=True)`
- **Input:** Provide `complex_ac`.
- **Output:** Returns a mapping with `complexAc`, `name`, `systematicName`, `synonyms` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.intact_get_complex_details(
    complex_ac='CPX-915',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "complexAc": "<string>",
  "name": "<string>",
  "systematicName": "<string>",
  "synonyms": [
    "<string>"
  ],
  "functions": [
    "<string>"
  ],
  "properties": [
    "<string>"
  ],
  "ac": "<string>",
  "species": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_get_complex_details-parameters)

## 17) Retrieve an interaction network centered on a protein {#dnakit.databases.intact_get_interaction_network}

- **Purpose:** Retrieve an interaction network centered on a protein.
- **Method:** Use `IntAct` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.intact_get_interaction_network(*, identifier=None, uniprot_id=None, protein_id=None, gene_symbol=None, gene_name=None, depth=None, limit=None, size=None, format=None, protein_name=None, protein=None, gene=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list for “Retrieve an interaction network centered on a protein”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.intact_get_interaction_network(
    identifier='P04637',
    depth=1,
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_get_interaction_network-parameters)

## 18) Retrieve interactions involving a protein or gene {#dnakit.databases.intact_get_interactions}

- **Purpose:** Retrieve interactions involving a protein or gene.
- **Method:** Use `IntAct` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.intact_get_interactions(*, identifier, format=None, size=None, verbose=True)`
- **Input:** Provide `identifier`.
- **Output:** Returns a list for “Retrieve interactions involving a protein or gene”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.intact_get_interactions(
    identifier='P04637',
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_get_interactions-parameters)

## 19) Search IntAct molecular interactions by keyword {#dnakit.databases.intact_search_interactions}

- **Purpose:** Search IntAct molecular interactions by keyword.
- **Method:** Use `IntAct` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.intact_search_interactions(*, query, first=None, max=None, format=None, species=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list for “Search IntAct molecular interactions by keyword”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.intact_search_interactions(
    query='P04637',
    max=10,
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.intact_search_interactions-parameters)

## 20) Retrieve protein domain and family entries {#dnakit.databases.inter_pro_get_entries_for_protein}

- **Purpose:** Retrieve protein domain and family entries.
- **Method:** Use `InterProEntry` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.inter_pro_get_entries_for_protein(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a mapping with `protein_accession`, `total_entries`, `entries` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.inter_pro_get_entries_for_protein(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "protein_accession": "<string>",
  "total_entries": "<integer>",
  "entries": [
    {
      "accession": "<string>",
      "name": "<string>",
      "type": "<string>",
      "source_database": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.inter_pro_get_entries_for_protein-parameters)

## 21) Retrieve protein domain and functional annotations {#dnakit.databases.inter_pro_get_protein_domains}

- **Purpose:** Retrieve protein domain and functional annotations.
- **Method:** Use `InterPro` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.inter_pro_get_protein_domains(*, protein_id=None, accession=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `metadata`, `proteins` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.inter_pro_get_protein_domains(
    protein_id='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "metadata": {
      "accession": "<string>",
      "name": "<string>",
      "source_database": "<string>",
      "type": "<string>",
      "integrated": "<null>",
      "member_databases": "<object>",
      "go_terms": "<array>"
    },
    "proteins": [
      "<object>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.inter_pro_get_protein_domains-parameters)

## 22) Search KEGG genes by keyword and organism {#dnakit.databases.kegg_find_genes}

- **Purpose:** Search KEGG genes by keyword and organism.
- **Method:** Use `KEGGFindGenes` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.kegg_find_genes(*, keyword=None, organism=None, query=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `gene_id`, `description` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.kegg_find_genes(
    keyword='insulin',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "gene_id": "<string>",
    "description": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_find_genes-parameters)

## 23) Retrieve KEGG gene sequences, functions and pathways {#dnakit.databases.kegg_get_gene_info}

- **Purpose:** Retrieve KEGG gene sequences, functions and pathways.
- **Method:** Use `KEGGGetGeneInfo` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.kegg_get_gene_info(*, gene_id, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `status`, `data`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.kegg_get_gene_info(
    gene_id='hsa:348',
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
    "gene_id": "<string>",
    "raw_data": "<string>",
    "line_count": "<integer>"
  },
  "url": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_gene_info-parameters)

## 24) Retrieve KEGG pathways containing a gene {#dnakit.databases.kegg_get_gene_pathways}

- **Purpose:** Retrieve KEGG pathways containing a gene.
- **Method:** Use `KEGGExt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.kegg_get_gene_pathways(*, gene_id=None, gene_symbol=None, gene=None, organism=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `gene_id`, `pathway_count`, `pathways` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.kegg_get_gene_pathways(
    gene_id='hsa:7157',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_id": "<string>",
  "pathway_count": "<integer>",
  "pathways": [
    {
      "pathway_id": "<string>",
      "pathway_name": "<string>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_gene_pathways-parameters)

## 25) List genes belonging to a KEGG pathway {#dnakit.databases.kegg_get_pathway_genes}

- **Purpose:** List genes belonging to a KEGG pathway.
- **Method:** Use `KEGGExt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.kegg_get_pathway_genes(*, pathway_id, verbose=True)`
- **Input:** Provide `pathway_id`.
- **Output:** Returns a mapping with `pathway_id`, `pathway_name`, `gene_count`, `genes` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.kegg_get_pathway_genes(
    pathway_id='hsa04115',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "pathway_id": "<string>",
  "pathway_name": "<string>",
  "gene_count": "<integer>",
  "genes": [
    "<string>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_pathway_genes-parameters)

## 26) Retrieve KEGG pathway genes, compounds and reactions {#dnakit.databases.kegg_get_pathway_info}

- **Purpose:** Retrieve KEGG pathway genes, compounds and reactions.
- **Method:** Use `KEGGGetPathwayInfo` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.kegg_get_pathway_info(*, pathway_id, verbose=True)`
- **Input:** Provide `pathway_id`.
- **Output:** Returns a mapping with `status`, `data`, `url` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.kegg_get_pathway_info(
    pathway_id='hsa00010',
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
    "pathway_id": "<string>",
    "raw_data": "<string>",
    "line_count": "<integer>"
  },
  "url": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_get_pathway_info-parameters)

## 27) Search KEGG pathways by keyword {#dnakit.databases.kegg_search_pathway}

- **Purpose:** Search KEGG pathways by keyword.
- **Method:** Use `KEGGSearchPathway` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.kegg_search_pathway(*, keyword, verbose=True)`
- **Input:** Provide `keyword`.
- **Output:** Returns a list of records; each record includes `pathway_id`, `description` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.kegg_search_pathway(
    keyword='diabetes',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "pathway_id": "<string>",
    "description": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.kegg_search_pathway-parameters)

## 28) Retrieve an NCBI protein record summary by identifier {#dnakit.databases.ncbiprotein_get_summary}

- **Purpose:** Retrieve an NCBI protein record summary by identifier.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.ncbiprotein_get_summary(*, id, verbose=True)`
- **Input:** Provide `id`.
- **Output:** Returns a mapping with `result` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.ncbiprotein_get_summary(
    id='119395750',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "result": {
    "uids": [
      "<string>"
    ]
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.ncbiprotein_get_summary-parameters)

## 29) Retrieve cross-species orthologs from PANTHER {#dnakit.databases.panther_ortholog}

- **Purpose:** Retrieve cross-species orthologs from PANTHER.
- **Method:** Use `PANTHER` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.panther_ortholog(*, gene_id, organism=None, target_organism=None, ortholog_type=None, verbose=True)`
- **Input:** Provide `gene_id`.
- **Output:** Returns a mapping with `gene_id`, `source_organism`, `target_organism`, `ortholog_type` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.panther_ortholog(
    gene_id='P04637',
    organism=9606,
    target_organism=10090,
    ortholog_type='LDO',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_id": "<string>",
  "source_organism": "<integer>",
  "target_organism": "<integer>",
  "ortholog_type": "<string>",
  "mapping": {
    "source_gene": "<string>",
    "target_gene": "<string>",
    "target_gene_symbol": "<string>",
    "ortholog_type": "<string>",
    "persistent_id": "<string>",
    "target_persistent_id": "<string>"
  },
  "mappings": [
    {
      "source_gene": "<string>",
      "target_gene": "<string>",
      "target_gene_symbol": "<string>",
      "ortholog_type": "<string>",
      "persistent_id": "<string>",
      "target_persistent_id": "<string>"
    }
  ],
  "total_mappings": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.panther_ortholog-parameters)

## 30) Retrieve PDB coordinate file text {#dnakit.databases.pdb_get_structure}

- **Purpose:** Retrieve PDB coordinate file text.
- **Method:** Use `API` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pdb_get_structure(*, pdb_id, format='cif', max_bytes=10000000, verbose=True)`
- **Input:** Provide `pdb_id`.
- **Output:** Mapping with pdb_id, format, data (coordinate file text), url and sha256; does not write to disk.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pdb_get_structure(
    pdb_id='1CRN',
    format='cif',
    verbose=False,
)
print(result)
```

- **Example output:** Mapping with pdb_id, format, data (coordinate file text), url and sha256; does not write to disk.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "pdb_id": "<str>",
  "format": "<str>",
  "data": "<coordinate-file text>",
  "url": "<str>",
  "sha256": "<str>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.pdb_get_structure-parameters)

## 31) Map PDB chains to UniProt residue intervals {#dnakit.databases.pdbe_get_uniprot_mappings}

- **Purpose:** Map PDB chains to UniProt residue intervals.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pdbe_get_uniprot_mappings(*, pdb_id=None, uniprot_id=None, limit=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping keyed by record; each record includes `UniProt` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pdbe_get_uniprot_mappings(
    pdb_id='4hhb',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "<record key>": {
    "UniProt": {
      "<record key>": "<object>"
    }
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.pdbe_get_uniprot_mappings-parameters)

## 32) Rank protein structures by coverage and resolution {#dnakit.databases.pdbe_sifts_get_best_structures}

- **Purpose:** Rank protein structures by coverage and resolution.
- **Method:** Use `PDBeSIFTS` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pdbe_sifts_get_best_structures(*, uniprot_accession, limit=None, offset=None, verbose=True)`
- **Input:** Provide `uniprot_accession`.
- **Output:** Returns a mapping with `uniprot_accession`, `structures`, `total_structures`, `returned_structures` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pdbe_sifts_get_best_structures(
    uniprot_accession='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "uniprot_accession": "<string>",
  "structures": [
    {
      "pdb_id": "<string>",
      "chain_id": "<string>",
      "uniprot_start": "<integer>",
      "uniprot_end": "<integer>",
      "resolution": "<number>",
      "experimental_method": "<string>",
      "coverage": "<number>",
      "tax_id": "<integer>"
    }
  ],
  "total_structures": "<integer>",
  "returned_structures": "<integer>",
  "offset": "<integer>",
  "note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.pdbe_sifts_get_best_structures-parameters)

## 33) Retrieve Pfam protein domains and residue positions {#dnakit.databases.pfam_get_protein_annotations}

- **Purpose:** Retrieve Pfam protein domains and residue positions.
- **Method:** Use `Pfam` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.pfam_get_protein_annotations(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a mapping with `accession`, `protein_length`, `domain_count`, `domains` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.pfam_get_protein_annotations(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "accession": "<string>",
  "protein_length": "<integer>",
  "domain_count": "<integer>",
  "domains": [
    {
      "pfam_accession": "<string>",
      "name": "<string>",
      "short_name": "<string>",
      "type": "<string>",
      "integrated_interpro": "<string>",
      "start": "<integer>",
      "end": "<integer>",
      "score": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.pfam_get_protein_annotations-parameters)

## 34) Retrieve plant pathway reactions and annotations {#dnakit.databases.plant_reactome_get_pathway}

- **Purpose:** Retrieve plant pathway reactions and annotations.
- **Method:** Use `PlantReactome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.plant_reactome_get_pathway(*, pathway_id, verbose=True)`
- **Input:** Provide `pathway_id`.
- **Output:** Returns a mapping with `stId`, `displayName`, `speciesName`, `isInDisease` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.plant_reactome_get_pathway(
    pathway_id='R-OSA-1119312',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "stId": "<string>",
  "displayName": "<string>",
  "speciesName": "<string>",
  "isInDisease": "<boolean>",
  "isInferred": "<boolean>",
  "hasEvent": [
    "<value>"
  ],
  "compartment": [
    "<value>"
  ],
  "goBiologicalProcess": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.plant_reactome_get_pathway-parameters)

## 35) List species covered by Plant Reactome {#dnakit.databases.plant_reactome_list_species}

- **Purpose:** List species covered by Plant Reactome.
- **Method:** Use `PlantReactome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.plant_reactome_list_species(*, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `dbId`, `displayName`, `taxId`, `abbreviation` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.plant_reactome_list_species(
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "dbId": "<integer>",
    "displayName": "<string>",
    "taxId": "<string>",
    "abbreviation": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.plant_reactome_list_species-parameters)

## 36) Search plant metabolic and regulatory pathways {#dnakit.databases.plant_reactome_search_pathways}

- **Purpose:** Search plant metabolic and regulatory pathways.
- **Method:** Use `PlantReactome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.plant_reactome_search_pathways(*, query, species=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a list of records; each record includes `stId`, `name`, `species`, `exact_type` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.plant_reactome_search_pathways(
    query='photosynthesis',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "stId": "<string>",
    "name": "<string>",
    "species": "<string>",
    "exact_type": "<string>",
    "compartment_names": [
      "<value>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.plant_reactome_search_pathways-parameters)

## 37) Retrieve protein sequence features and modified residues {#dnakit.databases.proteins_api_get_features}

- **Purpose:** Retrieve protein sequence features and modified residues.
- **Method:** Use `ProteinsAPI` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.proteins_api_get_features(*, accession, format=None, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a list for “Retrieve protein sequence features and modified residues”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.proteins_api_get_features(
    accession='P05067',
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.proteins_api_get_features-parameters)

## 38) Search proteins by gene, name or accession {#dnakit.databases.proteins_api_search}

- **Purpose:** Search proteins by gene, name or accession.
- **Method:** Use `ProteinsAPI` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.proteins_api_search(*, query=None, size=None, offset=None, format=None, gene_symbol=None, gene=None, reviewed=None, organism=None, taxid=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list for “Search proteins by gene, name or accession”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.proteins_api_search(
    query='BRCA1',
    size=5,
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.proteins_api_search-parameters)

## 39) Retrieve GO annotations and evidence for a gene product {#dnakit.databases.quick_go_annotations_by_gene}

- **Purpose:** Retrieve GO annotations and evidence for a gene product.
- **Method:** Use `QuickGO` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.quick_go_annotations_by_gene(*, gene_product_id, aspect=None, taxon_id=None, limit=None, verbose=True)`
- **Input:** Provide `gene_product_id`.
- **Output:** Returns a list of records; each record includes `go_id`, `go_name`, `go_aspect`, `gene_product_id` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.quick_go_annotations_by_gene(
    gene_product_id='UniProtKB:P04637',
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "go_id": "<string>",
    "go_name": "<string>",
    "go_aspect": "<string>",
    "gene_product_id": "<string>",
    "symbol": "<string>",
    "qualifier": "<string>",
    "go_evidence": "<string>",
    "assigned_by": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.quick_go_annotations_by_gene-parameters)

## 40) Retrieve PDB experimental methods and structure metadata {#dnakit.databases.rcsbdata_get_entry}

- **Purpose:** Retrieve PDB experimental methods and structure metadata.
- **Method:** Use `RCSBData` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.rcsbdata_get_entry(*, pdb_id=None, entry_id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Select at least one argument group: `pdb_id` / `entry_id`.
- **Output:** Returns a mapping with `pdb_id`, `title`, `method`, `resolution` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.rcsbdata_get_entry(
    pdb_id='4HHB',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "pdb_id": "<string>",
  "title": "<string>",
  "method": "<string>",
  "resolution": "<number>",
  "deposit_date": "<string>",
  "release_date": "<string>",
  "revision_date": "<string>",
  "polymer_entity_count": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.rcsbdata_get_entry-parameters)

## 41) Search Reactome pathways, reactions and entities {#dnakit.databases.reactome_content_search}

- **Purpose:** Search Reactome pathways, reactions and entities.
- **Method:** Use `ReactomeContent` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.reactome_content_search(*, query, species=None, types=None, rows=None, start=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `query`, `species`, `count`, `total_results` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.reactome_content_search(
    query='apoptosis',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "query": "<string>",
  "species": "<string>",
  "count": "<integer>",
  "total_results": "<integer>",
  "start": "<integer>",
  "has_more": "<boolean>",
  "results": [
    {
      "type": "<string>",
      "stId": "<string>",
      "name": "<string>",
      "species": "<array>",
      "is_disease": "<boolean>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_content_search-parameters)

## 42) Retrieve molecules and complexes participating in a reaction {#dnakit.databases.reactome_get_participants}

- **Purpose:** Retrieve molecules and complexes participating in a reaction.
- **Method:** Use `Reactome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.reactome_get_participants(*, stId, verbose=True)`
- **Input:** Provide `stId`.
- **Output:** Returns a list of records; each record includes `peDbId`, `displayName`, `schemaClass`, `refEntities` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.reactome_get_participants(
    stId='R-HSA-111289',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "peDbId": "<integer>",
    "displayName": "<string>",
    "schemaClass": "<string>",
    "refEntities": [
      "<object>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_get_participants-parameters)

## 43) Retrieve Reactome pathway details {#dnakit.databases.reactome_get_pathway}

- **Purpose:** Retrieve Reactome pathway details.
- **Method:** Use `Reactome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.reactome_get_pathway(*, stId, verbose=True)`
- **Input:** Provide `stId`.
- **Output:** Returns a mapping with `dbId`, `displayName`, `stId`, `stIdVersion` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.reactome_get_pathway(
    stId='R-HSA-73817',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "dbId": "<integer>",
  "displayName": "<string>",
  "stId": "<string>",
  "stIdVersion": "<string>",
  "isInDisease": "<boolean>",
  "isInferred": "<boolean>",
  "maxDepth": "<integer>",
  "name": [
    "<string>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_get_pathway-parameters)

## 44) List reactions and subpathways within a pathway {#dnakit.databases.reactome_get_pathway_reactions}

- **Purpose:** List reactions and subpathways within a pathway.
- **Method:** Use `Reactome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.reactome_get_pathway_reactions(*, stId, verbose=True)`
- **Input:** Provide `stId`.
- **Output:** Returns a list for “List reactions and subpathways within a pathway”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.reactome_get_pathway_reactions(
    stId='R-HSA-73817',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  "<integer>"
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_get_pathway_reactions-parameters)

## 45) Map a UniProt protein to Reactome pathways {#dnakit.databases.reactome_map_uniprot_to_pathways}

- **Purpose:** Map a UniProt protein to Reactome pathways.
- **Method:** Use `Reactome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.reactome_map_uniprot_to_pathways(*, uniprot_id=None, id=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `dbId`, `displayName`, `stId`, `stIdVersion` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.reactome_map_uniprot_to_pathways(
    uniprot_id='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "dbId": "<integer>",
    "displayName": "<string>",
    "stId": "<string>",
    "stIdVersion": "<string>",
    "isInDisease": "<boolean>",
    "isInferred": "<boolean>",
    "maxDepth": "<integer>",
    "name": [
      "<string>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.reactome_map_uniprot_to_pathways-parameters)

## 46) Retrieve STRING protein functional annotations {#dnakit.databases.string_get_functional_annotations}

- **Purpose:** Retrieve STRING protein functional annotations.
- **Method:** Use `STRINGExt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.string_get_functional_annotations(*, identifiers, species=None, category=None, verbose=True)`
- **Input:** Provide `identifiers`.
- **Output:** Returns a mapping with `identifiers`, `species`, `total_annotations`, `category_summary` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.string_get_functional_annotations(
    identifiers='TP53',
    species=9606,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "identifiers": "<string>",
  "species": "<integer>",
  "total_annotations": "<integer>",
  "category_summary": "<object>",
  "annotations": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_get_functional_annotations-parameters)

## 47) Retrieve protein interaction partners and confidence scores {#dnakit.databases.string_get_interaction_partners}

- **Purpose:** Retrieve protein interaction partners and confidence scores.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.string_get_interaction_partners(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `stringId_A`, `stringId_B`, `preferredName_A`, `preferredName_B` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.string_get_interaction_partners(
    identifiers='TP53',
    species=9606,
    limit=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "stringId_A": "<string>",
    "stringId_B": "<string>",
    "preferredName_A": "<string>",
    "preferredName_B": "<string>",
    "ncbiTaxonId": "<string>",
    "score": "<number>",
    "nscore": "<number>",
    "fscore": "<number>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_get_interaction_partners-parameters)

## 48) Retrieve an interaction network for proteins and their neighbors {#dnakit.databases.string_get_network}

- **Purpose:** Retrieve an interaction network for proteins and their neighbors.
- **Method:** Use `Base` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.string_get_network(*, identifiers=None, species=None, limit=None, required_score=None, protein_ids=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a list of records; each record includes `stringId_A`, `stringId_B`, `preferredName_A`, `preferredName_B` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.string_get_network(
    identifiers='TP53',
    species=9606,
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
    "stringId_A": "<string>",
    "stringId_B": "<string>",
    "preferredName_A": "<string>",
    "preferredName_B": "<string>",
    "ncbiTaxonId": "<string>",
    "score": "<number>",
    "nscore": "<number>",
    "fscore": "<number>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_get_network-parameters)

## 49) Test protein interaction enrichment in a protein set {#dnakit.databases.string_ppi_enrichment}

- **Purpose:** Test protein interaction enrichment in a protein set.
- **Method:** Use `STRING` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.string_ppi_enrichment(*, protein_ids, species=None, confidence_score=None, verbose=True)`
- **Input:** Provide `protein_ids`.
- **Output:** Returns a mapping with `number_of_nodes`, `number_of_edges`, `expected_number_of_edges`, `p_value` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.string_ppi_enrichment(
    protein_ids=['TP53', 'MDM2', 'ATM', 'CHEK2', 'BRCA1', 'BRCA2'],
    species=9606,
    confidence_score=0.4,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "number_of_nodes": "<integer>",
  "number_of_edges": "<integer>",
  "expected_number_of_edges": "<number>",
  "p_value": "<number>",
  "avg_node_degree": "<number>",
  "local_clustering_coefficient": "<number>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.string_ppi_enrichment-parameters)

## 50) Annotate residue interfaces, pockets and solvent exposure {#dnakit.databases.structure_annotate_per_residue}

- **Purpose:** Annotate residue interfaces, pockets and solvent exposure.
- **Method:** Use `StructureAnnotation` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.structure_annotate_per_residue(*, operation=None, pdb_id=None, pdb_content=None, target_chain=None, partner_chains=None, ligand_resnames=None, distance_cutoff=None, core_rsa_cutoff=None, include_secondary_structure=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults.
- **Output:** Returns a mapping with `status`, `pdb_id`, `target_chain`, `partner_chains` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.structure_annotate_per_residue(
    operation='annotate_per_residue',
    pdb_id='1CRN',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "pdb_id": "<string>",
  "target_chain": "<string>",
  "partner_chains": [
    "<string>"
  ],
  "ligand_resnames": [
    "<string>"
  ],
  "distance_cutoff": "<number>",
  "core_rsa_cutoff": "<number>",
  "n_residues": "<integer>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.structure_annotate_per_residue-parameters)

## 51) Extract UniProt protein variants and annotations {#dnakit.databases.uni_prot_get_disease_variants_by_accession}

- **Purpose:** Extract UniProt protein variants and annotations.
- **Method:** Use `UniProt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_prot_get_disease_variants_by_accession(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a list of records; each record includes `type`, `location`, `description`, `featureId` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_prot_get_disease_variants_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "type": "<string>",
    "location": {
      "start": "<object>",
      "end": "<object>"
    },
    "description": "<string>",
    "featureId": "<string>",
    "alternativeSequence": {
      "originalSequence": "<string>",
      "alternativeSequences": "<array>"
    },
    "featureCrossReferences": [
      "<object>"
    ],
    "evidences": [
      "<object>"
    ]
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_disease_variants_by_accession-parameters)

## 52) Retrieve a complete UniProt entry by accession {#dnakit.databases.uni_prot_get_entry_by_accession}

- **Purpose:** Retrieve a complete UniProt entry by accession.
- **Method:** Use `UniProt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_prot_get_entry_by_accession(*, accession, compact=None, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a mapping with `entryType`, `primaryAccession`, `uniProtkbId`, `protein_name` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_prot_get_entry_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "entryType": "<string>",
  "primaryAccession": "<string>",
  "uniProtkbId": "<string>",
  "protein_name": "<string>",
  "gene_names": [
    "<string>"
  ],
  "organism": {
    "scientificName": "<string>",
    "commonName": "<string>",
    "taxonId": "<integer>",
    "lineage": [
      "<string>"
    ]
  },
  "sequence": {
    "value": "<string>",
    "length": "<integer>",
    "molWeight": "<integer>",
    "crc64": "<string>",
    "md5": "<string>"
  },
  "annotationScore": "<number>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_entry_by_accession-parameters)

## 53) Extract protein domains, active sites and sequence features {#dnakit.databases.uni_prot_get_features_by_accession}

- **Purpose:** Extract protein domains, active sites and sequence features.
- **Method:** Use `UniProt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_prot_get_features_by_accession(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a list of records; each record includes `type`, `location`, `description` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_prot_get_features_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "type": "<string>",
    "location": "<object>",
    "description": "<string>"
  }
]
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_features_by_accession-parameters)

## 54) Extract UniProt protein function annotations {#dnakit.databases.uni_prot_get_function_by_accession}

- **Purpose:** Extract UniProt protein function annotations.
- **Method:** Use `UniProt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_prot_get_function_by_accession(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a list for “Extract UniProt protein function annotations”; the item shape is shown below.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_prot_get_function_by_accession(
    accession='P04637',
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

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_function_by_accession-parameters)

## 55) Extract a canonical UniProt protein sequence {#dnakit.databases.uni_prot_get_sequence_by_accession}

- **Purpose:** Extract a canonical UniProt protein sequence.
- **Method:** Use `UniProt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_prot_get_sequence_by_accession(*, accession, verbose=True)`
- **Input:** Provide `accession`.
- **Output:** Returns a string.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_prot_get_sequence_by_accession(
    accession='P04637',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
"<string>"
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_get_sequence_by_accession-parameters)

## 56) Search UniProt proteins with structured queries {#dnakit.databases.uni_prot_search}

- **Purpose:** Search UniProt proteins with structured queries.
- **Method:** Use `UniProt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.uni_prot_search(*, query, organism=None, limit=None, min_length=None, max_length=None, fields=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `total_results`, `returned`, `results`, `normalization_note` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.uni_prot_search(
    query='gene:TP53',
    organism='homo sapiens',
    limit=5,
    min_length=0,
    max_length=10000,
    fields=[],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "total_results": "<integer>",
  "returned": "<integer>",
  "results": [
    {
      "accession": "<string>",
      "id": "<string>",
      "protein_name": "<string>",
      "gene_names": "<array>",
      "organism": "<string>",
      "length": "<integer>"
    }
  ],
  "normalization_note": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.uni_prot_search-parameters)

## 57) Search WikiPathways by pathway title {#dnakit.databases.wiki_pathways_search}

- **Purpose:** Search WikiPathways by pathway title.
- **Method:** Use `WikiPathwaysSearch` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.wiki_pathways_search(*, query, organism=None, limit=None, verbose=True)`
- **Input:** Provide `query`.
- **Output:** Returns a mapping with `status`, `data`, `metadata` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.wiki_pathways_search(
    query='apoptosis',
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
  "status": "success",
  "data": [
    {
      "wpid": "<str>",
      "uri": "<str>",
      "title": "<str>",
      "organism": "<str>"
    }
  ],
  "metadata": {
    "query": "<str>",
    "count": "<int>",
    "backend": "WikiPathways SPARQL"
  }
}
```

[Full parameters and usage notes](../../faq/api/22_database_query_pathways.md#dnakit.databases.wiki_pathways_search-parameters)

<!-- END GENERATED API: pathways -->
