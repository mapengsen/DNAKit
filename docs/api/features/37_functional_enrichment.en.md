# Gene Set and Pathway Enrichment

See the [FAQ](../../faq/tutorial_integrations.md#optional-bioinformatics) for dependencies and return conventions.

<!-- BEGIN GENERATED API: enrichment -->

<span id="enrichment"></span>

## 1) Run gene set enrichment against an Enrichr library {#dnakit.databases.enrichr_enrich}

- **Purpose:** Run gene set enrichment against an Enrichr library.
- **Method:** Use `EnrichrExt` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.enrichr_enrich(*, operation=None, gene_list, library=None, top_n=None, species=None, verbose=True)`
- **Input:** Provide `gene_list`.
- **Output:** Returns enriched terms, p-values, z-scores, combined scores and overlapping genes.
- **Example code:**

```python
import dnakit

result = dnakit.databases.enrichr_enrich(
    gene_list=['TP53', 'BRCA1', 'EGFR', 'KRAS', 'MYC', 'BRAF', 'PTEN'],
    library='GO_Biological_Process_2023',
    top_n=5,
    verbose=False,
)
print(result)
```

- **Example output:** Returns enriched terms, p-values, z-scores, combined scores and overlapping genes.

[Full parameters and usage notes](../../faq/api/37_functional_enrichment.md#dnakit.databases.enrichr_enrich-parameters)

## 2) Retrieve connections between genes and enriched pathways {#dnakit.databases.enrichr_gene_enrichment_analysis}

- **Purpose:** Retrieve connections between genes and enriched pathways.
- **Method:** Use `Enrichr` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.enrichr_gene_enrichment_analysis(*, gene_list, libs=None, verbose=True)`
- **Input:** Provide `gene_list`.
- **Output:** Returns a mapping with `connected_paths`, `connections` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.enrichr_gene_enrichment_analysis(
    gene_list=['TP53', 'BRCA1', 'EGFR'],
    libs=['KEGG_2021_Human'],
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "connected_paths": "<object>",
  "connections": "<object>"
}
```

[Full parameters and usage notes](../../faq/api/37_functional_enrichment.md#dnakit.databases.enrichr_gene_enrichment_analysis-parameters)

## 3) Run PANTHER gene set overrepresentation analysis {#dnakit.databases.panther_enrichment}

- **Purpose:** Run PANTHER gene set overrepresentation analysis.
- **Method:** Use `PANTHER` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.panther_enrichment(*, gene_list, organism=None, annotation_dataset=None, verbose=True)`
- **Input:** Provide `gene_list`.
- **Output:** Returns a mapping with `gene_list`, `organism`, `annotation_dataset`, `result_count` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.panther_enrichment(
    gene_list='TP53,BRCA1,EGFR,KRAS,PTEN',
    organism=9606,
    annotation_dataset='GO:0008150',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "gene_list": "<string>",
  "organism": "<integer>",
  "annotation_dataset": "<string>",
  "result_count": "<integer>",
  "enriched_terms": [
    {
      "term_id": "<string>",
      "term_label": "<string>",
      "number_in_list": "<integer>",
      "number_in_reference": "<integer>",
      "expected": "<number>",
      "fold_enrichment": "<number>",
      "pvalue": "<number>",
      "fdr": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/37_functional_enrichment.md#dnakit.databases.panther_enrichment-parameters)

## 4) Run Reactome pathway enrichment analysis {#dnakit.databases.reactome_analysis_pathway_enrichment}

- **Purpose:** Run Reactome pathway enrichment analysis.
- **Method:** Use `ReactomeAnalysis` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.reactome_analysis_pathway_enrichment(*, identifiers, page_size=None, include_disease=None, projection=None, verbose=True)`
- **Input:** Provide `identifiers`.
- **Output:** Returns a mapping with `token`, `analysis_type`, `projection`, `identifiers_not_found` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.reactome_analysis_pathway_enrichment(
    identifiers='TP53\nBRCA1\nCDH1\nEGFR\nKRAS',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "token": "<string>",
  "analysis_type": "<string>",
  "projection": "<boolean>",
  "identifiers_not_found": "<integer>",
  "pathways_found": "<integer>",
  "pathways": [
    {
      "pathway_id": "<string>",
      "name": "<string>",
      "species": "<string>",
      "is_disease": "<boolean>",
      "is_lowest_level": "<boolean>",
      "entities_found": "<integer>",
      "entities_total": "<integer>",
      "entities_coverage": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/37_functional_enrichment.md#dnakit.databases.reactome_analysis_pathway_enrichment-parameters)

## 5) Analyze functional and pathway enrichment for a protein set {#dnakit.databases.string_functional_enrichment}

- **Purpose:** Analyze functional and pathway enrichment for a protein set.
- **Method:** Use `STRING` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.string_functional_enrichment(*, protein_ids, species=None, category=None, verbose=True)`
- **Input:** Provide `protein_ids`.
- **Output:** Returns a list of records; each record includes `category`, `term`, `number_of_genes`, `number_of_genes_in_background` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.string_functional_enrichment(
    protein_ids=['TP53', 'MDM2', 'ATM', 'CHEK2', 'BRCA1'],
    species=9606,
    category='Process',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
[
  {
    "category": "<string>",
    "term": "<string>",
    "number_of_genes": "<integer>",
    "number_of_genes_in_background": "<integer>",
    "ncbiTaxonId": "<integer>",
    "inputGenes": "<string>",
    "preferredNames": "<string>",
    "p_value": "<number>"
  }
]
```

[Full parameters and usage notes](../../faq/api/37_functional_enrichment.md#dnakit.databases.string_functional_enrichment-parameters)

## 6) Run STRING enrichment from newline-separated protein identifiers {#dnakit.databases.string_get_functional_enrichment}

- **Purpose:** Run STRING enrichment from newline-separated protein identifiers.
- **Method:** Use `API` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.string_get_functional_enrichment(*, identifiers, species=9606, category=None, verbose=True)`
- **Input:** Provide `identifiers`.
- **Output:** Returns STRING_functional_enrichment results; see [string_functional_enrichment](37_functional_enrichment.md#dnakit.databases.string_functional_enrichment).
- **Example code:**

```python
import dnakit

result = dnakit.databases.string_get_functional_enrichment(
    identifiers='TP53\nBRCA1\nBRCA2',
    species=9606,
    verbose=False,
)
print(result)
```

- **Example output:** Returns STRING_functional_enrichment results; see [string_functional_enrichment](37_functional_enrichment.md#dnakit.databases.string_functional_enrichment).

[Full parameters and usage notes](../../faq/api/37_functional_enrichment.md#dnakit.databases.string_get_functional_enrichment-parameters)

<!-- END GENERATED API: enrichment -->
