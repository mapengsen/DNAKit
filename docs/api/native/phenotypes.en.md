# Diseases, phenotypes and ontologies

This page documents **45 newly added public functions**. [Complete index](index.md) · [Conventions and examples](usage.md)

Each function lists its current Python signature, parameters, returns and registered execution requirements. Provider descriptions retain source wording. In generated database wrappers, None normally omits a parameter; actual required fields and conditional constraints come from the input schema. Follow each operation's coordinate and assembly conventions.

## Function index

- [`dnakit.databases.dis_ge_net_get_disease_genes`](#dnakit.databases.dis_ge_net_get_disease_genes)
- [`dnakit.databases.dis_ge_net_get_gda`](#dnakit.databases.dis_ge_net_get_gda)
- [`dnakit.databases.dis_ge_net_get_vda`](#dnakit.databases.dis_ge_net_get_vda)
- [`dnakit.databases.dis_ge_net_search_disease`](#dnakit.databases.dis_ge_net_search_disease)
- [`dnakit.databases.dis_ge_net_search_gene`](#dnakit.databases.dis_ge_net_search_gene)
- [`dnakit.databases.gen_cc_get_classifications`](#dnakit.databases.gen_cc_get_classifications)
- [`dnakit.databases.gen_cc_search_disease`](#dnakit.databases.gen_cc_search_disease)
- [`dnakit.databases.gen_cc_search_gene`](#dnakit.databases.gen_cc_search_gene)
- [`dnakit.databases.get_joint_associated_diseases_by_hpo_id_list`](#dnakit.databases.get_joint_associated_diseases_by_hpo_id_list)
- [`dnakit.databases.hpo_get_diseases_by_phenotype`](#dnakit.databases.hpo_get_diseases_by_phenotype)
- [`dnakit.databases.hpo_get_genes_by_phenotype`](#dnakit.databases.hpo_get_genes_by_phenotype)
- [`dnakit.databases.hpo_get_term`](#dnakit.databases.hpo_get_term)
- [`dnakit.databases.hpo_search_terms`](#dnakit.databases.hpo_search_terms)
- [`dnakit.databases.marrvel_get_gene`](#dnakit.databases.marrvel_get_gene)
- [`dnakit.databases.marrvel_get_omim_phenotypes`](#dnakit.databases.marrvel_get_omim_phenotypes)
- [`dnakit.databases.monarch_get_gene_diseases`](#dnakit.databases.monarch_get_gene_diseases)
- [`dnakit.databases.monarch_get_gene_phenotypes`](#dnakit.databases.monarch_get_gene_phenotypes)
- [`dnakit.databases.monarch_search`](#dnakit.databases.monarch_search)
- [`dnakit.databases.monarch_search_gene`](#dnakit.databases.monarch_search_gene)
- [`dnakit.databases.monarch_v3_get_associations`](#dnakit.databases.monarch_v3_get_associations)
- [`dnakit.databases.monarch_v3_get_entity`](#dnakit.databases.monarch_v3_get_entity)
- [`dnakit.databases.monarch_v3_get_histopheno`](#dnakit.databases.monarch_v3_get_histopheno)
- [`dnakit.databases.monarch_v3_get_mappings`](#dnakit.databases.monarch_v3_get_mappings)
- [`dnakit.databases.monarch_v3_phenotype_similarity_search`](#dnakit.databases.monarch_v3_phenotype_similarity_search)
- [`dnakit.databases.monarch_v3_search`](#dnakit.databases.monarch_v3_search)
- [`dnakit.databases.ols_get_term_ancestors`](#dnakit.databases.ols_get_term_ancestors)
- [`dnakit.databases.ols_get_term_children`](#dnakit.databases.ols_get_term_children)
- [`dnakit.databases.ols_get_term_info`](#dnakit.databases.ols_get_term_info)
- [`dnakit.databases.ols_search_terms`](#dnakit.databases.ols_search_terms)
- [`dnakit.databases.omim_get_clinical_synopsis`](#dnakit.databases.omim_get_clinical_synopsis)
- [`dnakit.databases.omim_get_entry`](#dnakit.databases.omim_get_entry)
- [`dnakit.databases.omim_get_gene_map`](#dnakit.databases.omim_get_gene_map)
- [`dnakit.databases.omim_search`](#dnakit.databases.omim_search)
- [`dnakit.databases.orphanet_get_classification`](#dnakit.databases.orphanet_get_classification)
- [`dnakit.databases.orphanet_get_disease`](#dnakit.databases.orphanet_get_disease)
- [`dnakit.databases.orphanet_get_disease_genes`](#dnakit.databases.orphanet_get_disease_genes)
- [`dnakit.databases.orphanet_get_epidemiology`](#dnakit.databases.orphanet_get_epidemiology)
- [`dnakit.databases.orphanet_get_gene_diseases`](#dnakit.databases.orphanet_get_gene_diseases)
- [`dnakit.databases.orphanet_get_genes`](#dnakit.databases.orphanet_get_genes)
- [`dnakit.databases.orphanet_get_icd_mapping`](#dnakit.databases.orphanet_get_icd_mapping)
- [`dnakit.databases.orphanet_get_natural_history`](#dnakit.databases.orphanet_get_natural_history)
- [`dnakit.databases.orphanet_get_phenotypes`](#dnakit.databases.orphanet_get_phenotypes)
- [`dnakit.databases.orphanet_search_by_name`](#dnakit.databases.orphanet_search_by_name)
- [`dnakit.databases.orphanet_search_diseases`](#dnakit.databases.orphanet_search_diseases)
- [`dnakit.databases.osl_get_efo_id_by_disease_name`](#dnakit.databases.osl_get_efo_id_by_disease_name)

## `dnakit.databases.dis_ge_net_get_disease_genes` {#dnakit.databases.dis_ge_net_get_disease_genes}

Get all genes associated with a disease from DisGeNET. Returns ranked gene list with association scores and evidence counts. Filter by minimum score to get high-confidence associations. Useful for disease gene prioritization.

```text
dnakit.databases.dis_ge_net_get_disease_genes(*, operation=None, disease, min_score=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_disease_genes) {&quot;const&quot;: &quot;get_disease_genes&quot;, &quot;default&quot;: &quot;get_disease_genes&quot;} |
| `disease` | string | Required by input schema | Disease ID (UMLS CUI) or disease name |
| `min_score` | number | None | Minimum association score (0-1). Default: no filter {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1} |
| `limit` | integer | None | Maximum genes to return (default: 50) {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_disease_genes",
          "description": "Operation type (fixed: get_disease_genes)",
          "default": "get_disease_genes"
        },
        "disease": {
          "type": "string",
          "description": "Disease ID (UMLS CUI) or disease name"
        },
        "min_score": {
          "type": "number",
          "description": "Minimum association score (0-1). Default: no filter",
          "minimum": 0,
          "maximum": 1
        },
        "limit": {
          "type": "integer",
          "description": "Maximum genes to return (default: 50)",
          "default": 50
        }
      },
      "required": [
        "disease"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "string"
            },
            "genes": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "symbol": {
                    "type": "string"
                  },
                  "score": {
                    "type": "number"
                  },
                  "evidence_count": {
                    "type": "integer"
                  }
                }
              }
            },
            "gene_count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DisGeNET_get_disease_genes`。

**Provider credentials:** `DISGENET_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dis_ge_net_get_disease_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dis_ge_net_get_gda` {#dnakit.databases.dis_ge_net_get_gda}

Get gene-disease associations with filtering. Filter by data source (CURATED, LITERATURE, ANIMAL_MODELS) and minimum score. Returns comprehensive GDA data with evidence levels. Requires DISGENET_API_KEY.

```text
dnakit.databases.dis_ge_net_get_gda(*, operation=None, gene=None, disease=None, source=None, min_score=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_gda) {&quot;const&quot;: &quot;get_gda&quot;, &quot;default&quot;: &quot;get_gda&quot;} |
| `gene` | string | None | Gene symbol (optional if disease provided) |
| `disease` | string | None | Disease ID/UMLS CUI (optional if gene provided) |
| `source` | string | None | Data source filter: CURATED, ANIMAL_MODELS, LITERATURE, etc. {&quot;enum&quot;: &#91;&quot;CURATED&quot;, &quot;ANIMAL_MODELS&quot;, &quot;LITERATURE&quot;, &quot;ALL&quot;&#93;} |
| `min_score` | number | None | Minimum GDA score (0-1). Higher = stronger evidence. {&quot;minimum&quot;: 0, &quot;maximum&quot;: 1} |
| `limit` | integer | None | Maximum results (default: 25) {&quot;default&quot;: 25} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_gda",
          "description": "Operation type (fixed: get_gda)",
          "default": "get_gda"
        },
        "gene": {
          "type": "string",
          "description": "Gene symbol (optional if disease provided)"
        },
        "disease": {
          "type": "string",
          "description": "Disease ID/UMLS CUI (optional if gene provided)"
        },
        "source": {
          "type": "string",
          "description": "Data source filter: CURATED, ANIMAL_MODELS, LITERATURE, etc.",
          "enum": [
            "CURATED",
            "ANIMAL_MODELS",
            "LITERATURE",
            "ALL"
          ]
        },
        "min_score": {
          "type": "number",
          "description": "Minimum GDA score (0-1). Higher = stronger evidence.",
          "minimum": 0,
          "maximum": 1
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results (default: 25)",
          "default": 25
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "associations": {
              "type": "array"
            },
            "count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DisGeNET_get_gda`。

**Provider credentials:** `DISGENET_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dis_ge_net_get_gda", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dis_ge_net_get_vda` {#dnakit.databases.dis_ge_net_get_vda}

Get variant-disease associations from DisGeNET. Query by rsID or gene symbol to find disease-associated variants. Returns variant IDs, diseases, and evidence scores. Useful for pharmacogenomics and variant interpretation.

```text
dnakit.databases.dis_ge_net_get_vda(*, operation=None, variant=None, gene=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_vda) {&quot;const&quot;: &quot;get_vda&quot;, &quot;default&quot;: &quot;get_vda&quot;} |
| `variant` | string | None | Variant rsID (e.g., rs1234) - optional if gene provided |
| `gene` | string | None | Gene symbol to get all variants - optional if variant provided |
| `limit` | integer | None | Maximum results (default: 25) {&quot;default&quot;: 25} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_vda",
          "description": "Operation type (fixed: get_vda)",
          "default": "get_vda"
        },
        "variant": {
          "type": "string",
          "description": "Variant rsID (e.g., rs1234) - optional if gene provided"
        },
        "gene": {
          "type": "string",
          "description": "Gene symbol to get all variants - optional if variant provided"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results (default: 25)",
          "default": 25
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "associations": {
              "type": "array"
            },
            "count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DisGeNET_get_vda`。

**Provider credentials:** `DISGENET_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dis_ge_net_get_vda", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dis_ge_net_search_disease` {#dnakit.databases.dis_ge_net_search_disease}

Search DisGeNET for genes associated with a disease. Requires a UMLS CUI, NOT a free-text disease name (e.g. 'C0006142' for breast cancer, not 'breast cancer') -- resolve a disease name to its CUI first via umls_search_concepts. Returns associated genes with scores and evidence. Requires DISGENET_API_KEY.

```text
dnakit.databases.dis_ge_net_search_disease(*, operation=None, disease, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_disease) {&quot;const&quot;: &quot;search_disease&quot;, &quot;default&quot;: &quot;search_disease&quot;} |
| `disease` | string | Required by input schema | UMLS CUI for the disease, e.g. &#x27;C0006142&#x27; (breast cancer). Free-text disease names are NOT accepted -- resolve one via umls_search_concepts first. |
| `limit` | integer | None | Maximum associations to return (default: 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search_disease",
          "description": "Operation type (fixed: search_disease)",
          "default": "search_disease"
        },
        "disease": {
          "type": "string",
          "description": "UMLS CUI for the disease, e.g. 'C0006142' (breast cancer). Free-text disease names are NOT accepted -- resolve one via umls_search_concepts first."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum associations to return (default: 10)",
          "default": 10
        }
      },
      "required": [
        "disease"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "string"
            },
            "associations": {
              "type": "array"
            },
            "count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DisGeNET_search_disease`。

**Provider credentials:** `DISGENET_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dis_ge_net_search_disease", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.dis_ge_net_search_gene` {#dnakit.databases.dis_ge_net_search_gene}

Search DisGeNET for diseases associated with a gene. DisGeNET is the largest public database of gene-disease associations. Returns diseases, association scores, and evidence counts. Use the 'gene' parameter (e.g., 'BRCA1', 'TP53') — note: parameter name is 'gene', not 'gene_symbol'. Requires DISGENET_API_KEY env var (register free at disgenet.org).

```text
dnakit.databases.dis_ge_net_search_gene(*, operation=None, gene, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_gene) {&quot;const&quot;: &quot;search_gene&quot;, &quot;default&quot;: &quot;search_gene&quot;} |
| `gene` | string | Required by input schema | Gene symbol (e.g., BRCA1, TP53, EGFR) |
| `limit` | integer | None | Maximum associations to return (default: 10) {&quot;default&quot;: 10} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search_gene",
          "description": "Operation type (fixed: search_gene)",
          "default": "search_gene"
        },
        "gene": {
          "type": "string",
          "description": "Gene symbol (e.g., BRCA1, TP53, EGFR)"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum associations to return (default: 10)",
          "default": 10
        }
      },
      "required": [
        "gene"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "gene": {
              "type": "string"
            },
            "associations": {
              "type": "array",
              "description": "Gene-disease associations"
            },
            "count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `DisGeNET_search_gene`。

**Provider credentials:** `DISGENET_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("dis_ge_net_search_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gen_cc_get_classifications` {#dnakit.databases.gen_cc_get_classifications}

Get summary statistics of gene-disease validity classifications from GenCC. Returns counts per classification level (Definitive, Strong, Moderate, Limited, etc.). Optionally filter by gene_symbol (e.g., 'GBA'), disease name substring, or submitter organization. Without filters, returns global stats across all 40,000+ curated submissions.

```text
dnakit.databases.gen_cc_get_classifications(*, operation=None, gene_symbol=None, disease=None, submitter=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_classifications) {&quot;const&quot;: &quot;get_classifications&quot;, &quot;default&quot;: &quot;get_classifications&quot;} |
| `gene_symbol` | string | None | Optional HGNC gene symbol to filter classification stats (e.g., &#x27;GBA&#x27;, &#x27;BRCA2&#x27;) |
| `disease` | string | None | Optional disease name substring to filter classification stats (e.g., &#x27;Gaucher&#x27;, &#x27;breast cancer&#x27;) |
| `submitter` | string | None | Optional filter by submitting organization name (e.g., &#x27;ClinGen&#x27;, &#x27;Ambry&#x27;, &#x27;Genomics England&#x27;) {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_classifications",
          "description": "Operation type (fixed: get_classifications)",
          "default": "get_classifications"
        },
        "gene_symbol": {
          "type": "string",
          "description": "Optional HGNC gene symbol to filter classification stats (e.g., 'GBA', 'BRCA2')"
        },
        "disease": {
          "type": "string",
          "description": "Optional disease name substring to filter classification stats (e.g., 'Gaucher', 'breast cancer')"
        },
        "submitter": {
          "type": "string",
          "description": "Optional filter by submitting organization name (e.g., 'ClinGen', 'Ambry', 'Genomics England')",
          "default": ""
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "classifications": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "classification": {
                    "type": "string"
                  },
                  "count": {
                    "type": "integer"
                  },
                  "rank": {
                    "type": "integer"
                  }
                }
              }
            },
            "total_submissions": {
              "type": "integer"
            },
            "unique_genes": {
              "type": "integer"
            },
            "unique_diseases": {
              "type": "integer"
            },
            "top_submitters": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "name": {
                    "type": "string"
                  },
                  "count": {
                    "type": "integer"
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `GenCC_get_classifications`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("gen_cc_get_classifications", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gen_cc_search_disease` {#dnakit.databases.gen_cc_search_disease}

Find genes with validity evidence for a disease from GenCC (Gene Curation Coalition). Search by disease name or MONDO/OMIM identifier. Returns gene-disease validity classifications from multiple expert curators. Useful for identifying established gene-disease relationships for clinical interpretation.

```text
dnakit.databases.gen_cc_search_disease(*, operation=None, disease=None, disease_name=None, classification=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_disease) {&quot;const&quot;: &quot;search_disease&quot;, &quot;default&quot;: &quot;search_disease&quot;} |
| `disease` | string | None | Disease name or ID to search (e.g., &#x27;Marfan syndrome&#x27;, &#x27;MONDO:0007947&#x27;, &#x27;breast cancer&#x27;) |
| `disease_name` | string | None | Alias for disease. Disease name or ID to search. |
| `classification` | string | None | Optional filter by classification level (e.g., Definitive, Strong) {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search_disease",
          "description": "Operation type (fixed: search_disease)",
          "default": "search_disease"
        },
        "disease": {
          "type": "string",
          "description": "Disease name or ID to search (e.g., 'Marfan syndrome', 'MONDO:0007947', 'breast cancer')"
        },
        "disease_name": {
          "type": "string",
          "description": "Alias for disease. Disease name or ID to search."
        },
        "classification": {
          "type": "string",
          "description": "Optional filter by classification level (e.g., Definitive, Strong)",
          "default": ""
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "disease": {
              "type": "string"
            },
            "submissions": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "gene_symbol": {
                    "type": "string"
                  },
                  "gene_curie": {
                    "type": "string"
                  },
                  "disease_title": {
                    "type": "string"
                  },
                  "disease_curie": {
                    "type": "string"
                  },
                  "classification": {
                    "type": "string"
                  },
                  "mode_of_inheritance": {
                    "type": "string"
                  },
                  "submitter": {
                    "type": "string"
                  },
                  "submitted_date": {
                    "type": "string"
                  }
                }
              }
            },
            "submission_count": {
              "type": "integer"
            },
            "unique_genes": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `GenCC_search_disease`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("gen_cc_search_disease", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.gen_cc_search_gene` {#dnakit.databases.gen_cc_search_gene}

Get gene-disease validity classifications for a gene from GenCC (Gene Curation Coalition). Returns curated evidence levels (Definitive, Strong, Moderate, Limited, Disputed, Refuted) from multiple expert panels including ClinGen, Ambry Genetics, Genomics England PanelApp, and Orphanet. Essential for assessing clinical validity of gene-disease associations in genetic testing.

```text
dnakit.databases.gen_cc_search_gene(*, operation=None, gene_symbol, classification=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_gene) {&quot;const&quot;: &quot;search_gene&quot;, &quot;default&quot;: &quot;search_gene&quot;} |
| `gene_symbol` | string | Required by input schema | HGNC gene symbol (e.g., BRCA2, TP53, FBN1, CFTR) |
| `classification` | string | None | Optional filter by classification level (e.g., Definitive, Strong, Moderate, Limited) {&quot;default&quot;: &quot;&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search_gene",
          "description": "Operation type (fixed: search_gene)",
          "default": "search_gene"
        },
        "gene_symbol": {
          "type": "string",
          "description": "HGNC gene symbol (e.g., BRCA2, TP53, FBN1, CFTR)"
        },
        "classification": {
          "type": "string",
          "description": "Optional filter by classification level (e.g., Definitive, Strong, Moderate, Limited)",
          "default": ""
        }
      },
      "required": [
        "gene_symbol"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "gene_symbol": {
              "type": "string"
            },
            "submissions": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "gene_symbol": {
                    "type": "string"
                  },
                  "gene_curie": {
                    "type": "string"
                  },
                  "disease_title": {
                    "type": "string"
                  },
                  "disease_curie": {
                    "type": "string"
                  },
                  "classification": {
                    "type": "string"
                  },
                  "mode_of_inheritance": {
                    "type": "string"
                  },
                  "submitter": {
                    "type": "string"
                  },
                  "submitted_date": {
                    "type": "string"
                  }
                }
              }
            },
            "submission_count": {
              "type": "integer"
            },
            "unique_diseases": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `GenCC_search_gene`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("gen_cc_search_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.get_joint_associated_diseases_by_hpo_id_list` {#dnakit.databases.get_joint_associated_diseases_by_hpo_id_list}

Retrieve diseases associated with a list of phenotypes or symptoms by a list of HPO IDs.

```text
dnakit.databases.get_joint_associated_diseases_by_hpo_id_list(*, HPO_ID_list, limit=None, offset=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `HPO_ID_list` | array | Required by input schema | List of phenotypes or symptoms {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;, &quot;description&quot;: &quot;The HPO ID of the phenotype or symptom.&quot;}} |
| `limit` | integer | None | Number of entries to fetch. |
| `offset` | integer | None | Number of initial entries to skip. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "HPO_ID_list": {
          "type": "array",
          "description": "List of phenotypes or symptoms",
          "items": {
            "type": "string",
            "description": "The HPO ID of the phenotype or symptom."
          }
        },
        "limit": {
          "type": "integer",
          "description": "Number of entries to fetch."
        },
        "offset": {
          "type": "integer",
          "description": "Number of initial entries to skip."
        }
      },
      "required": [
        "HPO_ID_list"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of disease names associated with the intersection of ALL provided HPO phenotype IDs. Returned when every input HPO ID matched at least one disease."
        },
        {
          "type": "object",
          "description": "Returned whenever the plain list would hide something: one or more input HPO IDs matched zero diseases, a phenotype's associations were truncated, or the intersection was larger than 'limit'. 'warning' is present only in the first two cases; 'total_matched' is always the full intersection size before 'limit'/'offset' were applied.",
          "properties": {
            "diseases": {
              "type": "array",
              "items": {
                "type": "string"
              },
              "description": "The page of the intersection selected by 'limit' and 'offset', in the first phenotype's upstream relevance order."
            },
            "total_matched": {
              "type": "integer",
              "description": "Number of diseases carrying every usable input phenotype, before 'limit' and 'offset' were applied. Use this as the denominator; len(diseases) is only the page size."
            },
            "warning": {
              "type": "string",
              "description": "Present only when a phenotype contributed no associations or its associations were truncated. Absent when the only reason for the object form is that the intersection exceeded 'limit'."
            }
          },
          "required": [
            "diseases",
            "total_matched"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `get_joint_associated_diseases_by_HPO_ID_list`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("get_joint_associated_diseases_by_hpo_id_list", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpo_get_diseases_by_phenotype` {#dnakit.databases.hpo_get_diseases_by_phenotype}

Get the diseases that present with an HPO phenotype term, from the JAX HPO network annotations. Returns disease IDs (ORPHA/OMIM/DECIPHER), names, and the cross-referenced MONDO ID. Example: HP:0001250 (Seizure) is a feature of thousands of diseases. Use to build a differential diagnosis from a set of observed phenotypes. A phenotype term can be linked to thousands of diseases but 'limit' caps a single response at 500 -- check metadata.has_more in the response and pass 'offset' to page through the rest; a well-known disease for a common phenotype may not appear in the first page.

```text
dnakit.databases.hpo_get_diseases_by_phenotype(*, term_id=None, limit=None, offset=None, hp_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `term_id` | string | None | HPO term identifier (e.g., &#x27;HP:0001250&#x27;). Must start with &#x27;HP:&#x27; followed by 7 digits. |
| `limit` | integer &#124; null | None | Maximum number of results to return (1-500, default 50). Phenotype terms can be linked to thousands of diseases. |
| `offset` | integer &#124; null | None | Number of results to skip before returning &#x27;limit&#x27; results (default 0). Use with the response&#x27;s metadata.has_more/metadata.total to page through phenotypes with more diseases than fit in one call. |
| `hp_id` | string | None | DNAKit optional alias for term_id |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "term_id": {
          "type": "string",
          "description": "HPO term identifier (e.g., 'HP:0001250'). Must start with 'HP:' followed by 7 digits."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of results to return (1-500, default 50). Phenotype terms can be linked to thousands of diseases."
        },
        "offset": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Number of results to skip before returning 'limit' results (default 0). Use with the response's metadata.has_more/metadata.total to page through phenotypes with more diseases than fit in one call."
        },
        "hp_id": {
          "type": "string",
          "description": "DNAKit optional alias for term_id"
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "diseases": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "id": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Identifier (e.g. 'NCBIGene:5649' for genes, 'ORPHA:..'/'OMIM:..' for diseases)"
                  },
                  "name": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Name"
                  },
                  "mondo_id": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Cross-referenced MONDO disease ID"
                  }
                }
              }
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

??? info "DNAKit argument extension rules"

    ```json
    {
      "aliases": {
        "hp_id": "term_id"
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `HPO_get_diseases_by_phenotype`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("hpo_get_diseases_by_phenotype", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpo_get_genes_by_phenotype` {#dnakit.databases.hpo_get_genes_by_phenotype}

Get the genes associated with an HPO phenotype term (genes in which variants cause this phenotype), from the JAX HPO network annotations. Returns NCBI gene IDs and symbols. Example: HP:0001250 (Seizure) is linked to genes such as RELN, SCN1A, and ~2000 others. Use to go from an observed patient phenotype to candidate disease genes for variant prioritization. A phenotype term can be linked to thousands of genes but 'limit' caps a single response at 500 -- check metadata.has_more in the response and pass 'offset' to page through the rest; a well-known gene for a common phenotype may not appear in the first page.

```text
dnakit.databases.hpo_get_genes_by_phenotype(*, term_id, limit=None, offset=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `term_id` | string | Required by input schema | HPO term identifier (e.g., &#x27;HP:0001250&#x27;). Must start with &#x27;HP:&#x27; followed by 7 digits. |
| `limit` | integer &#124; null | None | Maximum number of results to return (1-500, default 50). Phenotype terms can be linked to thousands of genes. |
| `offset` | integer &#124; null | None | Number of results to skip before returning &#x27;limit&#x27; results (default 0). Use with the response&#x27;s metadata.has_more/metadata.total to page through phenotypes with more genes than fit in one call. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "term_id": {
          "type": "string",
          "description": "HPO term identifier (e.g., 'HP:0001250'). Must start with 'HP:' followed by 7 digits."
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of results to return (1-500, default 50). Phenotype terms can be linked to thousands of genes."
        },
        "offset": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Number of results to skip before returning 'limit' results (default 0). Use with the response's metadata.has_more/metadata.total to page through phenotypes with more genes than fit in one call."
        }
      },
      "required": [
        "term_id"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "genes": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "id": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Identifier (e.g. 'NCBIGene:5649' for genes, 'ORPHA:..'/'OMIM:..' for diseases)"
                  },
                  "name": {
                    "type": [
                      "string",
                      "null"
                    ],
                    "description": "Name"
                  }
                }
              }
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `HPO_get_genes_by_phenotype`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("hpo_get_genes_by_phenotype", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpo_get_term` {#dnakit.databases.hpo_get_term}

Get detailed information about a Human Phenotype Ontology (HPO) term by its ID. HPO provides standardized descriptions of phenotypic abnormalities in human disease with over 18,000 terms. Returns the term name, definition, synonyms, cross-references to UMLS/SNOMED, descendant count, and translations. Example: HP:0001250 (Seizure) returns definition, 346 descendants, synonyms including 'Epileptic seizure' and 'Epilepsy', and cross-references to UMLS and SNOMED-CT.

```text
dnakit.databases.hpo_get_term(*, term_id, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `term_id` | string | Required by input schema | HPO term identifier. Must start with &#x27;HP:&#x27; followed by 7 digits. Examples: &#x27;HP:0001250&#x27; (Seizure), &#x27;HP:0001263&#x27; (Global developmental delay), &#x27;HP:0000252&#x27; (Microcephaly), &#x27;HP:0001249&#x27; (Intellectual disability). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "term_id": {
          "type": "string",
          "description": "HPO term identifier. Must start with 'HP:' followed by 7 digits. Examples: 'HP:0001250' (Seizure), 'HP:0001263' (Global developmental delay), 'HP:0000252' (Microcephaly), 'HP:0001249' (Intellectual disability)."
        }
      },
      "required": [
        "term_id"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "id": {
              "type": [
                "string",
                "null"
              ],
              "description": "HPO term ID (e.g., HP:0001250)"
            },
            "name": {
              "type": [
                "string",
                "null"
              ],
              "description": "Term name"
            },
            "definition": {
              "type": [
                "string",
                "null"
              ],
              "description": "Full term definition"
            },
            "comment": {
              "type": [
                "string",
                "null"
              ],
              "description": "Additional comments or notes"
            },
            "synonyms": {
              "type": "array",
              "items": {
                "type": "string"
              },
              "description": "Alternative names for the term"
            },
            "descendant_count": {
              "type": [
                "integer",
                "null"
              ],
              "description": "Number of more specific child terms"
            },
            "xrefs": {
              "type": "array",
              "items": {
                "type": "string"
              },
              "description": "Cross-references (UMLS, SNOMED-CT, etc.)"
            },
            "translations": {
              "type": "array",
              "items": {
                "type": "object"
              },
              "description": "Translations in other languages"
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `HPO_get_term`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("hpo_get_term", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.hpo_search_terms` {#dnakit.databases.hpo_search_terms}

Search for Human Phenotype Ontology (HPO) terms by keyword. HPO contains over 18,000 standardized phenotype terms used in clinical genetics and rare disease research. Returns matching terms with definitions, synonym counts, and descendant counts. Useful for finding the correct HPO ID for a clinical phenotype before querying disease or gene associations. Example: searching 'seizure' returns terms like 'Bilateral tonic-clonic seizure' (HP:0007334), 'Generalized-onset seizure' (HP:0002197), 'Focal-onset seizure' (HP:0007359).

```text
dnakit.databases.hpo_search_terms(*, query, max_results=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search keyword for phenotype terms. Examples: &#x27;seizure&#x27;, &#x27;intellectual disability&#x27;, &#x27;cardiomyopathy&#x27;, &#x27;scoliosis&#x27;, &#x27;microcephaly&#x27;. |
| `max_results` | integer &#124; null | None | Maximum number of results to return (default 10, max 50). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search keyword for phenotype terms. Examples: 'seizure', 'intellectual disability', 'cardiomyopathy', 'scoliosis', 'microcephaly'."
        },
        "max_results": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Maximum number of results to return (default 10, max 50)."
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "id": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "HPO term ID"
              },
              "name": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Term name"
              },
              "definition": {
                "type": [
                  "string",
                  "null"
                ],
                "description": "Term definition"
              },
              "descendant_count": {
                "type": [
                  "integer",
                  "null"
                ],
                "description": "Number of child terms"
              },
              "synonyms": {
                "type": "array",
                "items": {
                  "type": "string"
                },
                "description": "First 5 synonyms"
              }
            }
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `HPO_search_terms`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("hpo_search_terms", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.marrvel_get_gene` {#dnakit.databases.marrvel_get_gene}

Get aggregated identity and annotation for a human gene by symbol from MARRVEL: cross-references (HGNC, Ensembl, Entrez, UniProt), chromosome/cytogenetic location, gene type, aliases, previous symbols, and the Entrez gene summary. Use as a one-call gene-identity resolver in rare-disease / Mendelian variant triage workflows. Does NOT include an OMIM ID -- MARRVEL's own gene-lookup endpoint returns an unreliable placeholder value for it (confirmed live), use MARRVEL_get_omim_phenotypes for OMIM associations instead. MARRVEL indexes only current approved HGNC symbols, so a retired symbol or alias (e.g. 'RP20'/'LCA2' for RPE65) is resolved once via HGNC and re-queried; when that happens the substitution is disclosed in metadata.resolved_symbol and metadata.resolved_note, and metadata.query_symbol keeps the symbol you asked for. No API key required.

```text
dnakit.databases.marrvel_get_gene(*, symbol=None, gene_symbol=None, gene=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string &#124; null | None | HGNC gene symbol, e.g. &#x27;CFTR&#x27;, &#x27;BRCA1&#x27;. |
| `gene_symbol` | string &#124; null | None | Alias for symbol. |
| `gene` | string &#124; null | None | Alias for symbol. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "symbol": {
          "type": [
            "string",
            "null"
          ],
          "description": "HGNC gene symbol, e.g. 'CFTR', 'BRCA1'."
        },
        "gene_symbol": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for symbol."
        },
        "gene": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for symbol."
        }
      },
      "required": [],
      "anyOf": [
        {
          "required": [
            "symbol"
          ]
        },
        {
          "required": [
            "gene_symbol"
          ]
        },
        {
          "required": [
            "gene"
          ]
        }
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "object"
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `MARRVEL_get_gene`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("marrvel_get_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.marrvel_get_omim_phenotypes` {#dnakit.databases.marrvel_get_omim_phenotypes}

Get OMIM phenotype/disease associations for a human gene by symbol via MARRVEL: each association's phenotype name, gene and phenotype MIM numbers, mode of inheritance, and phenotypic series. Use to enumerate the Mendelian diseases linked to a gene and their inheritance patterns. MARRVEL indexes only current approved HGNC symbols, so a retired symbol or alias (e.g. 'RP20'/'LCA2' for RPE65) is resolved once via HGNC and re-queried; when that happens the substitution is disclosed in metadata.resolved_symbol and metadata.resolved_note. No API key required.

```text
dnakit.databases.marrvel_get_omim_phenotypes(*, symbol=None, gene_symbol=None, gene=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `symbol` | string &#124; null | None | HGNC gene symbol, e.g. &#x27;CFTR&#x27;, &#x27;BRCA1&#x27;. |
| `gene_symbol` | string &#124; null | None | Alias for symbol. |
| `gene` | string &#124; null | None | Alias for symbol. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "symbol": {
          "type": [
            "string",
            "null"
          ],
          "description": "HGNC gene symbol, e.g. 'CFTR', 'BRCA1'."
        },
        "gene_symbol": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for symbol."
        },
        "gene": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alias for symbol."
        }
      },
      "required": [],
      "anyOf": [
        {
          "required": [
            "symbol"
          ]
        },
        {
          "required": [
            "gene_symbol"
          ]
        },
        {
          "required": [
            "gene"
          ]
        }
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "oneOf": [
        {
          "type": "array",
          "items": {
            "type": "object"
          }
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ]
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `MARRVEL_get_omim_phenotypes`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("marrvel_get_omim_phenotypes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_get_gene_diseases` {#dnakit.databases.monarch_get_gene_diseases}

Get diseases associated with a gene from Monarch Initiative. Input is a gene CURIE (e.g., 'HGNC:11998' for TP53). Returns disease associations with evidence. Bridges gene→disease across species for hypothesis generation.

```text
dnakit.databases.monarch_get_gene_diseases(*, subject, limit=None, offset=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `subject` | string | Required by input schema | Gene CURIE identifier (e.g., &#x27;HGNC:11998&#x27; for TP53). Use Monarch_search_gene to find this. |
| `limit` | integer | None | Maximum disease associations to return (default: 100) {&quot;default&quot;: 100} |
| `offset` | integer | None | Number of initial entries to skip. {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "subject": {
          "type": "string",
          "description": "Gene CURIE identifier (e.g., 'HGNC:11998' for TP53). Use Monarch_search_gene to find this."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum disease associations to return (default: 100)",
          "default": 100
        },
        "offset": {
          "type": "integer",
          "description": "Number of initial entries to skip.",
          "default": 0
        }
      },
      "required": [
        "subject"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "items": {
          "type": "array",
          "description": "List of gene-disease associations"
        },
        "total": {
          "type": "integer"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Monarch_get_gene_diseases`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_get_gene_diseases", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_get_gene_phenotypes` {#dnakit.databases.monarch_get_gene_phenotypes}

Get phenotypes associated with a gene from Monarch Initiative (cross-species). Input is a gene CURIE (e.g., 'HGNC:11998' for TP53). Returns human and model organism phenotypes (HP/MP terms). Use Monarch_search_gene first to find the CURIE.

```text
dnakit.databases.monarch_get_gene_phenotypes(*, subject, limit=None, offset=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `subject` | string | Required by input schema | Gene CURIE identifier (e.g., &#x27;HGNC:11998&#x27; for TP53, &#x27;HGNC:1100&#x27; for BRCA1). Use Monarch_search_gene to find this. |
| `limit` | integer | None | Maximum phenotype associations to return (default: 100) {&quot;default&quot;: 100} |
| `offset` | integer | None | Number of initial entries to skip. {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "subject": {
          "type": "string",
          "description": "Gene CURIE identifier (e.g., 'HGNC:11998' for TP53, 'HGNC:1100' for BRCA1). Use Monarch_search_gene to find this."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum phenotype associations to return (default: 100)",
          "default": 100
        },
        "offset": {
          "type": "integer",
          "description": "Number of initial entries to skip.",
          "default": 0
        }
      },
      "required": [
        "subject"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "items": {
          "type": "array",
          "description": "List of gene-phenotype associations with phenotype IDs and labels"
        },
        "total": {
          "type": "integer"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Monarch_get_gene_phenotypes`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_get_gene_phenotypes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_search` {#dnakit.databases.monarch_search}

Search the Monarch Initiative knowledge graph for biomedical entities by name or keyword. Returns genes, diseases, phenotypes, and other entities matching the query, with CURIE identifiers and categories. Filter by category (biolink:Gene, biolink:Disease, biolink:PhenotypicFeature) to narrow results. Use returned CURIEs with Monarch_get_entity and Monarch_get_associations for detailed lookups.

```text
dnakit.databases.monarch_search(*, query, category=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Search query: gene names, disease names, phenotype terms. Examples: &#x27;BRCA1&#x27;, &#x27;Alzheimer&#x27;, &#x27;seizure&#x27;, &#x27;insulin resistance&#x27;. |
| `category` | string | None | Filter by entity category. Options: &#x27;biolink:Gene&#x27;, &#x27;biolink:Disease&#x27;, &#x27;biolink:PhenotypicFeature&#x27;. Default: all categories. |
| `limit` | integer | None | Maximum results (default: 10, max: 50). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query: gene names, disease names, phenotype terms. Examples: 'BRCA1', 'Alzheimer', 'seizure', 'insulin resistance'."
        },
        "category": {
          "type": "string",
          "description": "Filter by entity category. Options: 'biolink:Gene', 'biolink:Disease', 'biolink:PhenotypicFeature'. Default: all categories."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results (default: 10, max: 50)."
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "id": {
            "type": [
              "string",
              "null"
            ]
          },
          "name": {
            "type": [
              "string",
              "null"
            ]
          },
          "category": {
            "type": [
              "string",
              "null"
            ]
          },
          "symbol": {
            "type": [
              "string",
              "null"
            ]
          },
          "description": {
            "type": [
              "string",
              "null"
            ]
          },
          "taxon": {
            "type": [
              "string",
              "null"
            ]
          },
          "taxon_label": {
            "type": [
              "string",
              "null"
            ]
          },
          "deprecated": {
            "type": "boolean",
            "description": "True if this ontology term is obsolete. Obsolete terms remain queryable but carry no annotations."
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Monarch_search`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_search_gene` {#dnakit.databases.monarch_search_gene}

Search Monarch Initiative knowledge graph for a gene by symbol or name. Returns HGNC/NCBIGene CURIEs, cross-species orthologs, and gene metadata. Use this first to get the gene CURIE needed for Monarch_get_gene_phenotypes and Monarch_get_gene_diseases.

```text
dnakit.databases.monarch_search_gene(*, query, limit=None, offset=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | Required by input schema | Gene symbol or name (e.g., &#x27;TP53&#x27;, &#x27;BRCA1&#x27;, &#x27;WDR7&#x27;) |
| `limit` | integer | None | Number of results to return (default: 10) {&quot;default&quot;: 10} |
| `offset` | integer | None | Number of initial entries to skip. {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Gene symbol or name (e.g., 'TP53', 'BRCA1', 'WDR7')"
        },
        "limit": {
          "type": "integer",
          "description": "Number of results to return (default: 10)",
          "default": 10
        },
        "offset": {
          "type": "integer",
          "description": "Number of initial entries to skip.",
          "default": 0
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "items": {
          "type": "array",
          "description": "List of matching gene entities with IDs, names, and categories"
        },
        "total": {
          "type": "integer"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Monarch_search_gene`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_search_gene", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_v3_get_associations` {#dnakit.databases.monarch_v3_get_associations}

Query cross-species associations between biomedical entities in the Monarch Initiative knowledge graph. Find gene-phenotype, disease-phenotype, gene-disease, or gene-pathway associations. Uses Biolink model categories to specify association type. Integrates data from OMIM, ClinVar, HPO, MGI, ZFIN, FlyBase, WormBase. Associations are DIRECTED: every category has a fixed subject side and object side, so anchor your entity on the correct side. Supply 'subject', 'object', or both (at least one is required). Example: HGNC:11998 (TP53) as subject with category GeneToPhenotypicFeatureAssociation returns phenotypes like hepatocellular carcinoma; MONDO:0006559 as OBJECT with category CausalGeneToDiseaseAssociation returns the genes that cause that disease (NCSTN, PSENEN, PSEN1).

```text
dnakit.databases.monarch_v3_get_associations(*, subject=None, object=None, category, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `subject` | string | None | Subject-side entity CURIE, i.e. the entity on the LEFT of the association (the arrow&#x27;s source). For gene-&gt;disease categories (CausalGeneToDiseaseAssociation, CorrelatedGeneToDiseaseAssociation) the subject is the GENE; for gene-&gt;phenotype and gene-&gt;pathway the subject is the gene; for DiseaseToPhenotypicFeatureAssociation the subject is the disease. Examples: &#x27;HGNC:11998&#x27; (TP53), &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;HP:0001250&#x27; (seizure). Optional, but at least one of &#x27;subject&#x27; / &#x27;object&#x27; must be supplied. |
| `object` | string | None | Object-side entity CURIE, i.e. the entity on the RIGHT of the association (the arrow&#x27;s target). Use this to query the association backwards. For gene-&gt;disease categories the object is the DISEASE, so &#x27;which genes cause this disease?&#x27; is object=&#x27;MONDO:0006559&#x27; with category=&#x27;biolink:CausalGeneToDiseaseAssociation&#x27;; for gene/disease-&gt;phenotype categories the object is the HP phenotype term. Optional, but at least one of &#x27;subject&#x27; / &#x27;object&#x27; must be supplied. Supplying both narrows to associations linking exactly that pair. |
| `category` | string | Required by input schema | Biolink association category, written as subject-&gt;object. Options: &#x27;biolink:GeneToPhenotypicFeatureAssociation&#x27; (gene-&gt;phenotypes), &#x27;biolink:DiseaseToPhenotypicFeatureAssociation&#x27; (disease-&gt;phenotypes), &#x27;biolink:CorrelatedGeneToDiseaseAssociation&#x27; (gene-&gt;diseases, correlated), &#x27;biolink:CausalGeneToDiseaseAssociation&#x27; (gene-&gt;diseases, causal), &#x27;biolink:GeneToPathwayAssociation&#x27; (gene-&gt;pathways), &#x27;biolink:VariantToDiseaseAssociation&#x27; (variant-&gt;diseases), &#x27;biolink:GeneToExpressionSiteAssociation&#x27; (gene-&gt;expression sites/tissues). The left-hand side goes in &#x27;subject&#x27; and the right-hand side goes in &#x27;object&#x27;. |
| `limit` | integer | None | Maximum number of associations to return (default: 20, max: 200). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "subject": {
          "type": "string",
          "description": "Subject-side entity CURIE, i.e. the entity on the LEFT of the association (the arrow's source). For gene->disease categories (CausalGeneToDiseaseAssociation, CorrelatedGeneToDiseaseAssociation) the subject is the GENE; for gene->phenotype and gene->pathway the subject is the gene; for DiseaseToPhenotypicFeatureAssociation the subject is the disease. Examples: 'HGNC:11998' (TP53), 'MONDO:0005148' (type 2 diabetes), 'HP:0001250' (seizure). Optional, but at least one of 'subject' / 'object' must be supplied."
        },
        "object": {
          "type": "string",
          "description": "Object-side entity CURIE, i.e. the entity on the RIGHT of the association (the arrow's target). Use this to query the association backwards. For gene->disease categories the object is the DISEASE, so 'which genes cause this disease?' is object='MONDO:0006559' with category='biolink:CausalGeneToDiseaseAssociation'; for gene/disease->phenotype categories the object is the HP phenotype term. Optional, but at least one of 'subject' / 'object' must be supplied. Supplying both narrows to associations linking exactly that pair."
        },
        "category": {
          "type": "string",
          "description": "Biolink association category, written as subject->object. Options: 'biolink:GeneToPhenotypicFeatureAssociation' (gene->phenotypes), 'biolink:DiseaseToPhenotypicFeatureAssociation' (disease->phenotypes), 'biolink:CorrelatedGeneToDiseaseAssociation' (gene->diseases, correlated), 'biolink:CausalGeneToDiseaseAssociation' (gene->diseases, causal), 'biolink:GeneToPathwayAssociation' (gene->pathways), 'biolink:VariantToDiseaseAssociation' (variant->diseases), 'biolink:GeneToExpressionSiteAssociation' (gene->expression sites/tissues). The left-hand side goes in 'subject' and the right-hand side goes in 'object'."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of associations to return (default: 20, max: 200)."
        }
      },
      "required": [
        "category"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array"
    }
    ```

**Execution and provenance:** `network`；upstream name `MonarchV3_get_associations`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_v3_get_associations", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_v3_get_entity` {#dnakit.databases.monarch_v3_get_entity}

Look up detailed information about any biomedical entity in the Monarch Initiative knowledge graph by its CURIE identifier. Supports genes (HGNC), diseases (MONDO, OMIM), phenotypes (HP), variants, and more. Returns the entity name, category, description, synonyms, cross-references, and species. Example: HGNC:11998 returns TP53 gene with cross-refs to ENSEMBL, OMIM; MONDO:0005148 returns type 2 diabetes mellitus with detailed description.

```text
dnakit.databases.monarch_v3_get_entity(*, entity_id, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entity_id` | string | Required by input schema | Entity CURIE identifier. Examples: &#x27;HGNC:11998&#x27; (TP53), &#x27;HGNC:1100&#x27; (BRCA1), &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;HP:0001250&#x27; (seizure), &#x27;OMIM:191170&#x27; (Li-Fraumeni syndrome). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "entity_id": {
          "type": "string",
          "description": "Entity CURIE identifier. Examples: 'HGNC:11998' (TP53), 'HGNC:1100' (BRCA1), 'MONDO:0005148' (type 2 diabetes), 'HP:0001250' (seizure), 'OMIM:191170' (Li-Fraumeni syndrome)."
        }
      },
      "required": [
        "entity_id"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object"
    }
    ```

**Execution and provenance:** `network`；upstream name `MonarchV3_get_entity`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_v3_get_entity", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_v3_get_histopheno` {#dnakit.databases.monarch_v3_get_histopheno}

Get the phenotype profile of a disease grouped by body system from the Monarch Initiative. Returns counts of phenotypic features (HPO terms) organized by anatomical system (e.g., nervous_system, skeletal_system, eye, cardiovascular_system). Only systems with at least one phenotype are returned, sorted by count descending. Useful for understanding which organ systems are most affected by a disease. Example: MONDO:0007947 (Marfan syndrome) returns high counts for head_neck (58), musculature (58), skeletal_system (56), eye (31).

```text
dnakit.databases.monarch_v3_get_histopheno(*, entity_id, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entity_id` | string | Required by input schema | Disease CURIE identifier. Examples: &#x27;MONDO:0007947&#x27; (Marfan syndrome), &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;MONDO:0004975&#x27; (Alzheimer disease), &#x27;MONDO:0005015&#x27; (diabetes mellitus). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "entity_id": {
          "type": "string",
          "description": "Disease CURIE identifier. Examples: 'MONDO:0007947' (Marfan syndrome), 'MONDO:0005148' (type 2 diabetes), 'MONDO:0004975' (Alzheimer disease), 'MONDO:0005015' (diabetes mellitus)."
        }
      },
      "required": [
        "entity_id"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array"
    }
    ```

**Execution and provenance:** `network`；upstream name `MonarchV3_get_histopheno`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_v3_get_histopheno", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_v3_get_mappings` {#dnakit.databases.monarch_v3_get_mappings}

Get cross-ontology mappings for a biomedical entity from the Monarch Initiative. Returns equivalence mappings between ontology identifiers, such as MONDO to OMIM, ICD10, SNOMED CT, NCIT, MESH, DOID, and Orphanet. Useful for translating disease identifiers across coding systems (e.g., for clinical data integration, ICD-10 billing codes, or linking to OMIM entries). Example: MONDO:0005148 (type 2 diabetes) maps to ICD10CM:E11, DOID:9352, OMIM:125853, MESH:D003924.

```text
dnakit.databases.monarch_v3_get_mappings(*, entity_id, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `entity_id` | string | Required by input schema | Entity CURIE identifier to find mappings for. Examples: &#x27;MONDO:0005148&#x27; (type 2 diabetes), &#x27;MONDO:0007947&#x27; (Marfan syndrome), &#x27;HP:0001250&#x27; (seizure). |
| `limit` | integer | None | Maximum number of mappings to return (default: 20, max: 100). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "entity_id": {
          "type": "string",
          "description": "Entity CURIE identifier to find mappings for. Examples: 'MONDO:0005148' (type 2 diabetes), 'MONDO:0007947' (Marfan syndrome), 'HP:0001250' (seizure)."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of mappings to return (default: 20, max: 100)."
        }
      },
      "required": [
        "entity_id"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array"
    }
    ```

**Execution and provenance:** `network`；upstream name `MonarchV3_get_mappings`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_v3_get_mappings", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_v3_phenotype_similarity_search` {#dnakit.databases.monarch_v3_phenotype_similarity_search}

Find diseases that best match a set of phenotypes (HPO terms) using semantic similarity from the Monarch Initiative. Given a patient phenotype profile as a list of HPO IDs, returns diseases ranked by phenotype overlap score. Uses information content-based similarity (AncestorInformationContent metric). Essential for rare disease differential diagnosis: input observed phenotypes, get ranked candidate diseases. Example: [HP:0001250 (seizure), HP:0001249 (intellectual disability), HP:0002353 (EEG abnormality)] returns diseases like intellectual disability AR 31, Dravet syndrome ranked by semantic similarity score.

```text
dnakit.databases.monarch_v3_phenotype_similarity_search(*, phenotypes, group=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `phenotypes` | array | Required by input schema | List of HPO term CURIEs representing observed phenotypes. Examples: &#91;&#x27;HP:0001250&#x27;, &#x27;HP:0001249&#x27;, &#x27;HP:0002353&#x27;&#93; for seizure + intellectual disability + EEG abnormality. {&quot;items&quot;: {&quot;type&quot;: &quot;string&quot;}} |
| `group` | string | None | Entity group to search against. Default: &#x27;Human Diseases&#x27;. Options: &#x27;Human Diseases&#x27;, &#x27;Mouse&#x27;, &#x27;Rat&#x27;, &#x27;Zebrafish&#x27;, &#x27;C. elegans&#x27;, &#x27;Drosophila&#x27;. |
| `limit` | integer | None | Maximum number of matching diseases to return (default: 10, max: 50). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "phenotypes": {
          "type": "array",
          "items": {
            "type": "string"
          },
          "description": "List of HPO term CURIEs representing observed phenotypes. Examples: ['HP:0001250', 'HP:0001249', 'HP:0002353'] for seizure + intellectual disability + EEG abnormality."
        },
        "group": {
          "type": "string",
          "description": "Entity group to search against. Default: 'Human Diseases'. Options: 'Human Diseases', 'Mouse', 'Rat', 'Zebrafish', 'C. elegans', 'Drosophila'."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of matching diseases to return (default: 10, max: 50)."
        }
      },
      "required": [
        "phenotypes"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array"
    }
    ```

**Execution and provenance:** `network`；upstream name `MonarchV3_phenotype_similarity_search`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_v3_phenotype_similarity_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.monarch_v3_search` {#dnakit.databases.monarch_v3_search}

Search the Monarch Initiative knowledge graph for biomedical entities by name or keyword. Returns genes, diseases, phenotypes, and other entities matching the query, with their CURIE identifiers and categories. Useful for finding entity IDs to use with MonarchV3_get_entity and MonarchV3_get_associations. Example: searching 'BRCA1' returns HGNC:1100 (gene), 'Alzheimer' returns MONDO disease entries.

```text
dnakit.databases.monarch_v3_search(*, query=None, limit=None, category=None, q=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `query` | string | None | Search query. Can be gene names, disease names, phenotype terms, etc. Examples: &#x27;BRCA1&#x27;, &#x27;breast cancer&#x27;, &#x27;seizure&#x27;, &#x27;insulin resistance&#x27;. |
| `limit` | integer | None | Maximum results to return (default: 10, max: 50). |
| `category` | string | None | Filter by entity category. Options: &#x27;biolink:Gene&#x27;, &#x27;biolink:Disease&#x27;, &#x27;biolink:PhenotypicFeature&#x27;, &#x27;biolink:NamedThing&#x27;. Default: all categories. |
| `q` | string | None | Alias for query parameter. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "query": {
          "type": "string",
          "description": "Search query. Can be gene names, disease names, phenotype terms, etc. Examples: 'BRCA1', 'breast cancer', 'seizure', 'insulin resistance'."
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results to return (default: 10, max: 50)."
        },
        "category": {
          "type": "string",
          "description": "Filter by entity category. Options: 'biolink:Gene', 'biolink:Disease', 'biolink:PhenotypicFeature', 'biolink:NamedThing'. Default: all categories."
        },
        "q": {
          "type": "string",
          "description": "Alias for query parameter."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "array"
    }
    ```

**Execution and provenance:** `network`；upstream name `MonarchV3_search`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("monarch_v3_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ols_get_term_ancestors` {#dnakit.databases.ols_get_term_ancestors}

Get ancestor terms of a specific term in an ontology. Use `term_id` with CURIE format (e.g., 'HP:0001928') — the ontology is auto-inferred from the prefix. Or supply `term_iri` + `ontology` explicitly.

```text
dnakit.databases.ols_get_term_ancestors(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (get_term_ancestors) {&quot;const&quot;: &quot;get_term_ancestors&quot;, &quot;default&quot;: &quot;get_term_ancestors&quot;} |
| `term_iri` | string | None | The full IRI of the term (e.g., http://www.ebi.ac.uk/efo/EFO_0000408) |
| `term_id` | string | None | Short-form CURIE (e.g., &#x27;HP:0001928&#x27;, &#x27;GO:0008150&#x27;). Ontology is auto-inferred from prefix. |
| `ontology` | string | None | Ontology ID (e.g., &#x27;hp&#x27;, &#x27;go&#x27;, &#x27;efo&#x27;). Auto-inferred from term_id prefix when not supplied. |
| `include_obsolete` | boolean | None | Include obsolete terms (default: false) {&quot;default&quot;: false} |
| `size` | integer | None | Number of results to return (default: 20) {&quot;default&quot;: 20} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_term_ancestors",
          "description": "The operation to perform (get_term_ancestors)",
          "default": "get_term_ancestors"
        },
        "term_iri": {
          "type": "string",
          "description": "The full IRI of the term (e.g., http://www.ebi.ac.uk/efo/EFO_0000408)"
        },
        "term_id": {
          "type": "string",
          "description": "Short-form CURIE (e.g., 'HP:0001928', 'GO:0008150'). Ontology is auto-inferred from prefix."
        },
        "ontology": {
          "type": "string",
          "description": "Ontology ID (e.g., 'hp', 'go', 'efo'). Auto-inferred from term_id prefix when not supplied."
        },
        "include_obsolete": {
          "type": "boolean",
          "description": "Include obsolete terms (default: false)",
          "default": false
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return (default: 20)",
          "default": 20
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "page": {
          "type": "integer"
        },
        "numElements": {
          "type": "integer"
        },
        "totalPages": {
          "type": "integer"
        },
        "totalElements": {
          "type": "integer"
        },
        "elements": {
          "type": "array"
        },
        "term_iri": {
          "type": "string"
        },
        "ontology": {
          "type": "string"
        },
        "filters": {
          "type": "object",
          "properties": {
            "include_obsolete": {
              "type": "boolean"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ols_get_term_ancestors`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ols_get_term_ancestors", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ols_get_term_children` {#dnakit.databases.ols_get_term_children}

Get child terms of a specific term in an ontology. Use `term_id` with CURIE format (e.g., 'HP:0001928') — the ontology is auto-inferred from the prefix. Or supply `term_iri` + `ontology` explicitly.

```text
dnakit.databases.ols_get_term_children(*, operation=None, term_iri=None, term_id=None, ontology=None, include_obsolete=None, size=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (get_term_children) {&quot;const&quot;: &quot;get_term_children&quot;, &quot;default&quot;: &quot;get_term_children&quot;} |
| `term_iri` | string | None | The full IRI of the term (e.g., http://www.ebi.ac.uk/efo/EFO_0000408) |
| `term_id` | string | None | Short-form CURIE (e.g., &#x27;HP:0001928&#x27;, &#x27;GO:0008150&#x27;). Ontology is auto-inferred from prefix. |
| `ontology` | string | None | Ontology ID (e.g., &#x27;hp&#x27;, &#x27;go&#x27;, &#x27;efo&#x27;). Auto-inferred from term_id prefix when not supplied. |
| `include_obsolete` | boolean | None | Include obsolete terms (default: false) {&quot;default&quot;: false} |
| `size` | integer | None | Number of results to return (default: 20) {&quot;default&quot;: 20} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_term_children",
          "description": "The operation to perform (get_term_children)",
          "default": "get_term_children"
        },
        "term_iri": {
          "type": "string",
          "description": "The full IRI of the term (e.g., http://www.ebi.ac.uk/efo/EFO_0000408)"
        },
        "term_id": {
          "type": "string",
          "description": "Short-form CURIE (e.g., 'HP:0001928', 'GO:0008150'). Ontology is auto-inferred from prefix."
        },
        "ontology": {
          "type": "string",
          "description": "Ontology ID (e.g., 'hp', 'go', 'efo'). Auto-inferred from term_id prefix when not supplied."
        },
        "include_obsolete": {
          "type": "boolean",
          "description": "Include obsolete terms (default: false)",
          "default": false
        },
        "size": {
          "type": "integer",
          "description": "Number of results to return (default: 20)",
          "default": 20
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "page": {
          "type": "integer"
        },
        "numElements": {
          "type": "integer"
        },
        "totalPages": {
          "type": "integer"
        },
        "totalElements": {
          "type": "integer"
        },
        "elements": {
          "type": "array"
        },
        "term_iri": {
          "type": "string"
        },
        "ontology": {
          "type": "string"
        },
        "filters": {
          "type": "object",
          "properties": {
            "include_obsolete": {
              "type": "boolean"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ols_get_term_children`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ols_get_term_children", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ols_get_term_info` {#dnakit.databases.ols_get_term_info}

Get detailed information about a specific term in OLS

```text
dnakit.databases.ols_get_term_info(*, operation=None, id=None, term_id=None, term_iri=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (get_term_info) {&quot;const&quot;: &quot;get_term_info&quot;, &quot;default&quot;: &quot;get_term_info&quot;} |
| `id` | string | None | Ontology term ID or IRI (e.g. HP:0001903, http://purl.obolibrary.org/obo/HP_0001903) |
| `term_id` | string | None | Alias for id — ontology term ID (e.g. HP:0001903) |
| `term_iri` | string | None | Alias for id. Full IRI or short-form ID (e.g., GO:0008150, EFO:0003971). Same as id. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_term_info",
          "description": "The operation to perform (get_term_info)",
          "default": "get_term_info"
        },
        "id": {
          "type": "string",
          "description": "Ontology term ID or IRI (e.g. HP:0001903, http://purl.obolibrary.org/obo/HP_0001903)"
        },
        "term_id": {
          "type": "string",
          "description": "Alias for id — ontology term ID (e.g. HP:0001903)"
        },
        "term_iri": {
          "type": "string",
          "description": "Alias for id. Full IRI or short-form ID (e.g., GO:0008150, EFO:0003971). Same as id."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "iri": {
          "type": "string"
        },
        "ontologyName": {
          "type": "string"
        },
        "shortForm": {
          "type": "string"
        },
        "label": {
          "type": "string"
        },
        "oboId": {
          "type": "string"
        },
        "isObsolete": {
          "type": "boolean"
        },
        "description": {
          "type": "array"
        },
        "synonyms": {
          "type": "array"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ols_get_term_info`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ols_get_term_info", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.ols_search_terms` {#dnakit.databases.ols_search_terms}

Search for terms in OLS (Ontology Lookup Service)

```text
dnakit.databases.ols_search_terms(*, operation=None, query, rows=None, size=None, ontology=None, exact_match=None, include_obsolete=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | The operation to perform (search_terms) {&quot;const&quot;: &quot;search_terms&quot;, &quot;default&quot;: &quot;search_terms&quot;} |
| `query` | string | Required by input schema | The search query for terms |
| `rows` | integer | None | Number of results to return (default: 10). Alias: size. {&quot;default&quot;: 10} |
| `size` | integer | None | Alias for rows. Number of results to return (default: 10). {&quot;default&quot;: 10} |
| `ontology` | string | None | Filter by specific ontology (optional) |
| `exact_match` | boolean | None | Restrict results to exact matches (default: false). The query must equal a term&#x27;s label or one of its exact synonyms in full - substring and word-overlap hits are excluded, so &#x27;fibroblast&#x27; returns CL:0000057 but not &#x27;skin fibroblast&#x27;. If the query is instead an identifier (a CURIE such as &#x27;CL:0000084&#x27;, the OBO underscore form &#x27;CL_0000084&#x27;, or a full IRI), the exact match is applied to the obo_id/short_form/iri fields, which returns that identifier in every ontology that declares it - pass `ontology` to narrow it to one. {&quot;default&quot;: false} |
| `include_obsolete` | boolean | None | Include obsolete terms (default: false) {&quot;default&quot;: false} |
| `limit` | integer &#124; null | None | Alias for rows. Maximum number of results to return (default: 10, max: 100). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search_terms",
          "description": "The operation to perform (search_terms)",
          "default": "search_terms"
        },
        "query": {
          "type": "string",
          "description": "The search query for terms"
        },
        "rows": {
          "type": "integer",
          "description": "Number of results to return (default: 10). Alias: size.",
          "default": 10
        },
        "size": {
          "type": "integer",
          "description": "Alias for rows. Number of results to return (default: 10).",
          "default": 10
        },
        "ontology": {
          "type": "string",
          "description": "Filter by specific ontology (optional)"
        },
        "exact_match": {
          "type": "boolean",
          "description": "Restrict results to exact matches (default: false). The query must equal a term's label or one of its exact synonyms in full - substring and word-overlap hits are excluded, so 'fibroblast' returns CL:0000057 but not 'skin fibroblast'. If the query is instead an identifier (a CURIE such as 'CL:0000084', the OBO underscore form 'CL_0000084', or a full IRI), the exact match is applied to the obo_id/short_form/iri fields, which returns that identifier in every ontology that declares it - pass `ontology` to narrow it to one.",
          "default": false
        },
        "include_obsolete": {
          "type": "boolean",
          "description": "Include obsolete terms (default: false)",
          "default": false
        },
        "limit": {
          "type": [
            "integer",
            "null"
          ],
          "description": "Alias for rows. Maximum number of results to return (default: 10, max: 100)."
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "response": {
          "type": "object",
          "properties": {
            "docs": {
              "type": "array"
            },
            "numFound": {
              "type": "integer"
            },
            "start": {
              "type": "integer"
            }
          }
        },
        "responseHeader": {
          "type": "object",
          "properties": {
            "QTime": {
              "type": "integer"
            },
            "status": {
              "type": "integer"
            }
          }
        },
        "facet_counts": {
          "type": "object",
          "properties": {
            "facet_fields": {
              "type": "object",
              "properties": {
                "ontologyPreferredPrefix": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "isDefiningOntology": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "ontologyId": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "ontologyIri": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "isObsolete": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "type": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                }
              }
            }
          }
        },
        "query": {
          "type": "string"
        },
        "filters": {
          "type": "object",
          "properties": {
            "ontology": {
              "type": "string"
            },
            "exact_match": {
              "type": "boolean"
            },
            "exact_match_fields": {
              "type": "string",
              "description": "Comma-separated OLS fields the exact match was applied to. Present only when exact_match is true."
            },
            "include_obsolete": {
              "type": "boolean"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `ols_search_terms`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("ols_search_terms", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.omim_get_clinical_synopsis` {#dnakit.databases.omim_get_clinical_synopsis}

Get clinical synopsis (phenotype features) for an OMIM phenotype entry. Returns structured clinical features organized by body system (e.g., neurologic, cardiovascular). Essential for phenotype-based disease matching. Requires OMIM_API_KEY.

```text
dnakit.databases.omim_get_clinical_synopsis(*, operation=None, mim_number, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_clinical_synopsis) {&quot;const&quot;: &quot;get_clinical_synopsis&quot;, &quot;default&quot;: &quot;get_clinical_synopsis&quot;} |
| `mim_number` | string | Required by input schema | OMIM MIM number for phenotype (e.g., 219700 for cystic fibrosis) |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_clinical_synopsis",
          "description": "Operation type (fixed: get_clinical_synopsis)",
          "default": "get_clinical_synopsis"
        },
        "mim_number": {
          "type": "string",
          "description": "OMIM MIM number for phenotype (e.g., 219700 for cystic fibrosis)"
        }
      },
      "required": [
        "mim_number"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "mimNumber": {
              "type": "integer"
            },
            "inheritance": {
              "type": "string",
              "description": "Inheritance pattern"
            },
            "neurologicCentralNervousSystem": {
              "type": "array",
              "description": "CNS features"
            },
            "cardiovascular": {
              "type": "array",
              "description": "Cardiovascular features"
            }
          },
          "description": "Clinical features organized by organ system"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OMIM_get_clinical_synopsis`。

**Provider credentials:** `OMIM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("omim_get_clinical_synopsis", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.omim_get_entry` {#dnakit.databases.omim_get_entry}

Get detailed OMIM entry by MIM number. Returns comprehensive information including text description, clinical features, inheritance, molecular genetics, and references. MIM numbers are 6-digit identifiers (e.g., 164730 for BRAF gene). Requires OMIM_API_KEY.

```text
dnakit.databases.omim_get_entry(*, operation=None, mim_number, include=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_entry) {&quot;const&quot;: &quot;get_entry&quot;, &quot;default&quot;: &quot;get_entry&quot;} |
| `mim_number` | string | Required by input schema | OMIM MIM number (e.g., 164730, 219700). Can include OMIM: prefix. |
| `include` | string | None | Data sections to include (default: text,clinicalSynopsis,geneMap) {&quot;default&quot;: &quot;text,clinicalSynopsis,geneMap&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_entry",
          "description": "Operation type (fixed: get_entry)",
          "default": "get_entry"
        },
        "mim_number": {
          "type": "string",
          "description": "OMIM MIM number (e.g., 164730, 219700). Can include OMIM: prefix."
        },
        "include": {
          "type": "string",
          "description": "Data sections to include (default: text,clinicalSynopsis,geneMap)",
          "default": "text,clinicalSynopsis,geneMap"
        }
      },
      "required": [
        "mim_number"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "mimNumber": {
              "type": "integer",
              "description": "OMIM MIM number"
            },
            "titles": {
              "type": "object",
              "description": "Entry titles"
            },
            "textSectionList": {
              "type": "array",
              "description": "Text sections with descriptions"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OMIM_get_entry`。

**Provider credentials:** `OMIM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("omim_get_entry", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.omim_get_gene_map` {#dnakit.databases.omim_get_gene_map}

Get OMIM gene-disease mapping information. Returns genes mapped to chromosomal locations with associated phenotypes and inheritance patterns. Query by MIM number or chromosome. Useful for finding all disease genes on a chromosome.

```text
dnakit.databases.omim_get_gene_map(*, operation=None, mim_number=None, chromosome=None, limit=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_gene_map) {&quot;const&quot;: &quot;get_gene_map&quot;, &quot;default&quot;: &quot;get_gene_map&quot;} |
| `mim_number` | string | None | OMIM MIM number for specific gene (optional if chromosome provided) |
| `chromosome` | string | None | Chromosome number (1-22, X, Y) to get all genes (optional if mim_number provided) |
| `limit` | integer | None | Maximum results for chromosome query (default: 50, max: 100) {&quot;default&quot;: 50} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_gene_map",
          "description": "Operation type (fixed: get_gene_map)",
          "default": "get_gene_map"
        },
        "mim_number": {
          "type": "string",
          "description": "OMIM MIM number for specific gene (optional if chromosome provided)"
        },
        "chromosome": {
          "type": "string",
          "description": "Chromosome number (1-22, X, Y) to get all genes (optional if mim_number provided)"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results for chromosome query (default: 50, max: 100)",
          "default": 50
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "gene_maps": {
              "type": "array",
              "description": "List of gene-disease mappings"
            },
            "count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OMIM_get_gene_map`。

**Provider credentials:** `OMIM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("omim_get_gene_map", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.omim_search` {#dnakit.databases.omim_search}

Search OMIM database for Mendelian diseases and genes. OMIM is the authoritative source for human genetic disorders. Search by gene name, disease name, or phenotype keywords. Returns MIM numbers, titles, and entry types. Requires OMIM_API_KEY env var (register at omim.org/api).

```text
dnakit.databases.omim_search(*, operation=None, query, limit=None, start=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search) {&quot;const&quot;: &quot;search&quot;, &quot;default&quot;: &quot;search&quot;} |
| `query` | string | Required by input schema | Search query - gene name (BRCA1), disease (cystic fibrosis), or phenotype (intellectual disability) |
| `limit` | integer | None | Maximum results (default: 10, max: 20) {&quot;default&quot;: 10, &quot;minimum&quot;: 1, &quot;maximum&quot;: 20} |
| `start` | integer | None | Offset for pagination (default: 0) {&quot;default&quot;: 0} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search",
          "description": "Operation type (fixed: search)",
          "default": "search"
        },
        "query": {
          "type": "string",
          "description": "Search query - gene name (BRCA1), disease (cystic fibrosis), or phenotype (intellectual disability)"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum results (default: 10, max: 20)",
          "default": 10,
          "minimum": 1,
          "maximum": 20
        },
        "start": {
          "type": "integer",
          "description": "Offset for pagination (default: 0)",
          "default": 0
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string",
          "enum": [
            "success",
            "error"
          ]
        },
        "data": {
          "type": "object",
          "properties": {
            "total_results": {
              "type": "integer",
              "description": "Total matching entries"
            },
            "entries": {
              "type": "array",
              "description": "List of matching OMIM entries"
            }
          }
        },
        "error": {
          "type": "string"
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `OMIM_search`。

**Provider credentials:** `OMIM_API_KEY`

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("omim_search", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_classification` {#dnakit.databases.orphanet_get_classification}

Get disease classification hierarchy from Orphanet. Returns parent and child disease categories showing where a disease fits in the rare disease taxonomy. Useful for understanding disease relationships.

```text
dnakit.databases.orphanet_get_classification(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_classification) {&quot;const&quot;: &quot;get_classification&quot;, &quot;default&quot;: &quot;get_classification&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_classification",
          "description": "Operation type (fixed: get_classification)",
          "default": "get_classification"
        },
        "orpha_code": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Orphanet ORPHA code"
        },
        "lang": {
          "type": "string",
          "description": "Language code (default: en)",
          "default": "en"
        },
        "orphacode": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        },
        "orpha_id": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "orpha_code": {
              "type": "string"
            },
            "parent": {
              "type": [
                "object",
                "null"
              ],
              "description": "Preferential parent disease/category in the rare disease taxonomy tree, or null if this is a top-level category."
            },
            "children": {
              "type": "array",
              "description": "Child diseases/categories in the rare disease taxonomy tree. Empty for leaf diseases."
            },
            "member_of_classifications": {
              "type": "object",
              "description": "Which named Orphanet classification systems (e.g. 'Orphanet classification of rare respiratory diseases') this disease is filed under. This is classification-system MEMBERSHIP, not the parent/child taxonomy tree -- see parent/children for that."
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_classification`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_classification", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_disease` {#dnakit.databases.orphanet_get_disease}

Get detailed information about a rare disease by ORPHA code. Returns disease name, definition, synonyms, and classification. ORPHA codes are unique Orphanet identifiers (e.g., 558 for Marfan syndrome). No authentication required.

```text
dnakit.databases.orphanet_get_disease(*, operation=None, orpha_code=None, lang=None, orphacode=None, orpha_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_disease) {&quot;const&quot;: &quot;get_disease&quot;, &quot;default&quot;: &quot;get_disease&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558, 166024). Can include ORPHA: prefix. |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_disease",
          "description": "Operation type (fixed: get_disease)",
          "default": "get_disease"
        },
        "orpha_code": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Orphanet ORPHA code (e.g., 558, 166024). Can include ORPHA: prefix."
        },
        "lang": {
          "type": "string",
          "description": "Language code (default: en)",
          "default": "en"
        },
        "orphacode": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        },
        "orpha_id": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "ORPHAcode": {
              "type": "string",
              "description": "Orphanet identifier"
            },
            "Preferred term": {
              "type": "string",
              "description": "Official disease name"
            },
            "Definition": {
              "type": "string",
              "description": "Disease definition"
            },
            "Synonyms": {
              "type": "array",
              "description": "Alternative names"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_disease`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_disease", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_disease_genes` {#dnakit.databases.orphanet_get_disease_genes}

Retrieve Orphanet's reported gene associations for an ORPHA disorder ID.

```text
dnakit.databases.orphanet_get_disease_genes(*, orpha_code, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `orpha_code` | string | Required by input schema | Orphanet ORPHA disorder identifier. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the Orphanet_get_genes association result; see [orphanet_get_genes](phenotypes.md#dnakit.databases.orphanet_get_genes).

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "orpha_code": {
          "type": "string"
        }
      },
      "required": [
        "orpha_code"
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_disease_genes`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_disease_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_epidemiology` {#dnakit.databases.orphanet_get_epidemiology}

Get epidemiology data (prevalence, incidence) for a rare disease from Orphanet. Returns prevalence estimates by geographic region and type (point prevalence, birth prevalence, annual incidence). Data includes validation status and literature sources.

```text
dnakit.databases.orphanet_get_epidemiology(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_epidemiology) {&quot;const&quot;: &quot;get_epidemiology&quot;, &quot;default&quot;: &quot;get_epidemiology&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis) |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_epidemiology",
          "description": "Operation type (fixed: get_epidemiology)",
          "default": "get_epidemiology"
        },
        "orpha_code": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis)"
        },
        "orphacode": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        },
        "orpha_id": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "orpha_code": {
              "type": "string"
            },
            "preferred_term": {
              "type": "string"
            },
            "prevalences": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "type": {
                    "type": "string"
                  },
                  "class": {
                    "type": "string"
                  },
                  "geographic": {
                    "type": "string"
                  },
                  "mean_value": {
                    "type": "string"
                  },
                  "source": {
                    "type": "string"
                  },
                  "validation_status": {
                    "type": "string"
                  }
                }
              }
            },
            "prevalence_count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_epidemiology`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_epidemiology", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_gene_diseases` {#dnakit.databases.orphanet_get_gene_diseases}

Get rare diseases associated with a gene from Orphanet. Accepts gene symbols (e.g., 'FBN1', 'BRCA1', 'HTT') or gene name keywords (e.g., 'fibrillin', 'huntingtin', 'dystrophin'). Gene symbols are auto-resolved to full names. Returns disease ORPHA codes, gene symbols, association types, and genomic locus information.

```text
dnakit.databases.orphanet_get_gene_diseases(*, operation=None, gene_name=None, gene_symbol=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_gene_diseases) {&quot;const&quot;: &quot;get_gene_diseases&quot;, &quot;default&quot;: &quot;get_gene_diseases&quot;} |
| `gene_name` | string | None | Gene symbol or name keyword (e.g., &#x27;FBN1&#x27;, &#x27;BRCA1&#x27;, &#x27;huntingtin&#x27;). Use gene_symbol for standard HGNC symbols. Note: Orphanet uses current HGNC symbols (e.g., &#x27;GBA1&#x27; not &#x27;GBA&#x27;). |
| `gene_symbol` | string &#124; null | None | Alternative to gene_name: standard HGNC gene symbol (e.g., &#x27;GBA1&#x27;, &#x27;FBN1&#x27;). Note: Orphanet uses current HGNC symbols (e.g., &#x27;GBA1&#x27; not &#x27;GBA&#x27;). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_gene_diseases",
          "description": "Operation type (fixed: get_gene_diseases)",
          "default": "get_gene_diseases"
        },
        "gene_name": {
          "type": "string",
          "description": "Gene symbol or name keyword (e.g., 'FBN1', 'BRCA1', 'huntingtin'). Use gene_symbol for standard HGNC symbols. Note: Orphanet uses current HGNC symbols (e.g., 'GBA1' not 'GBA')."
        },
        "gene_symbol": {
          "type": [
            "string",
            "null"
          ],
          "description": "Alternative to gene_name: standard HGNC gene symbol (e.g., 'GBA1', 'FBN1'). Note: Orphanet uses current HGNC symbols (e.g., 'GBA1' not 'GBA')."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "gene_name": {
              "type": "string"
            },
            "diseases": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "orpha_code": {
                    "type": [
                      "string",
                      "integer"
                    ]
                  },
                  "preferred_term": {
                    "type": "string"
                  },
                  "disorder_group": {
                    "type": "string"
                  },
                  "typology": {
                    "type": "string"
                  },
                  "genes": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "symbol": {
                          "type": "string"
                        },
                        "name": {
                          "type": "string"
                        },
                        "association_type": {
                          "type": "string"
                        }
                      }
                    }
                  }
                }
              }
            },
            "disease_count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_gene_diseases`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_gene_diseases", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_genes` {#dnakit.databases.orphanet_get_genes}

Get genes associated with a rare disease from Orphanet. Returns gene symbols, names, and association types (causative, modifier, susceptibility). Essential for rare disease genetic diagnosis.

```text
dnakit.databases.orphanet_get_genes(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_genes) {&quot;const&quot;: &quot;get_genes&quot;, &quot;default&quot;: &quot;get_genes&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code for the disease |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_genes",
          "description": "Operation type (fixed: get_genes)",
          "default": "get_genes"
        },
        "orpha_code": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Orphanet ORPHA code for the disease"
        },
        "orphacode": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        },
        "orpha_id": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "orpha_code": {
              "type": "string"
            },
            "genes": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "Symbol": {
                    "type": "string",
                    "description": "Gene symbol"
                  },
                  "Name": {
                    "type": "string",
                    "description": "Gene name"
                  },
                  "AssociationType": {
                    "type": "string",
                    "description": "Type of gene-disease association"
                  }
                }
              }
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_genes`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_genes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_icd_mapping` {#dnakit.databases.orphanet_get_icd_mapping}

Get cross-references between Orphanet ORPHA codes and other medical coding systems (ICD-10, ICD-11, OMIM, SNOMED-CT). Essential for mapping rare disease codes across clinical and research nomenclatures. Returns mapping relationships and validation status.

```text
dnakit.databases.orphanet_get_icd_mapping(*, operation=None, orpha_code=None, coding_system=None, lang=None, orphacode=None, orpha_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_icd_mapping) {&quot;const&quot;: &quot;get_icd_mapping&quot;, &quot;default&quot;: &quot;get_icd_mapping&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis) |
| `coding_system` | string | None | Which coding system to retrieve: all, icd10, icd11, omim, snomed. Default: all {&quot;default&quot;: &quot;all&quot;} |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_icd_mapping",
          "description": "Operation type (fixed: get_icd_mapping)",
          "default": "get_icd_mapping"
        },
        "orpha_code": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis)"
        },
        "coding_system": {
          "type": "string",
          "description": "Which coding system to retrieve: all, icd10, icd11, omim, snomed. Default: all",
          "default": "all"
        },
        "lang": {
          "type": "string",
          "description": "Language code (default: en)",
          "default": "en"
        },
        "orphacode": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        },
        "orpha_id": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "orpha_code": {
              "type": "string"
            },
            "mappings": {
              "type": "object",
              "description": "Cross-reference mappings keyed by coding system (ICD10, ICD11, OMIM, SNOMED-CT)"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_icd_mapping`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_icd_mapping", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_natural_history` {#dnakit.databases.orphanet_get_natural_history}

Get natural history data for a rare disease from Orphanet. Returns average age of onset and type of inheritance. Useful for understanding disease progression patterns and genetic counseling.

```text
dnakit.databases.orphanet_get_natural_history(*, operation=None, orpha_code=None, orphacode=None, orpha_id=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_natural_history) {&quot;const&quot;: &quot;get_natural_history&quot;, &quot;default&quot;: &quot;get_natural_history&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 399 for Huntington) |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `orpha_id` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_natural_history",
          "description": "Operation type (fixed: get_natural_history)",
          "default": "get_natural_history"
        },
        "orpha_code": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Orphanet ORPHA code (e.g., 558 for Marfan, 399 for Huntington)"
        },
        "orphacode": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        },
        "orpha_id": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "orpha_code": {
              "type": "string"
            },
            "preferred_term": {
              "type": "string"
            },
            "average_age_of_onset": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "type_of_inheritance": {
              "type": "array",
              "items": {
                "type": "string"
              }
            },
            "disorder_group": {
              "type": "string"
            },
            "typology": {
              "type": "string"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_natural_history`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_natural_history", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_get_phenotypes` {#dnakit.databases.orphanet_get_phenotypes}

Get HPO (Human Phenotype Ontology) phenotypes associated with a rare disease from Orphanet. Returns HPO terms with their frequency (Very frequent, Frequent, Occasional, Very rare, Excluded) and whether they are diagnostic criteria. Essential for phenotype-driven rare disease diagnosis.

```text
dnakit.databases.orphanet_get_phenotypes(*, operation=None, orpha_code=None, orpha_id=None, orphacode=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: get_phenotypes) {&quot;const&quot;: &quot;get_phenotypes&quot;, &quot;default&quot;: &quot;get_phenotypes&quot;} |
| `orpha_code` | string &#124; integer | None | Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis) |
| `orpha_id` | string &#124; integer | None | Alias for orpha_code. Orphanet ORPHA code (e.g., 558 for Marfan) |
| `orphacode` | integer &#124; string | None | Alias for orpha_code (the Orphanet disease code). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "get_phenotypes",
          "description": "Operation type (fixed: get_phenotypes)",
          "default": "get_phenotypes"
        },
        "orpha_code": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Orphanet ORPHA code (e.g., 558 for Marfan, 586 for Cystic Fibrosis)"
        },
        "orpha_id": {
          "type": [
            "string",
            "integer"
          ],
          "description": "Alias for orpha_code. Orphanet ORPHA code (e.g., 558 for Marfan)"
        },
        "orphacode": {
          "type": [
            "integer",
            "string"
          ],
          "description": "Alias for orpha_code (the Orphanet disease code)."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "orpha_code": {
              "type": "string"
            },
            "preferred_term": {
              "type": "string"
            },
            "phenotypes": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "hpo_id": {
                    "type": "string"
                  },
                  "hpo_term": {
                    "type": "string"
                  },
                  "frequency": {
                    "type": "string"
                  },
                  "diagnostic_criteria": {
                    "type": [
                      "string",
                      "null"
                    ]
                  }
                }
              }
            },
            "phenotype_count": {
              "type": "integer"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_get_phenotypes`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_get_phenotypes", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_search_by_name` {#dnakit.databases.orphanet_search_by_name}

Search for rare diseases by exact or partial name match in Orphanet. More precise than keyword search. Returns ORPHA codes and disease names. Supports exact matching mode for disambiguation.

```text
dnakit.databases.orphanet_search_by_name(*, operation=None, name=None, exact=None, lang=None, limit=None, query=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_by_name) {&quot;const&quot;: &quot;search_by_name&quot;, &quot;default&quot;: &quot;search_by_name&quot;} |
| `name` | string | None | Disease name to search (e.g., &#x27;Marfan syndrome&#x27;) |
| `exact` | boolean | None | If true, match exact name only. Default: false {&quot;default&quot;: false} |
| `lang` | string | None | Language code (default: en) {&quot;default&quot;: &quot;en&quot;} |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 200) |
| `query` | string | None | Alias for name. Disease name to search for in Orphanet. |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search_by_name",
          "description": "Operation type (fixed: search_by_name)",
          "default": "search_by_name"
        },
        "name": {
          "type": "string",
          "description": "Disease name to search (e.g., 'Marfan syndrome')"
        },
        "exact": {
          "type": "boolean",
          "description": "If true, match exact name only. Default: false",
          "default": false
        },
        "lang": {
          "type": "string",
          "description": "Language code (default: en)",
          "default": "en"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 20, max: 200)"
        },
        "query": {
          "type": "string",
          "description": "Alias for name. Disease name to search for in Orphanet."
        }
      },
      "required": []
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string"
        },
        "data": {
          "type": "object",
          "properties": {
            "results": {
              "type": "array"
            },
            "count": {
              "type": "integer"
            },
            "exact_match": {
              "type": "boolean"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_search_by_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_search_by_name", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.orphanet_search_diseases` {#dnakit.databases.orphanet_search_diseases}

Search Orphanet for rare diseases by keyword. Orphanet is the reference portal for rare diseases and orphan drugs. Returns disease names, ORPHA codes, and synonyms. No authentication required. Supports multiple languages (en, fr, de, es, it, pt, pl, nl).

```text
dnakit.databases.orphanet_search_diseases(*, operation=None, query, limit=None, lang=None, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `operation` | string | None | Operation type (fixed: search_diseases) {&quot;const&quot;: &quot;search_diseases&quot;, &quot;default&quot;: &quot;search_diseases&quot;} |
| `query` | string | Required by input schema | Search query - disease name or keyword (e.g., &#x27;Marfan&#x27;, &#x27;muscular dystrophy&#x27;) |
| `limit` | integer | None | Maximum number of results to return (default: 20, max: 200). The API may return thousands of fuzzy matches; use this to control output size. {&quot;default&quot;: 20} |
| `lang` | string | None | Language code (en, fr, de, es, it, pt, pl, nl). Default: en {&quot;default&quot;: &quot;en&quot;} |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "operation": {
          "type": "string",
          "const": "search_diseases",
          "description": "Operation type (fixed: search_diseases)",
          "default": "search_diseases"
        },
        "query": {
          "type": "string",
          "description": "Search query - disease name or keyword (e.g., 'Marfan', 'muscular dystrophy')"
        },
        "limit": {
          "type": "integer",
          "description": "Maximum number of results to return (default: 20, max: 200). The API may return thousands of fuzzy matches; use this to control output size.",
          "default": 20
        },
        "lang": {
          "type": "string",
          "description": "Language code (en, fr, de, es, it, pt, pl, nl). Default: en",
          "default": "en"
        }
      },
      "required": [
        "query"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "type": "object",
      "properties": {
        "status": {
          "type": "string",
          "enum": [
            "success",
            "error"
          ]
        },
        "data": {
          "type": "object",
          "properties": {
            "results": {
              "type": "array",
              "description": "Matching rare diseases"
            },
            "count": {
              "type": "integer",
              "description": "Number of results returned"
            },
            "total_count": {
              "type": "integer",
              "description": "Total matches before truncation"
            },
            "query": {
              "type": "string"
            }
          }
        }
      }
    }
    ```

**Execution and provenance:** `network`；upstream name `Orphanet_search_diseases`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("orphanet_search_diseases", verbose=False)
print(info["input_schema"])
```

## `dnakit.databases.osl_get_efo_id_by_disease_name` {#dnakit.databases.osl_get_efo_id_by_disease_name}

Legacy helper to find an EFO term ID from a disease name using EMBL-EBI OLS search.

Returns the best match as `{efo_id, name}` (or `null` if none). Note: the tool name contains a legacy typo (`OSL_` vs `OLS_`). If you need richer ontology data (IRI, synonyms, children), use the newer `ols_*` tools in `efo_tools.json`.

```text
dnakit.databases.osl_get_efo_id_by_disease_name(*, disease, verbose=True)
```

| Parameter | Type | Required / default | Description and constraints |
| --- | --- | --- | --- |
| `disease` | string | Required by input schema | Disease name or free-text query (e.g., &#x27;diabetes mellitus&#x27;). |
| `verbose` | bool | True | Progress flag. |

**Returns**

Returns the provider result after DNAKit argument handling, error checks and size limits. The registered provider result schema below is a structural reference, not a live response. Data envelopes and pagination vary by operation.

??? info "Complete input schema, aliases and conditional constraints"

    ```json
    {
      "type": "object",
      "properties": {
        "disease": {
          "type": "string",
          "description": "Disease name or free-text query (e.g., 'diabetes mellitus')."
        }
      },
      "required": [
        "disease"
      ]
    }
    ```

??? info "Registered provider result schema"

    ```json
    {
      "description": "Either a best-match `{efo_id,name}` result, an error object, or null when no results are found.",
      "anyOf": [
        {
          "type": "object",
          "properties": {
            "efo_id": {
              "type": "string"
            },
            "name": {
              "type": "string"
            }
          },
          "required": [
            "efo_id",
            "name"
          ],
          "additionalProperties": true
        },
        {
          "type": "object",
          "properties": {
            "error": {
              "type": "string"
            },
            "details": {
              "type": "string"
            }
          },
          "required": [
            "error"
          ],
          "additionalProperties": true
        },
        {
          "type": "null"
        }
      ]
    }
    ```

**Execution and provenance:** `network`；upstream name `OSL_get_efo_id_by_disease_name`。

**Inspect parameters offline without calling the provider:**

```python
import dnakit
info = dnakit.databases.describe("osl_get_efo_id_by_disease_name", verbose=False)
print(info["input_schema"])
```

